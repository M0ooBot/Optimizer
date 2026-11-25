"""Verify wing.step (reimported) against the close-rounded-wingtip skill.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/flat/blunt/verify_cap.py

Writes readback.json next to this script. All geometry is read back from
wing.step: B-spline faces are evaluated exactly from their STEP poles and knots
(cross-checked against Open CASCADE's own D2 evaluation).
"""
import hashlib, json, os, sys, time
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize, brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import design as D  # noqa: E402  (only for the parent definition and constants)
PP = D.PP

from OCP.STEPControl import STEPControl_Reader
from OCP.TopExp import TopExp_Explorer, TopExp
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID, TopAbs_SHELL, TopAbs_EDGE
from OCP.TopoDS import TopoDS
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Plane
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds
from OCP.gp import gp_Pnt, gp_Vec
from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape

STEP = os.path.join(HERE, 'wing.step')
PARENT_STEP = os.path.join(D.PARENT_DIR, 'parent.step')
OUT = os.path.join(HERE, 'readback.json')
RULE_T, RULE_R = 0.8, 0.4


# --------------------------------------------------------------------------- helpers
def read(path):
    r = STEPControl_Reader(); assert r.ReadFile(path) == 1; r.TransferRoots(); return r.OneShape()


def explore(shape, kind, cast):
    out, ex = [], TopExp_Explorer(shape, kind)
    while ex.More():
        out.append(cast(ex.Current())); ex.Next()
    return out


class Surf:
    """Exact evaluation of a non-rational Geom_BSplineSurface from its poles and knots."""

    def __init__(self, bs, face=None):
        assert not bs.IsURational() and not bs.IsVRational()
        self.bs, self.face = bs, face
        self.ku, self.kv = bs.UDegree(), bs.VDegree()
        nu, nv = bs.NbUPoles(), bs.NbVPoles()
        self.P = np.array([[[bs.Pole(i, j).X(), bs.Pole(i, j).Y(), bs.Pole(i, j).Z()] for j in range(1, nv + 1)] for i in range(1, nu + 1)])
        self.tu = np.repeat([bs.UKnot(i) for i in range(1, bs.NbUKnots() + 1)], [bs.UMultiplicity(i) for i in range(1, bs.NbUKnots() + 1)])
        self.tv = np.repeat([bs.VKnot(i) for i in range(1, bs.NbVKnots() + 1)], [bs.VMultiplicity(i) for i in range(1, bs.NbVKnots() + 1)])
        self.nu, self.nv = nu, nv
        self.su = [BSpline(self.tu, self.P.reshape(nu, -1), self.ku)]
        zero = lambda spl, k, n: (spl.derivative(k) if k <= spl.k else (lambda x: np.zeros((len(np.atleast_1d(x)), n))))
        self.su += [zero(self.su[0], 1, nv * 3), zero(self.su[0], 2, nv * 3)]
        eye = np.eye(nv)
        self.sv = [BSpline(self.tv, eye, self.kv)]
        self.sv += [zero(self.sv[0], 1, nv), zero(self.sv[0], 2, nv)]
        self.urange = (self.tu[0], self.tu[-1]); self.vrange = (self.tv[0], self.tv[-1])

    def grid(self, u, v):
        """Derivatives on a grid: dict with keys P, Pu, Pv, Puu, Puv, Pvv, shape (len(u), len(v), 3)."""
        u = np.clip(np.asarray(u, float), *self.urange); v = np.clip(np.asarray(v, float), *self.vrange)
        A = [s(u).reshape(len(u), self.nv, 3) for s in self.su]
        Bv = [s(v) for s in self.sv]
        g = lambda a, b: np.einsum('ajc,bj->abc', A[a], Bv[b])
        return {'P': g(0, 0), 'Pu': g(1, 0), 'Pv': g(0, 1), 'Puu': g(2, 0), 'Puv': g(1, 1), 'Pvv': g(0, 2)}

    def at(self, u, v):
        """Derivatives at scattered points (arrays of equal length)."""
        u = np.atleast_1d(np.asarray(u, float)); v = np.atleast_1d(np.asarray(v, float))
        u = np.clip(u, *self.urange); v = np.clip(v, *self.vrange)
        A = [s(u).reshape(len(u), self.nv, 3) for s in self.su]
        Bv = [s(v) for s in self.sv]
        g = lambda a, b: np.einsum('ajc,aj->ac', A[a], Bv[b])
        return {'P': g(0, 0), 'Pu': g(1, 0), 'Pv': g(0, 1), 'Puu': g(2, 0), 'Puv': g(1, 1), 'Pvv': g(0, 2)}


def geometry(d, flip=1.0):
    """Unit normal, fundamental forms, principal curvatures, normal derivatives. flip orients n outward."""
    Pu, Pv = d['Pu'], d['Pv']
    n = np.cross(Pu, Pv); n = flip * n / np.linalg.norm(n, axis=-1, keepdims=True)
    E = (Pu * Pu).sum(-1); F = (Pu * Pv).sum(-1); G = (Pv * Pv).sum(-1)
    L = (d['Puu'] * n).sum(-1); M = (d['Puv'] * n).sum(-1); N = (d['Pvv'] * n).sum(-1)
    det = E * G - F * F
    H = (E * N - 2 * F * M + G * L) / (2 * det); K = (L * N - M * M) / det
    disc = np.sqrt(np.maximum(H * H - K, 0))
    k1, k2 = H + disc, H - disc
    nu_ = ((F * M - G * L)[..., None] * Pu + (F * L - E * M)[..., None] * Pv) / det[..., None]
    nv_ = ((F * N - G * M)[..., None] * Pu + (F * M - E * N)[..., None] * Pv) / det[..., None]
    return dict(n=n, E=E, F=F, G=G, L=L, M=M, N=N, det=det, k1=k1, k2=k2,
                kmax=np.maximum(np.abs(k1), np.abs(k2)), nu=nu_, nv=nv_)


def II(g, d, t):
    """Normal curvature II(t, t)/I(t, t) for 3D tangent vectors t (shape (..., 3))."""
    # express t in (Pu, Pv): solve [E F; F G][a b] = [t.Pu, t.Pv]
    tu = (t * d['Pu']).sum(-1); tv = (t * d['Pv']).sum(-1)
    a = (g['G'] * tu - g['F'] * tv) / g['det']; b = (g['E'] * tv - g['F'] * tu) / g['det']
    num = g['L'] * a * a + 2 * g['M'] * a * b + g['N'] * b * b
    den = g['E'] * a * a + 2 * g['F'] * a * b + g['G'] * b * b
    return num / den, a, b


def bilinear(g, a1, b1, a2, b2):
    return g['L'] * a1 * a2 + g['M'] * (a1 * b2 + a2 * b1) + g['N'] * b1 * b2


def thickness(X):
    X = np.asarray(X, float)
    t = PP.section_thickness(np.clip(X / D.CHORD, 0, 1), D.ATT)
    return np.where(np.isnan(t), D.H_TE * 2, t)


def required_radius(X):
    t = thickness(X)
    return np.where(t >= RULE_T, RULE_R, t / 2), t


# --------------------------------------------------------------------------- main
def main():
    T0 = time.time()
    out = {}
    raw = open(STEP, 'rb').read()
    out['step'] = {'path': os.path.relpath(STEP), 'sha256': hashlib.sha256(raw).hexdigest(), 'bytes': len(raw)}
    shape = read(STEP)
    solids = explore(shape, TopAbs_SOLID, TopoDS.Solid_s)
    shells = explore(shape, TopAbs_SHELL, TopoDS.Shell_s)
    fl = explore(shape, TopAbs_FACE, TopoDS.Face_s)
    fb = ShapeAnalysis_FreeBounds(shape, 1e-7, False, False)
    from OCP.TopoDS import TopoDS_Compound
    def nedges(s):
        return len(explore(s, TopAbs_EDGE, TopoDS.Edge_s)) if s is not None and not s.IsNull() else 0
    props = GProp_GProps(); BRepGProp.VolumeProperties_s(shape, props, 1e-9)
    sprops = GProp_GProps(); BRepGProp.SurfaceProperties_s(shape, sprops, 1e-9)
    box = Bnd_Box(); BRepBndLib.AddOptimal_s(shape, box, False, False)
    bb = box.Get()
    an = BRepCheck_Analyzer(shape)
    from OCP.BRep import BRep_Tool as BT
    out['solid'] = {
        'solids': len(solids), 'shells': len(shells), 'faces': len(fl),
        'brepcheck_valid': bool(an.IsValid()),
        'free_boundary_edges': nedges(fb.GetClosedWires()) + nedges(fb.GetOpenWires()),
        'shell_closed_flag': bool(shells[0].Closed()) if shells else None,
        'volume_mm3': props.Mass(), 'surface_area_mm2': sprops.Mass(),
        'bounding_box_mm': [round(x, 9) for x in bb],
        'face_types': [str(BRepAdaptor_Surface(f).GetType()).split('.')[-1] for f in fl],
    }
    # max tolerance of vertices/edges
    from OCP.ShapeAnalysis import ShapeAnalysis_ShapeTolerance
    stol = ShapeAnalysis_ShapeTolerance()
    out['solid']['max_shape_tolerance_mm'] = stol.Tolerance(shape, 1)

    # ---- classify faces
    caps, parent_skin, te, closures = {}, {}, None, {}
    for f in fl:
        ad = BRepAdaptor_Surface(f)
        if ad.GetType() == GeomAbs_Plane:
            c = ad.Plane().Location(); bx = Bnd_Box(); BRepBndLib.Add_s(f, bx); y = 0.5 * (bx.Get()[1] + bx.Get()[4])
            closures['right' if y > 0 else 'left'] = f
            continue
        bs = ad.BSpline()
        if abs(bs.VKnot(1) + 1) < 1e-12 and abs(bs.VKnot(bs.NbVKnots()) - 1) < 1e-12 and bs.UDegree() == 5:
            S = Surf(bs, f); y = S.P[..., 1].mean()
            caps['right' if y > 0 else 'left'] = S
        elif bs.UDegree() == 5:
            S = Surf(bs, f); z = S.P[..., 2].mean()
            parent_skin['upper' if z > 0 else 'lower'] = S
        else:
            te = Surf(bs, f)
    assert set(caps) == {'right', 'left'} and set(parent_skin) == {'upper', 'lower'} and te is not None and len(closures) == 2
    # orientation: outward normal flip per cap (check at the lip, mid chord)
    flips = {}
    for side, S in caps.items():
        d = S.at([0.5], [0.0]); n = np.cross(d['Pu'], d['Pv'])[0]
        flips[side] = 1.0 if n[1] * (1 if side == 'right' else -1) > 0 else -1.0
    pflip = {}
    for k, S in parent_skin.items():
        d = S.at([0.5], [0.0]); n = np.cross(d['Pu'], d['Pv'])[0]
        pflip[k] = 1.0 if n[2] * (1 if k == 'upper' else -1) > 0 else -1.0

    # cross-check exact evaluation with OCC D2
    ad = BRepAdaptor_Surface(caps['right'].face)
    errs = []
    for (u, v) in [(0.3, 0.2), (0.77, -0.6), (0.999, 0.95), (0.05, 0.0)]:
        p, d1u, d1v, d2u, d2v, d2uv = gp_Pnt(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec()
        ad.D2(u, v, p, d1u, d1v, d2u, d2v, d2uv)
        dd = caps['right'].at([u], [v])
        occ = [p, d1u, d1v, d2u, d2uv, d2v]
        mine = [dd['P'][0], dd['Pu'][0], dd['Pv'][0], dd['Puu'][0], dd['Puv'][0], dd['Pvv'][0]]
        errs.append(max(np.abs(np.array([o.X(), o.Y(), o.Z()]) - m).max() / max(1, np.abs(m).max()) for o, m in zip(occ, mine)))
    out['evaluation_crosscheck_max_relative_difference_vs_OCC_D2'] = float(max(errs))

    # ---- parent preserved
    pstep = read(PARENT_STEP)
    pfaces = {}
    for f in explore(pstep, TopAbs_FACE, TopoDS.Face_s):
        a = BRepAdaptor_Surface(f)
        if a.GetType() == GeomAbs_BSplineSurface and a.BSpline().UDegree() == 5:
            S = Surf(a.BSpline()); pfaces['upper' if S.P[..., 2].mean() > 0 else 'lower'] = S
        elif a.GetType() == GeomAbs_BSplineSurface:
            pfaces['te'] = Surf(a.BSpline())
    sg = np.unique(np.concatenate([np.linspace(0, 1, 801), np.geomspace(1e-6, 0.02, 50), 1 - np.geomspace(1e-6, 0.02, 50)]))
    yg = np.linspace(-350, 350, 57)
    pres = {}
    for k in ('upper', 'lower'):
        S = parent_skin[k]
        P = S.grid(sg, yg)['P']
        ref = PP.skin(sg[:, None] * np.ones_like(yg)[None], yg[None] * np.ones_like(sg)[:, None], +1 if k == 'upper' else -1)
        pres[k + '_vs_parent_py_max_mm'] = float(np.linalg.norm(P - ref, axis=-1).max())
        pres[k + '_poles_vs_parent_step_max_mm'] = float(np.abs(S.P - pfaces[k].P).max())
    tt = np.linspace(0, 1, 41)
    Pte = te.grid(yg, tt)['P']; ref = PP.te_face(tt[None] * np.ones_like(yg)[:, None], yg[:, None] * np.ones_like(tt)[None])
    pres['te_face_vs_parent_py_max_mm'] = float(np.linalg.norm(Pte - ref, axis=-1).max())
    pres['te_face_poles_vs_parent_step_max_mm'] = float(np.abs(te.P - pfaces['te'].P).max())
    out['parent_preserved'] = pres

    # ---- locked dimensions and outline
    lock = {}
    lock['chord_X_range_mm'] = [bb[0], bb[3]]
    lock['span_Y_range_mm'] = [bb[1], bb[4]]
    lock['projected_span_mm'] = bb[4] - bb[1]
    lock['thickness_Z_range_mm'] = [bb[2], bb[5]]
    # attachment: the cap's cut-edge rows lie at |Y| = 350
    lock['cap_cut_edge_max_abs_Y_minus_350_mm'] = float(max(np.abs(np.abs(caps[s].P[:, [0, -1], 1]) - 350).max() for s in caps))
    # outline: the lip (w = 0) against the locked tip outline (distance to the curve, refined)
    pp = np.unique(np.concatenate([np.linspace(0, 1, 400001), np.geomspace(1e-9, 1e-2, 40000), 1 - np.geomspace(1e-9, 1e-2, 40000)]))
    from scipy.spatial import cKDTree
    oc = PP.tip_outline(pp)[:, :2]; tree = cKDTree(oc)
    sl = np.unique(np.concatenate([np.linspace(0, 1, 40001), 1 - np.geomspace(1e-10, 0.1, 8000), np.geomspace(1e-10, 0.05, 2000)]))
    od = {}
    for side, S in caps.items():
        lip = S.at(sl, np.zeros_like(sl))['P']
        q = lip[:, :2] * np.array([1, 1 if side == 'right' else -1])
        dist, idx = tree.query(q)
        # distance to the polyline through the outline samples (segments either side of the nearest sample)
        def segd(k0, k1):
            a = oc[np.clip(k0, 0, len(oc) - 1)]; b = oc[np.clip(k1, 0, len(oc) - 1)]
            ab = b - a; t_ = np.clip(((q - a) * ab).sum(-1) / np.maximum((ab * ab).sum(-1), 1e-300), 0, 1)
            return np.linalg.norm(q - (a + t_[:, None] * ab), axis=-1)
        dist = np.minimum(segd(idx - 1, idx), segd(idx, idx + 1))
        # refine: exact distance to the curve near the nearest sample
        def refine(k):
            p0 = pp[idx[k]]; lo, hi = max(0, p0 - 2e-4), min(1, p0 + 2e-4)
            f = lambda p: np.sum((PP.tip_outline(np.array([p]))[0, :2] - q[k]) ** 2)
            from scipy.optimize import minimize_scalar
            r = minimize_scalar(f, bounds=(lo, hi), method='bounded', options={'xatol': 1e-14})
            return np.sqrt(r.fun), r.x
        top = np.argsort(-dist)[:40]
        best = max((refine(k) + (k,) for k in top), key=lambda x: x[0])
        mid = sl < 0.99
        kz = np.argmax(np.abs(lip[:, 2]))
        od[side] = {'lip_max_distance_to_locked_outline_mm': float(best[0]), 'at_s': float(sl[best[2]]),
                    'at_XY': lip[best[2], :2].tolist(),
                    'lip_max_distance_for_s_below_0.99_mm': float(dist[mid].max()),
                    'lip_max_abs_Z_mm': float(np.abs(lip[:, 2]).max()),
                    'tip_max_abs_Y_mm': float(np.abs(lip[:, 1]).max()),
                    'tip_at_X_mm': float(lip[np.argmax(np.abs(lip[:, 1])), 0])}
        # where does the lip leave the outline by more than 1e-5 mm?
        far = sl[dist > 1e-5]
        od[side]['s_range_where_distance_exceeds_1e-5_mm'] = [float(far.min()), float(far.max())] if len(far) else None
        if len(far):
            Xf = S.at(np.array([far.min()]), np.array([0.0]))['P'][0]
            od[side]['that_range_starts_at_XY'] = Xf[:2].tolist()
    lock['outline'] = od
    # closure face: half-round end of the flat trailing-edge face
    clo = {}
    for side, S in caps.items():
        w = np.linspace(-1, 1, 20001)
        d = S.at(np.ones_like(w), w)
        Pc = d['P']
        # apex radius in the face plane (planar curve in Y-Z)
        Yw, Zw = d['Pv'][:, 1], d['Pv'][:, 2]; Yww, Zww = d['Pvv'][:, 1], d['Pvv'][:, 2]
        kap = np.abs(Yw * Zww - Zw * Yww) / (Yw ** 2 + Zw ** 2) ** 1.5
        i0 = np.argmin(np.abs(w))
        # crease angle between cap and the X = 120 plane along the closure curve
        g = geometry(d, flips[side])
        crease = np.degrees(np.arccos(np.clip(g['n'][:, 0], -1, 1)))
        clo[side] = {'closure_curve_X_range_mm': [float(Pc[:, 0].min()), float(Pc[:, 0].max())],
                     'reach_beyond_attachment_mm': float(np.abs(np.abs(Pc[:, 1]) - 350).max()),
                     'height_mm': float(Pc[:, 2].max() - Pc[:, 2].min()),
                     'apex_radius_in_face_plane_mm': float(1 / kap[i0]),
                     'min_radius_of_closure_curve_in_face_plane_mm': float(1 / kap.max()),
                     'crease_angle_cap_to_TE_plane_deg': {'at_corners': [float(crease[0]), float(crease[-1])],
                                                          'at_apex': float(crease[i0]),
                                                          'min': float(crease.min()), 'max': float(crease.max())}}
    # the parent's own crease between skin and TE face, for reference
    d = parent_skin['upper'].at([1.0], [350.0]); g = geometry(d, pflip['upper'])
    clo['parent_skin_to_TE_face_crease_deg'] = float(np.degrees(np.arccos(np.clip(g['n'][0, 0], -1, 1))))
    lock['trailing_edge_closure'] = clo
    out['locked_dimensions'] = lock

    # ---- joins: cap to parent skin along the cut edges
    sj = np.unique(np.concatenate([np.linspace(1e-5, 1, 4001), np.geomspace(1e-6, 0.02, 200), 1 - np.geomspace(1e-7, 0.02, 200)]))
    joins = {}
    for side, S in caps.items():
        sgn = 1 if side == 'right' else -1
        for skin, wv in (('upper', 1.0), ('lower', -1.0)):
            dc = S.at(sj, np.full_like(sj, wv)); gc = geometry(dc, flips[side])
            dp = parent_skin[skin].at(sj, np.full_like(sj, sgn * 350.0)); gp = geometry(dp, pflip[skin])
            sep = np.linalg.norm(dc['P'] - dp['P'], axis=-1)
            ang = np.degrees(np.arccos(np.clip((gc['n'] * gp['n']).sum(-1), -1, 1)))
            # common tangent frame: t1 along the edge (parent chord direction), t2 = n x t1
            t1 = dp['Pu'] / np.linalg.norm(dp['Pu'], axis=-1, keepdims=True)
            t2 = np.cross(gp['n'], t1)
            res = {}
            for nm, (g_, d_) in (('cap', (gc, dc)), ('parent', (gp, dp))):
                _, a1, b1 = II(g_, d_, t1); _, a2, b2 = II(g_, d_, t2)
                # normalise so that a P_u + b P_v = t
                res[nm] = (bilinear(g_, a1, b1, a1, b1), bilinear(g_, a1, b1, a2, b2), bilinear(g_, a2, b2, a2, b2))
            dk = [np.abs(res['cap'][i] - res['parent'][i]) for i in range(3)]
            joins[f'{side}_{skin}'] = {
                'samples': len(sj), 's_range': [float(sj.min()), float(sj.max())],
                'nearest_sample_to_collapsed_LE_point_mm': float(np.linalg.norm(dp['P'][0] - np.array([0, sgn * 350, 0]))),
                'max_separation_mm': float(sep.max()),
                'max_tangent_plane_angle_deg': float(ang.max()), 'at_s': float(sj[np.argmax(ang)]),
                'max_abs_diff_normal_curvature_chordwise_per_mm': float(dk[0].max()),
                'max_abs_diff_twist_curvature_per_mm': float(dk[1].max()),
                'max_abs_diff_normal_curvature_spanwise_per_mm': float(dk[2].max()),
                'parent_curvature_range_chordwise_per_mm': [float(res['parent'][0].min()), float(res['parent'][0].max())],
                'cap_spanwise_curvature_max_abs_per_mm': float(np.abs(res['cap'][2]).max()),
            }
            # excluding the last 1e-4 of s near the LE degenerate point
            m = sj > 1e-3
            joins[f'{side}_{skin}']['max_tangent_plane_angle_deg_s_above_1e-3'] = float(ang[m].max())
    out['joins'] = joins

    # ---- curvature map over each cap face
    s_map = np.unique(np.concatenate([np.geomspace(1e-5, 0.02, 160), np.linspace(0.02, 0.98, 1201), 1 - np.geomspace(1e-7, 0.02, 220), [1.0]]))
    w_map = np.unique(np.concatenate([np.linspace(-1, 1, 801), 1 - np.geomspace(1e-6, 0.02, 60), -1 + np.geomspace(1e-6, 0.02, 60)]))
    cmap = {}
    for side, S in caps.items():
        d = S.grid(s_map, w_map); g = geometry(d, flips[side])
        r = 1 / g['kmax']
        X = d['P'][..., 0]
        req, t = required_radius(X)
        ratio = r / req
        # local minima of ratio in 3x3 neighbourhoods, refine the 25 smallest
        from scipy.ndimage import minimum_filter
        lm = (ratio == minimum_filter(ratio, size=3, mode='nearest'))
        cand = np.argwhere(lm)
        cand = cand[np.argsort(ratio[lm])][:25]

        def f(x):
            dd = S.at([x[0]], [x[1]]); gg = geometry(dd, flips[side])
            rr = 1 / gg['kmax'][0]; q, _ = required_radius(dd['P'][0, 0])
            return rr / q
        refined = []
        for (i, j) in cand:
            x0 = [s_map[i], w_map[j]]
            lo_s = s_map[max(i - 2, 0)]; hi_s = s_map[min(i + 2, len(s_map) - 1)]
            lo_w = w_map[max(j - 2, 0)]; hi_w = w_map[min(j + 2, len(w_map) - 1)]
            rr = minimize(f, x0, method='L-BFGS-B', bounds=[(lo_s, hi_s), (lo_w, hi_w)], options={'ftol': 1e-14, 'gtol': 1e-12})
            dd = S.at([rr.x[0]], [rr.x[1]]); gg = geometry(dd, flips[side])
            q, tq = required_radius(dd['P'][0, 0])
            refined.append({'s': float(rr.x[0]), 'w': float(rr.x[1]), 'XYZ': dd['P'][0].round(7).tolist(),
                            'radius_mm': float(1 / gg['kmax'][0]), 'attachment_thickness_mm': float(tq),
                            'required_mm': float(q), 'ratio': float(rr.fun)})
        refined.sort(key=lambda z: z['ratio'])
        ib = np.unravel_index(np.argmin(r), r.shape)
        thick = t >= RULE_T
        it = np.unravel_index(np.argmin(np.where(thick, r, np.inf)), r.shape)
        ithin = np.unravel_index(np.argmin(np.where(~thick, r, np.inf)), r.shape)
        # concave anywhere? (both principal curvatures should bend the same way on a convex cap)
        concave = (g['k1'] > 1e-9) & (g['k2'] > 1e-9)  # outward-normal convention: convex => k <= 0
        saddle = (g['k1'] * g['k2'] < -1e-12)
        cmap[side] = {
            'grid': [len(s_map), len(w_map)], 'points': int(r.size),
            'grid_s_min': float(s_map.min()), 'grid_reach_to_collapsed_LE_point_mm': float(np.linalg.norm(d['P'][0] - np.array([0, (1 if side == 'right' else -1) * 350.0, 0]), axis=-1).max()),
            'smallest_radius_grid_mm': float(r[ib]), 'smallest_radius_at': {'s': float(s_map[ib[0]]), 'w': float(w_map[ib[1]]), 'XYZ': d['P'][ib].round(6).tolist(), 'attachment_thickness_mm': float(t[ib])},
            'smallest_radius_where_thickness_ge_0.8_mm': float(r[it]),
            'its_location': {'s': float(s_map[it[0]]), 'w': float(w_map[it[1]]), 'XYZ': d['P'][it].round(6).tolist(), 'attachment_thickness_mm': float(t[it])},
            'smallest_radius_where_thickness_lt_0.8_mm': float(r[ithin]),
            'its_location_thin': {'s': float(s_map[ithin[0]]), 'w': float(w_map[ithin[1]]), 'XYZ': d['P'][ithin].round(7).tolist(), 'attachment_thickness_mm': float(t[ithin]), 'required_mm': float(req[ithin])},
            'min_ratio_radius_over_required_grid': float(ratio.min()),
            'refined_local_minima_of_ratio': refined[:8],
            'points_failing_rule': int((ratio < 1).sum()),
            'saddle_points_fraction': float(saddle.mean()),
            'saddle_region_s_range': [float(s_map[np.any(saddle, 1)].min()), float(s_map[np.any(saddle, 1)].max())] if saddle.any() else None,
            'concave_points': int(concave.sum()),
        }
    out['curvature_map'] = cmap
    # parent's own smallest radius (context): chordwise curvature of the parent section
    dp = parent_skin['upper'].grid(sg, [350.0]); gp = geometry(dp, pflip['upper'])
    rp = 1 / gp['kmax'][:, 0]
    out['parent_context'] = {'parent_skin_min_radius_mm': float(rp.min()), 'at_X_mm': float(dp['P'][np.argmin(rp), 0, 0]),
                             'parent_skin_min_radius_for_X_above_119_mm': float(rp[dp['P'][:, 0, 0] > 119].min())}

    # ---- sections
    out['sections'] = sections(caps, flips)

    # ---- zebra
    out['zebra'] = zebra(caps, flips, parent_skin, pflip)
    out['elapsed_s'] = time.time() - T0
    json.dump(out, open(OUT, 'w'), indent=1, default=float)
    print(json.dumps(summary(out), indent=1, default=float))


def trace_stats(P, kn, kplane, w_or_l, X, lip_index, label):
    """Statistics of a section trace: signed curvature reversals, min radius, lip radius."""
    sign = np.sign(kn)
    rev = np.where(np.diff(sign[np.abs(kn) > 1e-9]) != 0)[0]
    r_plane = 1 / np.maximum(np.abs(kplane), 1e-300)
    r_norm = 1 / np.maximum(np.abs(kn), 1e-300)
    i = np.argmin(r_plane); j = np.argmin(r_norm)
    return {'label': label, 'points': len(kn), 'reversals': int(len(rev)),
            'reversal_locations_XYZ': [P[np.abs(kn) > 1e-9][k].round(4).tolist() for k in rev[:6]],
            'min_trace_radius_mm': float(r_plane[i]), 'min_trace_radius_at_XYZ': P[i].round(5).tolist(),
            'min_normal_curvature_radius_mm': float(r_norm[j]), 'min_normal_radius_at_XYZ': P[j].round(5).tolist(),
            'lip_trace_radius_mm': float(r_plane[lip_index]), 'lip_normal_radius_mm': float(r_norm[lip_index]),
            'lip_XYZ': P[lip_index].round(5).tolist()}


def sections(caps, flips):
    res = {'note': ('trace radius = curvature radius of the planar section curve; normal radius = radius of the '
                    'surface normal curvature along the trace (equal when the plane contains the surface normal). '
                    'Signed: positive = convex (bending toward the solid).'), 'constant_X': [], 'normal_to_outline': [], 'constant_Y': []}
    S = caps['right']; fl = flips['right']
    # (a) constant-X planes: each is an s-iso curve of the cap (X depends on s only)
    s_list = [1e-4, 1e-3, 0.01, 0.03, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97,
              0.98, 0.99, 0.995, 0.999, 0.9999, 1.0]
    w = np.concatenate([np.linspace(1, -1, 8001)])
    for s in s_list:
        d = S.at(np.full_like(w, s), w); g = geometry(d, fl)
        P = d['P']; X = P[0, 0]
        t = float(thickness(X)); req = 0.4 if t >= 0.8 else t / 2
        Yw, Zw, Yww, Zww = d['Pv'][:, 1], d['Pv'][:, 2], d['Pvv'][:, 1], d['Pvv'][:, 2]
        kpl = (Yw * Zww - Zw * Yww) / (Yw ** 2 + Zw ** 2) ** 1.5   # +: CCW in (Y, Z) with w increasing
        kpl = -kpl  # trace runs w = +1 -> -1 on the right tip: clockwise is convex; make convex positive
        tvec = d['Pv'] / np.linalg.norm(d['Pv'], axis=-1, keepdims=True)
        kn, _, _ = II(g, d, tvec)
        kn = -kn  # outward normal: convex surfaces have negative II; make convex positive
        st = trace_stats(P, kn, kpl, w, X, np.argmin(np.abs(w)), f'X={X:.6f} (s={s})')
        # principal minimum radius along the trace (the map, restricted to this section)
        rmin = 1 / g['kmax']; k = np.argmin(rmin)
        # monotonic progression: is |k| increasing from each attachment to the lip?
        up = np.abs(kn[w >= 0])  # w from 1 down to 0
        prog_up = float(np.max(np.maximum.accumulate(up) - up))
        st.update({'X_mm': float(X), 'attachment_thickness_mm': t, 'required_radius_mm': req,
                   'principal_min_radius_on_section_mm': float(rmin[k]), 'principal_min_at_XYZ': P[k].round(5).tolist(),
                   'normal_curvature_max_drop_before_lip_per_mm': prog_up,
                   'normal_curvature_peak_before_lip': bool(np.argmax(up) < len(up) - 5 and up.max() > 1.05 * up[-1]),
                   'reach_mm': float(P[:, 1].max() - 350)})
        res['constant_X'].append(st)
    # (b) planes normal to the outline at the lip (plane contains Z and the outline normal)
    for p in [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98]:
        oj = PP.tip_outline_jet(np.array([p]))
        O = oj['point'][0]; T = oj['d_p'][0] / np.linalg.norm(oj['d_p'][0])
        st = plane_trace(S, fl, O, T, f'normal to outline at p={p}')
        if st: res['normal_to_outline'].append(st)
    # (c) constant-Y planes (chordwise cuts across the cap, the normal planes near the LE/TE ends)
    for Y in [350.05, 350.2, 350.5, 352, 355, 360, 370, 380, 390, 395, 398, 399.5]:
        st = plane_trace(S, fl, np.array([60.0, Y, 0.0]), np.array([0.0, 1.0, 0.0]), f'Y={Y}')
        if st: res['constant_Y'].append(st)
    return res


def plane_trace(S, fl, O, m, label, ns=1600, nw=900):
    """Trace the cap's intersection with the plane through O with normal m; curvature from surface derivatives."""
    import contourpy
    sg = np.unique(np.concatenate([np.geomspace(1e-6, 0.02, 200), np.linspace(0.02, 0.98, ns), 1 - np.geomspace(1e-8, 0.02, 300)]))
    wg = np.linspace(-1, 1, nw)
    P = S.grid(sg, wg)['P']
    f = ((P - O) * m).sum(-1)
    gen = contourpy.contour_generator(wg, sg, f)
    lines = gen.lines(0.0)
    if not lines:
        return None
    line = max(lines, key=len)
    wv, sv = line[:, 0], line[:, 1]
    # Newton-project each point onto the plane along the surface gradient of f
    for _ in range(6):
        d = S.at(sv, wv)
        fv = ((d['P'] - O) * m).sum(-1); fs = (d['Pu'] * m).sum(-1); fw = (d['Pv'] * m).sum(-1)
        gg = fs * fs + fw * fw
        sv = np.clip(sv - fv * fs / gg, 0, 1); wv = np.clip(wv - fv * fw / gg, -1, 1)
    d = S.at(sv, wv); g = geometry(d, fl)
    resid = np.abs(((d['P'] - O) * m).sum(-1)).max()
    tvec = np.cross(g['n'], m); tvec /= np.linalg.norm(tvec, axis=-1, keepdims=True)
    kn, _, _ = II(g, d, tvec); kn = -kn
    nu = np.cross(m, tvec)  # in-plane unit normal of the trace
    cosphi = np.abs((nu * g['n']).sum(-1))
    kpl = kn / np.maximum(cosphi, 1e-300)
    lip = np.argmin(np.abs(wv))
    X = d['P'][:, 0]
    st = trace_stats(d['P'], kn, kpl, wv, X, lip, label)
    req, t = required_radius(X)
    rmin = 1 / g['kmax']
    st.update({'plane_residual_mm': float(resid), 'w_range': [float(wv.min()), float(wv.max())], 's_range': [float(sv.min()), float(sv.max())],
               'min_ratio_trace_radius_over_required': float((1 / np.abs(kpl) / req).min()),
               'min_ratio_normal_radius_over_required': float((1 / np.abs(kn) / req).min()),
               'principal_min_radius_on_section_mm': float(rmin.min()),
               'attachment_thickness_at_tightest_trace_point_mm': float(t[np.argmin(1 / np.abs(kpl))])})
    return st


def zebra(caps, flips, parent_skin, pflip):
    """Numerical zebra across the cap-to-skin joins from CAD normals and derivatives."""
    period_deg = 4.0
    dtheta = np.radians(period_deg)
    sj = np.unique(np.concatenate([np.linspace(1e-4, 1, 6001), 1 - np.geomspace(1e-7, 1e-2, 300)]))
    views = []
    for el in (90, 60, 35, 15, -15, -35, -60, -90):
        for az in (0, 45, 90, 135, 180, 225, 270, 315):
            if abs(el) == 90 and az:
                continue
            e, a = np.radians(el), np.radians(az)
            views.append(np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)]))  # direction toward the eye
    stripe_angles = np.radians(np.arange(0, 180, 22.5))
    worst = {'phase_mismatch_stripes': 0, 'gradient_mismatch_relative': 0, 'crossing_direction_mismatch_deg': 0}
    where = {}
    ncross = 0
    for side, S in caps.items():
        sgn = 1 if side == 'right' else -1
        for skin, wv in (('upper', 1.0), ('lower', -1.0)):
            dc = S.at(sj, np.full_like(sj, wv)); gc = geometry(dc, flips[side])
            dp = parent_skin[skin].at(sj, np.full_like(sj, sgn * 350.0)); gp = geometry(dp, pflip[skin])
            for V in views:
                Vm = V * np.array([1, sgn, 1])  # mirrored lighting for the left tip
                eye = -Vm  # incident direction (from eye to surface)
                if (Vm * gp['n'][len(sj) // 2]).sum() <= 0:
                    continue  # join seen from behind
                up = np.array([0, 0, 1.0]) if abs(Vm[2]) < 0.9 else np.array([1.0, 0, 0])
                b1 = np.cross(Vm, up); b1 /= np.linalg.norm(b1); b2 = np.cross(Vm, b1)
                for ang in stripe_angles:
                    a = np.cos(ang) * b1 + np.sin(ang) * b2
                    a = a * np.array([1, sgn, 1])
                    vals = {}
                    for nm, (g, d) in (('cap', (gc, dc)), ('parent', (gp, dp))):
                        n = g['n']; vn = (eye * n).sum(-1)
                        r = eye - 2 * vn[:, None] * n
                        ra = np.clip((r * a).sum(-1), -1, 1)
                        phi = np.arcsin(ra) / dtheta
                        # dphi/du = a . dr/du / (sqrt(1 - ra^2) dtheta), dr = -2[(eye.dn) n + (eye.n) dn]
                        def dr(dn):
                            return -2 * ((eye * dn).sum(-1)[:, None] * n + vn[:, None] * dn)
                        den = np.sqrt(np.maximum(1 - ra * ra, 1e-300)) * dtheta
                        pu = (dr(g['nu']) * a).sum(-1) / den; pv = (dr(g['nv']) * a).sum(-1) / den
                        # 3D gradient in the tangent plane
                        ga = (g['G'] * pu - g['F'] * pv) / g['det']; gb = (g['E'] * pv - g['F'] * pu) / g['det']
                        grad = ga[:, None] * d['Pu'] + gb[:, None] * d['Pv']
                        vals[nm] = (phi, grad, n)
                    phc, grc, nc = vals['cap']; php, grp, npn = vals['parent']
                    dph = np.abs(phc - php)
                    gm = np.linalg.norm(grc - grp, axis=-1) / np.maximum(np.linalg.norm(grp, axis=-1), 1e-9)
                    # actual light/dark crossings on the join (phase crosses an integer on the parent side)
                    fl_ = np.floor(php)
                    cr = np.where(np.diff(fl_) != 0)[0]
                    if len(cr):
                        dc_ = np.cross(nc[cr], grc[cr]); dp_ = np.cross(npn[cr], grp[cr])
                        ok = (np.linalg.norm(dc_, axis=-1) > 1e-12) & (np.linalg.norm(dp_, axis=-1) > 1e-12)
                        angd = np.degrees(np.arctan2(np.linalg.norm(np.cross(dc_, dp_), axis=-1), np.abs((dc_ * dp_).sum(-1))))[ok]
                        ncross += int(ok.sum())
                        if len(angd) and angd.max() > worst['crossing_direction_mismatch_deg']:
                            k = cr[ok][np.argmax(angd)]
                            worst['crossing_direction_mismatch_deg'] = float(angd.max())
                            where['crossing_direction_mismatch'] = {'join': f'{side}_{skin}', 's': float(sj[k]), 'XYZ': dp['P'][k].round(5).tolist(), 'view': Vm.round(3).tolist()}
                    k = np.argmax(dph)
                    if dph[k] > worst['phase_mismatch_stripes']:
                        worst['phase_mismatch_stripes'] = float(dph[k]); where['phase'] = {'join': f'{side}_{skin}', 's': float(sj[k]), 'XYZ': dp['P'][k].round(5).tolist()}
                    k = np.argmax(gm)
                    if gm[k] > worst['gradient_mismatch_relative']:
                        worst['gradient_mismatch_relative'] = float(gm[k]); where['gradient'] = {'join': f'{side}_{skin}', 's': float(sj[k]), 'XYZ': dp['P'][k].round(5).tolist(), 'view': Vm.round(3).tolist()}
    return {'stripe_period_deg_of_reflected_ray': period_deg, 'views': len(views), 'stripe_angles': len(stripe_angles),
            'joins': 4, 'samples_per_join': len(sj), 's_range': [float(sj.min()), float(sj.max())],
            'excluded': 's < 1e-4 at the collapsed leading-edge point of the cap (its wrap shrinks to the point; derivatives in w vanish there)',
            'light_dark_crossings_checked': ncross, 'worst': worst, 'where': where}


def summary(o):
    s = {'step': o['step'], 'solid': {k: o['solid'][k] for k in ('solids', 'faces', 'brepcheck_valid', 'free_boundary_edges', 'volume_mm3', 'max_shape_tolerance_mm')},
         'parent_preserved': o['parent_preserved'],
         'span': o['locked_dimensions']['projected_span_mm'], 'outline': o['locked_dimensions']['outline'],
         'closure': o['locked_dimensions']['trailing_edge_closure'],
         'joins_worst_angle_deg': max(j['max_tangent_plane_angle_deg'] for j in o['joins'].values()),
         'joins_worst_sep_mm': max(j['max_separation_mm'] for j in o['joins'].values()),
         'joins_worst_curv_diff': max(max(j['max_abs_diff_normal_curvature_chordwise_per_mm'], j['max_abs_diff_twist_curvature_per_mm'], j['max_abs_diff_normal_curvature_spanwise_per_mm']) for j in o['joins'].values()),
         'map': {k: {kk: v[kk] for kk in ('smallest_radius_grid_mm', 'smallest_radius_at', 'smallest_radius_where_thickness_ge_0.8_mm', 'its_location', 'smallest_radius_where_thickness_lt_0.8_mm', 'its_location_thin', 'min_ratio_radius_over_required_grid', 'points_failing_rule', 'saddle_region_s_range', 'concave_points')} | {'best_refined': v['refined_local_minima_of_ratio'][:2]} for k, v in o['curvature_map'].items()},
         'parent_context': o['parent_context'],
         'zebra': o['zebra']['worst'], 'zebra_where': o['zebra']['where'],
         'sections_constant_X': [(x['label'], round(x['min_trace_radius_mm'], 4), round(x['min_normal_curvature_radius_mm'], 4), round(x['lip_normal_radius_mm'], 4), round(x['principal_min_radius_on_section_mm'], 4), x['reversals'], x['normal_curvature_peak_before_lip'], round(x['attachment_thickness_mm'], 4)) for x in o['sections']['constant_X']],
         'sections_normal': [(x['label'], round(x['min_trace_radius_mm'], 4), round(x['lip_trace_radius_mm'], 4), round(x['min_normal_curvature_radius_mm'], 4), x['reversals'], round(x['min_ratio_trace_radius_over_required'], 3), round(x['principal_min_radius_on_section_mm'], 4)) for x in o['sections']['normal_to_outline']],
         'sections_Y': [(x['label'], round(x['min_trace_radius_mm'], 4), round(x['lip_trace_radius_mm'], 4), round(x['min_normal_curvature_radius_mm'], 4), x['reversals'], round(x['min_ratio_normal_radius_over_required'], 3), round(x['principal_min_radius_on_section_mm'], 4)) for x in o['sections']['constant_Y']],
         'joins': {k: {kk: v[kk] for kk in ('max_separation_mm', 'max_tangent_plane_angle_deg', 'at_s', 'max_abs_diff_normal_curvature_chordwise_per_mm', 'max_abs_diff_twist_curvature_per_mm', 'max_abs_diff_normal_curvature_spanwise_per_mm', 'nearest_sample_to_collapsed_LE_point_mm')} for k, v in o['joins'].items()},
         'grid_reach_LE': {k: v['grid_reach_to_collapsed_LE_point_mm'] for k, v in o['curvature_map'].items()},
         'elapsed_s': o['elapsed_s']}
    return s


if __name__ == '__main__':
    main()
