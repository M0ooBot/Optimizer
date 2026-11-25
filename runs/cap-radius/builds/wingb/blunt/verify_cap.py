"""Verify wing.step (Wing B blunt, rounded tip cap) on the reimported STEP.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/wingb/blunt/verify_cap.py
Writes readback.json and zebra renders (PNG) next to this script.

Independent of the build: it uses only wing.step, the parent definition
(parent.py / parent.step / spec.json) and the skill's printable rule.
"""
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
PARENT_DIR = ROOT / "runs" / "cap-radius" / "parents" / "wingb" / "blunt"
sys.path.insert(0, str(PARENT_DIR))
import parent as P  # noqa: E402

from OCP.BRep import BRep_Tool  # noqa: E402
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface  # noqa: E402
from OCP.BRepCheck import BRepCheck_Analyzer  # noqa: E402
from OCP.BRepGProp import BRepGProp  # noqa: E402
from OCP.BRepTools import BRepTools  # noqa: E402
from OCP.Bnd import Bnd_Box  # noqa: E402
from OCP.BRepBndLib import BRepBndLib  # noqa: E402
from OCP.GProp import GProp_GProps  # noqa: E402
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Plane  # noqa: E402
from OCP.IFSelect import IFSelect_RetDone  # noqa: E402
from OCP.STEPControl import STEPControl_Reader  # noqa: E402
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_REVERSED, TopAbs_SHELL, TopAbs_SOLID  # noqa: E402
from OCP.TopExp import TopExp, TopExp_Explorer  # noqa: E402
from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape  # noqa: E402
from OCP.TopoDS import TopoDS  # noqa: E402
from OCP.BRep import BRep_Tool as BT  # noqa: E402,F811

STEP = HERE / "wing.step"
SPEC = json.loads((PARENT_DIR / "spec.json").read_text())
ATT = 540.0
NOZZLE, THIN = 0.4, 0.8
TE_ZONE = 2.0      # mm of span next to the attachment where the flat TE face fades out (by design)
NOSE_ZONE = 1.0    # mm of span next to the attachment where the parent's nose micro-spot fades out


# ----------------------------------------------------------------------------
# printable rule, judged against the attachment section's thickness at the same
# chord position (spec.json table; chord fraction from the section frame origin)
# ----------------------------------------------------------------------------
_tab = np.array([[r["chord_fraction"], r["thickness_mm"]] for r in SPEC["attachment_thickness"]["table"]])
_tab = _tab[np.argsort(_tab[:, 0])]


def chord_position(x, z):
    s = np.maximum(-np.asarray(x, float), 0.0)
    f = s / P.PROJECTED_SPAN
    z0 = -(35.0 * f + 45.0 * f * f)
    c = 200.0 - 125.0 * f
    return np.clip((z0 - np.asarray(z, float)) / c, 0.0, 1.0)


def att_thickness(xi):
    return np.interp(xi, _tab[:, 0], _tab[:, 1], left=0.0, right=_tab[-1, 1])


def rule_radius(xi):
    t = att_thickness(xi)
    return np.where(t >= THIN, NOZZLE, t / 2.0)


# ----------------------------------------------------------------------------
# B-spline faces evaluated with scipy from the reimported poles/knots
# ----------------------------------------------------------------------------
class Face:
    def __init__(self, face, label):
        self.face, self.label = face, label
        self.rev = face.Orientation() == TopAbs_REVERSED
        s = BRep_Tool.Surface_s(face)
        ad = BRepAdaptor_Surface(face)
        self.kind = ad.GetType()
        self.bounds = BRepTools.UVBounds_s(face)
        self.surf = s
        if self.kind == GeomAbs_BSplineSurface:
            b = ad.BSpline()
            self.bs = b

            def knots(n, k, m):
                out = []
                for i in range(1, n + 1):
                    out += [k(i)] * m(i)
                return np.array(out)

            self.tu = knots(b.NbUKnots(), b.UKnot, b.UMultiplicity)
            self.tv = knots(b.NbVKnots(), b.VKnot, b.VMultiplicity)
            self.du, self.dv = b.UDegree(), b.VDegree()
            nu, nv = b.NbUPoles(), b.NbVPoles()
            self.poles = np.array([[b.Pole(i, j).Coord() for j in range(1, nv + 1)] for i in range(1, nu + 1)])
            self.w = np.array([[b.Weight(i, j) for j in range(1, nv + 1)] for i in range(1, nu + 1)])

    def _ev(self, u, v, a, b):
        nu, nv, _ = self.poles.shape
        if a > self.du or b > self.dv:
            return np.zeros((len(u), len(v), 4))
        hp = np.concatenate([self.poles * self.w[..., None], self.w[..., None]], -1)  # homogeneous
        bv = BSpline(self.tv, np.transpose(hp, (1, 0, 2)).reshape(nv, -1), self.dv)
        Gv = (bv.derivative(b) if b else bv)(v).reshape(len(v), nu, 4)
        bu = BSpline(self.tu, np.transpose(Gv, (1, 0, 2)).reshape(nu, -1), self.du)
        return (bu.derivative(a) if a else bu)(u).reshape(len(u), len(v), 4)

    def derivs(self, u, v):
        """Grid evaluation: point and first/second derivatives (rational quotient rule)."""
        u, v = np.atleast_1d(u).astype(float), np.atleast_1d(v).astype(float)
        A = {(a, b): self._ev(u, v, a, b) for a, b in ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2))}
        w = {k: x[..., 3:] for k, x in A.items()}
        c = {k: x[..., :3] for k, x in A.items()}
        S = c[0, 0] / w[0, 0]
        Su = (c[1, 0] - w[1, 0] * S) / w[0, 0]
        Sv = (c[0, 1] - w[0, 1] * S) / w[0, 0]
        Suu = (c[2, 0] - 2 * w[1, 0] * Su - w[2, 0] * S) / w[0, 0]
        Svv = (c[0, 2] - 2 * w[0, 1] * Sv - w[0, 2] * S) / w[0, 0]
        Suv = (c[1, 1] - w[1, 0] * Sv - w[0, 1] * Su - w[1, 1] * S) / w[0, 0]
        return S, Su, Sv, Suu, Suv, Svv

    def points(self, uv, chunk=400):
        """Positions only at scattered (u, v)."""
        uv = np.asarray(uv, float)
        out = []
        for a in range(0, len(uv), chunk):
            c = uv[a:a + chunk]
            h = self._ev(c[:, 0], c[:, 1], 0, 0)
            idx = np.arange(len(c))
            h = h[idx, idx]
            out.append(h[:, :3] / h[:, 3:])
        return np.concatenate(out) if out else np.zeros((0, 3))

    def geometry(self, u, v):
        S, Su, Sv, Suu, Suv, Svv = self.derivs(u, v)
        n = np.cross(Su, Sv)
        nn = np.linalg.norm(n, axis=-1, keepdims=True)
        n = n / np.where(nn > 0, nn, 1)
        if self.rev:
            n = -n
        E, F, Gm = (Su * Su).sum(-1), (Su * Sv).sum(-1), (Sv * Sv).sum(-1)
        L, M, N = (Suu * n).sum(-1), (Suv * n).sum(-1), (Svv * n).sum(-1)
        den = E * Gm - F * F
        K = (L * N - M * M) / den
        H = (E * N - 2 * F * M + Gm * L) / (2 * den)
        disc = np.sqrt(np.maximum(H * H - K, 0))
        return {"S": S, "Su": Su, "Sv": Sv, "n": n, "k1": H + disc, "k2": H - disc,
                "I": (E, F, Gm), "II": (L, M, N), "nn": nn[..., 0]}


def load():
    r = STEPControl_Reader()
    assert r.ReadFile(str(STEP)) == IFSelect_RetDone
    r.TransferRoots()
    return r.OneShape()


def explore(shape, kind):
    out, ex = [], TopExp_Explorer(shape, kind)
    while ex.More():
        out.append(ex.Current())
        ex.Next()
    return out


def bbox(shape):
    b = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape, b, False, False)
    return list(b.Get())


def area_centroid(face):
    pr = GProp_GProps()
    BRepGProp.SurfaceProperties_s(face, pr)
    c = pr.CentreOfMass()
    return pr.Mass(), np.array(c.Coord())


def label_faces(shape):
    faces = []
    for f in explore(shape, TopAbs_FACE):
        f = TopoDS.Face_s(f)
        a, c = area_centroid(f)
        ad = BRepAdaptor_Surface(f)
        bb = bbox(f)
        faces.append((f, a, c, ad.GetType(), bb))
    labels = {}
    for f, a, c, kind, bb in faces:
        if kind == GeomAbs_Plane:
            lab = "root"
        elif bb[0] > -ATT - 1e-3:          # entirely inboard of the attachment
            lab = "parent_skin" if a > 1000 else "parent_te"
        elif a < 5.0:
            lab = "cap_strip"
        else:
            lab = None
        labels[id(f)] = lab
    caps = [(f, c) for f, a, c, kind, bb in faces if labels[id(f)] is None]
    caps.sort(key=lambda fc: -fc[1][1])
    out = {}
    for f, a, c, kind, bb in faces:
        lab = labels[id(f)]
        if lab is None:
            lab = "cap_upper" if f.IsSame(caps[0][0]) else "cap_lower"
        out[lab] = Face(f, lab)
    return out


# ----------------------------------------------------------------------------
# outline (locked plan perimeter)
# ----------------------------------------------------------------------------
_OP = np.linspace(0.0, 1.0, 40001)
_OXZ = P.tip_outline(_OP)


def outline_distance(x, z):
    """Signed plan distance to the locked outline (positive = inside the cap) and foot p."""
    from scipy.spatial import cKDTree
    tree = getattr(outline_distance, "_tree", None)
    if tree is None:
        tree = outline_distance._tree = cKDTree(_OXZ)
    q = np.stack([np.ravel(x), np.ravel(z)], -1)
    _, i = tree.query(q)
    best = np.full(len(q), np.inf)
    foot = np.zeros(len(q))
    sgn = np.zeros(len(q))
    for di in (-1, 0):
        a = np.clip(i + di, 0, len(_OP) - 2)
        A, B = _OXZ[a], _OXZ[a + 1]
        t = np.clip(((q - A) * (B - A)).sum(-1) / ((B - A) ** 2).sum(-1), 0, 1)
        C = A + t[:, None] * (B - A)
        dd = np.linalg.norm(q - C, axis=-1)
        tan = (B - A) / np.linalg.norm(B - A, axis=-1, keepdims=True)
        nrm = np.stack([-tan[:, 1], tan[:, 0]], -1)       # interior side
        s = np.sign(((q - C) * nrm).sum(-1))
        better = dd < best
        best = np.where(better, dd, best)
        foot = np.where(better, _OP[a] + t * (_OP[a + 1] - _OP[a]), foot)
        sgn = np.where(better, s, sgn)
    return (sgn * best).reshape(np.shape(x)), foot.reshape(np.shape(x))


# ----------------------------------------------------------------------------
# checks
# ----------------------------------------------------------------------------
def topology(shape):
    solids = explore(shape, TopAbs_SOLID)
    shells = explore(shape, TopAbs_SHELL)
    faces = explore(shape, TopAbs_FACE)
    m = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(shape, TopAbs_EDGE, TopAbs_FACE, m)
    not_two, degenerate, tol = 0, 0, 0.0
    for i in range(1, m.Extent() + 1):
        e = TopoDS.Edge_s(m.FindKey(i))
        tol = max(tol, BRep_Tool.Tolerance_s(e))
        if BRep_Tool.Degenerated_s(e):
            degenerate += 1
            continue
        if m.FindFromIndex(i).Size() != 2:
            not_two += 1
    pr = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, pr, 1e-9, True)
    return {"solids": len(solids), "shells": len(shells), "faces": len(faces),
            "edges_not_shared_by_two_faces": not_two, "degenerated_edges": degenerate,
            "closed": not_two == 0, "brepcheck_valid": bool(BRepCheck_Analyzer(shape).IsValid()),
            "max_edge_tolerance_mm": tol, "volume_mm3": pr.Mass(), "bounding_box_mm": bbox(shape)}


def parent_preserved(F):
    r = STEPControl_Reader()
    assert r.ReadFile(str(P.STEP_PATH)) == IFSelect_RetDone
    r.TransferRoots()
    orig = {}
    for f in explore(r.OneShape(), TopAbs_FACE):
        f = TopoDS.Face_s(f)
        a, c = area_centroid(f)
        kind = BRepAdaptor_Surface(f).GetType()
        if kind == GeomAbs_Plane:
            orig["root" if abs(c[0]) < 1e-6 else "attachment"] = (f, a, c)
        else:
            orig["parent_skin" if a > 1000 else "parent_te"] = (f, a, c)
    out = {}
    for lab in ("parent_skin", "parent_te"):
        o = Face(orig[lab][0], lab)
        n = F[lab]
        same_shape = o.poles.shape == n.poles.shape and len(o.tu) == len(n.tu) and len(o.tv) == len(n.tv)
        out[lab] = {"same_pole_and_knot_counts": bool(same_shape),
                    "max_pole_difference_mm": float(np.abs(o.poles - n.poles).max()) if same_shape else None,
                    "max_knot_difference": float(max(np.abs(o.tu - n.tu).max(), np.abs(o.tv - n.tv).max())) if same_shape else None,
                    "max_weight_difference": float(np.abs(o.w - n.w).max()) if same_shape else None,
                    "uv_bounds_original": list(BRepTools.UVBounds_s(orig[lab][0])),
                    "uv_bounds_readback": list(n.bounds),
                    "area_original_mm2": orig[lab][1], "area_readback_mm2": area_centroid(n.face)[0]}
    a, c = area_centroid(F["root"].face)
    out["root"] = {"area_original_mm2": orig["root"][1], "area_readback_mm2": a,
                   "centroid_difference_mm": float(np.linalg.norm(c - orig["root"][2]))}
    return out


def edge_pairs(shape, F):
    m = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(shape, TopAbs_EDGE, TopAbs_FACE, m)
    lab_of = lambda f: next(k for k, v in F.items() if v.face.IsSame(f))
    pairs = {}
    for i in range(1, m.Extent() + 1):
        e = TopoDS.Edge_s(m.FindKey(i))
        if BRep_Tool.Degenerated_s(e):
            continue
        lst = m.FindFromIndex(i)
        fs = [TopoDS.Face_s(f) for f in lst]
        if len(fs) != 2:
            continue
        key = tuple(sorted((lab_of(fs[0]), lab_of(fs[1]))))
        pairs.setdefault(key, []).append(e)
    return pairs


def _uv_dir(geo, k, e):
    """Express 3D tangent direction e (..., 3) in (du, dv) at sample k."""
    Su, Sv = geo["Su"][k], geo["Sv"][k]
    A = np.stack([Su, Sv], -1)
    sol, *_ = np.linalg.lstsq(A, e, rcond=None)
    return sol


def normal_curvature(geo, k, e):
    a, b = _uv_dir(geo, k, e)
    E, Fm, Gm = (x[k] for x in geo["I"])
    L, M, N = (x[k] for x in geo["II"])
    return (L * a * a + 2 * M * a * b + N * b * b) / (E * a * a + 2 * Fm * a * b + Gm * b * b)


def sample_edge(edge, face_a, face_b, n=400):
    from OCP.BRepAdaptor import BRepAdaptor_Curve2d
    from OCP.gp import gp_Pnt, gp_Vec
    c3 = BRepAdaptor_Curve(edge)
    t0, t1 = c3.FirstParameter(), c3.LastParameter()
    ts = np.linspace(t0, t1, n)
    pts, tans = [], []
    for t in ts:
        p, v = gp_Pnt(), gp_Vec()
        c3.D1(t, p, v)
        pts.append(p.Coord())
        tans.append(v.Coord())
    uvs = []
    for f in (face_a, face_b):
        c2 = BRepAdaptor_Curve2d(edge, f.face)
        s0, s1 = c2.FirstParameter(), c2.LastParameter()
        # the 2D and 3D parametrisations of one edge share their range in OCCT
        uvs.append(np.array([c2.Value(s0 + (t - t0) / (t1 - t0) * (s1 - s0)).Coord() for t in ts]))
    tans = np.array(tans)
    return np.array(pts), tans / np.maximum(np.linalg.norm(tans, axis=-1, keepdims=True), 1e-300), uvs


def point_geometry(face, uv, chunk=250):
    """Geometry at scattered (u, v) points (chunks of the grid evaluator, diagonal kept)."""
    uv = np.asarray(uv, float)
    parts = []
    for a in range(0, len(uv), chunk):
        c = uv[a:a + chunk]
        g = face.geometry(c[:, 0], c[:, 1])
        idx = np.arange(len(c))
        parts.append({k: (tuple(x[idx, idx] for x in v) if isinstance(v, tuple) else v[idx, idx])
                      for k, v in g.items()})
    out = {}
    for k in parts[0]:
        if isinstance(parts[0][k], tuple):
            out[k] = tuple(np.concatenate([p[k][j] for p in parts]) for j in range(3))
        else:
            out[k] = np.concatenate([p[k] for p in parts])
    return out



# ----------------------------------------------------------------------------
# joins: separation, tangent-plane angle, curvature in a common frame, zebra
# ----------------------------------------------------------------------------
VIEWS = {"top": (0, -1, 0), "bottom": (0, 1, 0), "front_high": (0, -0.6, -1), "rear_high": (0, -0.6, 1),
         "outboard_high": (1, -0.7, 0), "oblique": (0.6, -0.6, -0.5)}
AXES = {"X": (1, 0, 0), "Z": (0, 0, 1), "XZ": (1, 0, 1), "XmZ": (1, 0, -1)}
STRIPES = 8          # light/dark pairs per half turn of the reflected ray
FLAT_BAND_MM = 20.0
RENDER_STRIPES = 24  # finer stripes for the pictures  # stripes wider than this count as flat when comparing gradients/directions


def zebra_phase(n, view, axis):
    V = np.asarray(view, float)
    V = V / np.linalg.norm(V)
    A = np.asarray(axis, float)
    A = A / np.linalg.norm(A)
    B = np.cross(A, [0.3, 0.9, 0.1])
    B /= np.linalg.norm(B)
    C = np.cross(A, B)
    R = V - 2 * (n @ V)[..., None] * n
    return np.arctan2(R @ C, R @ B)


def zebra_join(fa, fb, uva, uvb, pts, tan, h=1e-6):
    """Phase mismatch, gradient mismatch and stripe-direction change at the light/dark
    crossings along one join, over all views and stripe axes."""
    def offsets(face, uv):
        u0, u1, v0, v1 = face.bounds
        hu, hv = h * (u1 - u0), h * (v1 - v0)
        res = []
        for du, dv in ((0, 0), (hu, 0), (-hu, 0), (0, hv), (0, -hv)):
            q = np.stack([np.clip(uv[:, 0] + du, u0, u1), np.clip(uv[:, 1] + dv, v0, v1)], -1)
            res.append((point_geometry(face, q), q))
        return res
    oa, ob = offsets(fa, uva), offsets(fb, uvb)
    worst = {"phase_mismatch_stripes": 0.0, "gradient_mismatch_rel": 0.0, "direction_change_deg": 0.0,
             "crossings": 0}
    where = {}
    for vn, V in VIEWS.items():
        for an, Ax in AXES.items():
            res = []
            for off in (oa, ob):
                th = [zebra_phase(g["n"], V, Ax) for g, _ in off]
                g0 = off[0][0]
                du = np.unwrap(np.stack([th[1], th[2]]), axis=0)
                dv = np.unwrap(np.stack([th[3], th[4]]), axis=0)
                tu = (du[0] - du[1]) / np.maximum(off[1][1][:, 0] - off[2][1][:, 0], 1e-300)
                tv = (dv[0] - dv[1]) / np.maximum(off[3][1][:, 1] - off[4][1][:, 1], 1e-300)
                E, Fm, Gm = g0["I"]
                det = E * Gm - Fm * Fm
                # 3D gradient of the phase on the surface
                a = (Gm * tu - Fm * tv) / det
                b = (E * tv - Fm * tu) / det
                grad = a[:, None] * g0["Su"] + b[:, None] * g0["Sv"]
                iso = np.cross(g0["n"], grad)
                iso /= np.maximum(np.linalg.norm(iso, axis=-1, keepdims=True), 1e-300)
                res.append((th[0], grad, iso))
            dth = np.angle(np.exp(1j * (res[0][0] - res[1][0])))
            ph = np.abs(dth) * STRIPES / np.pi
            floor = (np.pi / STRIPES) / FLAT_BAND_MM
            ga_, gb_ = np.linalg.norm(res[0][1], axis=-1), np.linalg.norm(res[1][1], axis=-1)
            gm = np.linalg.norm(res[0][1] - res[1][1], axis=-1) / np.maximum(0.5 * (ga_ + gb_), floor)
            band = np.sign(np.sin(STRIPES * res[0][0]))
            cross = np.where((band[1:] != band[:-1]) & (ga_[:-1] > floor) & (gb_[:-1] > floor))[0]
            dang = 0.0
            for c in cross:
                ca = np.degrees(np.arccos(np.clip(abs((res[0][2][c] * res[1][2][c]).sum()), 0, 1)))
                if ca > dang:
                    dang = ca
                    if ca > worst["direction_change_deg"]:
                        where["direction_change"] = {"view": vn, "axis": an, "point_mm": pts[c].tolist()}
            worst["crossings"] += len(cross)
            if ph.max() > worst["phase_mismatch_stripes"]:
                worst["phase_mismatch_stripes"] = float(ph.max())
                where["phase"] = {"view": vn, "axis": an, "point_mm": pts[np.argmax(ph)].tolist()}
            if gm.max() > worst["gradient_mismatch_rel"]:
                worst["gradient_mismatch_rel"] = float(gm.max())
                where["gradient"] = {"view": vn, "axis": an, "point_mm": pts[np.argmax(gm)].tolist()}
            worst["direction_change_deg"] = max(worst["direction_change_deg"], dang)
    worst["where"] = where
    return worst


def joins(shape, F, exclude_mm=0.05):
    pairs = edge_pairs(shape, F)
    corners = [np.array(P.skin(float(P.nose_u(ATT)), ATT)), np.array(P.te_face(0.0, ATT)), np.array(P.te_face(1.0, ATT))]
    out = {}
    for key, edges in sorted(pairs.items()):
        if "root" in key or not any(k.startswith("cap") for k in key):
            continue
        fa, fb = F[key[0]], F[key[1]]
        rows = []
        P3, T3, UA, UB = [], [], [], []
        for e in edges:
            pts, tan, (uva, uvb) = sample_edge(e, fa, fb, 600)
            P3.append(pts); T3.append(tan); UA.append(uva); UB.append(uvb)
        pts, tan, uva, uvb = (np.concatenate(x) for x in (P3, T3, UA, UB))
        dmin = np.min([np.linalg.norm(pts - c, axis=-1) for c in corners], axis=0)
        keep = dmin > exclude_mm
        pts, tan, uva, uvb = pts[keep], tan[keep], uva[keep], uvb[keep]
        ga, gb = point_geometry(fa, uva), point_geometry(fb, uvb)
        fin = np.isfinite(ga["k1"]) & np.isfinite(gb["k1"]) & (ga["nn"] > 1e-9) & (gb["nn"] > 1e-9)
        if not fin.all():
            pts, tan, uva, uvb = pts[fin], tan[fin], uva[fin], uvb[fin]
            ga, gb = point_geometry(fa, uva), point_geometry(fb, uvb)
        sep = np.linalg.norm(ga["S"] - gb["S"], axis=-1)
        on_edge = np.maximum(np.linalg.norm(ga["S"] - pts, axis=-1), np.linalg.norm(gb["S"] - pts, axis=-1))
        ang = np.degrees(np.arccos(np.clip((ga["n"] * gb["n"]).sum(-1), -1, 1)))
        e1 = tan
        e2 = np.cross(ga["n"], e1)
        e2 /= np.linalg.norm(e2, axis=-1, keepdims=True)
        ka1 = np.array([normal_curvature(ga, k, e1[k]) for k in range(len(pts))])
        kb1 = np.array([normal_curvature(gb, k, e1[k]) for k in range(len(pts))])
        ka2 = np.array([normal_curvature(ga, k, e2[k]) for k in range(len(pts))])
        kb2 = np.array([normal_curvature(gb, k, e2[k]) for k in range(len(pts))])
        dk = np.abs(ka2 - kb2)
        i = int(np.argmax(dk))
        out["|".join(key)] = {
            "samples": int(len(pts)), "excluded_within_mm_of_collapsed_corners": exclude_mm,
            "max_separation_mm": float(sep.max()), "max_distance_to_shared_edge_mm": float(on_edge.max()),
            "max_tangent_plane_angle_deg": float(ang.max()),
            "tangent_angle_at": pts[int(np.argmax(ang))].tolist(),
            "max_cross_edge_curvature_difference_per_mm": float(dk.max()),
            "cross_edge_curvatures_there_per_mm": [float(ka2[i]), float(kb2[i])],
            "cross_edge_curvature_difference_at": pts[i].tolist(),
            "median_cross_edge_curvature_difference_per_mm": float(np.median(dk)),
            "max_along_edge_curvature_difference_per_mm": float(np.abs(ka1 - kb1).max()),
            "zebra": zebra_join(fa, fb, uva, uvb, pts, tan),
        }
    return out



# ----------------------------------------------------------------------------
# sections
# ----------------------------------------------------------------------------
CAP = ("cap_upper", "cap_lower", "cap_strip")


def _grid(face, nu, nv):
    u0, u1, v0, v1 = face.bounds
    a = np.linspace(0, 1, nu)
    u = u0 + (u1 - u0) * np.unique(np.concatenate([a, 0.5 * (1 - np.cos(np.pi * a))]))
    b = np.linspace(0, 1, nv)
    v = v0 + (v1 - v0) * np.unique(np.concatenate([b, 1 - (1 - b) ** 3, b ** 2]))
    return u, v


def plane_points(face, P0, nrm, nu=500, nv=160, it=30):
    u, v = _grid(face, nu, nv)
    S = face.derivs(u, v)[0]
    f = (S - P0) @ nrm
    segs = []
    iu, jv = np.where(np.sign(f[1:, :]) != np.sign(f[:-1, :]))
    for i, j in zip(iu, jv):
        segs.append(((u[i], v[j]), (u[i + 1], v[j]), f[i, j]))
    iu, jv = np.where(np.sign(f[:, 1:]) != np.sign(f[:, :-1]))
    for i, j in zip(iu, jv):
        segs.append(((u[i], v[j]), (u[i], v[j + 1]), f[i, j]))
    if not segs:
        return None
    A = np.array([s[0] for s in segs])
    B = np.array([s[1] for s in segs])
    fa = (face.points(A) - P0) @ nrm
    fb = (face.points(B) - P0) @ nrm
    for _ in range(it):                      # Illinois regula falsi along each grid segment
        w = fa / np.where(fa - fb != 0, fa - fb, 1e-300)
        M = A + np.clip(w, 0, 1)[:, None] * (B - A)
        fm = (face.points(M) - P0) @ nrm
        left = np.sign(fm) == np.sign(fa)
        A = np.where(left[:, None], M, A)
        B = np.where(left[:, None], B, M)
        fa_new = np.where(left, fm, fa)
        fb_new = np.where(left, fb * 0.5, fm)
        fa_new = np.where(~left, fa * 0.5, fa_new)
        fa, fb = fa_new, fb_new
        if np.max(np.abs(fm)) < 1e-11:
            break
    uv = M
    g = point_geometry(face, uv)
    n = g["n"]
    T = np.cross(nrm, n)
    T /= np.maximum(np.linalg.norm(T, axis=-1, keepdims=True), 1e-300)
    kn = np.array([normal_curvature(g, k, T[k]) for k in range(len(uv))])
    nin = np.linalg.norm(n - (n @ nrm)[:, None] * nrm, axis=-1)
    conv = -kn / np.maximum(nin, 1e-12)       # > 0: bending away from the outward normal
    return {"S": g["S"], "conv": conv, "nin": nin, "uv": uv, "kn": kn}


def section_report(F, P0, nrm, kind, lip_dir=None, lip_point=None):
    pts = []
    for lab in CAP:
        r = plane_points(F[lab], P0, nrm)
        if r is not None:
            r["face"] = np.array([lab] * len(r["S"]))
            pts.append(r)
    if not pts:
        return None
    S = np.concatenate([r["S"] for r in pts])
    conv = np.concatenate([r["conv"] for r in pts])
    nin = np.concatenate([r["nin"] for r in pts])
    kn = np.concatenate([r["kn"] for r in pts])
    good = nin > 0.2
    S, conv, kn, nin = S[good], conv[good], kn[good], nin[good]
    R_plane = 1.0 / np.maximum(np.abs(conv), 1e-12)        # radius of the planar section curve
    R = 1.0 / np.maximum(np.abs(kn), 1e-12)                # surface (normal-section) radius along the trace
    xi = chord_position(S[:, 0], S[:, 2])
    rule = rule_radius(xi)
    d = -S[:, 0] - ATT
    te_zone = (d < TE_ZONE) & (xi > 0.98)
    nose_zone = (d < NOSE_ZONE) & (xi < 0.01)
    ok = ~(te_zone | nose_zone)
    rep = {"kind": kind, "plane_point_mm": np.round(P0, 4).tolist(), "plane_normal": np.round(nrm, 6).tolist(),
           "points": int(len(S)),
           "radius_definition": "surface normal-curvature radius along the trace (the in-plane radius of the "
                                "planar curve is also given; it is smaller where the plane meets the surface obliquely)"}
    if not ok.any():
        rep["note"] = "all points inside the trailing-edge or nose transition zones"
        return rep
    ratio = R / np.maximum(rule, 1e-9)
    i = np.where(ok)[0][np.argmin(R[ok])]
    j = np.where(ok)[0][np.argmin(ratio[ok])]
    ip = np.where(ok)[0][np.argmin(R_plane[ok])]
    rep.update({"min_radius_mm": float(R[i]), "min_radius_at_mm": np.round(S[i], 4).tolist(),
                "min_in_plane_radius_mm": float(R_plane[ip]), "min_in_plane_radius_at_mm": np.round(S[ip], 4).tolist(),
                "plane_obliquity_there": float(nin[ip]),
                "min_radius_chord_position": float(xi[i]), "attachment_thickness_there_mm": float(att_thickness(xi[i])),
                "rule_radius_there_mm": float(rule[i]),
                "min_radius_over_rule": float(ratio[j]), "worst_ratio_at_mm": np.round(S[j], 4).tolist(),
                "worst_ratio_radius_mm": float(R[j]), "worst_ratio_rule_mm": float(rule[j])})
    conc = conv < -1e-3
    rep["reversals"] = {"concave_points": int(conc.sum())}
    if conc.any():
        k = np.argmax(np.where(conc, -conv, -np.inf))
        rep["reversals"].update({"tightest_concave_radius_mm": float(1 / -conv[k]),
                                 "at_mm": np.round(S[k], 4).tolist(),
                                 "x_range_mm": [float(S[conc, 0].min()), float(S[conc, 0].max())]})
    if lip_dir is not None:
        k = np.argmax(S @ lip_dir)
        rep["lip"] = {"point_mm": np.round(S[k], 4).tolist(), "radius_mm": float(R[k]),
                      "in_plane_radius_mm": float(R_plane[k]),
                      "chord_position": float(xi[k]), "attachment_thickness_mm": float(att_thickness(xi[k])),
                      "rule_radius_mm": float(rule[k])}
        if lip_point is not None:
            rep["lip"]["plan_distance_to_outline_point_mm"] = float(np.hypot(S[k, 0] - lip_point[0], S[k, 2] - lip_point[1]))
    if (~ok).any():
        rep["transition_zone_min_radius_mm"] = float(R[~ok].min())
    return rep


def sections(F):
    out = []
    dp = 1e-5
    for p in [0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6,
              0.65, 0.7, 0.75, 0.8, 0.85, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 0.99, 0.995]:
        a, b = P.tip_outline(np.array([p - dp, p + dp]))
        o = P.tip_outline(np.array([p]))[0]
        t = (b - a) / np.linalg.norm(b - a)
        nrm = np.array([t[0], 0.0, t[1]])
        out_dir = np.array([t[1], 0.0, -t[0]])       # outward plan normal
        P0 = np.array([o[0], 50.0, o[1]])
        r = section_report(F, P0, nrm, f"perpendicular to the outline at p={p}", out_dir, o)
        if r:
            out.append(r)
    for d in [0.5, 1, 2, 3, 5, 10, 20, 30, 40, 50, 55, 58, 59.5]:
        P0 = np.array([-ATT - d, 50.0, -110.0])
        r = section_report(F, P0, np.array([1.0, 0, 0]), f"span station X={-ATT - d}")
        if r:
            for nm, direc in (("leading", np.array([0, 0, 1.0])), ("trailing", np.array([0, 0, -1.0]))):
                rr = section_report(F, P0, np.array([1.0, 0, 0]), "", direc)
                if rr and "lip" in rr:
                    r[f"{nm}_edge_point"] = rr["lip"]
            out.append(r)
    for xi in [0.0005, 0.002, 0.01, 0.05, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.97, 0.98, 0.99, 0.995]:
        f = ATT / P.PROJECTED_SPAN
        z = -(35 * f + 45 * f * f) - xi * (200 - 125 * f)
        P0 = np.array([-560.0, 50.0, z])
        r = section_report(F, P0, np.array([0, 0, 1.0]), f"chord station xi={xi} (plane Z={z:.4f}, parallel to the span)",
                           np.array([-1.0, 0, 0]))
        if r:
            out.append(r)
    return out


# ----------------------------------------------------------------------------
# principal curvature map with refinement
# ----------------------------------------------------------------------------
def curvature_map(F, nu=1400, nv=420, refine=60, log=print):
    res = {}
    allpts = []
    for lab in CAP:
        f = F[lab]
        u, v = _grid(f, nu, nv if lab != "cap_strip" else 40)
        g = f.geometry(u, v)
        S = g["S"]
        R = 1.0 / np.maximum(np.maximum(np.abs(g["k1"]), np.abs(g["k2"])), 1e-12)
        xi = chord_position(S[..., 0], S[..., 2])
        rule = rule_radius(xi)
        d = -S[..., 0] - ATT
        te_zone = (d < TE_ZONE) & (xi > 0.98)
        nose_zone = (d < NOSE_ZONE) & (xi < 0.01)
        u0, u1 = f.bounds[0], f.bounds[1]
        collapsed = (g["nn"] < 1e-9) | ((u[:, None] - u0) < 1e-6 * (u1 - u0)) | ((u1 - u[:, None]) < 1e-6 * (u1 - u0))
        collapsed = np.broadcast_to(collapsed, R.shape) | ~np.isfinite(R)
        R = np.where(np.isfinite(R), R, np.inf)
        ok = ~(te_zone | nose_zone | collapsed)
        ratio = R / np.maximum(rule, 1e-9)
        allpts.append((lab, u, v, S, R, ratio, ok, te_zone, nose_zone, collapsed, xi, rule))
        res[lab] = {"grid": [int(len(u)), int(len(v))]}
        if ok.any():
            res[lab]["min_radius_mm_outside_zones"] = float(R[ok].min())
            res[lab]["min_ratio_outside_zones"] = float(ratio[ok].min())
        else:
            res[lab]["note"] = "entirely inside the trailing-edge transition zone"
        if te_zone.any():
            res[lab]["te_zone_min_radius_mm"] = float(R[te_zone & ~collapsed].min())
            bad = te_zone & (ratio < 1)
            if bad.any():
                res[lab]["te_zone_below_rule_span_mm"] = [float(d[bad].min()), float(d[bad].max())]
        if nose_zone.any():
            nz = nose_zone & ~collapsed
            res[lab]["nose_zone_min_radius_mm"] = float(R[nz].min())
            nose_pt = np.array(P.skin(float(P.nose_u(ATT)), ATT))
            dist = np.linalg.norm(S - nose_pt, axis=-1)
            res[lab]["nose_zone_min_radius_by_distance_from_nose_point_mm"] = {
                f"{a}-{b}": float(R[nz & (dist >= a) & (dist < b)].min()) if (nz & (dist >= a) & (dist < b)).any() else None
                for a, b in ((0, 0.001), (0.001, 0.01), (0.01, 0.1), (0.1, 1.0), (1.0, 2.0))}
            bad = nz & (ratio < 1)
            if bad.any():
                res[lab]["nose_zone_below_rule_extent"] = {"span_mm": [float(d[bad].min()), float(d[bad].max())],
                                                           "max_distance_from_nose_point_mm": float(dist[bad].max())}
        # how close the grid gets to the collapsed edges
        uu = np.broadcast_to(u[:, None], R.shape)
        near0 = uu[~collapsed].min() - u0
        near1 = u1 - uu[~collapsed].max()
        res[lab]["closest_grid_u_to_collapsed_ends"] = [float(near0), float(near1)]
    # refinement of the worst local minima (outside the zones)
    cands = []
    for lab, u, v, S, R, ratio, ok, *_ in allpts:
        rr = np.where(ok, ratio, np.inf)
        flat = np.argsort(rr, axis=None)[: refine * 20]
        picked = []
        for k in flat:
            i, j = np.unravel_index(k, rr.shape)
            if all(abs(i - a) > 6 or abs(j - b) > 6 for a, b in picked):
                picked.append((i, j))
            if len(picked) >= refine:
                break
        for i, j in picked:
            cands.append((rr[i, j], lab, u[i], v[j]))
    cands.sort(key=lambda c: c[0])
    best = None
    for r0, lab, u0_, v0_ in cands[:refine]:
        f = F[lab]
        ub = f.bounds

        def obj(x):
            uu = np.clip(x[0], ub[0], ub[1])
            vv = np.clip(x[1], ub[2], ub[3])
            g = point_geometry(f, np.array([[uu, vv]]))
            S = g["S"][0]
            R = 1.0 / max(abs(g["k1"][0]), abs(g["k2"][0]), 1e-12)
            xi = chord_position(S[0], S[2])
            d = -S[0] - ATT
            if (d < TE_ZONE and xi > 0.98) or (d < NOSE_ZONE and xi < 0.01):
                return 1e3
            return R / max(float(rule_radius(xi)), 1e-9)
        scale = [1e-4 * (ub[1] - ub[0]), 1e-3 * (ub[3] - ub[2])]
        m = minimize(obj, [u0_, v0_], method="Nelder-Mead",
                     options={"xatol": 1e-9, "fatol": 1e-7, "maxiter": 300,
                              "initial_simplex": [[u0_, v0_], [u0_ + scale[0], v0_], [u0_, v0_ + scale[1]]]})
        if best is None or m.fun < best[0]:
            best = (m.fun, lab, m.x)
    fun, lab, x = best
    f = F[lab]
    g = point_geometry(f, np.array([[np.clip(x[0], f.bounds[0], f.bounds[1]), np.clip(x[1], f.bounds[2], f.bounds[3])]]))
    S = g["S"][0]
    R = 1.0 / max(abs(g["k1"][0]), abs(g["k2"][0]))
    xi = float(chord_position(S[0], S[2]))
    res["refined_worst"] = {"face": lab, "uv": [float(x[0]), float(x[1])], "point_mm": np.round(S, 5).tolist(),
                            "radius_mm": float(R), "chord_position": xi,
                            "attachment_thickness_mm": float(att_thickness(xi)),
                            "rule_radius_mm": float(rule_radius(xi)), "radius_over_rule": float(fun),
                            "span_beyond_attachment_mm": float(-S[0] - ATT)}
    # smallest absolute radius outside the zones
    k = min(((R[ok].min(), lab, np.unravel_index(np.argmin(np.where(ok, R, np.inf)), R.shape), S, xi_, rule)
             for lab, u, v, S, R, ratio, ok, te, nz, col, xi_, rule in allpts if ok.any()), key=lambda a: a[0])
    i, j = k[2]
    res["smallest_radius_outside_zones"] = {"face": k[1], "radius_mm": float(k[0]), "point_mm": np.round(k[3][i, j], 4).tolist(),
                                            "chord_position": float(k[4][i, j]), "rule_radius_mm": float(k[5][i, j])}
    return res, allpts



# ----------------------------------------------------------------------------
# zebra renders of the reimported model (edge lines hidden, fixed camera/light)
# ----------------------------------------------------------------------------
def render(F, path, view, up, axis, xlim=None, size=(1600, 1000), title="", dens=1.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    V = np.asarray(view, float)
    V /= np.linalg.norm(V)
    right = np.cross(V, up)
    right /= np.linalg.norm(right)
    upv = np.cross(right, V)
    S_all, C_all = [], []
    for lab, f in F.items():
        if lab == "root":
            continue
        u0, u1, v0, v1 = f.bounds
        if lab == "parent_skin":
            nu, nv = int(1600 * dens), int(700 * dens)
            if xlim is not None:
                vv = np.linspace(v0, v1, 4001)
                xs = f.derivs(np.array([0.5 * (u0 + u1)]), vv)[0][0, :, 0]
                v0 = vv[np.argmax(xs <= xlim[1] + 5)] if (xs <= xlim[1] + 5).any() else v0
        elif lab == "parent_te":
            nu, nv = 8, int(700 * dens)
        else:
            nu, nv = int(1500 * dens), int(500 * dens)
        a = np.linspace(0, 1, nu)
        u = u0 + (u1 - u0) * np.unique(np.concatenate([a, 0.5 * (1 - np.cos(np.pi * a))]))
        v = np.linspace(v0, v1, nv)
        if lab.startswith("cap"):
            b = np.linspace(0, 1, nv)
            v = v0 + (v1 - v0) * np.unique(np.concatenate([b, 1 - (1 - b) ** 2, 1 - (1 - b) ** 4]))
        Sl, nl = [], []
        for k in range(0, len(u), 300):
            g = f.geometry(u[k:k + 300], v)
            Sl.append(g["S"].reshape(-1, 3))
            nl.append(g["n"].reshape(-1, 3))
        S, n = np.concatenate(Sl), np.concatenate(nl)
        th = zebra_phase(n, V, axis)
        stripe = (np.sin(RENDER_STRIPES * th) > 0).astype(float)
        shade = 0.55 + 0.45 * np.clip(-(n @ V), 0, 1)
        col = (0.12 + 0.82 * stripe) * shade
        if xlim is not None:
            m = (S[:, 0] >= xlim[0]) & (S[:, 0] <= xlim[1])
            S, col = S[m], col[m]
        S_all.append(S)
        C_all.append(col)
    S = np.concatenate(S_all)
    col = np.concatenate(C_all)
    x, y, depth = S @ right, S @ upv, S @ V
    W, H = size
    sx = (x.max() - x.min()) / (W * 0.94)
    sy = (y.max() - y.min()) / (H * 0.94)
    sc = max(sx, sy)
    ix = ((x - 0.5 * (x.max() + x.min())) / sc + W / 2).astype(int)
    iy = (H / 2 - (y - 0.5 * (y.max() + y.min())) / sc).astype(int)
    img = np.ones((H, W))
    order = np.argsort(-depth)
    for dx in (0, 1):
        for dy in (0, 1):
            jx, jy = np.clip(ix + dx, 0, W - 1), np.clip(iy + dy, 0, H - 1)
            img[jy[order], jx[order]] = col[order]
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    ax.axis("off")
    ax.text(10, 24, title + f"   ({sc * 1000:.0f} um/px)", fontsize=11, family="monospace", color="black",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8))
    fig.savefig(path)
    plt.close(fig)


# ----------------------------------------------------------------------------
def main():
    t0 = time.time()
    log = lambda *a: print(f"[{time.time() - t0:6.1f}s]", *a, flush=True)
    data = STEP.read_bytes()
    shape = load()
    F = label_faces(shape)
    out = {"step_file": STEP.name, "step_sha256": hashlib.sha256(data).hexdigest(), "step_bytes": len(data),
           "faces": {k: {"type": "plane" if v.kind == GeomAbs_Plane else "bspline",
                         "area_mm2": area_centroid(v.face)[0],
                         "poles": list(v.poles.shape[:2]) if v.kind == GeomAbs_BSplineSurface else None,
                         "degrees": [v.du, v.dv] if v.kind == GeomAbs_BSplineSurface else None}
                     for k, v in F.items()}}
    out["topology"] = topology(shape)
    log("topology", out["topology"])
    out["parent_preserved"] = parent_preserved(F)
    log("parent preserved")
    bb = out["topology"]["bounding_box_mm"]
    # locked dimensions and outline
    cap_pts = []
    for lab in CAP:
        f = F[lab]
        u, v = _grid(f, 900, 200)
        cap_pts.append(f.derivs(u, v)[0].reshape(-1, 3))
    cp = np.concatenate(cap_pts)
    sd, foot = outline_distance(cp[:, 0], cp[:, 2])
    pairs = edge_pairs(shape, F)
    lip = []
    for key in (("cap_lower", "cap_upper"), ("cap_strip", "cap_upper"), ("cap_lower", "cap_strip")):
        for e in pairs.get(key, []):
            c = BRepAdaptor_Curve(e)
            for t in np.linspace(c.FirstParameter(), c.LastParameter(), 3000):
                lip.append(c.Value(t).Coord())
    lip = np.array(lip)
    ld, lfoot = outline_distance(lip[:, 0], lip[:, 2])
    out["locked"] = {
        "projected_half_span_mm": -bb[0], "attachment_station_mm": ATT,
        "join_edges_max_abs_x_deviation_from_attachment_mm": None,
        "cap_points_outside_outline_max_mm": float(max(0.0, -sd.min())),
        "perimeter_edge_to_outline_max_mm": float(np.abs(ld).max()),
        "perimeter_edge_outline_coverage_p": [float(lfoot.min()), float(lfoot.max())],
        "perimeter_elevation_mm": {f"p={p}": float(np.interp(p, np.sort(lfoot), lip[np.argsort(lfoot), 1]))
                                   for p in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)},
    }
    jx = []
    for key, edges in pairs.items():
        if ("parent_skin" in key or "parent_te" in key) and any(k.startswith("cap") for k in key):
            for e in edges:
                c = BRepAdaptor_Curve(e)
                jx += [c.Value(t).X() for t in np.linspace(c.FirstParameter(), c.LastParameter(), 200)]
    out["locked"]["join_edges_max_abs_x_deviation_from_attachment_mm"] = float(np.abs(np.array(jx) + ATT).max())
    log("locked", out["locked"])
    out["joins"] = joins(shape, F)
    log("joins done")
    out["sections"] = sections(F)
    log("sections done")
    cm, _ = curvature_map(F, log=log)
    out["curvature_map"] = cm
    log("map", cm["refined_worst"])
    out["rule"] = {"minimum_radius_mm": NOZZLE, "thin_threshold_mm": THIN,
                   "judged_against": "attachment thickness at the point's chord position (spec.json table, chord "
                                     "fraction from the section frame origin at the point's span)",
                   "zones_reported_separately": {
                       "te_zone": f"span beyond attachment < {TE_ZONE} mm and chord position > 0.98: the parent's "
                                  "sharp TE corners (flat 0.4 mm face) fade into the half-round",
                       "nose_zone": f"span beyond attachment < {NOSE_ZONE} mm and chord position < 0.01: the parent's "
                                    "own nose micro-spot fades out",
                       "collapsed": "degenerate patch corners (nose point, TE corners) where derivatives vanish"}}
    for nm, view, up, axis, xlim, title in (
            ("zebra_full_top.png", (0, -1, 0), (0, 0, 1), (1, 0, 0), None, "full half-wing, top, stripes along X"),
            ("zebra_tip_top.png", (0, -1, 0), (0, 0, 1), (1, 0, 0), (-605, -500), "tip, top, stripes along X"),
            ("zebra_tip_oblique.png", (0.5, -0.6, -0.6), (0, 1, 0), (0, 0, 1), (-605, -500), "tip, oblique from front-outboard, stripes along Z"),
            ("zebra_tip_bottom.png", (0, 1, 0), (0, 0, 1), (1, 0, 1), (-605, -500), "tip, bottom, stripes along X+Z"),
            ("zebra_te_corner.png", (0.35, -0.6, 0.7), (0, 1, 0), (1, 0, 0), (-560, -530), "trailing-edge corner, rear-outboard"),
            ("zebra_le_corner.png", (0.3, -0.4, -0.85), (0, 1, 0), (1, 0, 0), (-560, -530), "leading-edge corner, front-outboard")):
        render(F, HERE / nm, view, up, axis, xlim, title=title, dens=1.0 if xlim is None else 1.6)
    out["renders"] = [p.name for p in sorted(HERE.glob("zebra_*.png"))]
    log("renders done")
    (HERE / "readback.json").write_text(json.dumps(out, indent=1, default=float))
    log("wrote readback.json")


if __name__ == "__main__":
    main()
