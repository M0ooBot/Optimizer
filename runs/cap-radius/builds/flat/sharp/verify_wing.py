"""Verify runs/cap-radius/builds/flat/sharp/wing.step on the reimported STEP; writes readback.json and PNGs.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/flat/sharp/verify_wing.py

Everything measured here comes from the reimported STEP (OCCT reader). Dense
grids evaluate the reimported B-spline faces with NumPy from the poles and knots
read back from the file (cross-checked against OCCT's GeomLProp on random points);
refinements and the final reported extremes use OCCT directly. The parent
definition (parent.py) is used only as the reference the build is compared with.
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
PARENT_DIR = ROOT / "runs/cap-radius/parents/flat/sharp"
sys.path.insert(0, str(PARENT_DIR))
import parent  # noqa: E402

from OCP.BRep import BRep_Tool  # noqa: E402
from OCP.BRepBndLib import BRepBndLib  # noqa: E402
from OCP.Bnd import Bnd_Box  # noqa: E402
from OCP.BRepCheck import BRepCheck_Analyzer  # noqa: E402
from OCP.BRepGProp import BRepGProp  # noqa: E402
from OCP.GProp import GProp_GProps  # noqa: E402
from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf  # noqa: E402
from OCP.GeomAdaptor import GeomAdaptor_Surface  # noqa: E402
from OCP.GeomLProp import GeomLProp_SLProps  # noqa: E402
from OCP.IFSelect import IFSelect_RetDone  # noqa: E402
from OCP.STEPControl import STEPControl_Reader  # noqa: E402
from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds  # noqa: E402
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_REVERSED, TopAbs_SHELL, TopAbs_SOLID  # noqa: E402
from OCP.TopExp import TopExp_Explorer  # noqa: E402
from OCP.TopoDS import TopoDS  # noqa: E402
from OCP.gp import gp_Pnt  # noqa: E402

STEP = HERE / "wing.step"
RULE_MIN_R = 0.4
RULE_T = 0.8
GRAD_FLOOR = 1e-3          # stripes per mm: below this (stripe spacing > 1 m) a stripe direction is not visible


# ---------------------------------------------------------------- reading
def read_step(path):
    reader = STEPControl_Reader()
    if reader.ReadFile(str(path)) != IFSelect_RetDone:
        raise RuntimeError(f"cannot read {path}")
    reader.TransferRoots()
    return reader.OneShape()


def explore(shape, kind):
    out = []
    exp = TopExp_Explorer(shape, kind)
    while exp.More():
        out.append(exp.Current())
        exp.Next()
    return out


class Face:
    """A reimported B-spline face: OCCT handle plus a NumPy evaluator of the same poles and knots."""

    def __init__(self, face):
        self.face = TopoDS.Face_s(face)
        self.surf = BRep_Tool.Surface_s(self.face)
        ad = GeomAdaptor_Surface(self.surf)
        self.type = str(ad.GetType()).split(".")[-1]
        self.sign = -1.0 if self.face.Orientation() == TopAbs_REVERSED else 1.0
        if self.type != "GeomAbs_BSplineSurface":
            return
        bs = ad.BSpline()
        self.bs = bs
        if bs.IsURational() or bs.IsVRational():
            raise RuntimeError("rational surface not expected")
        self.P = np.array([[[bs.Pole(i, j).X(), bs.Pole(i, j).Y(), bs.Pole(i, j).Z()]
                            for j in range(1, bs.NbVPoles() + 1)] for i in range(1, bs.NbUPoles() + 1)])
        self.tu = np.repeat([bs.UKnot(i) for i in range(1, bs.NbUKnots() + 1)],
                            [bs.UMultiplicity(i) for i in range(1, bs.NbUKnots() + 1)])
        self.tv = np.repeat([bs.VKnot(i) for i in range(1, bs.NbVKnots() + 1)],
                            [bs.VMultiplicity(i) for i in range(1, bs.NbVKnots() + 1)])
        self.du, self.dv = bs.UDegree(), bs.VDegree()
        self.urange = (self.tu[0], self.tu[-1])
        self.vrange = (self.tv[0], self.tv[-1])
        self.Bu = [BSpline(self.tu, np.eye(self.P.shape[0]), self.du).derivative(d) if d else
                   BSpline(self.tu, np.eye(self.P.shape[0]), self.du) for d in range(3)]
        self.Bv = [BSpline(self.tv, np.eye(self.P.shape[1]), self.dv).derivative(d) if d else
                   BSpline(self.tv, np.eye(self.P.shape[1]), self.dv) for d in range(3)]

    def grid(self, u, v):
        """Point and derivatives on the tensor grid u x v: dict of (nu, nv, 3) arrays."""
        Nu = [np.nan_to_num(b(np.clip(u, *self.urange))) for b in self.Bu]
        Nv = [np.nan_to_num(b(np.clip(v, *self.vrange))) for b in self.Bv]
        ev = lambda a, b: np.einsum("ia,jb,abk->ijk", Nu[a], Nv[b], self.P, optimize=True)
        return {"p": ev(0, 0), "su": ev(1, 0), "sv": ev(0, 1), "suu": ev(2, 0), "suv": ev(1, 1), "svv": ev(0, 2)}

    def pts(self, u, v):
        """Point and derivatives at paired parameters (arrays of equal length)."""
        Nu = [np.nan_to_num(b(np.clip(u, *self.urange))) for b in self.Bu]
        Nv = [np.nan_to_num(b(np.clip(v, *self.vrange))) for b in self.Bv]
        ev = lambda a, b: np.einsum("na,nb,abk->nk", Nu[a], Nv[b], self.P, optimize=True)
        return {"p": ev(0, 0), "su": ev(1, 0), "sv": ev(0, 1), "suu": ev(2, 0), "suv": ev(1, 1), "svv": ev(0, 2)}

    def occ(self, u, v):
        """OCCT point, outward normal, principal curvatures (signed for the outward normal)."""
        pr = GeomLProp_SLProps(self.surf, float(u), float(v), 2, 1e-12)
        p = pr.Value()
        if not pr.IsNormalDefined():
            return np.array([p.X(), p.Y(), p.Z()]), None, None
        n = pr.Normal()
        n = self.sign * np.array([n.X(), n.Y(), n.Z()])
        k = np.array([pr.MaxCurvature(), pr.MinCurvature()]) * self.sign
        return np.array([p.X(), p.Y(), p.Z()]), n, k

    def value(self, u, v):
        p = self.surf.Value(float(u), float(v))
        return np.array([p.X(), p.Y(), p.Z()])

    def project(self, xyz):
        pr = GeomAPI_ProjectPointOnSurf(gp_Pnt(*map(float, xyz)), self.surf)
        if pr.NbPoints() == 0:
            return None
        u, v = pr.LowerDistanceParameters()
        return u, v, pr.LowerDistance()


def geometry(d, sign):
    """Outward unit normal, shape-operator data and principal curvatures from derivative arrays."""
    su, sv = d["su"], d["sv"]
    nraw = np.cross(su, sv) * sign
    nn = np.linalg.norm(nraw, axis=-1)
    n = nraw / nn[..., None]
    E = (su * su).sum(-1); F = (su * sv).sum(-1); G = (sv * sv).sum(-1)
    L = (d["suu"] * n).sum(-1); M = (d["suv"] * n).sum(-1); N = (d["svv"] * n).sum(-1)
    det = E * G - F * F
    H = (E * N - 2 * F * M + G * L) / (2 * det)
    K = (L * N - M * M) / det
    disc = np.sqrt(np.maximum(H * H - K, 0))
    k1, k2 = H + disc, H - disc        # curvature positive where the surface bends away from the outward normal? see note
    return {"n": n, "E": E, "F": F, "G": G, "L": L, "M": M, "N": N, "k1": k1, "k2": k2, "area_el": nn}


def required_radius(X):
    xi = np.clip(np.asarray(X, float) / parent.CHORD, 0, 1)
    t = parent.section_thickness(xi, parent.ATTACHMENT)
    t = np.nan_to_num(t, nan=0.0)
    return np.where(t >= RULE_T, RULE_MIN_R, t / 2), t


def volume(shape):
    pr = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, pr, 1e-9)
    return pr.Mass()


# ---------------------------------------------------------------- checks
def check_topology(shape):
    solids = explore(shape, TopAbs_SOLID)
    shells = explore(shape, TopAbs_SHELL)
    fcs = explore(shape, TopAbs_FACE)
    fb = ShapeAnalysis_FreeBounds(shape, 1e-7, False, False)
    free = len(explore(fb.GetOpenWires(), TopAbs_EDGE)) + len(explore(fb.GetClosedWires(), TopAbs_EDGE))
    box = Bnd_Box(); box.SetGap(0.0)
    BRepBndLib.AddOptimal_s(shape, box, False, False)
    bb = box.Get()
    degenerate = sum(1 for e in explore(shape, TopAbs_EDGE) if BRep_Tool.Degenerated_s(TopoDS.Edge_s(e)))
    return {"solids": len(solids), "shells": len(shells), "faces": len(fcs),
            "face_types": [Face(f).type for f in fcs],
            "brepcheck_valid": bool(BRepCheck_Analyzer(shape).IsValid()),
            "free_boundary_edges": free, "shell_closed": free == 0,
            "degenerated_edges": degenerate,
            "volume_mm3": volume(shape), "bounding_box_mm": [float(x) for x in bb]}


def classify(fcs):
    """Upper skin, lower skin, right cap, left cap among the reimported faces."""
    out = {}
    for f in fcs:
        g = f.grid(np.array([0.5 * sum(f.urange)]), np.linspace(*f.vrange, 5))["p"][0]
        ymax = g[:, 1].max(); ymin = g[:, 1].min(); zm = g[:, 2].mean()
        if ymax - ymin > 600:
            out["upper_skin" if zm > 0 else "lower_skin"] = f
        else:
            out["right_cap" if ymin > 0 else "left_cap"] = f
    return out


def check_skin_preserved(F):
    """Reimported skins against parent.py (exact) on a dense grid; also against parent.step faces."""
    s = np.linspace(0, 1, 241); y = np.linspace(-350, 350, 141)
    res = {}
    pshape = read_step(PARENT_DIR / "parent.step")
    pfaces = [Face(f) for f in explore(pshape, TopAbs_FACE)]
    pf = {("upper_skin" if f.grid(np.array([0.5]), np.array([0.0]))["p"][0, 0, 2] > 0 else "lower_skin"): f
          for f in pfaces if f.type == "GeomAbs_BSplineSurface"}
    for name, side in (("upper_skin", +1), ("lower_skin", -1)):
        f = F[name]
        ref = parent.skin(s[:, None], y[None, :], side)
        got = f.grid(s, y)["p"]          # skin parameters are native u and Y (checked below by projection)
        dev_py = np.linalg.norm(got - ref, axis=-1).max()
        same_poles = bool(np.array_equal(f.P, pf[name].P)) and bool(np.array_equal(f.tu, pf[name].tu)) \
            and bool(np.array_equal(f.tv, pf[name].tv))
        pole_dev = float(np.abs(f.P - pf[name].P).max())
        # parameterisation-free check: project parent.py points onto the reimported face
        rng = np.random.default_rng(1)
        proj = 0.0
        for _ in range(300):
            ss, yy = rng.uniform(0, 1), rng.uniform(-350, 350)
            q = parent.skin(ss, yy, side)
            proj = max(proj, f.project(q)[2])
        res[name] = {"max_dev_to_parent_py_grid_mm": float(dev_py), "max_projection_dist_parent_py_mm": float(proj),
                     "poles_knots_identical_to_parent_step": same_poles, "max_pole_diff_to_parent_step_mm": pole_dev,
                     "parameter_range_v": list(map(float, f.vrange))}
    return res


def extents(F):
    """Exact coordinate extremes of the solid from its faces (dense grid, then Nelder-Mead on the CAD surface)."""
    out = {}
    for ax, name in enumerate("XYZ"):
        lo, hi = np.inf, -np.inf
        for f in F.values():
            u = np.linspace(*f.urange, 801); v = np.linspace(*f.vrange, 801)
            g = f.grid(u, v)["p"][..., ax]
            for sgn in (+1, -1):
                i, j = np.unravel_index(np.argmax(sgn * g), g.shape)
                r = minimize(lambda x: -sgn * f.value(np.clip(x[0], *f.urange), np.clip(x[1], *f.vrange))[ax],
                             [u[i], v[j]], method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-14})
                val = f.value(np.clip(r.x[0], *f.urange), np.clip(r.x[1], *f.vrange))[ax]
                lo, hi = (min(lo, val), hi) if sgn < 0 else (lo, max(hi, val))
        out[name] = [float(lo), float(hi)]
    return out


def check_outline(F):
    """The lip is the locked outline: projection distance, Z, horizontal normal; and nothing pokes outside."""
    p = np.linspace(0, 1, 2001)
    res = {}
    for name, side in (("right_cap", +1), ("left_cap", -1)):
        f = F[name]
        ol = parent.tip_outline(p, side)
        dmax, nz, ndot = 0.0, 0.0, 0.0
        tang = parent.tip_outline_jet(p, side)["d_p"]
        for i in range(1, len(p) - 1, 4):
            u, v, d = f.project(ol[i])
            dmax = max(dmax, d)
            _, n, _ = f.occ(u, v)
            t = tang[i] / np.linalg.norm(tang[i])
            nz = max(nz, abs(n[2])); ndot = max(ndot, abs(n @ t))
        # outside check on a dense grid: |Y| beyond the outline at the same X
        g = f.grid(np.linspace(0, 1, 801), np.linspace(*f.vrange, 801))["p"].reshape(-1, 3)
        xi = np.clip(g[:, 0] / parent.CHORD, 0, 1)
        pp = np.empty_like(xi)
        lo, hi = np.zeros_like(xi), np.ones_like(xi)
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            xm = mid**3 * (10 - 15 * mid + 6 * mid * mid)
            lo = np.where(xm < xi, mid, lo); hi = np.where(xm < xi, hi, mid)
        pp = 0.5 * (lo + hi)
        y_ol = parent.ATTACHMENT + 200 * pp * (1 - pp)
        poke = (np.abs(g[:, 1]) - y_ol).max()
        inboard = (parent.ATTACHMENT - np.abs(g[:, 1])).max()
        res[name] = {"max_outline_to_cap_distance_mm": float(dmax), "max_abs_normal_z_on_lip": float(nz),
                     "max_abs_normal_dot_outline_tangent": float(ndot),
                     "max_cap_beyond_outline_mm": float(poke), "max_cap_inboard_of_attachment_mm": float(inboard),
                     "samples": int(len(range(1, len(p) - 1, 4)))}
    return res


def frame_curvatures(geo_i, d_i, e1, e2):
    """Normal curvatures k(e1), k(e2) and twist for unit tangent vectors e1, e2 (common frame)."""
    su, sv = d_i["su"], d_i["sv"]
    E, Fv, G = geo_i["E"], geo_i["F"], geo_i["G"]
    det = E * G - Fv * Fv

    def comps(e):
        a = (e * su).sum(-1); b = (e * sv).sum(-1)
        cu = (G * a - Fv * b) / det; cv = (E * b - Fv * a) / det
        return cu, cv

    def II(x, y):
        return geo_i["L"] * x[0] * y[0] + geo_i["M"] * (x[0] * y[1] + x[1] * y[0]) + geo_i["N"] * x[1] * y[1]

    c1, c2 = comps(e1), comps(e2)
    return II(c1, c1), II(c2, c2), II(c1, c2)


def check_joins(F, n=1201):
    """Separation, tangent-plane angle and curvature in a common frame along the four wing-to-cap joins."""
    u = np.concatenate([np.geomspace(1e-4, 0.02, 60), np.linspace(0.02, 0.98, n), 1 - np.geomspace(0.02, 1e-4, 60)])
    out = {}
    joins = {}
    for cap, side in (("right_cap", +1), ("left_cap", -1)):
        for skin, zs in (("upper_skin", +1), ("lower_skin", -1)):
            fs, fc = F[skin], F[cap]
            jet = parent.skin_jet(u, side * parent.ATTACHMENT, zs)
            q = jet["point"]
            # locate the join on both reimported faces by projection
            ps, pc = [], []
            for x in q:
                ps.append(fs.project(x)); pc.append(fc.project(x))
            ps, pc = np.array(ps), np.array(pc)
            ds, dc = fs.pts(ps[:, 0], ps[:, 1]), fc.pts(pc[:, 0], pc[:, 1])
            gs, gc = geometry(ds, fs.sign), geometry(dc, fc.sign)
            sep = np.linalg.norm(ds["p"] - dc["p"], axis=-1)
            ang = np.degrees(np.arccos(np.clip((gs["n"] * gc["n"]).sum(-1), -1, 1)))
            e1 = jet["d_chord"] / np.linalg.norm(jet["d_chord"], axis=-1)[:, None]
            nrm = gs["n"]
            e1 = e1 - (e1 * nrm).sum(-1)[:, None] * nrm
            e1 /= np.linalg.norm(e1, axis=-1)[:, None]
            e2 = np.cross(nrm, e1)
            ks = frame_curvatures(gs, ds, e1, e2)
            kc = frame_curvatures(gc, dc, e1, e2)
            dk = [np.abs(a - b) for a, b in zip(ks, kc)]
            scale = np.maximum(np.abs(ks[0]), 1e-3)
            key = f"{cap}/{skin}"
            joins[key] = {"u": u, "skin": (fs, ps, ds, gs), "cap": (fc, pc, dc, gc), "e1": e1, "e2": e2}
            inner = (u > 1e-3) & (u < 1 - 1e-3)
            out[key] = {
                "samples": int(len(u)), "u_range": [float(u.min()), float(u.max())],
                "max_separation_mm": float(sep.max()),
                "max_projection_distance_mm": float(max(ps[:, 2].max(), pc[:, 2].max())),
                "max_tangent_plane_angle_deg": float(np.nanmax(ang)),
                "max_tangent_plane_angle_deg_at_u": float(u[np.nanargmax(ang)]),
                "max_abs_dk_chordwise_1_per_mm": float(np.nanmax(dk[0])),
                "max_rel_dk_chordwise": float(np.nanmax(dk[0] / scale)),
                "max_abs_dk_spanwise_1_per_mm": float(np.nanmax(dk[1])),
                "max_abs_dk_twist_1_per_mm": float(np.nanmax(dk[2])),
                "max_abs_dk_spanwise_1_per_mm_u_in_1e-3_to_1-1e-3": float(np.nanmax(dk[1][inner])),
                "skin_spanwise_curvature_max_abs": float(np.nanmax(np.abs(ks[1]))),
                "cap_spanwise_curvature_max_abs": float(np.nanmax(np.abs(kc[1]))),
                "cap_param_v_on_join": [float(pc[:, 1].min()), float(pc[:, 1].max())],
                "cap_param_u_minus_native_u_max": float(np.abs(pc[:, 0] - u).max()),
            }
    return out, joins


def trace_x(f, X0, v):
    """Points of the plane section X = X0 on face f, one per v (Newton in u)."""
    ug = np.linspace(0, 1, 4001)
    xg = f.grid(ug, v)["p"][..., 0]
    u = ug[np.argmin(np.abs(xg - X0), axis=0)]
    for _ in range(60):
        d = f.pts(u, v)
        step = (d["p"][:, 0] - X0) / np.where(np.abs(d["su"][:, 0]) > 1e-300, d["su"][:, 0], 1e-300)
        u = np.clip(u - step, 0, 1)
        if np.abs(step).max() < 1e-15:
            break
    return u, f.pts(u, v)


def section_curvature(f, d, T):
    """Plane-section curvature by Meusnier from CAD derivatives, for a plane with unit normal T.

    Returns kappa_n (normal curvature of the surface along the trace, signed: > 0 convex),
    kappa_plane (curvature of the trace inside its plane, same sign) and n.m, the cosine
    between the surface normal and the trace's in-plane normal (1 = normal section, -> 0 = grazing).
    """
    g = geometry(d, f.sign)
    n = g["n"]
    t = np.cross(T, n); t /= np.linalg.norm(t, axis=-1)[:, None]
    m = np.cross(T, t)
    su, sv = d["su"], d["sv"]
    E, F_, G = g["E"], g["F"], g["G"]
    det = E * G - F_ * F_
    a = (t * su).sum(-1); b = (t * sv).sum(-1)
    cu = (G * a - F_ * b) / det; cv = (E * b - F_ * a) / det
    II = g["L"] * cu * cu + 2 * g["M"] * cu * cv + g["N"] * cv * cv
    kn = -II                      # outward normal: a convex surface bends away from it
    nm = np.abs((n * m).sum(-1))
    return kn, kn / np.maximum(nm, 1e-300), nm


def check_sections(F):
    """Chord-normal sections X = const (the wrap guides and between them), skin -> lip -> skin."""
    fr = np.array([1.7e-5, 1e-4, 1e-3, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
                   0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999,
                   0.99975, 0.9999, 0.99999])
    out = []
    TX = np.array([1.0, 0, 0])
    for cap, side in (("right_cap", +1), ("left_cap", -1)):
        f = F[cap]
        v = np.linspace(0, 1, 6001)
        for xi in fr:
            X0 = xi * parent.CHORD
            u, d = trace_x(f, X0, v)
            kn, kp, nm = section_curvature(f, d, TX)
            i = int(np.argmax(np.abs(kp)))
            vv = np.linspace(v[max(i - 3, 0)], v[min(i + 3, len(v) - 1)], 2001)
            u2, d2 = trace_x(f, X0, vv)
            kn2, kp2, nm2 = section_curvature(f, d2, TX)
            j = int(np.argmax(np.abs(kp2)))
            jn = int(np.argmax(np.abs(kn)))
            i_lip = int(np.argmax(np.abs(d["p"][:, 1])))
            p_at = d2["p"][j]
            g_at = geometry(f.pts(np.array([u2[j]]), np.array([vv[j]])), f.sign)
            Rprin = 1 / max(abs(g_at["k1"][0]), abs(g_at["k2"][0]))
            t = float(np.nan_to_num(parent.section_thickness(xi, parent.ATTACHMENT), nan=0.0))
            req = RULE_MIN_R if t >= RULE_T else t / 2
            sig = np.where(np.abs(kp) > 1e-6 * np.abs(kp).max(), np.sign(kp), 0)
            sig = sig[sig != 0]
            out.append({"cap": cap, "chord_fraction": float(xi), "X_mm": float(X0),
                        "attachment_thickness_mm": t, "required_min_radius_mm": req,
                        "in_plane_smallest_radius_mm": float(1 / abs(kp2[j])),
                        "in_plane_smallest_radius_at_YZ": [float(p_at[1]), float(p_at[2])],
                        "plane_vs_surface_cos_at_that_point": float(nm2[j]),
                        "principal_radius_at_that_point_mm": float(Rprin),
                        "in_plane_lip_radius_mm": float(1 / abs(kp[i_lip])), "lip_YZ": [float(d["p"][i_lip, 1]), float(d["p"][i_lip, 2])],
                        "lip_plane_vs_surface_cos": float(nm[i_lip]),
                        "normal_curvature_smallest_radius_along_trace_mm": float(1 / abs(kn[jn])),
                        "min_signed_curvature_1_per_mm": float(kp.min()), "bending_reversals": int((np.diff(sig) != 0).sum()),
                        "max_X_drift_mm": float(np.abs(d["p"][:, 0] - X0).max()),
                        "pass_in_plane": bool(1 / abs(kp2[j]) >= req * (1 - 1e-9))})
    return out


def trace_plane(f, O, T, window):
    """Section of face f by the plane through O with unit normal T; the branch through O within 'window' mm."""
    import contourpy
    s = np.linspace(0, 1, 1201)
    ug = 0.5 - 0.5 * np.cos(np.pi * s)
    vg = np.linspace(*f.vrange, 1201)
    G = ((f.grid(ug, vg)["p"] - O) @ T)
    lines = contourpy.contour_generator(x=vg, y=ug, z=G).lines(0.0)
    best, bi, bd = None, None, np.inf
    for ln in lines:
        P = f.pts(ln[:, 1], ln[:, 0])["p"]
        dist = np.linalg.norm(P - O, axis=-1)
        k = int(dist.argmin())
        if dist[k] < bd:
            best, bi, bd = ln, k, dist[k]
    u, v = best[:, 1].copy(), best[:, 0].copy()
    # densify x8 by linear interpolation in parameter space, then project onto the plane by Newton
    tt = np.linspace(0, len(u) - 1, 8 * (len(u) - 1) + 1)
    u = np.interp(tt, np.arange(len(u)), u); v = np.interp(tt, np.arange(len(v)), v)
    for _ in range(8):
        d = f.pts(u, v)
        g = (d["p"] - O) @ T
        gu, gv = d["su"] @ T, d["sv"] @ T
        nn = gu * gu + gv * gv
        u = np.clip(u - g * gu / nn, *f.urange); v = np.clip(v - g * gv / nn, *f.vrange)
    d = f.pts(u, v)
    plan = np.linalg.norm((d["p"] - O)[:, :2], axis=-1)
    k = int(np.argmin(np.linalg.norm(d["p"] - O, axis=-1)))
    lo = k
    while lo > 0 and plan[lo - 1] <= window:
        lo -= 1
    hi = k
    while hi < len(u) - 1 and plan[hi + 1] <= window:
        hi += 1
    sl = slice(lo, hi + 1)
    return u[sl], v[sl], {kk: vv_[sl] for kk, vv_ in d.items()}, float(np.abs((d["p"] - O) @ T).max())


def check_cross_lip_sections(F):
    """Normal sections across the lip: vertical planes through outline points, normal to the outline."""
    ps = np.array([0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95,
                   0.98, 0.99, 0.995, 0.998])
    out = []
    for cap, side in (("right_cap", +1), ("left_cap", -1)):
        f = F[cap]
        for p in ps:
            jet = parent.tip_outline_jet(p, side)
            O = jet["point"]; T = jet["d_p"] / np.linalg.norm(jet["d_p"])
            window = max(200 * p * (1 - p), 3.0)
            u, v, d, resid = trace_plane(f, O, T, window)
            kn, kp, nm = section_curvature(f, d, T)
            g = geometry(d, f.sign)
            Rp = 1 / np.maximum(np.abs(g["k1"]), np.abs(g["k2"]))
            req, t = required_radius(d["p"][:, 0])
            j = int(np.argmax(np.abs(kp)))
            jr = int(np.argmin(np.abs(kp) ** -1 / req))
            ilip = int(np.argmin(np.linalg.norm(d["p"] - O, axis=-1)))
            sig = np.where(np.abs(kp) > 1e-6 * np.abs(kp).max(), np.sign(kp), 0)
            sig = sig[sig != 0]
            tO = float(np.nan_to_num(parent.section_thickness(O[0] / parent.CHORD, parent.ATTACHMENT), nan=0.0))
            out.append({"cap": cap, "outline_p": float(p), "lip_point_mm": O.tolist(), "window_plan_mm": window,
                        "trace_points": int(len(u)), "max_off_plane_mm": resid,
                        "min_plane_vs_surface_cos": float(nm.min()),
                        "smallest_radius_mm": float(1 / abs(kp[j])), "smallest_radius_at_XYZ": d["p"][j].tolist(),
                        "attachment_thickness_there_mm": float(t[j]), "required_there_mm": float(req[j]),
                        "worst_radius_over_required": float(1 / abs(kp[jr]) / req[jr]),
                        "lip_radius_mm": float(1 / abs(kp[ilip])), "attachment_thickness_at_lip_chord_mm": tO,
                        "lip_radius_over_half_thickness": float(1 / abs(kp[ilip]) / (tO / 2)) if tO > 0 else None,
                        "principal_min_radius_on_trace_mm": float(Rp.min()),
                        "min_signed_curvature_1_per_mm": float(kp.min()), "bending_reversals": int((np.diff(sig) != 0).sum()),
                        "pass": bool(1 / abs(kp[jr]) >= req[jr] * (1 - 1e-9))})
    return out


def principal_map(F, nu=1400, nv=1401):
    """Larger principal curvature over each cap on a dense grid, refined at local extremes, judged by the rule."""
    res = {}
    for cap in ("right_cap", "left_cap"):
        f = F[cap]
        s = np.linspace(0, 1, nu)
        u = 0.5 - 0.5 * np.cos(np.pi * s)                 # cluster towards both collapsed ends
        u = u[(u > 2e-4) & (u < 1 - 2e-4)]
        v = np.linspace(*f.vrange, nv)
        worst_ratio = np.full((len(u), len(v)), np.inf)
        kmax = np.zeros((len(u), len(v)))
        X = np.zeros((len(u), len(v)))
        for a in range(0, len(u), 100):
            d = f.grid(u[a:a + 100], v)
            g = geometry(d, f.sign)
            km = np.maximum(np.abs(g["k1"]), np.abs(g["k2"]))
            kmax[a:a + 100] = km
            X[a:a + 100] = d["p"][..., 0]
        req, t = required_radius(X)
        Rm = 1 / np.maximum(kmax, 1e-300)
        ratio = Rm / req
        # cross-check the NumPy evaluator against OCCT at random grid points
        rng = np.random.default_rng(7)
        diff = 0.0
        for _ in range(200):
            i, j = rng.integers(len(u)), rng.integers(len(v))
            _, n, k = f.occ(u[i], v[j])
            if k is not None:
                diff = max(diff, abs(np.abs(k).max() - kmax[i, j]) / max(kmax[i, j], 1e-9))
        # local minima of the ratio and of the radius, refined with OCCT
        cand = []
        for arr in (ratio, Rm):
            flat = np.argsort(arr, axis=None)[:4000]
            for idx in flat:
                i, j = np.unravel_index(idx, arr.shape)
                if all(abs(i - a) > 20 or abs(j - b) > 20 for a, b in cand):
                    cand.append((i, j))
                if len(cand) > 12:
                    break

        def occ_ratio(x):
            uu = np.clip(x[0], u[0], u[-1]); vv = np.clip(x[1], *f.vrange)
            p, n, k = f.occ(uu, vv)
            if k is None:
                return np.inf, None
            req_, t_ = required_radius(p[0])
            return float(1 / np.abs(k).max() / req_), (p, k, float(t_), float(req_), uu, vv)

        refined = []
        for i, j in cand:
            r = minimize(lambda x: occ_ratio(x)[0], [u[i], v[j]], method="Nelder-Mead",
                         options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 400})
            val, info = occ_ratio(r.x)
            p, k, t_, req_, uu, vv = info
            refined.append({"u": float(uu), "v": float(vv), "XYZ_mm": p.tolist(), "radius_mm": float(1 / np.abs(k).max()),
                            "principal_curvatures_1_per_mm": k.tolist(), "attachment_thickness_mm": t_,
                            "required_mm": req_, "radius_over_required": val})
        refined.sort(key=lambda r: r["radius_over_required"])
        iR = np.unravel_index(np.argmin(Rm), Rm.shape)
        res[cap] = {"grid": [int(len(u)), int(len(v))], "u_closest_to_collapsed_ends": [float(u[0]), float(1 - u[-1])],
                    "X_closest_to_ends_mm": [float(X[0].min()), float(parent.CHORD - X[-1].max())],
                    "grid_min_radius_mm": float(Rm.min()), "grid_min_radius_at_XYZ": f.grid(u[[iR[0]]], v[[iR[1]]])["p"][0, 0].tolist(),
                    "grid_min_radius_over_required": float(ratio.min()),
                    "numpy_vs_occt_max_rel_curvature_diff": float(diff),
                    "refined_extremes": refined[:6]}
    return res


# ---------------------------------------------------------------- zebra
def stripe_phase(n, V, axis, n_stripes):
    """Reflection-line phase (in stripes) for view direction V and stripe axis 'axis'."""
    R = V - 2 * (n @ V)[..., None] * n
    a = axis / np.linalg.norm(axis)
    e1 = np.cross(a, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    return np.arctan2(R @ e2, R @ e1) * n_stripes / np.pi, R, (e1, e2)


def phase_gradient(n, dn, V, axis, n_stripes):
    """Derivative of the phase along a surface direction with normal derivative dn."""
    a = axis / np.linalg.norm(axis)
    e1 = np.cross(a, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    nV = (n * V).sum(-1)
    R = V - 2 * nV[:, None] * n
    dR = -2 * ((dn * V).sum(-1)[:, None] * n + nV[:, None] * dn)
    x, y = R @ e1, R @ e2
    dx, dy = dR @ e1, dR @ e2
    return (x * dy - y * dx) / (x * x + y * y) * n_stripes / np.pi


def dn_along(g, d, e):
    """Derivative of the unit normal along unit tangent e (Weingarten)."""
    su, sv = d["su"], d["sv"]
    E, F, G = g["E"], g["F"], g["G"]
    det = E * G - F * F
    a = (e * su).sum(-1); b = (e * sv).sum(-1)
    cu = (G * a - F * b) / det; cv = (E * b - F * a) / det
    # dn/du = (FM - GL)/det su + (FL - EM)/det sv, dn/dv = (FN - GM)/det su + (FM - EN)/det sv
    L, M, N = g["L"], g["M"], g["N"]
    dnu = ((F * M - G * L) / det)[:, None] * su + ((F * L - E * M) / det)[:, None] * sv
    dnv = ((F * N - G * M) / det)[:, None] * su + ((F * M - E * N) / det)[:, None] * sv
    return cu[:, None] * dnu + cv[:, None] * dnv


def check_zebra(joins, n_stripes=24):
    views = {"top": [0, 0, -1], "bottom": [0, 0, 1], "from_front": [1, 0, 0], "from_behind": [-1, 0, 0],
             "oblique_front_above": [0.6, 0.0, -0.8], "oblique_outboard_above": [0.3, 0.45, -0.84],
             "oblique_rear_below": [-0.5, 0.3, 0.81]}
    axes = {"X": [1, 0, 0], "Y": [0, 1, 0], "Z": [0, 0, 1], "XY45": [1, 1, 0], "XZ45": [1, 0, 1]}
    worst = {"phase_mismatch_stripes": (0, None), "gradient_mismatch_rel": (0, None), "crossing_direction_change_deg": (0, None),
             "gradient_mismatch_stripes_per_mm": (0, None)}
    per = []
    crossings_total = 0
    for key, J in joins.items():
        side = 1 if key.startswith("right") else -1
        u = J["u"]
        m = (u > 1e-3) & (u < 1 - 1e-3)
        fs, ps, ds, gs = J["skin"]
        fc, pc, dc, gc = J["cap"]
        e1, e2 = J["e1"][m], J["e2"][m]
        sl = lambda dct: {k: v[m] for k, v in dct.items()}
        dsm, dcm, gsm, gcm = sl(ds), sl(dc), sl(gs), sl(gc)
        for vn, V in views.items():
            V = np.array(V, float); V[1] *= side; V /= np.linalg.norm(V)       # mirrored camera
            for an, A in axes.items():
                A = np.array(A, float); A[1] *= side                            # mirrored stripe lighting
                front = ((gsm["n"] @ V) < -0.05) & ((gcm["n"] @ V) < -0.05)     # both faces seen from the front
                if front.sum() < 10:
                    continue
                ph_s = stripe_phase(gsm["n"], V, A, n_stripes)[0]
                ph_c = stripe_phase(gcm["n"], V, A, n_stripes)[0]
                dph = np.abs((ph_s - ph_c + n_stripes) % (2 * n_stripes) - n_stripes)
                grads = {}
                for nm, g, d in (("s", gsm, dsm), ("c", gcm, dcm)):
                    g1 = phase_gradient(g["n"], dn_along(g, d, e1), V, A, n_stripes)
                    g2 = phase_gradient(g["n"], dn_along(g, d, e2), V, A, n_stripes)
                    grads[nm] = np.stack([g1, g2], -1)
                gs_n = np.linalg.norm(grads["s"], axis=-1); gc_n = np.linalg.norm(grads["c"], axis=-1)
                gmax = np.maximum(gs_n, gc_n)
                floor = max(1e-2 * np.percentile(gmax[front], 95), GRAD_FLOOR)
                gabs = np.linalg.norm(grads["s"] - grads["c"], axis=-1)
                gmis = np.where(front, gabs / np.maximum(gmax, floor), 0.0)
                gabs = np.where(front, gabs, 0.0)
                dph = np.where(front, dph, 0.0)
                # light/dark crossings: phase passing through an integer along the join
                fl = np.floor(ph_s)
                cross = np.nonzero((np.diff(fl) != 0) & front[:-1] & front[1:])[0]
                crossings_total += len(cross)
                dirchg = 0.0
                at = None
                for i in cross:
                    if min(gs_n[i], gc_n[i]) < floor:
                        continue
                    a1 = np.arctan2(grads["s"][i, 1], grads["s"][i, 0])
                    a2 = np.arctan2(grads["c"][i, 1], grads["c"][i, 0])
                    dd = abs((np.degrees(a1 - a2) + 180) % 360 - 180)
                    if dd > dirchg:
                        dirchg, at = dd, float(u[m][i])
                rec = {"join": key, "view": vn, "axis": an, "max_phase_mismatch_stripes": float(dph.max()),
                       "max_gradient_mismatch_rel": float(gmis.max()), "max_gradient_mismatch_stripes_per_mm": float(gabs.max()), "crossings": int(len(cross)),
                       "max_crossing_direction_change_deg": float(dirchg), "at_u": at}
                per.append(rec)
                for k_, val, arr in (("phase_mismatch_stripes", dph.max(), dph), ("gradient_mismatch_rel", gmis.max(), gmis),
                                     ("gradient_mismatch_stripes_per_mm", gabs.max(), gabs),
                                     ("crossing_direction_change_deg", dirchg, None)):
                    if val > worst[k_][0]:
                        worst[k_] = (float(val), {"join": key, "view": vn, "axis": an,
                                                  "u": float(u[m][int(np.argmax(arr))]) if arr is not None else at})
    evaluated = len(per)
    return {"combinations_evaluated": evaluated, "n_stripes_per_180deg": n_stripes, "views": list(views), "stripe_axes": list(axes),
            "sampled": "each join at native u in [1e-3, 1 - 1e-3] (collapsed LE/TE corners excluded), only where both faces face the camera (n.V < -0.05); gradient ratios use a floor of max(1% of the 95th-percentile gradient on that join, 1e-3 stripes/mm); crossing directions are compared only where both gradients exceed that floor",
            "light_dark_crossings_checked": crossings_total,
            "worst": {k: {"value": v[0], "where": v[1]} for k, v in worst.items()}}, per


# ---------------------------------------------------------------- rendering
def render(F, view, up, center, width_mm, px, out, axis, n_stripes=24, faces=None, title="", skin_y=None, max_samples=5000):
    """Orthographic z-buffer splat of CAD-normal zebra stripes (no edge lines)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    V = np.array(view, float); V /= np.linalg.norm(V)
    upv = np.array(up, float); upv -= (upv @ V) * V; upv /= np.linalg.norm(upv)
    right = np.cross(V, upv)
    W, H = px
    scale = W / width_mm
    zbuf = np.full((H, W), np.inf); img = np.full((H, W), np.nan)
    def to_screen(p):
        p = p - center
        return (p @ right) * scale + W / 2, H / 2 - (p @ upv) * scale, p @ V

    for name in faces or F:
        f = F[name]
        vr = f.vrange
        if "skin" in name and skin_y is not None:
            vr = (max(vr[0], skin_y[0]), min(vr[1], skin_y[1]))
            if vr[0] >= vr[1]:
                continue
        # parameter window that lands in the image, then a density that keeps samples < 0.6 px apart
        uc, vc = np.linspace(*f.urange, 401), np.linspace(*vr, 401)
        x, y, _ = to_screen(f.grid(uc, vc)["p"])
        inside = (x > -20) & (x < W + 20) & (y > -20) & (y < H + 20)
        if not inside.any():
            continue
        iu, iv = np.nonzero(inside)
        ur = (uc[max(iu.min() - 1, 0)], uc[min(iu.max() + 1, 400)])
        vr = (vc[max(iv.min() - 1, 0)], vc[min(iv.max() + 1, 400)])
        uc, vc = np.linspace(*ur, 201), np.linspace(*vr, 201)
        x, y, _ = to_screen(f.grid(uc, vc)["p"])
        su = np.hypot(np.diff(x, axis=0), np.diff(y, axis=0)).max()
        sv = np.hypot(np.diff(x, axis=1), np.diff(y, axis=1)).max()
        nu = int(min(max(200 * su / 0.6, 50), max_samples))
        nv = int(min(max(200 * sv / 0.6, 50), max_samples))
        uu_all = np.linspace(*ur, nu)
        v = np.linspace(*vr, nv)
        for a in range(0, nu, 64):
            u = uu_all[a:a + 64]
            d = f.grid(u, v)
            g = geometry(d, f.sign)
            n = g["n"].reshape(-1, 3)
            x, y, z = to_screen(d["p"].reshape(-1, 3))
            vis = (n @ V) < 0
            ok = vis & (x >= 0) & (x < W) & (y >= 0) & (y < H) & np.isfinite(n).all(1)
            xi, yi, z, n = x[ok].astype(int), y[ok].astype(int), z[ok], n[ok]
            ph = stripe_phase(n, V, np.array(axis, float), n_stripes)[0]
            val = (np.floor(ph) % 2) * 0.8 + 0.1 + 0.1 * np.abs(n @ V)
            order = np.lexsort((z, yi, xi))
            xi, yi, z, val = xi[order], yi[order], z[order], val[order]
            first = np.ones(len(xi), bool); first[1:] = (np.diff(xi) != 0) | (np.diff(yi) != 0)
            xi, yi, z, val = xi[first], yi[first], z[first], val[first]
            closer = z < zbuf[yi, xi]
            zbuf[yi[closer], xi[closer]] = z[closer]; img[yi[closer], xi[closer]] = val[closer]
    fig = plt.figure(figsize=(W / 100, H / 100 + 0.3), dpi=100)
    ax = fig.add_axes([0, 0, 1, H / (H + 30)])
    ax.imshow(img, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    ax.set_axis_off()
    fig.text(0.01, 0.995, title, va="top", fontsize=8)
    fig.savefig(out, dpi=100, facecolor="white")
    plt.close(fig)


def renders_only():
    shape = read_step(STEP)
    F = classify([Face(f) for f in explore(shape, TopAbs_FACE)])
    return F


def main():
    data = STEP.read_bytes()
    shape = read_step(STEP)
    rb = {"file": "wing.step", "step_sha256": hashlib.sha256(data).hexdigest(), "step_bytes": len(data)}
    rb["reimport"] = check_topology(shape)
    F = classify([Face(f) for f in explore(shape, TopAbs_FACE)])
    rb["faces_found"] = sorted(F)
    rb["parent_skin_preserved"] = check_skin_preserved(F)
    ex = extents(F)
    rb["locked_dimensions"] = {
        "note": "extremes evaluated on the CAD faces (OCCT's optimal bounding box is only approximate on the left tip)",
        "chord_X_extent_mm": ex["X"], "finished_span_mm": ex["Y"][1] - ex["Y"][0],
        "tip_Y_mm": ex["Y"], "thickness_Z_extent_mm": ex["Z"],
        "attachment_Y_of_joins_mm": [float(min(F["left_cap"].P[:, 0, 1].min(), F["left_cap"].P[:, -1, 1].min())),
                                      float(max(F["right_cap"].P[:, 0, 1].max(), F["right_cap"].P[:, -1, 1].max()))],
        "parent_volume_mm3": 864000.0, "cap_volume_each_mm3": (rb["reimport"]["volume_mm3"] - 864000.0) / 2}
    rb["outline"] = check_outline(F)
    rb["joins"], joins = check_joins(F)
    rb["sections_chord_normal"] = check_sections(F)
    rb["sections_across_lip"] = check_cross_lip_sections(F)
    rb["principal_curvature_map"] = principal_map(F)
    rb["zebra"], per = check_zebra(joins)
    (HERE / "zebra_join_samples.json").write_text(json.dumps(per, indent=0))
    # summary of the printable radius rule
    secs = rb["sections_chord_normal"]
    xl = rb["sections_across_lip"]
    worst_sec = min(secs, key=lambda s: s["in_plane_smallest_radius_mm"] / s["required_min_radius_mm"])
    worst_xl = min(xl, key=lambda s: s["worst_radius_over_required"])
    pm = [e for c in rb["principal_curvature_map"].values() for e in c["refined_extremes"]]
    worst_pm = min(pm, key=lambda e: e["radius_over_required"])
    smallest_pm = min(pm, key=lambda e: e["radius_mm"])
    rb["printable_radius_rule"] = {
        "principal_smallest_radius": smallest_pm, "principal_worst_vs_rule": worst_pm,
        "across_lip_sections_worst_vs_rule": worst_xl,
        "across_lip_sections_failing": [s for s in xl if not s["pass"]],
        "chord_normal_sections_worst_in_plane": worst_sec,
        "chord_normal_sections_below_rule_in_plane": [
            {k: s[k] for k in ("cap", "chord_fraction", "attachment_thickness_mm", "in_plane_smallest_radius_mm",
                               "required_min_radius_mm", "plane_vs_surface_cos_at_that_point", "principal_radius_at_that_point_mm")}
            for s in secs if not s["pass_in_plane"]]}
    (HERE / "readback.json").write_text(json.dumps(rb, indent=1))

    render(F, [0, 0, -1], [-1, 0, 0], np.array([60, 0, 0]), 840, (1680, 300), HERE / "zebra_top_full.png", [0, 1, 0],
           title="wing.step reimported, top view, centred camera, stripes about the span axis (mirror-symmetric lighting)")
    for side, nm in ((1, "right"), (-1, "left")):
        sy = (300, 351) if side > 0 else (-351, -300)
        render(F, [0, 0, -1], [-1, 0, 0], np.array([60, side * 370, 0]), 150, (900, 840),
               HERE / f"zebra_top_{nm}_tip.png", [0, 1, 0], title=f"{nm} tip, top view, stripes about the span axis", skin_y=sy)
        render(F, [0.35, side * 0.55, -0.76], [0, 0, 1], np.array([60, side * 380, 0]), 150, (900, 700),
               HERE / f"zebra_oblique_{nm}_tip.png", [0, side * 1, 0.4],
               title=f"{nm} tip, oblique from outboard/front/above, mirrored camera and stripes", skin_y=sy)
    render(F, [1, 0, 0], [0, 0, 1], np.array([0, 375, 0]), 60, (1200, 400), HERE / "zebra_front_right_tip.png", [0, 1, 0],
           title="right tip seen from ahead (looking aft), stripes about the span axis", skin_y=(330, 351))
    render(F, [-1, 0, 0], [0, 0, 1], np.array([120, 375, 0]), 60, (1200, 400), HERE / "zebra_rear_right_tip.png", [0, 1, 0],
           title="right tip seen from behind (looking forward), stripes about the span axis", skin_y=(330, 351))
    print(json.dumps({"sha256": rb["step_sha256"], "reimport": rb["reimport"],
                      "rule": rb["printable_radius_rule"], "zebra": rb["zebra"]["worst"]}, indent=1))


if __name__ == "__main__":
    main()
