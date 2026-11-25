"""Verify wing.step (fin5 / sharp) on the reimported STEP and write readback.json.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/fin5/sharp/verify_cap.py

Every measurement below is taken on the geometry read back from wing.step.
Curvature uses the convention "convex positive" with outward normals.
"""
import hashlib
import json
import pathlib
import sys

import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cap_design as cd  # noqa: E402

par = cd.par
from OCP.STEPControl import STEPControl_Reader
from OCP.IFSelect import IFSelect_RetDone
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID, TopAbs_SHELL, TopAbs_EDGE, TopAbs_REVERSED
from OCP.TopoDS import TopoDS
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface, BRepAdaptor_Curve
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Plane
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds, ShapeAnalysis_Surface
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.gp import gp_Pln, gp_Pnt, gp_Dir, gp_Vec
from OCP.GeomLProp import GeomLProp_SLProps
from OCP.TopTools import TopTools_IndexedMapOfShape
from OCP.TopExp import TopExp

STEP = HERE / "wing.step"
OUT = HERE / "readback.json"
RULE_MIN_R = 0.4
RULE_T = 0.8
V_POLE_CUT = 1 - 1e-4   # closest approach to the collapsed crown pole for derivative-based measures


def rnd(x, n=6):
    if isinstance(x, dict):
        return {k: rnd(v, n) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, n) for v in x]
    if isinstance(x, (float, np.floating)):
        return float(f"{float(x):.{n}g}")
    if isinstance(x, np.integer):
        return int(x)
    return x


# ------------------------------------------------------------------ surfaces
class Surf:
    """Tensor-product B-spline evaluator built from a reimported face (non-rational)."""

    def __init__(self, face):
        self.face = face
        self.geom = BRep_Tool.Surface_s(face)
        a = BRepAdaptor_Surface(face)
        s = a.BSpline()
        assert not s.IsURational() and not s.IsVRational()
        self.bs = s
        self.p, self.q = s.UDegree(), s.VDegree()
        self.tu = np.repeat([s.UKnot(i) for i in range(1, s.NbUKnots() + 1)],
                            [s.UMultiplicity(i) for i in range(1, s.NbUKnots() + 1)])
        self.tv = np.repeat([s.VKnot(i) for i in range(1, s.NbVKnots() + 1)],
                            [s.VMultiplicity(i) for i in range(1, s.NbVKnots() + 1)])
        self.P = np.array([[[s.Pole(i, j).X(), s.Pole(i, j).Y(), s.Pole(i, j).Z()]
                            for j in range(1, s.NbVPoles() + 1)] for i in range(1, s.NbUPoles() + 1)])
        self.sign = -1.0 if face.Orientation() == TopAbs_REVERSED else 1.0
        self.u0, self.u1 = self.tu[0], self.tu[-1]
        self.v0, self.v1 = self.tv[0], self.tv[-1]
        nu, nv = self.P.shape[:2]
        self._bu = [BSpline(self.tu, np.eye(nu), self.p).derivative(k) if k else BSpline(self.tu, np.eye(nu), self.p) for k in range(3)]
        self._bv = [BSpline(self.tv, np.eye(nv), self.q).derivative(k) if k else BSpline(self.tv, np.eye(nv), self.q) for k in range(3)]

    def jet(self, u, v):
        """Point and derivatives on the product grid u x v (arrays) -> dict of (len(u), len(v), 3)."""
        u = np.atleast_1d(np.asarray(u, float))
        v = np.atleast_1d(np.asarray(v, float))
        Bu = [b(u) for b in self._bu]
        Bv = [b(v) for b in self._bv]
        def ev(i, j):
            T = np.tensordot(Bu[i], self.P, axes=(1, 0))          # (a, nv, 3)
            return np.transpose(np.tensordot(T, Bv[j], axes=(1, 1)), (0, 2, 1))  # (a, b, 3)
        return {"P": ev(0, 0), "Su": ev(1, 0), "Sv": ev(0, 1), "Suu": ev(2, 0), "Suv": ev(1, 1), "Svv": ev(0, 2)}

    def jet_pts(self, u, v):
        """Same for paired points (u[i], v[i]) -> (n, 3)."""
        u = np.atleast_1d(np.asarray(u, float))
        v = np.atleast_1d(np.asarray(v, float))
        Bu = [b(u) for b in self._bu]
        Bv = [b(v) for b in self._bv]
        ev = lambda i, j: np.einsum("ijk,ni,nj->nk", self.P, Bu[i], Bv[j], optimize=True)
        return {"P": ev(0, 0), "Su": ev(1, 0), "Sv": ev(0, 1), "Suu": ev(2, 0), "Suv": ev(1, 1), "Svv": ev(0, 2)}


def geometry(j, sign):
    """Outward normal, fundamental forms, principal curvatures (convex positive), normal derivatives."""
    Su, Sv = j["Su"], j["Sv"]
    n = np.cross(Su, Sv)
    n = sign * n / np.linalg.norm(n, axis=-1, keepdims=True)
    E, F, G = (Su * Su).sum(-1), (Su * Sv).sum(-1), (Sv * Sv).sum(-1)
    L, M, N = -(j["Suu"] * n).sum(-1), -(j["Suv"] * n).sum(-1), -(j["Svv"] * n).sum(-1)
    det = E * G - F * F
    K = (L * N - M * M) / det
    H = (E * N - 2 * F * M + G * L) / (2 * det)
    d = np.sqrt(np.clip(H * H - K, 0, None))
    k1, k2 = H + d, H - d
    # n_u = a Su + b Sv with I [a b]^T = [L M]^T ; n_v from [M N]
    au, bu = (G * L - F * M) / det, (E * M - F * L) / det
    av, bv = (G * M - F * N) / det, (E * N - F * M) / det
    nu = au[..., None] * Su + bu[..., None] * Sv
    nv = av[..., None] * Su + bv[..., None] * Sv
    return {"n": n, "E": E, "F": F, "G": G, "L": L, "M": M, "N": N, "k1": k1, "k2": k2,
            "kmax": np.maximum(np.abs(k1), np.abs(k2)), "nu": nu, "nv": nv}


def normal_curvature(j, g, t):
    """Normal curvature (convex positive) along unit tangent t (..., 3)."""
    Su, Sv = j["Su"], j["Sv"]
    det = g["E"] * g["G"] - g["F"] ** 2
    tu, tv = (t * Su).sum(-1), (t * Sv).sum(-1)
    al = (g["G"] * tu - g["F"] * tv) / det
    be = (g["E"] * tv - g["F"] * tu) / det
    return (g["L"] * al * al + 2 * g["M"] * al * be + g["N"] * be * be) / (
        g["E"] * al * al + 2 * g["F"] * al * be + g["G"] * be * be)


def shape_matrix(j, g, e1, e2):
    """Second fundamental form in the orthonormal tangent frame (e1, e2): [[II11, II12], [II12, II22]]."""
    a = normal_curvature(j, g, e1)
    b = normal_curvature(j, g, e2)
    d = (e1 + e2) / np.sqrt(2)
    c = normal_curvature(j, g, d)
    off = c - 0.5 * (a + b)
    return np.stack([np.stack([a, off], -1), np.stack([off, b], -1)], -2)


# ------------------------------------------------------------------ helpers
def thickness_at(X, Z):
    """Attachment thickness at the sweep-following chord fraction of a point (X, Z)."""
    le = 65.0 * Z / 240.0
    c = 160.0 - 75.0 * Z / 240.0
    xi = np.clip((X - le) / c, 0.0, 1.0)
    return par.section_thickness(xi, 240.0), xi


def required_radius(t):
    return np.where(t >= RULE_T, RULE_MIN_R, t / 2)


def read():
    r = STEPControl_Reader()
    assert r.ReadFile(str(STEP)) == IFSelect_RetDone
    r.TransferRoots()
    return r.OneShape()


def count(shape, kind):
    m = TopTools_IndexedMapOfShape()
    TopExp.MapShapes_s(shape, kind, m)
    return m.Extent()


import time
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def main():
    data = STEP.read_bytes()
    out = {"file": str(STEP.relative_to(cd.ROOT)), "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
    shape = read()
    # ---------------------------------------------------------- topology / validity
    fb = ShapeAnalysis_FreeBounds(shape)
    ex = TopExp_Explorer(shape, TopAbs_SOLID)
    solid = TopoDS.Solid_s(ex.Current())
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(solid, props, 1e-9)
    sprops = GProp_GProps()
    BRepGProp.SurfaceProperties_s(solid, sprops, 1e-9)
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(solid, box, False, False)
    bb = box.Get()
    faces = []
    exf = TopExp_Explorer(shape, TopAbs_FACE)
    while exf.More():
        faces.append(TopoDS.Face_s(exf.Current()))
        exf.Next()
    types = [str(BRepAdaptor_Surface(f).GetType()).split(".")[-1] for f in faces]
    out["reimport"] = {
        "solids": count(shape, TopAbs_SOLID), "shells": count(shape, TopAbs_SHELL), "faces": len(faces),
        "face_types": types, "brepcheck_valid": BRepCheck_Analyzer(shape).IsValid(),
        "free_edges": count(fb.GetOpenWires(), TopAbs_EDGE) + count(fb.GetClosedWires(), TopAbs_EDGE),
        "volume_mm3": props.Mass(), "area_mm2": sprops.Mass(),
        "bounding_box_mm": list(bb),
    }

    # ---------------------------------------------------------- identify faces
    skins, caps, root = {}, {}, None
    for f in faces:
        a = BRepAdaptor_Surface(f)
        if a.GetType() == GeomAbs_Plane:
            root = f
            continue
        s = Surf(f)
        mid = s.jet([0.5 * (s.u0 + s.u1)], [0.5 * (s.v0 + s.v1)])["P"][0, 0]
        side = "upper" if mid[1] > 0 else "lower"
        (skins if s.v1 > 100 else caps)[side] = s
    assert set(skins) == set(caps) == {"upper", "lower"} and root is not None

    # ---------------------------------------------------------- parent preserved
    log("parent")
    import importlib  # noqa
    sides = {"upper": +1, "lower": -1}
    ss = np.linspace(0, 1, 401)
    zz = np.linspace(0, 240, 121)
    dev_py = 0.0
    for side, sg in sides.items():
        P = skins[side].jet(ss, zz)["P"]
        S, Zg = np.meshgrid(ss, zz, indexing="ij")
        ref = par.skin(S, Zg, sg)
        dev_py = max(dev_py, float(np.max(np.linalg.norm(P - ref, axis=-1))))
    rp = STEPControl_Reader()
    rp.ReadFile(str(cd.PARENT_DIR / "parent.step"))
    rp.TransferRoots()
    pshape = rp.OneShape()
    pole_dev = 0.0
    exp = TopExp_Explorer(pshape, TopAbs_FACE)
    while exp.More():
        f = TopoDS.Face_s(exp.Current())
        if BRepAdaptor_Surface(f).GetType() == GeomAbs_BSplineSurface:
            s = Surf(f)
            side = "upper" if s.P[16, 2, 1] > 0 else "lower"
            pole_dev = max(pole_dev, float(np.max(np.abs(s.P - skins[side].P))))
            assert np.array_equal(s.tu, skins[side].tu) and np.array_equal(s.tv, skins[side].tv)
        exp.Next()
    ra = BRepAdaptor_Surface(root)
    out["parent_preserved"] = {
        "skin_vs_parent_py_max_mm": dev_py,
        "skin_poles_vs_parent_step_max_mm": pole_dev,
        "skin_knots_identical_to_parent_step": True,
        "root_face_plane_z_mm": ra.Plane().Location().Z(),
        "sample_grid": "401 chord parameters x 121 stations per skin",
    }

    # ---------------------------------------------------------- locked dimensions
    log("locked")
    crown = caps["upper"].jet(np.linspace(0, 1, 101), [1.0])["P"][:, 0]
    vv = np.unique(np.concatenate([np.linspace(0, 1, 4001), 1 - np.geomspace(1e-4, 1e-2, 400)]))
    outline_dev, outline_y = 0.0, 0.0
    for side in caps:
        for u, branch in ((caps[side].u0, -1), (caps[side].u1, +1)):
            P = caps[side].jet([u], vv)["P"][0]
            h = np.clip((P[:, 2] - 240.0) / 5.0, 0, 1)
            ref = par.tip_outline(h, branch)      # outline point at the same height
            ok = 1 - h >= 2e-8                     # closer to the crown, Z round-off (~5e-14 mm) dominates X(h)
            outline_dev = max(outline_dev, float(np.max(np.abs(P[ok, 0] - ref[ok, 0]))))
            outline_y = max(outline_y, float(np.max(np.abs(P[:, 1]))))
    samples = json.loads((cd.PARENT_DIR / "spec.json").read_text())["locked_outline"]["samples"]
    spec_dev = 0.0
    for key, u_attr in (("leading_edge_branch_XYZ", "u0"), ("trailing_edge_branch_XYZ", "u1")):
        for (x, y, z) in samples[key]:
            vq = 1 - np.sqrt(max(0.0, 1 - (z - 240.0) / 5.0))
            for side in caps:
                P = caps[side].jet([getattr(caps[side], u_attr)], [vq])["P"][0, 0]
                spec_dev = max(spec_dev, float(np.linalg.norm(P - [x, y, z])))
    att_z = 0.0
    for side in caps:
        P = caps[side].jet(np.linspace(0, 1, 2001), [0.0])["P"][:, 0]
        att_z = max(att_z, float(np.max(np.abs(P[:, 2] - 240.0))))
    out["locked_dimensions"] = {
        "finished_height_mm": bb[5],
        "crown_point_mm": crown[0].tolist(),
        "crown_row_spread_mm": float(np.max(np.linalg.norm(crown - crown[0], axis=-1))),
        "attachment_edge_max_abs_z_minus_240_mm": att_z,
        "root_chord_extent_x_mm": [bb[0], bb[3]],
        "outline_max_deviation_mm": outline_dev,
        "outline_max_abs_y_mm": outline_y,
        "outline_spec_samples_max_deviation_mm": spec_dev,
        "note": "outline deviation = max |X - X_outline(h)| at the same height Z (h = (Z - 240)/5), 4,400 heights per branch up to 1 - h = 2e-8 (Z within 0.1 micrometre of the crown; above that the comparison is limited by round-off in Z), both perimeter branches of both cap faces; the normal distance to the outline is never larger. The crown point itself is compared separately.",
        "crown_vs_outline_crown_mm": float(np.linalg.norm(crown[0] - par.tip_outline(1.0, -1))),
    }

    # ---------------------------------------------------------- joins
    log("joins")
    def join_metrics(sa, ua, va, sb, ub, vb, ref_tangent=None):
        ja, jb = sa.jet_pts(ua, va), sb.jet_pts(ub, vb)
        ga, gb = geometry(ja, sa.sign), geometry(jb, sb.sign)
        sep = np.linalg.norm(ja["P"] - jb["P"], axis=-1)
        ang = np.degrees(np.arccos(np.clip((ga["n"] * gb["n"]).sum(-1), -1, 1)))
        # common tangent frame: e1 along the join (from side a), e2 = n x e1
        e1 = ja["Su"] if ref_tangent == "u" else ja["Sv"]
        e1 = e1 - (e1 * ga["n"]).sum(-1, keepdims=True) * ga["n"]
        e1 /= np.linalg.norm(e1, axis=-1, keepdims=True)
        e2 = np.cross(ga["n"], e1)
        IIa = shape_matrix(ja, ga, e1, e2)
        IIb = shape_matrix(jb, gb, e1, e2)
        dII = np.abs(IIa - IIb)
        scale = np.maximum(np.abs(IIa).max(axis=(-1, -2)), 1e-12)
        return {"P": ja["P"], "sep": sep, "ang": ang, "dII": dII.max(axis=(-1, -2)), "rel": dII.max(axis=(-1, -2)) / scale,
                "IIa": IIa, "IIb": IIb, "ja": ja, "jb": jb, "ga": ga, "gb": gb}

    uj = np.unique(np.concatenate([np.linspace(0, 1, 2001), np.geomspace(1e-7, 0.01, 200), 1 - np.geomspace(1e-7, 0.01, 200)]))
    joins = {}
    proj_sep = 0.0
    for side in ("upper", "lower"):
        m = join_metrics(skins[side], uj, np.full_like(uj, 240.0), caps[side], uj, np.zeros_like(uj), ref_tangent="u")
        # independent check of correspondence: project skin points onto the cap surface
        sa = ShapeAnalysis_Surface(caps[side].geom)
        for k in range(0, len(uj), 40):
            p = gp_Pnt(*map(float, m["P"][k]))
            uv = sa.ValueOfUV(p, 1e-12)
            proj_sep = max(proj_sep, sa.Gap())
        i = int(np.argmax(m["ang"]))
        joins[f"{side}_skin_to_cap"] = {
            "points": len(uj), "max_separation_mm": float(m["sep"].max()), "max_projection_gap_mm": proj_sep,
            "max_tangent_plane_angle_deg": float(m["ang"].max()), "at_xyz": m["P"][i].tolist(),
            "max_second_fundamental_form_difference_per_mm": float(m["dII"].max()),
            "max_relative_curvature_difference": float(m["rel"].max()),
            "at_xyz_curvature": m["P"][int(np.argmax(m["dII"]))].tolist(),
        }
        joins[f"_{side}"] = m
    vp = np.unique(np.concatenate([np.linspace(0, 0.999, 2000), 1 - np.geomspace(1 - V_POLE_CUT, 1e-3, 100)]))
    for name, u in (("perimeter_leading_edge_branch", 0.0), ("perimeter_trailing_edge_branch", 1.0)):
        m = join_metrics(caps["upper"], np.full_like(vp, u), vp, caps["lower"], np.full_like(vp, u), vp, ref_tangent="v")
        i = int(np.argmax(m["ang"]))
        joins[name] = {
            "points": len(vp), "v_max": float(vp.max()),
            "closest_sample_to_crown_mm": float(np.linalg.norm(m["P"][-1] - caps["upper"].jet([0.0], [1.0])["P"][0, 0])),
            "max_separation_mm": float(m["sep"].max()),
            "max_tangent_plane_angle_deg": float(m["ang"].max()), "at_xyz": m["P"][i].tolist(),
            "max_second_fundamental_form_difference_per_mm": float(m["dII"].max()),
            "max_relative_curvature_difference": float(m["rel"].max()),
        }
        joins[f"_{name}"] = m
    out["joins"] = {k: v for k, v in joins.items() if not k.startswith("_")}
    out["joins"]["note"] = ("skin and cap are evaluated at the same chord parameter (the cap reuses the parent's u parameter); "
                            "projection gap is an independent check. Curvature: second fundamental form of both faces in a common "
                            "orthonormal tangent frame (e1 along the join). The perimeter is the shared edge of the +Y and -Y cap "
                            "patches; it ends at the collapsed crown pole and is sampled up to v = 1 - 1e-4 (7 micrometres from the crown); "
                            "closer than that the degenerate parametrization loses the second derivatives to round-off (crown curvature from the limit is reported under curvature_map).")

    # ---------------------------------------------------------- curvature map
    log("cmap")
    th = np.concatenate([np.geomspace(1e-6, 0.05, 80), np.linspace(0.05, np.pi - 0.05, 700), np.pi - np.geomspace(1e-6, 0.05, 80)[::-1]])
    ug = np.unique(np.concatenate([[0.0, 1.0], th / np.pi]))
    vg = np.unique(np.concatenate([np.linspace(0, 1 - 1e-3, 700), 1 - np.geomspace(1 - V_POLE_CUT, 1e-3, 60)]))
    cmap = {}
    for side in ("upper", "lower"):
        s = caps[side]
        j = s.jet(ug, vg)
        g = geometry(j, s.sign)
        R = 1 / g["kmax"]
        t, xi = thickness_at(j["P"][..., 0], j["P"][..., 2])
        req = required_radius(t)
        margin = np.divide(R, req, out=np.full_like(R, np.inf), where=req > 0)
        # refine the tightest spot in each chordwise band (local minima of R)
        def Rfun(x):
            u, v = np.clip(x, [0, 0], [1, V_POLE_CUT])
            jj = s.jet_pts([u], [v])
            return float(1 / geometry(jj, s.sign)["kmax"][0])
        refined = []
        bands = [(0, 0.02), (0.02, 0.1), (0.1, 0.3), (0.3, 0.7), (0.7, 0.9), (0.9, 0.98), (0.98, 1.0)]
        for lo, hi in bands:
            sel = (ug >= lo) & (ug <= hi)
            sub = np.where(sel[:, None], R, np.inf)
            a, b = np.unravel_index(np.argmin(sub), sub.shape)
            r = minimize(Rfun, [ug[a], vg[b]], method="Nelder-Mead", options={"xatol": 1e-11, "fatol": 1e-13, "maxiter": 4000})
            u, v = np.clip(r.x, [lo, 0], [hi, V_POLE_CUT])
            Rr = Rfun([u, v])
            if Rr > R[a, b]:
                u, v, Rr = ug[a], vg[b], float(R[a, b])
            jj = s.jet_pts([u], [v])
            Pp = jj["P"][0]
            tt, xx = thickness_at(Pp[0], Pp[2])
            refined.append({"u_band": [lo, hi], "R_mm": Rr, "u": float(u), "v": float(v), "xyz": Pp.tolist(),
                            "chord_fraction_sweep": float(xx), "attachment_thickness_mm": float(tt),
                            "required_mm": float(required_radius(tt))})
        i = np.unravel_index(np.argmin(margin), margin.shape)
        k = np.unravel_index(np.argmin(R), R.shape)
        cmap[side] = {
            "grid": f"{len(ug)} u x {len(vg)} v (u dense toward both ends, v up to 1 - 1e-4)",
            "smallest_radius_grid_mm": float(R[k]), "at_xyz": j["P"][k].tolist(),
            "tightest_per_chordwise_band": refined,
            "smallest_radius_refined_mm": min(d["R_mm"] for d in refined),
            "smallest_margin_R_over_required": float(margin[i]), "margin_at_xyz": j["P"][i].tolist(),
            "margin_attachment_thickness_mm": float(t[i]), "margin_radius_mm": float(R[i]),
            "points_failing_rule": int(np.sum(margin < 1.0)),
            "closest_grid_point_to_crown_mm": float(np.min(np.linalg.norm(j["P"][:, -1] - crown[0], axis=-1))),
            "crown_principal_radii_limit_mm": sorted((1 / np.abs(np.array([g["k1"][0, -1], g["k2"][0, -1]]))).tolist()),
        }
        cmap[f"_{side}"] = (j, g, R)
    # cross-check the numpy evaluator against OCCT's GeomLProp on random samples
    rng = np.random.default_rng(1)
    worst = 0.0
    for side in ("upper", "lower"):
        s = caps[side]
        for _ in range(300):
            u, v = rng.uniform(0, 1), rng.uniform(0, 0.999)
            pr = GeomLProp_SLProps(s.geom, u, v, 2, 1e-12)
            kk = max(abs(pr.MaxCurvature()), abs(pr.MinCurvature()))
            mine = geometry(s.jet_pts([u], [v]), s.sign)["kmax"][0]
            worst = max(worst, abs(kk - mine) / max(kk, 1e-9))
    # parent skins near the tip, for context
    pk = []
    for side in ("upper", "lower"):
        s = skins[side]
        j = s.jet(ug, np.linspace(200, 240, 41))
        pk.append(float(np.min(1 / geometry(j, s.sign)["kmax"])))
    out["curvature_map"] = {k: v for k, v in cmap.items() if not k.startswith("_")}
    out["curvature_map"]["evaluator_vs_occt_geomlprop_max_relative_difference"] = worst
    out["curvature_map"]["parent_skin_smallest_radius_mm_z200_240"] = min(pk)
    out["curvature_map"]["rule"] = ("required radius = 0.4 mm where the attachment thickness at the point's sweep-following chord "
                                    "fraction (X - LE(Z))/c(Z) is >= 0.8 mm, else half that thickness (full half-round)")
    allR = min(cmap["upper"]["smallest_radius_refined_mm"], cmap["lower"]["smallest_radius_refined_mm"])
    best = min(cmap["upper"]["tightest_per_chordwise_band"] + cmap["lower"]["tightest_per_chordwise_band"], key=lambda d: d["R_mm"])
    out["curvature_map"]["smallest_cap_radius_mm"] = allR

    # ---------------------------------------------------------- sections
    log("sections")
    out["sections"] = sections(skins, caps)

    # ---------------------------------------------------------- zebra
    log("zebra")
    out["zebra"] = zebra(joins)

    # ---------------------------------------------------------- pictures
    log("pictures")
    out["pictures"] = render(skins, caps)

    # ---------------------------------------------------------- summary
    sec_min = min(s["smallest_surface_radius_mm"] for s in out["sections"]["traces"])
    out["summary"] = {
        "closed_valid_single_solid": bool(out["reimport"]["solids"] == 1 and out["reimport"]["brepcheck_valid"] and out["reimport"]["free_edges"] == 0),
        "finished_height_mm": bb[5],
        "max_join_tangent_angle_deg": max(out["joins"][k]["max_tangent_plane_angle_deg"] for k in out["joins"] if k != "note"),
        "max_join_separation_mm": max(out["joins"][k]["max_separation_mm"] for k in out["joins"] if k != "note"),
        "max_join_relative_curvature_difference": max(out["joins"][k]["max_relative_curvature_difference"] for k in out["joins"] if k != "note"),
        "smallest_cap_principal_radius_mm": allR,
        "smallest_cap_radius_location": best,
        "smallest_rule_margin_R_over_required": min(cmap["upper"]["smallest_margin_R_over_required"], cmap["lower"]["smallest_margin_R_over_required"]),
        "printable_rule_met_everywhere": bool(cmap["upper"]["points_failing_rule"] == 0 and cmap["lower"]["points_failing_rule"] == 0 and allR >= RULE_MIN_R),
        "smallest_section_surface_radius_mm": sec_min,
        "zebra_max_phase_jump_periods": out["zebra"]["max_phase_jump_periods"],
        "zebra_max_stripe_direction_change_deg": out["zebra"]["max_stripe_direction_change_at_crossings_deg"],
        "zebra_max_gradient_relative_difference": out["zebra"]["max_gradient_relative_difference"],
    }
    OUT.write_text(json.dumps(rnd(out, 8), indent=1) + "\n")
    print(json.dumps(rnd(out["summary"], 6), indent=1))


# ====================================================================== sections
CROWN_EXCLUDE = 0.01   # mm; section samples closer than this to the collapsed crown pole are skipped


def plane_traces(face_surfs, origin, normal, zmin=235.0, step=0.01):
    """Exact samples of plane-face intersections with in-plane and normal curvature."""
    pln = gp_Pln(gp_Pnt(*map(float, origin)), gp_Dir(*map(float, normal)))
    npl = np.asarray(normal, float) / np.linalg.norm(normal)
    pieces = []
    for name, s in face_surfs.items():
        sec = BRepAlgoAPI_Section(s.face, pln, True)
        ex = TopExp_Explorer(sec.Shape(), TopAbs_EDGE)
        sa = ShapeAnalysis_Surface(s.geom)
        while ex.More():
            e = TopoDS.Edge_s(ex.Current())
            c = BRepAdaptor_Curve(e)
            t0, t1 = c.FirstParameter(), c.LastParameter()
            tt = np.linspace(t0, t1, 2001)
            P = np.array([[p.X(), p.Y(), p.Z()] for p in (c.Value(t) for t in tt)])
            inside = np.nonzero(P[:, 2] >= zmin - 0.05)[0]
            if len(inside) == 0:
                ex.Next()
                continue
            a, b = max(inside[0] - 1, 0), min(inside[-1] + 1, len(tt) - 1)
            length = np.sum(np.linalg.norm(np.diff(P[a:b + 1], axis=0), axis=1))
            n = max(int(length / step), 50)
            ts = np.linspace(tt[a], tt[b], n)
            X3 = np.array([[p.X(), p.Y(), p.Z()] for p in (c.Value(t) for t in ts)])
            coarse = np.unique(np.concatenate([np.arange(0, n, 25), [n - 1]]))
            uvc = []
            for k in coarse:
                q = sa.ValueOfUV(gp_Pnt(*map(float, X3[k])), 1e-10)
                uvc.append([q.X(), q.Y()])
            uvc = np.array(uvc)
            uv = np.stack([np.interp(np.arange(n), coarse, uvc[:, 0]), np.interp(np.arange(n), coarse, uvc[:, 1])], -1)
            for _ in range(5):   # closest-point Newton (Gauss-Newton) toward the sampled section points
                jj = s.jet_pts(uv[:, 0], uv[:, 1])
                r = X3 - jj["P"]
                E_ = (jj["Su"] ** 2).sum(-1); F_ = (jj["Su"] * jj["Sv"]).sum(-1); G_ = (jj["Sv"] ** 2).sum(-1)
                bu_ = (r * jj["Su"]).sum(-1); bv_ = (r * jj["Sv"]).sum(-1)
                det_ = E_ * G_ - F_ * F_
                uv[:, 0] += (G_ * bu_ - F_ * bv_) / det_
                uv[:, 1] += (E_ * bv_ - F_ * bu_) / det_
                uv[:, 0] = np.clip(uv[:, 0], s.u0, s.u1)
                uv[:, 1] = np.clip(uv[:, 1], s.v0, s.v1)
            # Newton-correct each sample onto the plane while staying on the surface
            for _ in range(4):
                jj = s.jet_pts(uv[:, 0], uv[:, 1])
                f = ((jj["P"] - origin) * npl).sum(-1)
                gu, gv = (jj["Su"] * npl).sum(-1), (jj["Sv"] * npl).sum(-1)
                g2 = gu * gu + gv * gv
                uv[:, 0] -= f * gu / g2
                uv[:, 1] -= f * gv / g2
                uv[:, 0] = np.clip(uv[:, 0], s.u0, s.u1)
                uv[:, 1] = np.clip(uv[:, 1], s.v0, s.v1 - 1e-9 if s.v1 <= 1 else s.v1)
            jj = s.jet_pts(uv[:, 0], uv[:, 1])
            g = geometry(jj, s.sign)
            ns = g["n"]
            t = np.cross(npl, ns)
            t /= np.linalg.norm(t, axis=-1, keepdims=True)
            Nc = np.cross(npl, t)
            Nc *= -np.sign((Nc * ns).sum(-1, keepdims=True))
            kn = normal_curvature(jj, g, t)
            graze = -(Nc * ns).sum(-1)
            k_in = kn / graze
            crown = np.array([cd.m_of(cd.TOP), 0.0, cd.TOP])
            keep = (jj["P"][:, 2] >= zmin) & (np.linalg.norm(jj["P"] - crown, axis=-1) > CROWN_EXCLUDE) & np.isfinite(g["kmax"])
            if keep.sum() > 2:
                pieces.append({"face": name, "P": jj["P"][keep], "k_plane": k_in[keep], "k_normal": kn[keep],
                               "kmax": g["kmax"][keep], "cos_plane_normal": graze[keep], "uv": uv[keep]})
            ex.Next()
    return pieces


def sections(skins, caps):
    faces = {"upper_skin": skins["upper"], "lower_skin": skins["lower"], "upper_cap": caps["upper"], "lower_cap": caps["lower"]}
    traces = []
    xis = [0.0005, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
           0.95, 0.98, 0.99, 0.995, 0.997, 0.998, 0.9985, 0.999, 0.9995]
    cut_defs = []
    for xi in xis:
        x0 = 65 + 85 * xi
        cut_defs.append(("chordwise X = const", xi, np.array([x0, 0, 240.0]), np.array([1.0, 0, 0])))
    for xi in xis:
        x0 = 65 + 85 * xi
        th = np.arccos(1 - 2 * xi)
        k = cd.DM + cd.DA * (-np.cos(th))          # dX/dZ of the parent's constant-chord-fraction ruling
        cut_defs.append(("sweep-following (parent ruling)", xi, np.array([x0, 0, 240.0]), np.array([1.0, 0, -k])))
    for hq in [0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.85, 0.95, 0.99]:
        for branch in (-1, +1):
            jt = par.tip_outline_jet(np.array([hq]), branch)
            p0, tg = jt["point"][0], jt["d_p"][0]
            cut_defs.append((f"normal to perimeter ({'LE' if branch < 0 else 'TE'} branch, h={hq})", None, p0, tg))
    for kind, xi, origin, normal in cut_defs:
        pieces = plane_traces(faces, origin, normal)
        cap_pieces = [p for p in pieces if "cap" in p["face"]]
        if not cap_pieces:
            continue
        allk = np.concatenate([p["k_plane"] for p in pieces])
        allP = np.concatenate([p["P"] for p in pieces])
        allkn = np.concatenate([p["k_normal"] for p in pieces])
        allkm = np.concatenate([p["kmax"] for p in pieces])
        allc = np.concatenate([p["cos_plane_normal"] for p in pieces])
        capk = np.concatenate([p["k_plane"] for p in cap_pieces])
        capP = np.concatenate([p["P"] for p in cap_pieces])
        capkn = np.concatenate([p["k_normal"] for p in cap_pieces])
        capkm = np.concatenate([p["kmax"] for p in cap_pieces])
        capc = np.concatenate([p["cos_plane_normal"] for p in cap_pieces])
        # reversals: sign changes of in-plane curvature along each piece (ignoring |k| < 1e-4 /mm)
        rev = []
        for p in pieces:
            k = p["k_plane"]
            sig = np.sign(np.where(np.abs(k) < 1e-4, 0, k))
            nz = np.nonzero(sig)[0]
            for a, b in zip(nz[:-1], nz[1:]):
                if sig[a] != sig[b]:
                    rev.append(p["P"][b].tolist())
        lip = [p for p in cap_pieces]
        lip_pts = []
        for p in lip:
            for end in (0, -1):
                if abs(p["P"][end, 1]) < 1e-5:
                    lip_pts.append((p["P"][end], 1 / p["k_plane"][end], 1 / p["k_normal"][end], p["cos_plane_normal"][end]))
        i = int(np.argmax(np.abs(capk)))
        im = int(np.argmax(capkm))
        if xi is not None:
            t = float(par.section_thickness(xi, 240.0))
        else:
            t = float(thickness_at(origin[0], origin[2])[0])
        req = float(required_radius(t))
        traces.append({
            "kind": kind, "chord_fraction": xi, "plane_origin": origin.tolist(), "plane_normal": (normal / np.linalg.norm(normal)).tolist(),
            "attachment_thickness_mm": t, "required_radius_mm": req,
            "points_on_cap": int(len(capk)), "points_total": int(len(allk)),
            "smallest_section_radius_mm": float(1 / np.abs(capk[i])), "at_xyz": capP[i].tolist(),
            "cos_plane_to_surface_normal_there": float(capc[i]),
            "smallest_surface_radius_mm": float(1 / capkm[im]), "surface_radius_at_xyz": capP[im].tolist(),
            "smallest_normal_curvature_radius_along_trace_mm": float(1 / np.max(np.abs(capkn))),
            "lip": [{"xyz": a.tolist(), "section_radius_mm": float(b), "normal_radius_mm": float(c), "cos_plane_to_surface_normal": float(d)} for a, b, c, d in lip_pts],
            "curvature_reversals": len(rev), "reversal_xyz": rev[:6],
            "meets_rule_section_radius": bool(1 / np.abs(capk[i]) >= req),
            "meets_rule_surface_radius": bool(1 / capkm[im] >= req),
        })
    worst_sec = min(traces, key=lambda d: d["smallest_section_radius_mm"] / d["required_radius_mm"])
    worst_surf = min(traces, key=lambda d: d["smallest_surface_radius_mm"] / d["required_radius_mm"])
    return {"traces": traces,
            "count": len(traces),
            "total_reversals": int(sum(t["curvature_reversals"] for t in traces)),
            "worst_section_radius_vs_rule": {k: worst_sec[k] for k in ("kind", "chord_fraction", "smallest_section_radius_mm", "at_xyz", "cos_plane_to_surface_normal_there", "attachment_thickness_mm", "required_radius_mm")},
            "worst_surface_radius_vs_rule": {k: worst_surf[k] for k in ("kind", "chord_fraction", "smallest_surface_radius_mm", "surface_radius_at_xyz", "attachment_thickness_mm", "required_radius_mm")},
            "note": ("section radius = curvature of the planar cut curve; where the plane is far from containing the surface normal "
                     "(cos_plane_to_surface_normal small: the plane grazes the surface) it exceeds the surface's own curvature by 1/cos. "
                     "Surface radius = 1/(larger principal curvature) at trace points. Traces start at Z = 235 on the parent skin; "
                     "samples within 0.01 mm of the collapsed crown pole are skipped."),
            "traces_failing_rule_on_section_radius": [(t["kind"], t["chord_fraction"], t["smallest_section_radius_mm"], t["cos_plane_to_surface_normal_there"], t["attachment_thickness_mm"]) for t in traces if not t["meets_rule_section_radius"]],
            "traces_failing_rule_on_surface_radius": [(t["kind"], t["chord_fraction"]) for t in traces if not t["meets_rule_surface_radius"]]}


# ====================================================================== zebra
def view_set():
    views = []
    for side in (+1, -1):
        for az in (-60, -30, 0, 30, 60):
            for el in (-30, 0, 30, 60):
                a, e = np.radians(az), np.radians(el)
                # camera on the +/-Y side looking toward the fin
                d = np.array([np.sin(a) * np.cos(e), side * np.cos(a) * np.cos(e), np.sin(e)])
                views.append((f"{'+' if side > 0 else '-'}Y az{az} el{el}", -d))
    return views


def stripe_axes(V):
    up = np.array([0, 0, 1.0]) if abs(V[2]) < 0.9 else np.array([1.0, 0, 0])
    e1 = np.cross(V, up)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(V, e1)
    return [(ang, np.cos(np.radians(ang)) * e1 + np.sin(np.radians(ang)) * e2) for ang in range(0, 180, 30)]


PERIOD = 0.1  # stripe period in units of (reflected ray . stripe axis)


def phase_and_gradient(j, g, V, A):
    n = g["n"]
    R = V - 2 * (V * n).sum(-1, keepdims=True) * n
    p = (R * A).sum(-1) / PERIOD
    dRA = lambda dn: -2 * ((V * dn).sum(-1) * (n * A).sum(-1) + (V * n).sum(-1) * (dn * A).sum(-1)) / PERIOD
    pu, pv = dRA(g["nu"]), dRA(g["nv"])
    E, F, G = g["E"], g["F"], g["G"]
    det = E * G - F * F
    a = (G * pu - F * pv) / det
    b = (E * pv - F * pu) / det
    grad = a[..., None] * j["Su"] + b[..., None] * j["Sv"]
    return p, grad, (V * n).sum(-1)


def zebra(joins):
    res = {"views": 0, "stripe_axes_per_view": 6, "period": PERIOD, "joins": {}}
    worst_phase, worst_dir, worst_grad = (0, None), (0, None), (0, None)
    for jname in ("upper", "lower", "perimeter_leading_edge_branch", "perimeter_trailing_edge_branch"):
        m = joins[f"_{jname}"]
        jr = {"crossings": 0, "max_phase_jump_periods": 0.0, "max_gradient_relative_difference": 0.0, "gradient_difference_where": None,
              "max_stripe_direction_change_at_crossings_deg": 0.0}
        for vname, V in view_set():
            for ang, A in stripe_axes(V):
                pa, ga_, va = phase_and_gradient(m["ja"], m["ga"], V, A)
                pb, gb_, vb = phase_and_gradient(m["jb"], m["gb"], V, A)
                vis = (va < -1e-3) & (vb < -1e-3)
                if not vis.any():
                    continue
                dp = np.abs(pa - pb)[vis]
                k = int(np.argmax(dp))
                if dp[k] > jr["max_phase_jump_periods"]:
                    jr["max_phase_jump_periods"] = float(dp[k])
                if dp[k] > worst_phase[0]:
                    worst_phase = (float(dp[k]), {"join": jname, "view": vname, "stripe_axis_deg": ang, "xyz": m["P"][vis][k].tolist()})
                gn = np.maximum(np.maximum(np.linalg.norm(ga_, axis=-1), np.linalg.norm(gb_, axis=-1)), 1e-3)
                gr = (np.linalg.norm(ga_ - gb_, axis=-1) / gn)[vis]
                k = int(np.argmax(gr))
                if gr[k] > jr["max_gradient_relative_difference"]:
                    jr["max_gradient_relative_difference"] = float(gr[k])
                    jr["gradient_difference_where"] = {"view": vname, "stripe_axis_deg": ang, "xyz": m["P"][vis][k].tolist()}
                if gr[k] > worst_grad[0]:
                    worst_grad = (float(gr[k]), {"join": jname, "view": vname, "stripe_axis_deg": ang, "xyz": m["P"][vis][k].tolist(),
                                                 "gradient_per_mm": float(gn[vis][k])})
                # light/dark crossings along the join (boundaries at half-integer multiples of the period)
                s = np.sin(2 * np.pi * pa)
                idx = np.nonzero((np.sign(s[:-1]) != np.sign(s[1:])) & vis[:-1] & vis[1:])[0]
                for i in idx:
                    w = s[i] / (s[i] - s[i + 1])
                    da = (1 - w) * ga_[i] + w * ga_[i + 1]
                    db = (1 - w) * gb_[i] + w * gb_[i + 1]
                    na = (1 - w) * m["ga"]["n"][i] + w * m["ga"]["n"][i + 1]
                    nb = (1 - w) * m["gb"]["n"][i] + w * m["gb"]["n"][i + 1]
                    sa_, sb_ = np.cross(na, da), np.cross(nb, db)
                    if np.linalg.norm(sa_) < 1e-9 or np.linalg.norm(sb_) < 1e-9:
                        continue
                    c = abs(np.dot(sa_, sb_)) / (np.linalg.norm(sa_) * np.linalg.norm(sb_))
                    angd = float(np.degrees(np.arccos(min(1.0, c))))
                    jr["crossings"] += 1
                    jr["max_stripe_direction_change_at_crossings_deg"] = max(jr["max_stripe_direction_change_at_crossings_deg"], angd)
                    if angd > worst_dir[0]:
                        worst_dir = (angd, {"join": jname, "view": vname, "stripe_axis_deg": ang,
                                            "xyz": ((1 - w) * m["P"][i] + w * m["P"][i + 1]).tolist()})
        res["joins"][jname] = jr
    res["views"] = len(view_set())
    res["max_phase_jump_periods"] = worst_phase[0]
    res["max_phase_jump_where"] = worst_phase[1]
    res["max_gradient_relative_difference"] = worst_grad[0]
    res["max_gradient_difference_where"] = worst_grad[1]
    res["max_stripe_direction_change_at_crossings_deg"] = worst_dir[0]
    res["max_stripe_direction_change_where"] = worst_dir[1]
    res["note"] = ("phase p = (R . A)/period with R the view ray reflected about the exact CAD normal and A a stripe axis normal to the "
                   "view; 40 views (+/-Y sides, azimuth -60..60, elevation -30..60) x 6 stripe axes; only points visible from the view "
                   "on both faces count. Gradient from the exact shape operator. The perimeter join is sampled up to v = 1 - 1e-6 "
                   "from the collapsed crown, where surface derivatives are undefined. Gradient difference is relative to the larger of the two "
                   "gradient magnitudes (floored at 1e-3 stripe periods per mm).")
    return res


# ====================================================================== pictures
def render(skins, caps):
    """Zebra pictures from exact CAD normals: dense parameter grids point-splatted with a z-buffer.

    No tessellation and no edge lines; camera, zoom and stripe lighting are fixed in code so
    revisions can be compared like for like.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    period = 0.15
    oblique = -np.array([0.35, 0.75, 0.56]) / np.linalg.norm([0.35, 0.75, 0.56])
    above = -np.array([0.0, 0.45, 0.89]) / np.linalg.norm([0.0, 0.45, 0.89])
    shots = [
        # file, view ray, window (screen units, None = fit), stripe axis angle, lowest Z drawn, width px
        ("zebra_fin_full.png", np.array([0.0, -1.0, 0.0]), None, 75.0, 0.0, 900),
        ("zebra_tip_side.png", np.array([0.0, -1.0, 0.0]), None, 75.0, 228.0, 1600),
        ("zebra_tip_oblique.png", oblique, None, 90.0, 230.0, 1600),
        ("zebra_tip_oblique_axis30.png", oblique, None, 30.0, 230.0, 1600),
        ("zebra_tip_above.png", above, None, 90.0, 236.0, 1600),
    ]
    pics = []
    for fname, V, window, sang, zmin, px in shots:
        up = np.array([0, 0, 1.0])
        e1 = np.cross(up, V)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(V, e1)
        A = np.cos(np.radians(sang)) * e1 + np.sin(np.radians(sang)) * e2
        full = zmin < 1
        pts, cols = [], []
        grids = [(skins["upper"], (zmin, 240.0), 3000, 3000 if full else 700),
                 (skins["lower"], (zmin, 240.0), 3000, 3000 if full else 700),
                 (caps["upper"], (0.0, 1.0), 3000, 700), (caps["lower"], (0.0, 1.0), 3000, 700)]
        for s, (va, vb), nu, nv in grids:
            uu = 0.5 - 0.5 * np.cos(np.linspace(0, np.pi, nu))
            vv = np.linspace(va, vb, nv) if vb > 1 else 1 - (1 - np.linspace(0, 1, nv)) ** 2
            for chunk in np.array_split(np.arange(nu), 6):
                j = s.jet(uu[chunk], vv)
                g = geometry(j, s.sign)
                n = g["n"]
                R = V - 2 * (V * n).sum(-1, keepdims=True) * n
                p = (R * A).sum(-1) / period
                vn = (V * n).sum(-1)
                light = (np.floor(p) % 2 == 0)
                c = np.where(light, 0.95, 0.10)          # flat zebra, no shading
                front = vn < 0
                pts.append(j["P"][front])
                cols.append(c[front])
        P = np.concatenate(pts)
        C = np.concatenate(cols)
        x, y, d = P @ e1, P @ e2, P @ V
        if window is None:
            window = (x.min() - 1, x.max() + 1, y.min() - 1, y.max() + 1)
        x0, x1, y0, y1 = window
        W = px
        Hh = max(int(px * (y1 - y0) / (x1 - x0)), 50)
        img = np.full((Hh, W), 1.0)
        ix = ((x - x0) / (x1 - x0) * (W - 1)).astype(int)
        iy = ((y1 - y) / (y1 - y0) * (Hh - 1)).astype(int)
        order = np.argsort(-d)                       # far first, near last (near overwrites)
        ix, iy, C = ix[order], iy[order], C[order]
        for dx in (0, 1):                            # 2 x 2 splat closes sampling gaps
            for dy in (0, 1):
                ok = (ix + dx < W) & (iy + dy < Hh) & (ix + dx >= 0) & (iy + dy >= 0)
                img[iy[ok] + dy, ix[ok] + dx] = C[ok]
        plt.imsave(HERE / fname, img, cmap="gray", vmin=0, vmax=1)
        pics.append({"file": fname, "view_ray": [round(float(q), 4) for q in V], "stripe_axis_deg_in_view": sang,
                     "stripe_period": period, "lowest_z_drawn_mm": zmin, "pixels": [W, Hh]})
    return {"pictures": pics, "note": "orthographic, exact CAD normals on dense parameter grids (no tessellation), "
                                      "no edge lines, back faces culled; camera and stripes fixed in verify_cap.py"}


if __name__ == "__main__":
    main()
