"""Verify wing.step (reimported) against the skill's checks; writes readback.json.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/fin40/sharp/verify_cap.py
"""
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import numpy as np
import contourpy
from scipy.optimize import minimize
from scipy.spatial import cKDTree

from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_SOLID, TopAbs_SHELL, TopAbs_WIRE, TopAbs_EDGE, TopAbs_IN, TopAbs_OUT
from OCP.TopoDS import TopoDS
from OCP.BRep import BRep_Tool
from OCP.GeomAdaptor import GeomAdaptor_Surface
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Plane
from OCP.GeomLProp import GeomLProp_SLProps
from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds
from OCP.TopTools import TopTools_HSequenceOfShape
from OCP.gp import gp_Pnt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
PARENT_DIR = ROOT / "runs/cap-radius/parents/fin40/sharp"
sys.path.insert(0, str(PARENT_DIR))
sys.path.insert(0, str(HERE))
import parent  # noqa: E402
import build_cap as design  # noqa: E402  (guide stations and perimeter map only)
from surfkit import (BSurf, read_step, faces_of, sha256, unit, differential, principal,  # noqa: E402
                     second_form, shape_operator, zebra_phase, zebra_gradient)

STEP = HERE / "wing.step"
ZA = parent.ATTACHMENT
R_MIN = 0.4                  # printable minimum radius (mm)
U_EXCL = 2e-4                # grid stays this far (in u) from the collapsed LE/TE ends of the cap


def att_thickness(u):
    """Parent attachment thickness at the chord fraction of cap chord station u."""
    return parent.section_thickness(design.xi_of(np.asarray(u, float)), ZA)


def required_radius(t):
    return np.where(t >= 0.8, R_MIN, 0.5 * t)


def r3(p):
    return [round(float(x), 4) for x in p]


# ---------------------------------------------------------------------------
def topology(shape):
    n = lambda kind: sum(1 for _ in iter_shapes(shape, kind))
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props, 1e-9)
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape, box, False, False)
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
    fb = ShapeAnalysis_FreeBounds(shape)
    free_closed = sum(1 for _ in iter_shapes(fb.GetClosedWires(), TopAbs_WIRE)) + \
        sum(1 for _ in iter_shapes(fb.GetOpenWires(), TopAbs_WIRE))
    shell = TopoDS.Shell_s(next(iter_shapes(shape, TopAbs_SHELL)))
    return {
        "solids": n(TopAbs_SOLID), "shells": n(TopAbs_SHELL),
        "faces": len(faces_of(shape)),
        "face_types": [str(GeomAdaptor_Surface(BRep_Tool.Surface_s(f)).GetType()).split(".")[-1] for f in faces_of(shape)],
        "brepcheck_valid": bool(BRepCheck_Analyzer(shape).IsValid()),
        "shell_closed": bool(BRep_Tool.IsClosed_s(shell)),
        "free_boundary_wires": free_closed,
        "volume_mm3": props.Mass(),
        "bounding_box_mm": [xmin, xmax, ymin, ymax, zmin, zmax],
    }


def iter_shapes(shape, kind):
    ex = TopExp_Explorer(shape, kind)
    while ex.More():
        yield ex.Current()
        ex.Next()


def classify(shape):
    skins, cap, planes = {}, None, []
    for f in faces_of(shape):
        ad = GeomAdaptor_Surface(BRep_Tool.Surface_s(f))
        if ad.GetType() == GeomAbs_BSplineSurface:
            s = BSurf(f)
            if abs(s.vrange[1] - ZA) < 1e-9 and s.nv == 5:
                skins[+1 if s.at(0.5, 120.0)[1] > 0 else -1] = s
            else:
                cap = s
        elif ad.GetType() == GeomAbs_Plane:
            planes.append(f)
    return skins, cap, planes


# ---------------------------------------------------------------------------
def evaluator_crosscheck(surfs, solid, rng):
    """NumPy evaluator vs OCC GeomLProp; outward normal vs solid classifier."""
    worst = {"point_mm": 0.0, "normal_deg": 0.0, "k1_rel": 0.0, "outward_probe_failures": 0}
    clf = BRepClass3d_SolidClassifier(solid)
    for s in surfs:
        for _ in range(60):
            u = rng.uniform(*s.urange) * 0.98 + 0.01 * (s.urange[1] - s.urange[0])
            v = s.vrange[0] + (s.vrange[1] - s.vrange[0]) * rng.uniform(0.02, 0.98)
            d = differential(s.jet(u, v), s.sign)
            k1, k2, _, _ = principal(d)
            pr = GeomLProp_SLProps(s.geom, u, v, 2, 1e-9)
            p = pr.Value()
            nn = pr.Normal()
            no = s.sign * np.array([nn.X(), nn.Y(), nn.Z()])
            worst["point_mm"] = max(worst["point_mm"], float(np.linalg.norm(d["P"] - [p.X(), p.Y(), p.Z()])))
            worst["normal_deg"] = max(worst["normal_deg"], float(np.degrees(np.arccos(np.clip(no @ d["n"], -1, 1)))))
            occ = max(abs(pr.MaxCurvature()), abs(pr.MinCurvature()))
            worst["k1_rel"] = max(worst["k1_rel"], float(abs(max(abs(k1), abs(k2)) - occ) / max(occ, 1e-9)))
            q = d["P"] + 0.02 * d["n"]
            clf.Perform(gp_Pnt(*q), 1e-9)
            if clf.State() != TopAbs_OUT:
                worst["outward_probe_failures"] += 1
            q = d["P"] - 0.02 * d["n"]
            clf.Perform(gp_Pnt(*q), 1e-9)
            if clf.State() != TopAbs_IN:
                worst["outward_probe_failures"] += 1
    return worst


def parent_preserved(skins):
    orig = read_step(PARENT_DIR / "parent.step")
    oskins = {}
    for f in faces_of(orig):
        if GeomAdaptor_Surface(BRep_Tool.Surface_s(f)).GetType() == GeomAbs_BSplineSurface:
            s = BSurf(f)
            oskins[+1 if s.at(0.5, 120.0)[1] > 0 else -1] = s
    out = {}
    for side in (+1, -1):
        a, b = skins[side], oskins[side]
        pole_diff = float(np.abs(a.C - b.C).max()) if a.C.shape == b.C.shape else None
        knots_same = bool(a.C.shape == b.C.shape and np.array_equal(a.tu, b.tu) and np.array_equal(a.tv, b.tv))
        s = np.linspace(0, 1, 201)
        z = np.linspace(0, ZA, 61)
        ref = parent.skin(s[:, None], z[None, :], side)
        got = a.grid(s, z)
        out[parent.SIDES[side]] = {"max_pole_difference_vs_parent_step_mm": pole_diff, "same_knots_and_degrees": knots_same,
                                   "max_deviation_vs_parent_py_mm": float(np.linalg.norm(got - ref, axis=-1).max()),
                                   "samples": int(ref.size // 3)}
    return out


# ---------------------------------------------------------------------------
def outline_distance(pts, tree, w):
    """Distance from points to the locked outline: dense KD-tree seed, then exact 1-D refinement."""
    from scipy.optimize import minimize_scalar
    dist, idx = tree.query(pts)
    for k in range(len(pts)):
        w0 = w[idx[k]]

        def f(t):
            x, z = design.perimeter(np.clip(t, -1, 1))
            return (x - pts[k, 0]) ** 2 + pts[k, 1] ** 2 + (z - pts[k, 2]) ** 2
        r = minimize_scalar(f, bounds=(max(-1, w0 - 2e-5), min(1, w0 + 2e-5)), method="bounded",
                            options={"xatol": 1e-14})
        dist[k] = min(dist[k], np.sqrt(r.fun))
    return dist


def locked_dimensions(cap, skins, topo):
    u = np.linspace(0, 1, 4001)
    lip = cap.grid(u, np.array([0.0]))[:, 0]
    # distance of every lip point to the outline (dense outline + local refinement)
    w = np.linspace(-1, 1, 200001)
    ox, oz = design.perimeter(w)
    tree = cKDTree(np.c_[ox, np.zeros_like(ox), oz])
    dist = outline_distance(lip, tree, w)
    # outline samples projected onto the cap face (OCC)
    hs = np.linspace(0, 1, 101)
    proj_max = 0.0
    for br in (-1, 1):
        pts = parent.tip_outline(hs, br)
        for p in pts:
            pr = GeomAPI_ProjectPointOnSurf(gp_Pnt(*p), cap.geom)
            proj_max = max(proj_max, pr.LowerDistance())
    # side-view containment: every cap point lies inside the locked outline
    uu = np.linspace(0, 1, 1201)
    vv = np.linspace(-1, 1, 801)
    S = cap.grid(uu, vv)
    X, Y, Z = S[..., 0], S[..., 1], S[..., 2]
    above_crown = float((Z - (ZA + parent.REACH)).max())
    hh = np.clip((Z - ZA) / parent.REACH, 0, 1)
    half = design.a_of(Z) * np.sqrt(1 - hh**3)
    excess = np.maximum(np.maximum(design.m_of(Z) - half - X, X - design.m_of(Z) - half), 0)
    excess = np.maximum(excess, np.maximum(Z - (ZA + parent.REACH), 0))
    flag = excess > 0
    outside = 0.0
    if flag.any():
        pts = np.c_[X[flag], np.zeros(flag.sum()), Z[flag]]
        outside = float(outline_distance(pts, tree, w).max())
    crown_i = np.unravel_index(np.argmax(Z), Z.shape)
    zmax = topo["bounding_box_mm"][5]
    return {
        "finished_height_mm": zmax,
        "finished_height_target_mm": ZA + parent.REACH,
        "crown_point_mm": r3(S[crown_i]),
        "crown_target_mm": r3(parent.tip_outline(1.0, 1)),
        "lip_max_abs_Y_mm": float(np.abs(lip[:, 1]).max()),
        "lip_to_outline_max_mm": float(dist.max()),
        "outline_to_cap_max_mm": float(proj_max),
        "cap_min_Z_mm": float(Z.min()),
        "cap_points_above_finished_height_mm": max(0.0, above_crown),
        "cap_points_outside_outline_in_side_view": int(flag.sum()),
        "cap_max_distance_outside_outline_mm": outside,
        "containment_grid": [len(uu), len(vv)],
        "attachment_station_mm": ZA,
        "root_chord_mm": float(skins[+1].at(1.0, 0.0)[0] - skins[+1].at(0.0, 0.0)[0]),
        "attachment_chord_mm": float(skins[+1].at(1.0, ZA)[0] - skins[+1].at(0.0, ZA)[0]),
    }


# ---------------------------------------------------------------------------
def joins(cap, skins):
    u = np.concatenate([np.linspace(1e-3, 0.02, 200), np.linspace(0.02, 0.98, 1601)[1:-1], np.linspace(0.98, 1 - 1e-3, 200)])
    res = {}
    for side, vcap in ((+1, 1.0), (-1, -1.0)):
        dc = differential(cap.jet(u, np.full_like(u, vcap)), cap.sign)
        dp = differential(skins[side].jet(u, np.full_like(u, ZA)), skins[side].sign)
        sep = np.linalg.norm(dc["P"] - dp["P"], axis=-1)
        ang = np.degrees(np.arccos(np.clip(np.sum(dc["n"] * dp["n"], -1), -1, 1)))
        e1 = unit(dp["Su"])
        n = dp["n"]
        e2 = np.cross(n, e1)                      # across the join
        e2 *= np.sign(np.sum(e2 * [0, 0, 1], -1))[:, None]   # pointing up (into the cap)
        Wc = shape_operator(dc, e1, e2)
        Wp = shape_operator(dp, e1, e2)
        dW = np.abs(Wc - Wp)
        kscale = np.abs(Wp).max(axis=(1, 2))
        i_ang = int(np.argmax(ang))
        i_w = int(np.argmax(dW.max(axis=(1, 2))))
        res[parent.SIDES[side]] = {
            "samples": len(u), "u_range": [float(u[0]), float(u[-1])],
            "max_separation_mm": float(sep.max()),
            "max_tangent_plane_angle_deg": float(ang.max()), "at_u": float(u[i_ang]), "at_mm": r3(dp["P"][i_ang]),
            "max_curvature_difference_per_mm": {
                "along_join_k11": float(dW[:, 0, 0].max()), "twist_k12": float(dW[:, 0, 1].max()),
                "across_join_k22": float(dW[:, 1, 1].max()),
                "worst_at_u": float(u[i_w]), "worst_at_mm": r3(dp["P"][i_w]),
                "parent_curvature_at_worst_per_mm": float(kscale[i_w]),
            },
            "max_relative_curvature_difference": float((dW.max(axis=(1, 2)) / np.maximum(kscale, 1e-9)).max()),
            "endpoint_note": "u in [0.001, 0.999]; the LE/TE ends of the join are the collapsed corners of the cap face",
            "endpoint_distance_mm": [float(np.linalg.norm(dp["P"][0] - [65, 0, ZA])), float(np.linalg.norm(dp["P"][-1] - [150, 0, ZA]))],
        }
    return res


# ---------------------------------------------------------------------------
def parent_k2_near_attachment(skins):
    """The parent skin is ruled along its spanwise lines (Gaussian curvature <= 0), so its
    smaller principal curvature is slightly negative; report it for comparison with the cap."""
    s = np.linspace(1e-3, 1 - 1e-3, 2001)
    out = {}
    for side in (+1, -1):
        d = differential(skins[side].jet(s, np.array([200.0, 220.0, ZA]), grid=True), skins[side].sign)
        k1, k2, K, H = principal(d)
        out[parent.SIDES[side]] = {"min_k2_per_mm_Z200_to_240": float(k2.min()),
                                   "at_mm": r3(d["P"][np.unravel_index(np.argmin(k2), k2.shape)])}
    return out


def curvature_map(cap):
    nu, nv = 1601, 1201
    s = np.linspace(0, 1, nu)
    us = U_EXCL + (1 - 2 * U_EXCL) * 0.5 * (s + 0.5 * (1 - np.cos(np.pi * s)))
    vs = np.linspace(-1, 1, nv)
    d = differential(cap.jet(us, vs, grid=True), cap.sign)
    k1, k2, K, H = principal(d)
    # fold check: parametric normal must agree with the outward reference of the thickness law
    P = d["P"]
    m = design.m_of(P[..., 2])
    ref = np.stack([2 * (P[..., 0] - m), 2 * P[..., 1] / design.K**2,
                    np.ones_like(m)], -1)        # Z component sign is irrelevant away from the crown
    # (use the X/Y components only; they never vanish together off the crown line)
    agree = np.sum(d["n"][..., :2] * ref[..., :2], -1)
    fold_bad = int(np.sum(agree < -1e-9))
    # local maxima of k1
    pad = np.pad(k1, 1, constant_values=-np.inf)
    is_max = np.ones_like(k1, bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                is_max &= k1 >= pad[1 + di:1 + di + nu, 1 + dj:1 + dj + nv]
    cand = np.argwhere(is_max)
    cand = cand[np.argsort(-k1[tuple(cand.T)])][:12]
    refined = []
    for i, j in cand:
        def f(x):
            dd = differential(cap.jet(np.array([x[0]]), np.array([x[1]])), cap.sign)
            return -principal(dd)[0][0]
        lo = max(us[0], us[max(i - 3, 0)]), min(us[-1], us[min(i + 3, nu - 1)])
        bv = max(-1, vs[max(j - 3, 0)]), min(1, vs[min(j + 3, nv - 1)])
        r = minimize(f, [us[i], vs[j]], method="L-BFGS-B", bounds=[lo, bv], options={"ftol": 1e-14, "gtol": 1e-12})
        kk = -r.fun
        uu, vv = r.x
        t = float(att_thickness(uu))
        req = float(required_radius(t))
        refined.append({"u": float(uu), "v": float(vv), "point_mm": r3(cap.at(uu, vv)), "radius_mm": float(1 / kk),
                        "chord_fraction": float(design.xi_of(uu)), "attachment_thickness_mm": t,
                        "required_radius_mm": req, "pass": bool(1 / kk >= req)})
    refined.sort(key=lambda e: e["radius_mm"])
    # rule check on every grid point
    t = att_thickness(us)[:, None] * np.ones_like(k1)
    req = required_radius(t)
    R = 1 / np.maximum(k1, 1e-12)
    viol = int(np.sum(R < req))
    thick = t >= 0.8
    i_t = np.unravel_index(np.argmin(np.where(thick, R, np.inf)), R.shape)
    i_n = np.unravel_index(np.argmin(np.where(~thick, R, np.inf)), R.shape) if (~thick).any() else None
    end_dist = [float(np.linalg.norm(cap.at(us[0], vs) - [65, 0, ZA], axis=-1).max()),
                float(np.linalg.norm(cap.at(us[-1], vs) - [150, 0, ZA], axis=-1).max())]
    return {
        "grid": [nu, nv], "u_range": [float(us[0]), float(us[-1])], "v_range": [-1.0, 1.0],
        "grid_reaches_within_mm_of_collapsed_ends": {"LE_point": end_dist[0], "TE_point": end_dist[1]},
        "smallest_radius_mm": float(R.min()),
        "smallest_radius_at_mm": r3(P[np.unravel_index(np.argmin(R), R.shape)]),
        "smallest_radius_refined": refined[0],
        "local_maxima_refined": refined,
        "thick_zone_min_radius_mm": float(R[i_t]), "thick_zone_min_at_mm": r3(P[i_t]),
        "thin_zone_min_radius_mm": None if i_n is None else float(R[i_n]),
        "thin_zone_min_at_mm": None if i_n is None else r3(P[i_n]),
        "thin_zone_attachment_thickness_mm": None if i_n is None else float(t[i_n]),
        "grid_points_violating_rule": viol,
        "min_k2_per_mm": float(k2.min()), "min_k2_at_uv": [float(us[np.unravel_index(np.argmin(k2), k2.shape)[0]]), float(vs[np.unravel_index(np.argmin(k2), k2.shape)[1]])],
        "min_k2_at_mm": r3(P[np.unravel_index(np.argmin(k2), k2.shape)]),
        "concave_points_k2_below_-1e-6": int(np.sum(k2 < -1e-6)),
        "concave_points_k2_below_-1e-4": int(np.sum(k2 < -1e-4)),
        "min_k2_away_from_corners_per_mm": float(k2[(us > 0.01) & (us < 0.99)].min()),
        "min_k2_away_from_corners_at_mm": r3(P[(us > 0.01) & (us < 0.99)][np.unravel_index(np.argmin(k2[(us > 0.01) & (us < 0.99)]), k2[(us > 0.01) & (us < 0.99)].shape)]),
        "min_parametric_area_element": float(d["area"].min()),
        "fold_points": fold_bad,
    }, (us, vs, d, k1, k2)


def design_law_check(cap):
    """Fit-quality recheck with CAD normals: angle between the STEP cap normal and the
    normal of the analytic thickness law the cap was fitted to, plus the smoothness of the
    zebra phase along dense lines (used to rule out stair-steps seen in the rendered images)."""
    us = np.linspace(0.005, 0.995, 1500)
    vs = np.linspace(-0.995, 0.995, 900)
    d = differential(cap.jet(us, vs, grid=True), cap.sign)
    X, Y, Z = d["P"][..., 0], d["P"][..., 1], d["P"][..., 2]
    h = (Z - ZA) / parent.REACH
    m, a = design.m_of(Z), design.a_of(Z)
    dm, da = 27.5 / ZA, -37.5 / ZA
    dq = 2 * a * da * (1 - h**3) + a * a * (-3 * h * h / parent.REACH)
    g = unit(np.stack([2 * (X - m), 2 * Y / design.K**2, -2 * (X - m) * dm - dq], -1))
    ang = np.degrees(np.arccos(np.clip(np.sum(g * d["n"], -1), -1, 1)))
    i = np.unravel_index(np.argmax(ang), ang.shape)
    V = unit(np.array([1, -1.2, -0.35]))
    up = unit(np.array([0, 0, 1.0]) - V[2] * V)
    dvec = np.cos(np.radians(20)) * up + np.sin(np.radians(20)) * np.cross(V, up)
    worst2 = 0.0
    for v in (0.3, 0.6, 0.85):
        uu = np.linspace(0.2, 0.8, 20001)
        dd = differential(cap.jet(uu, np.full_like(uu, v)), cap.sign)
        ph, _ = zebra_phase(dd["n"], V, dvec, np.radians(1.5))
        worst2 = max(worst2, float(np.abs(np.diff(ph, 2)).max()))
    return {"max_normal_angle_to_design_law_deg": float(ang.max()), "at_mm": r3(d["P"][i]),
            "p99_normal_angle_deg": float(np.percentile(ang, 99)),
            "zebra_phase_max_second_difference_on_lines": worst2,
            "note": "u in [0.005, 0.995]; phase lines at v = 0.3, 0.6, 0.85, u in [0.2, 0.8], 20001 samples, 1.5 deg stripes, cap_oblique_front camera"}


# ---------------------------------------------------------------------------
def plane_trace(surf, grid_uv, Sgrid, p0, nP, sign):
    us, vs = grid_uv
    f = np.sum((Sgrid - p0) * nP, -1)
    gen = contourpy.contour_generator(vs, us, f)
    lines = gen.lines(0.0)
    out = []
    for ln in lines:
        if len(ln) < 4:
            continue
        v, u = ln[:, 0].copy(), ln[:, 1].copy()
        for _ in range(6):
            j = surf.jet(u, v)
            fv = np.sum((j[(0, 0)] - p0) * nP, -1)
            gu = j[(1, 0)] @ nP
            gv = j[(0, 1)] @ nP
            g2 = gu * gu + gv * gv + 1e-300
            u = np.clip(u - fv * gu / g2, *surf.urange)
            v = np.clip(v - fv * gv / g2, *surf.vrange)
        d = differential(surf.jet(u, v), sign)
        t = unit(np.cross(d["n"], nP))
        nt = unit(d["n"] - np.sum(d["n"] * nP, -1, keepdims=True) * nP)
        kn = second_form(d, t, t)
        kc = kn / np.sum(d["n"] * nt, -1)
        k1, k2, _, _ = principal(d)
        out.append(dict(u=u, v=v, P=d["P"], kc=kc, kn=kn, k1=k1))
    return out


def section_planes():
    planes = []
    ug = [0.005, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.995]
    umid = [0.035, 0.075, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.925, 0.965]
    for kind, lst in (("guide", ug), ("between_guides", umid)):
        for u in lst:
            xa = 65 + 85 * design.xi_of(u)
            xl, zl = design.perimeter(design.w_of(u))
            d = np.array([xl - xa, 0, zl - ZA])
            nP = unit(np.cross(d, [0, 1, 0]))
            planes.append((f"{kind} plane u={u}", kind, np.array([xa, 0, ZA]), nP, u))
    for u in ug + umid:
        w = design.w_of(u)
        xl, zl = design.perimeter(w)
        e = 1e-6
        x2, z2 = design.perimeter(w + e)
        x1, z1 = design.perimeter(w - e)
        tng = unit(np.array([x2 - x1, 0, z2 - z1]))
        planes.append((f"lip-normal plane at L(u={u})", "lip_normal", np.array([xl, 0, zl]), tng, u))
    for x in [65.3, 66, 70, 80, 90, 100, 110, 120, 130, 140, 146, 149, 149.8]:
        planes.append((f"constant-X plane X={x}", "constant_X", np.array([x, 0, 0]), np.array([1.0, 0, 0]), None))
    for z in [240.5, 245, 250, 255, 260, 265, 270, 275, 278, 279.5]:
        planes.append((f"horizontal plane Z={z}", "horizontal", np.array([0, 0, z]), np.array([0, 0, 1.0]), None))
    return planes


def sections(cap, skins, cap_grid):
    us, vs, d = cap_grid[0], cap_grid[1], cap_grid[2]
    Scap = d["P"]
    zs = np.linspace(200, ZA, 81)
    ss = np.linspace(0, 1, 801)
    sk_grid = {side: (ss, zs, skins[side].grid(ss, zs)) for side in (+1, -1)}
    results = []
    for name, kind, p0, nP, u_nom in section_planes():
        cap_tr = plane_trace(cap, (us, vs), Scap, p0, nP, cap.sign)
        par_tr = []
        for side in (+1, -1):
            a, b, G = sk_grid[side]
            par_tr += plane_trace(skins[side], (a, b), G, p0, nP, skins[side].sign)
        if not cap_tr:
            continue
        u = np.concatenate([t["u"] for t in cap_tr])
        v = np.concatenate([t["v"] for t in cap_tr])
        P = np.concatenate([t["P"] for t in cap_tr])
        kc = np.concatenate([t["kc"] for t in cap_tr])
        kn = np.concatenate([t["kn"] for t in cap_tr])
        k1 = np.concatenate([t["k1"] for t in cap_tr])
        # reversals: sign changes of section curvature along each traced polyline (cap and parent parts)
        rev = []
        for t in cap_tr + par_tr:
            sgn = np.sign(np.where(np.abs(t["kc"]) < 1e-6, 0, t["kc"]))
            nz = np.nonzero(sgn)[0]
            for a_, b_ in zip(nz[:-1], nz[1:]):
                if sgn[a_] != sgn[b_]:
                    rev.append(r3(t["P"][b_]))
        i = int(np.argmax(kc))
        ti = float(att_thickness(u[i]))
        # lip crossing(s) on the cap: Y changes sign
        lips = []
        for t in cap_tr:
            y = t["P"][:, 1]
            for j in np.nonzero(np.sign(y[:-1]) * np.sign(y[1:]) < 0)[0]:
                a_ = abs(y[j]) / (abs(y[j]) + abs(y[j + 1]))
                kcl = (1 - a_) * t["kc"][j] + a_ * t["kc"][j + 1]
                knl = (1 - a_) * t["kn"][j] + a_ * t["kn"][j + 1]
                uu = (1 - a_) * t["u"][j] + a_ * t["u"][j + 1]
                lips.append({"at_mm": r3((1 - a_) * t["P"][j] + a_ * t["P"][j + 1]),
                             "section_radius_mm": float(1 / kcl), "normal_radius_across_lip_mm": float(1 / knl),
                             "chord_fraction": float(design.xi_of(uu)),
                             "attachment_thickness_mm": float(att_thickness(uu))})
        tmin = float(att_thickness(u[np.argmin(np.abs(v) - 2 * (np.abs(v) > 0.999))]))
        results.append({
            "plane": name, "family": kind, "nominal_u": u_nom, "points_on_cap": int(len(u)),
            "points_on_parent": int(sum(len(t["u"]) for t in par_tr)),
            "curvature_reversals": rev,
            "min_section_radius_on_cap_mm": float(1 / kc[i]), "at_mm": r3(P[i]),
            "at_chord_fraction": float(design.xi_of(u[i])), "attachment_thickness_there_mm": ti,
            "required_radius_mm": float(required_radius(ti)),
            "min_normal_radius_along_trace_mm": float(1 / kn.max()),
            "min_principal_radius_on_trace_mm": float(1 / k1.max()),
            "lip": lips,
            "max_section_obliquity_factor": float((kc / np.maximum(kn, 1e-12)).max()),
        })
    return results


# ---------------------------------------------------------------------------
def zebra(cap, skins):
    width = np.radians(5.0)
    views = {
        "side +Y": (0, -1, 0), "side -Y": (0, 1, 0),
        "front-left-up": (1, -1, -0.5), "rear-left-up": (-1, -1, -0.4), "above-left": (0.2, -0.6, -1),
        "front-right-up": (1, 1, -0.5), "rear-right-up": (-1, 1, -0.4), "above-right": (0.2, 0.6, -1),
        "low-front-left": (1, -0.7, 0.4), "low-rear-right": (-1, 0.7, 0.4),
    }
    angles = np.radians(np.arange(0, 180, 22.5))
    u = np.linspace(2e-3, 1 - 2e-3, 3001)
    worst = {"phase_mismatch_stripes": (0, None), "gradient_rel": (0, None), "direction_deg": (0, None)}
    ncross = 0
    for side, vc in ((+1, 1.0), (-1, -1.0)):
        dc = differential(cap.jet(u, np.full_like(u, vc)), cap.sign)
        dp = differential(skins[side].jet(u, np.full_like(u, ZA)), skins[side].sign)
        e1 = unit(dp["Su"])
        e2 = np.cross(dp["n"], e1)
        # cap normal frame for its own gradient
        e1c = unit(dc["Su"] - np.sum(dc["Su"] * dc["n"], -1, keepdims=True) * dc["n"])
        e2c = np.cross(dc["n"], e1c)
        Wc = shape_operator(dc, e1c, e2c)
        Wp = shape_operator(dp, e1, e2)
        for vname, V in views.items():
            V = unit(np.array(V, float))
            vis = (np.sum(dp["n"] * V, -1) < -0.05) & (np.sum(dc["n"] * V, -1) < -0.05)
            if not vis.any():
                continue
            up = unit(np.array([0, 0, 1.0]) - V[2] * V) if abs(V[2]) < 0.95 else unit(np.array([1.0, 0, 0]) - V[0] * V)
            side_ax = np.cross(V, up)
            for a in angles:
                dvec = np.cos(a) * up + np.sin(a) * side_ax
                pc, _ = zebra_phase(dc["n"], V, dvec, width)
                pp, _ = zebra_phase(dp["n"], V, dvec, width)
                gc = zebra_gradient(dc["n"], Wc, e1c, e2c, V, dvec, width)
                gp = zebra_gradient(dp["n"], Wp, e1, e2, V, dvec, width)
                dph = np.abs(pc - pp)
                gn = np.linalg.norm(gp, axis=-1)
                gabs = np.linalg.norm(gc - gp, axis=-1)
                grel = gabs / np.maximum(gn, 1e-9)
                worst_abs = float(np.max(np.where(vis, gabs, 0)))
                if worst_abs > worst.get("gradient_abs", (0,))[0]:
                    worst["gradient_abs"] = (worst_abs, None)
                tag = lambda i: {"view": vname, "stripe_angle_deg": float(np.degrees(a)), "skin": parent.SIDES[side],
                                 "u": float(u[i]), "at_mm": r3(dp["P"][i])}
                i = int(np.argmax(np.where(vis, dph, -1)))
                if dph[i] > worst["phase_mismatch_stripes"][0]:
                    worst["phase_mismatch_stripes"] = (float(dph[i]), tag(i))
                ok = vis & (gn > 1e-3)
                if ok.any():
                    i = int(np.argmax(np.where(ok, grel, -1)))
                    if grel[i] > worst["gradient_rel"][0]:
                        worst["gradient_rel"] = (float(grel[i]), tag(i))
                # light/dark crossings on the join: integer phase crossings along the edge
                fl = np.floor(pp)
                for j in np.nonzero((fl[1:] != fl[:-1]) & vis[1:] & vis[:-1])[0]:
                    ncross += 1
                    sc = unit(np.cross(dc["n"][j], gc[j]))
                    sp = unit(np.cross(dp["n"][j], gp[j]))
                    ang = float(np.degrees(np.arccos(np.clip(abs(sc @ sp), -1, 1))))
                    if ang > worst["direction_deg"][0]:
                        worst["direction_deg"] = (ang, tag(j))
    return {
        "stripe_model": "reflection of the view direction in the CAD normal; phase = asin(r . d) / 5 deg, light when floor(phase) is even",
        "views": list(views), "stripe_angles_deg": [float(x) for x in np.degrees(angles)],
        "join_samples_per_side": len(u), "u_range": [float(u[0]), float(u[-1])],
        "light_dark_crossings_checked": ncross,
        "max_phase_mismatch_stripes": worst["phase_mismatch_stripes"][0], "phase_mismatch_at": worst["phase_mismatch_stripes"][1],
        "max_absolute_gradient_mismatch_stripes_per_mm": worst.get("gradient_abs", (0,))[0],
        "max_relative_gradient_mismatch": worst["gradient_rel"][0], "gradient_mismatch_at": worst["gradient_rel"][1],
        "max_stripe_direction_change_deg": worst["direction_deg"][0], "direction_change_at": worst["direction_deg"][1],
    }


# ---------------------------------------------------------------------------
def main():
    shape = read_step(STEP)
    solid = TopoDS.Solid_s(next(iter_shapes(shape, TopAbs_SOLID)))
    topo = topology(shape)
    skins, cap, planes = classify(shape)
    rng = np.random.default_rng(1)
    out = {"step_file": "wing.step", "step_sha256": sha256(STEP), "step_bytes": STEP.stat().st_size,
           "reimport": topo,
           "faces": {"parent_skins": 2, "parent_root_plane": len(planes), "cap": {
               "poles": [cap.nu, cap.nv], "degrees": [cap.pu, cap.pv]}}}
    out["evaluator_crosscheck_vs_occ"] = evaluator_crosscheck([cap, skins[1], skins[-1]], solid, rng)
    out["parent_preserved"] = parent_preserved(skins)
    out["locked_dimensions"] = locked_dimensions(cap, skins, topo)
    out["joins"] = joins(cap, skins)
    cmap, grid = curvature_map(cap)
    out["principal_curvature_map"] = cmap
    out["parent_k2_near_attachment"] = parent_k2_near_attachment(skins)
    out["fit_vs_design_law"] = design_law_check(cap)
    out["sections"] = sections(cap, skins, grid)
    out["zebra"] = zebra(cap, skins)
    (HERE / "readback.json").write_text(json.dumps(out, indent=1))
    s = out["sections"]
    summ = {
        "valid": topo["brepcheck_valid"], "solids": topo["solids"], "volume": topo["volume_mm3"],
        "bbox": topo["bounding_box_mm"], "eval": out["evaluator_crosscheck_vs_occ"],
        "parent": out["parent_preserved"], "locked": out["locked_dimensions"],
        "joins": {k: {kk: vv for kk, vv in v.items() if kk.startswith("max")} for k, v in out["joins"].items()},
        "map": {k: v for k, v in cmap.items() if k != "local_maxima_refined"},
        "sections_min": min((x["min_section_radius_on_cap_mm"], x["plane"]) for x in s),
        "sections_min_normal": min((x["min_normal_radius_along_trace_mm"], x["plane"]) for x in s),
        "sections_reversals": sum(len(x["curvature_reversals"]) for x in s), "n_sections": len(s),
        "zebra": {k: v for k, v in out["zebra"].items() if k.startswith("max") or k.endswith("_at") or k == "light_dark_crossings_checked"},
    }
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
