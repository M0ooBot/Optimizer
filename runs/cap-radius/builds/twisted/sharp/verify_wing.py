"""Verify wing.step on the reimported STEP and write readback.json.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/twisted/sharp/verify_wing.py

Checks (skill "Verify the finished cap", steps 1-4):
  1  solid / closure / validity, parent skin preserved, locked span, attachment and outline,
     separation, tangent-plane angle and curvature in a common frame along all four joins
  2  sections through the cap (guide planes, planes between guides, rim-normal cuts) traced from the
     parent skin round the lip and back, with signed curvature, reversals, smallest radius, lip radius
  3  larger principal curvature of both cap faces on a dense grid with local refinement, judged
     against the printable rule at the guide's chord position
  4  numerical zebra (reflection-line phase and its surface gradient) across the four joins
All geometry is evaluated from the reimported STEP faces (CAD derivatives through OCP).
"""
import json
import pathlib
import sys
import time
import warnings

import numpy as np
from scipy.optimize import minimize

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cap_model as C  # noqa: E402  (parent.py loader, thickness table, rule)
import step_geom as SG  # noqa: E402

from OCP.TopAbs import TopAbs_SOLID, TopAbs_SHELL, TopAbs_FACE
from OCP.TopoDS import TopoDS
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRep import BRep_Tool
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds
from OCP.TopTools import TopTools_HSequenceOfShape
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.gp import gp_Pnt
from OCP.TopAbs import TopAbs_OUT, TopAbs_IN
from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf

warnings.filterwarnings("ignore", category=RuntimeWarning)
P = C.P
WING = HERE / "wing.step"
SIDES = ("right", "left")
MIRROR = np.array([1.0, -1.0, 1.0])


def chord_fraction(p):
    return float(C.chord_fraction(p))


def thickness(p):
    return float(C.attachment_thickness(C.chord_fraction(p)))


def rule(p):
    t = thickness(p)
    return (0.4 if t >= 0.8 else t / 2), t


# ----------------------------------------------------------------------------------------------
def check_solid(shape):
    solids = SG.explore(shape, TopAbs_SOLID)
    shells = SG.explore(shape, TopAbs_SHELL)
    faces = SG.explore(shape, TopAbs_FACE)
    solid = solids[0]
    fb = ShapeAnalysis_FreeBounds(solid, 1e-7, False, False)
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_EDGE
    free = 0
    ex = TopExp_Explorer(fb.GetClosedWires(), TopAbs_EDGE)
    while ex.More():
        free += 1; ex.Next()
    ex = TopExp_Explorer(fb.GetOpenWires(), TopAbs_EDGE)
    while ex.More():
        free += 1; ex.Next()
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(solid, props, 1e-9)
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(solid, box, False, False)
    xmin, ymin, zmin, xmax, ymax, zmax = box.Get()
    return solid, {
        "solids": len(solids), "shells": len(shells), "faces": len(faces),
        "brepcheck_valid": bool(BRepCheck_Analyzer(solid).IsValid()),
        "shell_closed": bool(BRep_Tool.IsClosed_s(shells[0])),
        "free_boundary_edges": free,
        "volume_mm3": props.Mass(),
        "bounding_box_mm": [xmin, xmax, ymin, ymax, zmin, zmax],
        "projected_span_mm": ymax - ymin,
    }


def check_outward(solid, faces):
    """Offset a few points along each face's outward normal and classify them."""
    bad = 0
    for key, F in faces.items():
        for u in (0.3, 0.6):
            v = 0.5 * (F.v0 + F.v1) if key[1] != "cap" else 0.4
            fr = F.frame(u, v)
            for s, want in ((+1e-3, TopAbs_OUT), (-1e-3, TopAbs_IN)):
                q = fr["X"] + s * fr["N"]
                cl = BRepClass3d_SolidClassifier(solid, gp_Pnt(*q), 1e-9)
                bad += cl.State() != want
    return bad


def check_parent(faces):
    pshape = SG.read_step(C.PARENT_DIR / "parent.step")
    pf = {}
    for f in SG.explore(pshape, TopAbs_FACE):
        from OCP.BRepAdaptor import BRepAdaptor_Surface
        from OCP.GeomAbs import GeomAbs_BSplineSurface
        if BRepAdaptor_Surface(TopoDS.Face_s(f)).GetType() != GeomAbs_BSplineSurface:
            continue
        F = SG.Face(f)
        mid = F.value(0.5, 400.0)
        key = ("right" if mid[1] > 0 else "left", "upper" if F.value(0.5, 0.0)[2] > 3 else "lower")
        pf[key] = F
    pole_diff = max(float(np.abs(pf[k].poles() - faces[k].poles()).max()) for k in pf)
    # sampled against the analytic parent definition
    s = np.linspace(0, 1, 201)
    y = np.linspace(0, 440, 23)
    dmax = 0.0
    for side in SIDES:
        for skin, sg in (("upper", +1), ("lower", -1)):
            F = faces[(side, skin)]
            for yy in y:
                ref = P.skin(s, yy if side == "right" else -yy, sg)
                got = np.array([F.value(u, yy) for u in s])
                dmax = max(dmax, float(np.linalg.norm(ref - got, axis=1).max()))
    return {"parent_face_poles_max_difference_to_parent_step_mm": pole_diff,
            "parent_skin_max_distance_to_parent_py_mm": dmax,
            "parent_skin_samples": int(4 * len(s) * len(y))}


def check_outline(faces):
    p = np.linspace(0, 1, 401)
    out = {}
    for side in SIDES:
        F = faces[(side, "cap")]
        ref = P.tip_outline(p, +1 if side == "right" else -1)
        rim = np.array([F.value(u, 1.0) for u in p])
        proj = []
        for q in ref[1:-1]:
            pr = GeomAPI_ProjectPointOnSurf(gp_Pnt(*q), F.surface)
            proj.append(pr.LowerDistance() if pr.NbPoints() else np.nan)
        # the rim is the outermost point of every guide (largest |Y| along w)
        guide_outermost = []
        for u in np.linspace(0.02, 0.98, 49):
            ww = np.linspace(0, 2, 2001)
            Y = np.abs([F.value(u, w)[1] for w in ww])
            guide_outermost.append(abs(ww[np.argmax(Y)] - 1.0))
        att = np.array([[F.value(u, 0.0)[1], F.value(u, 2.0)[1]] for u in p])
        out[side] = {"rim_param_max_distance_to_locked_outline_mm": float(np.linalg.norm(rim - ref, axis=1).max()),
                     "locked_outline_max_distance_to_cap_face_mm": float(np.nanmax(proj)),
                     "max_abs_Y_of_cap_mm": float(np.abs(rim[:, 1]).max()),
                     "guide_outermost_point_offset_from_rim_in_w": float(max(guide_outermost)),
                     "attachment_edge_abs_Y_range_mm": [float(np.abs(att).min()), float(np.abs(att).max())]}
    return out


# ----------------------------------------------------------------------------------------------
def shape_operator_dN(fr, a):
    """Derivative of the outward unit normal along tangent vector a (Weingarten)."""
    Su, Sv = fr["Su"], fr["Sv"]
    I = np.array([[fr["E"], fr["F"]], [fr["F"], fr["G"]]])
    II = np.array([[fr["L"], fr["M"]], [fr["M"], fr["Nn"]]])
    ca = np.linalg.solve(I, np.array([a @ Su, a @ Sv]))
    coef = -np.linalg.solve(I, II @ ca)
    return coef[0] * Su + coef[1] * Sv


VIEWS = {
    "top": (0, 0, -1), "bottom": (0, 0, 1), "front": (1, 0, 0), "rear": (-1, 0, 0),
    "outboard": (0, -1, 0), "oblique_front_top": (1, -1, -1), "oblique_rear_top": (-1, -1, -1),
    "oblique_front_bottom": (1, -1, 1), "oblique_rear_bottom": (-1, -1, 1),
}
STRIPE_ANGLES = (0, 30, 60, 90, 120, 150)
N_BANDS = 24  # stripe frequency: phase = pi * N_BANDS * (R . d)


def screen_basis(V):
    V = np.asarray(V, float); V = V / np.linalg.norm(V)
    up = np.array([0, 0, 1.0]) if abs(V[2]) < 0.9 else np.array([1.0, 0, 0])
    e1 = np.cross(up, V); e1 /= np.linalg.norm(e1)
    e2 = np.cross(V, e1)
    return V, e1, e2


def zebra_at(fr, V, d, tdir, cdir):
    N = fr["N"]
    R = V - 2 * (V @ N) * N
    phase = np.pi * N_BANDS * (R @ d)
    grads = []
    for a in (tdir, cdir):
        dN = shape_operator_dN(fr, a)
        dR = -2 * ((V @ dN) * N + (V @ N) * dN)
        grads.append(np.pi * N_BANDS * (dR @ d))
    g = grads[0] * tdir + grads[1] * cdir
    return phase, g


def check_joins(faces, n=4001, zebra=True):
    u = np.linspace(1e-3, 1 - 1e-3, n)
    results = {}
    for side in SIDES:
        cap = faces[(side, "cap")]
        mir = MIRROR if side == "left" else np.ones(3)
        for skin, w0 in (("upper", 0.0), ("lower", 2.0)):
            par = faces[(side, skin)]
            sep, ang = np.zeros(n), np.zeros(n)
            dII = np.zeros((n, 3)); kref = np.zeros(n)
            frames = []
            for i, uu in enumerate(u):
                fp = par.frame(uu, 440.0)
                fc = cap.frame(uu, w0)
                sep[i] = np.linalg.norm(fp["X"] - fc["X"])
                ang[i] = np.degrees(np.arccos(np.clip(fp["N"] @ fc["N"], -1, 1)))
                t = fp["Su"] / np.linalg.norm(fp["Su"])
                c = np.cross(fp["N"], t)
                kp = [SG.second_form_in(fp, t, t), SG.second_form_in(fp, c, c), SG.second_form_in(fp, t, c)]
                # cap evaluated in the same 3D frame projected into its own tangent plane
                tc = t - (t @ fc["N"]) * fc["N"]; tc /= np.linalg.norm(tc)
                cc = np.cross(fc["N"], tc)
                kc = [SG.second_form_in(fc, tc, tc), SG.second_form_in(fc, cc, cc), SG.second_form_in(fc, tc, cc)]
                dII[i] = np.array(kc) - np.array(kp)
                kref[i] = max(abs(np.array(kp)).max(), 1e-12)
                frames.append((fp, fc, t, c, tc, cc))
            key = f"{side}_{skin}"
            i_sep, i_ang = sep.argmax(), ang.argmax()
            rel = np.abs(dII).max(1) / kref
            res = {
                "samples": n, "u_range": [float(u[0]), float(u[-1])],
                "max_separation_mm": float(sep.max()), "at_chord_fraction": chord_fraction(u[i_sep]),
                "max_tangent_plane_angle_deg": float(ang.max()), "angle_at_chord_fraction": chord_fraction(u[i_ang]),
                "max_tangent_plane_angle_deg_for_chord_fraction_0.001_to_0.999": float(ang[(C.chord_fraction(u) > 1e-3) & (C.chord_fraction(u) < 0.999)].max()),
                "curvature_common_frame": {
                    "components": "II(t,t) along the edge, II(c,c) across the edge, II(t,c) twist; outward normal, 1/mm",
                    "max_abs_difference_per_mm": [float(x) for x in np.abs(dII).max(0)],
                    "max_difference_relative_to_local_parent_curvature": float(rel.max()),
                    "at_chord_fraction": chord_fraction(u[rel.argmax()]),
                    "max_relative_difference_chord_fraction_0.001_to_0.999": float(rel[(C.chord_fraction(u) > 1e-3) & (C.chord_fraction(u) < 0.999)].max()),
                },
            }
            if zebra:
                res["zebra"] = zebra_join(frames, u, mir)
            results[key] = res
    return results


def zebra_join(frames, u, mir):
    worst = {"phase_mismatch_periods": (0, None), "gradient_mismatch_relative": (0, None),
             "stripe_direction_change_at_crossings_deg": (0, None)}
    n_cross = 0
    for vname, V0 in VIEWS.items():
        V, e1, e2 = screen_basis(np.asarray(V0, float) * mir)
        for ang in STRIPE_ANGLES:
            a = np.radians(ang)
            d = np.cos(a) * e1 + np.sin(a) * e2
            if mir[1] < 0:
                d = d * 1.0  # screen basis is already built from the mirrored view
            ph_p, ph_c, gp_, gc_ = [], [], [], []
            for fp, fc, t, c, tc, cc in frames:
                p1, g1 = zebra_at(fp, V, d, t, c)
                p2, g2 = zebra_at(fc, V, d, tc, cc)
                ph_p.append(p1); ph_c.append(p2); gp_.append(g1); gc_.append(g2)
            ph_p, ph_c, gp_, gc_ = map(np.array, (ph_p, ph_c, gp_, gc_))
            dph = np.abs(ph_c - ph_p) / (2 * np.pi)
            gscale = np.maximum(np.linalg.norm(gp_, axis=1), np.linalg.norm(gc_, axis=1))
            gm = np.linalg.norm(gc_ - gp_, axis=1) / np.maximum(gscale, 1e-9)
            ok = gscale > 1e-3 * gscale.max()   # ignore points where the stripes are locally stationary
            i = int(np.argmax(dph))
            if dph[i] > worst["phase_mismatch_periods"][0]:
                worst["phase_mismatch_periods"] = (float(dph[i]), dict(view=vname, stripe_angle=ang, chord_fraction=chord_fraction(u[i])))
            if np.any(ok):
                i = int(np.argmax(np.where(ok, gm, 0)))
                if gm[i] > worst["gradient_mismatch_relative"][0]:
                    worst["gradient_mismatch_relative"] = (float(gm[i]), dict(view=vname, stripe_angle=ang, chord_fraction=chord_fraction(u[i])))
            # light/dark crossings: phase = k pi on the parent side of the join
            s = np.sin(ph_p)
            idx = np.where(np.sign(s[:-1]) * np.sign(s[1:]) < 0)[0]
            for i in idx:
                n_cross += 1
                # stripe direction = tangent perpendicular to the phase gradient
                def sdir(g, N):
                    v = np.cross(N, g); return v / max(np.linalg.norm(v), 1e-300)
                fp, fc = frames[i][0], frames[i][1]
                a1, a2 = sdir(gp_[i], fp["N"]), sdir(gc_[i], fc["N"])
                dang = np.degrees(np.arccos(np.clip(abs(a1 @ a2), 0, 1)))
                if np.linalg.norm(gp_[i]) > 1e-3 * gscale.max() and dang > worst["stripe_direction_change_at_crossings_deg"][0]:
                    worst["stripe_direction_change_at_crossings_deg"] = (float(dang), dict(view=vname, stripe_angle=ang, chord_fraction=chord_fraction(u[i])))
    out = {k: {"value": v[0], "where": v[1]} for k, v in worst.items()}
    out["light_dark_crossings_checked"] = n_cross
    out["views"] = list(VIEWS)
    out["stripe_angles_deg"] = list(STRIPE_ANGLES)
    out["bands_per_unit_reflection"] = N_BANDS
    return out


# ----------------------------------------------------------------------------------------------
def kmax_at(F, u, v):
    return SG.principal_from_frame(F.frame(u, v))[0]


def curvature_map(F, side):
    s = np.linspace(0, 1, 601)
    p = np.unique(np.concatenate([0.5 - 0.5 * np.cos(np.pi * s), np.geomspace(1e-4, 0.02, 60), 1 - np.geomspace(1e-4, 0.02, 60)]))
    p = p[(p >= 1e-4) & (p <= 1 - 1e-4)]
    t = np.linspace(0, 1, 401)
    wu = np.unique(np.concatenate([t, 1 - (1 - t) ** 2]))
    w = np.unique(np.concatenate([wu, 2 - wu]))
    K = np.zeros((len(p), len(w)))
    for i, pp in enumerate(p):
        for j, ww in enumerate(w):
            K[i, j] = kmax_at(F, pp, ww)
    req = np.array([rule(pp)[0] for pp in p])
    tt = np.array([rule(pp)[1] for pp in p])
    R = 1 / K
    ratio = R / req[:, None]
    # local minima of the radius ratio -> refine
    cand = []
    for i in range(len(p)):
        j = int(np.argmin(ratio[i]))
        lo, hi = max(i - 1, 0), min(i + 1, len(p) - 1)
        if ratio[i, j] <= ratio[lo].min() and ratio[i, j] <= ratio[hi].min():
            cand.append((ratio[i, j], i, j))
    cand.sort()
    refined = []
    for _, i, j in cand[:40]:
        x0 = np.array([p[i], w[j]])
        dp = 0.5 * (p[min(i + 1, len(p) - 1)] - p[max(i - 1, 0)])
        dw = 0.01

        def f(x):
            pp, ww = np.clip(x[0], 1e-4, 1 - 1e-4), np.clip(x[1], 0, 2)
            return (1 / kmax_at(F, pp, ww)) / rule(pp)[0]
        o = minimize(f, x0, method="Nelder-Mead",
                     options=dict(initial_simplex=[x0, x0 + [dp, 0], x0 + [0, dw]], xatol=1e-7, fatol=1e-9, maxiter=300))
        pp, ww = float(np.clip(o.x[0], 1e-4, 1 - 1e-4)), float(np.clip(o.x[1], 0, 2))
        r = 1 / kmax_at(F, pp, ww)
        req_, t_ = rule(pp)
        refined.append(dict(p=pp, w=ww, chord_fraction=chord_fraction(pp), radius_mm=r, attachment_thickness_mm=t_,
                            required_mm=req_, ratio=r / req_, point_mm=[float(x) for x in F.value(pp, ww)]))
    refined.sort(key=lambda d: d["ratio"])
    thick = tt >= 0.8
    i, j = np.unravel_index(np.argmin(np.where(thick[:, None], R, np.inf)), R.shape)
    thick_min = dict(radius_mm=float(R[i, j]), p=float(p[i]), w=float(w[j]), chord_fraction=chord_fraction(p[i]),
                     attachment_thickness_mm=float(tt[i]), point_mm=[float(x) for x in F.value(p[i], w[j])])
    # refine the thick-zone minimum
    def g(x):
        pp, ww = np.clip(x[0], 1e-4, 1 - 1e-4), np.clip(x[1], 0, 2)
        if thickness(pp) < 0.8:
            return 1e3
        return 1 / kmax_at(F, pp, ww)
    o = minimize(g, [p[i], w[j]], method="Nelder-Mead", options=dict(xatol=1e-8, fatol=1e-10, maxiter=400))
    thick_min_ref = dict(radius_mm=float(o.fun), p=float(o.x[0]), w=float(o.x[1]), chord_fraction=chord_fraction(o.x[0]),
                         attachment_thickness_mm=thickness(o.x[0]), point_mm=[float(x) for x in F.value(o.x[0], o.x[1])])
    thin_aft = (~thick) & (p > 0.5)
    thin_le = (~thick) & (p < 0.5)
    lip_j = int(np.argmin(np.abs(w - 1.0)))
    lip_ratio = R[:, lip_j] / (tt / 2)
    rmin_row = R.min(1)
    le_pt = F.value(0.0, 1.0); te_pt = F.value(1.0, 1.0)
    reach_le = max(np.linalg.norm(F.value(p[0], ww) - le_pt) for ww in w[::20])
    reach_te = max(np.linalg.norm(F.value(p[-1], ww) - te_pt) for ww in w[::20])
    by_p = []
    for pp_ in (0.0005, 0.005, 0.0133, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.88, 0.8925, 0.9, 0.95, 0.98, 0.995, 0.999):
        i = int(np.argmin(np.abs(p - pp_)))
        jj = int(np.argmin(R[i]))
        by_p.append(dict(p=float(p[i]), chord_fraction=chord_fraction(p[i]), attachment_thickness_mm=float(tt[i]),
                         required_mm=float(req[i]), min_radius_mm=float(R[i, jj]), at_w=float(w[jj]),
                         lip_radius_mm=float(R[i, lip_j]), ratio=float(R[i, jj] / req[i])))
    return {
        "grid": {"p_points": int(len(p)), "w_points": int(len(w)), "p_range": [float(p[0]), float(p[-1])],
                 "collapsed_endpoints_excluded": "p = 0 (leading-edge point) and p = 1 (trailing-edge point)",
                 "closest_grid_row_to_LE_point_mm": float(reach_le), "closest_grid_row_to_TE_point_mm": float(reach_te)},
        "chord_position_used": "chord fraction of the guide, sin^2(pi p / 2) (the rim point's chord fraction)",
        "thick_zone_min_radius_grid": thick_min,
        "thick_zone_min_radius_refined": thick_min_ref,
        "thick_zone_passes_0p4": bool(thick_min_ref["radius_mm"] >= 0.4 and thick_min["radius_mm"] >= 0.4),
        "thin_aft_zone": {"chord_fraction_range": [chord_fraction(p[thin_aft].min()), chord_fraction(p[thin_aft].max())],
                          "min_radius_over_half_thickness": float((rmin_row / (tt / 2))[thin_aft].min()),
                          "lip_radius_over_half_thickness_range": [float(lip_ratio[thin_aft].min()), float(lip_ratio[thin_aft].max())]},
        "thin_le_zone": {"chord_fraction_range": [chord_fraction(p[thin_le].min()), chord_fraction(p[thin_le].max())],
                         "min_radius_mm": float(rmin_row[thin_le].min()),
                         "max_half_thickness_mm": float((tt / 2)[thin_le].max())},
        "worst_refined_ratios": refined[:8],
        "by_chord_position": by_p,
    }


# ----------------------------------------------------------------------------------------------
def trace(F, O, m, uv0, direction, h_rel=0.03, h_max=0.5, max_len=1e9, stop=None, vrange=None):
    """March along {X : (X - O) . m = 0} on face F from uv0; returns list of (u, v)."""
    vlo, vhi = vrange if vrange else (F.v0, F.v1)

    def f_and_grad(uv):
        X, Su, Sv, *_ = F.d2(*uv)
        return (X - O) @ m, np.array([Su @ m, Sv @ m]), Su, Sv

    def correct(uv):
        for _ in range(20):
            f, g, _, _ = f_and_grad(uv)
            if abs(f) < 1e-12:
                break
            uv = uv - f * g / (g @ g)
        return uv
    uv = correct(np.array(uv0, float))
    pts = [uv.copy()]
    length = 0.0
    prev_t = None
    for _ in range(200000):
        f, g, Su, Sv = f_and_grad(uv)
        tau = np.array([-g[1], g[0]]) * direction
        dS = Su * tau[0] + Sv * tau[1]
        if prev_t is not None and dS @ prev_t < 0:
            tau, dS = -tau, -dS
        k = max(kmax_at(F, *uv), 1e-6)
        h = min(h_max, h_rel / k)
        step = h / np.linalg.norm(dS)
        new = correct(uv + tau * step)
        X0 = F.value(*uv)
        if new[0] < 0 or new[0] > 1 or new[1] < vlo or new[1] > vhi:
            # land on the boundary (bisection on the step)
            lo, hi = 0.0, step
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                c_ = correct(uv + tau * mid)
                if c_[0] < 0 or c_[0] > 1 or c_[1] < vlo or c_[1] > vhi:
                    hi = mid
                else:
                    lo = mid
            new = correct(uv + tau * lo)
            new[1] = min(max(new[1], vlo), vhi)
            pts.append(new)
            break
        prev_t = F.value(*new) - X0
        length += np.linalg.norm(prev_t)
        pts.append(new)
        uv = new
        if length > max_len or (stop and stop(uv)):
            break
    return pts


def section_curvature(F, uvs, m):
    rows = []
    for k in range(len(uvs)):
        a = uvs[max(k - 1, 0)]; b = uvs[min(k + 1, len(uvs) - 1)]
        fr = F.frame(*uvs[k])
        T = F.value(*b) - F.value(*a)
        T = T - (T @ fr["N"]) * fr["N"]
        T /= np.linalg.norm(T)
        bvec = np.cross(m, T)
        kn = SG.second_form_in(fr, T, T)
        denom = bvec @ fr["N"]
        kappa_convex = -kn / abs(denom) if abs(denom) > 1e-12 else np.nan
        rows.append((fr["X"], kappa_convex, kn))
    return rows


def run_section(faces, side, p0, kind):
    """Trace one section plane through the cap (and onto the parent skins for guide-type planes)."""
    cap = faces[(side, "cap")]
    rim = cap.value(p0, 1.0)
    if kind in ("guide", "between"):
        A, Bp = cap.value(p0, 0.0), cap.value(p0, 2.0)
        m = np.cross(rim - A, Bp - A); m /= np.linalg.norm(m)
        O = rim
    else:  # rim-normal cut: plane through the rim point normal to the rim tangent
        _, Su, *_ = cap.d2(p0, 1.0)
        m = Su / np.linalg.norm(Su); O = rim
    t_att = thickness(p0)
    max_len = max(6 * t_att, 3.0) if kind == "rim_normal" else 1e9
    stop = (lambda uv: uv[1] < 0.15 or uv[1] > 1.85) if kind == "rim_normal" else None
    seg_up = trace(cap, O, m, (p0, 1.0), +1, max_len=max_len, stop=stop)
    seg_lo = trace(cap, O, m, (p0, 1.0), -1, max_len=max_len, stop=stop)
    # orient: upper part reversed so the trace runs upper -> lip -> lower
    if np.mean([uv[1] for uv in seg_up]) > 1.0:
        seg_up, seg_lo = seg_lo, seg_up
    cap_uv = seg_up[::-1] + seg_lo[1:]
    rows_cap = section_curvature(cap, cap_uv, m)
    parent_rows = []
    if kind in ("guide", "between"):
        for skin, end in (("upper", cap_uv[0]), ("lower", cap_uv[-1])):
            if not (abs(end[1]) < 1e-6 or abs(end[1] - 2) < 1e-6):
                continue
            par = faces[(side, skin)]
            seg = trace(par, O, m, (end[0], 440.0), +1, max_len=15.0, vrange=(0.0, 440.0))
            if len(seg) > 2 and seg[2][1] > 440 - 1e-9:
                seg = trace(par, O, m, (end[0], 440.0), -1, max_len=15.0, vrange=(0.0, 440.0))
            parent_rows.append((skin, section_curvature(par, seg, m)))
    kap = np.array([r[1] for r in rows_cap])
    X = np.array([r[0] for r in rows_cap])
    rad = 1 / np.abs(kap)
    i_min = int(np.nanargmin(rad))
    lip_i = int(np.argmin([abs(uv[1] - 1.0) for uv in cap_uv]))
    sgn = np.sign(kap)
    big = np.abs(kap) > 1e-3 * np.nanmax(np.abs(kap))
    rev = [i for i in range(1, len(kap)) if big[i] and big[i - 1] and sgn[i] != sgn[i - 1]]
    # reversals are also checked with sign changes separated by small-curvature stretches
    signs_big = [(i, sgn[i]) for i in range(len(kap)) if big[i]]
    rev2 = [signs_big[k][0] for k in range(1, len(signs_big)) if signs_big[k][1] != signs_big[k - 1][1]]
    par_min = None
    for skin, rows in parent_rows:
        kk = np.array([r[1] for r in rows])
        if len(kk):
            par_min = min(par_min or 1e9, float(np.nanmin(1 / np.abs(kk))))
    req, t = rule(p0)
    _, Tr, *_ = cap.d2(p0, 1.0)
    Tr = Tr / np.linalg.norm(Tr)
    return {
        "rim_tangent_dot_plane_normal": float(abs(Tr @ m)),
        "kind": kind, "p": float(p0), "chord_fraction": chord_fraction(p0), "attachment_thickness_mm": t,
        "required_mm": req, "cap_points": len(cap_uv),
        "cap_trace_ends_w": [float(cap_uv[0][1]), float(cap_uv[-1][1])],
        "parent_trace_points": [len(r) for _, r in parent_rows],
        "min_radius_mm": float(rad[i_min]), "min_at_w": float(cap_uv[i_min][1]), "min_at_mm": [float(x) for x in X[i_min]],
        "lip_radius_mm": float(rad[lip_i]),
        "bending_reversals_on_cap": len(rev2),
        "reversal_w": [float(cap_uv[i][1]) for i in rev2][:10],
        "min_radius_on_parent_stretch_mm": par_min,
    }


def sections(faces):
    out = []
    guides = [0.005, 0.0133, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.88, 0.8925, 0.9, 0.95, 0.98]
    between = [0.035, 0.075, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.825, 0.865, 0.925, 0.965]
    for side in SIDES:
        for p0 in guides:
            out.append(dict(side=side, **run_section(faces, side, p0, "guide")))
        for p0 in between:
            out.append(dict(side=side, **run_section(faces, side, p0, "between")))
        for p0 in guides + between:
            out.append(dict(side=side, **run_section(faces, side, p0, "rim_normal")))
    return out


# ----------------------------------------------------------------------------------------------
def main():
    t0 = time.time()
    shape = SG.read_step(WING)
    solid, solid_info = check_solid(shape)
    faces = SG.wing_faces(shape)
    rb = {"file": "wing.step", "sha256": SG.sha256(WING), "bytes": WING.stat().st_size,
          "units": "mm", "axes": P.AXES, "solid": solid_info,
          "outward_normal_classifier_failures": check_outward(solid, faces),
          "face_summary": {f"{k[0]}_{k[1]}": {"u_poles": F.nu, "v_poles": F.nv, "v_range": [F.v0, F.v1]} for k, F in faces.items()}}
    rb["parent_preserved"] = check_parent(faces)
    rb["locked_outline_and_attachment"] = check_outline(faces)
    print("solid/parent/outline done", time.time() - t0, flush=True)
    rb["joins"] = check_joins(faces)
    print("joins done", time.time() - t0, flush=True)
    secs = sections(faces)
    rb["sections"] = secs
    def worst(kind):
        sel = [s for s in secs if s["kind"] == kind]
        thick = [s for s in sel if s["attachment_thickness_mm"] >= 0.8]
        thin = [s for s in sel if s["attachment_thickness_mm"] < 0.8]
        o = {"count": len(sel)}
        if thick:
            w = min(thick, key=lambda s: s["min_radius_mm"])
            o["thick_min_radius_mm"] = w["min_radius_mm"]; o["thick_min_at"] = dict(side=w["side"], p=w["p"], w=w["min_at_w"])
        if thin:
            w = min(thin, key=lambda s: s["min_radius_mm"] / (s["attachment_thickness_mm"] / 2))
            o["thin_min_radius_over_half_thickness"] = w["min_radius_mm"] / (w["attachment_thickness_mm"] / 2)
            o["thin_lip_radius_over_half_thickness_range"] = [min(s["lip_radius_mm"] / (s["attachment_thickness_mm"] / 2) for s in thin),
                                                              max(s["lip_radius_mm"] / (s["attachment_thickness_mm"] / 2) for s in thin)]
        o["sections_with_bending_reversal_on_cap"] = sum(1 for s in sel if s["bending_reversals_on_cap"])
        return o
    rb["sections_summary"] = {k: worst(k) for k in ("guide", "between", "rim_normal")}
    # guide-plane lip radius against the true (rim-normal) lip radius: an oblique plane through the
    # rim sees the lip radius scaled by |rim tangent . plane normal| (Meusnier), which is why guide
    # sections near the leading and trailing edges show small in-plane radii
    obl = []
    for s in secs:
        if s["kind"] != "guide":
            continue
        rn = [r for r in secs if r["kind"] == "rim_normal" and r["side"] == s["side"] and r["p"] == s["p"]][0]
        obl.append(dict(side=s["side"], p=s["p"], chord_fraction=s["chord_fraction"],
                        guide_plane_lip_radius_mm=s["lip_radius_mm"], rim_normal_lip_radius_mm=rn["lip_radius_mm"],
                        ratio=s["lip_radius_mm"] / rn["lip_radius_mm"], rim_tangent_dot_plane_normal=s["rim_tangent_dot_plane_normal"]))
    rb["sections_obliquity"] = obl
    print("sections done", time.time() - t0, flush=True)
    rb["curvature_map"] = {side: curvature_map(faces[(side, "cap")], side) for side in SIDES}
    print("map done", time.time() - t0, flush=True)
    rb["runtime_s"] = time.time() - t0
    (HERE / "readback.json").write_text(json.dumps(rb, indent=1))
    print(json.dumps({k: rb[k] for k in ("sha256", "solid", "parent_preserved", "locked_outline_and_attachment", "sections_summary")}, indent=1))


if __name__ == "__main__":
    main()
