"""Synthetic calibration: geometry whose curvature is known in closed form, measured with audit.py.

1. Primitives: a sphere (R = 3 mm) and a torus (R = 10, r = 2 mm), NURBS-converted, checked point by point
   against their principal curvatures and against the in-plane curvature of plane sections, signs included.
2. Whole models: the flat parents (sharp and blunt) with each attachment section's upper curve revolved
   180 degrees about the attachment chord line. Every constant-X cut of such a cap is a half circle whose
   radius is the local half thickness h = t/2, so the section audit must return min radius = lip radius = t/2
   (ratio 1) with no reversals; the cap's principal curvatures are 1/(h sqrt(1 + h_x^2)) round the circle and
   -h_xx/(1 + h_x^2)^1.5 along the meridian; the join curvature mismatch is the round-the-circle term; the lip
   lies at |Y| = 350 + h, so its offset from the locked outline is known; the volume adds pi * int h^2 dX.
   On the blunt parent the 0.4 mm TE face is closed by a flat half disc of radius 0.2 mm.

    python runs/cap-radius/audit/calibration/synthetic.py
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeSolid, BRepBuilderAPI_MakeWire,
                                BRepBuilderAPI_NurbsConvert, BRepBuilderAPI_Sewing)
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepPrimAPI import BRepPrimAPI_MakeRevol, BRepPrimAPI_MakeSphere, BRepPrimAPI_MakeTorus
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.ShapeFix import ShapeFix_Solid
from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_SHELL
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Ax1, gp_Dir, gp_Pln, gp_Pnt
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import audit as A  # noqa: E402

OUT = HERE / "synthetic"


# ---------------------------------------------------------------------------
# 1. Primitives
# ---------------------------------------------------------------------------

def primitive_checks():
    rows = {}
    # Sphere: both principal curvatures 1/3, convex.
    sphere = BRepBuilderAPI_NurbsConvert(BRepPrimAPI_MakeSphere(3.0).Shape(), True).Shape()
    faces, _ = A.shapes_of(sphere, TopAbs_FACE)
    fr = A.Face(0, TopoDS.Face_s(faces[0]))
    u0, u1, v0, v1 = fr.uv
    errs = []
    for u in np.linspace(u0, u1, 23)[1:-1]:
        for v in np.linspace(v0, v1, 23)[1:-1]:
            j = fr.jet(u, v)
            if j["N"] is not None:
                errs.append(np.max(np.abs(j["k"] - 1 / 3)))
    rows["sphere_principal_max_abs_error_per_mm"] = float(max(errs))
    rows["sphere_points"] = len(errs)
    # Torus about Z: tube curvature 1/r (convex), ring curvature cos(phi)/(R + r cos(phi)) (concave inside).
    R, r = 10.0, 2.0
    torus = BRepBuilderAPI_NurbsConvert(BRepPrimAPI_MakeTorus(R, r).Shape(), True).Shape()
    faces, _ = A.shapes_of(torus, TopAbs_FACE)
    fr = A.Face(0, TopoDS.Face_s(faces[0]))
    u0, u1, v0, v1 = fr.uv
    errs, signs = [], 0
    for u in np.linspace(u0, u1, 25)[1:-1]:
        for v in np.linspace(v0, v1, 25)[1:-1]:
            j = fr.jet(u, v)
            if j["N"] is None:
                continue
            P = j["P"]
            rho = math.hypot(P[0], P[1])
            cphi = (rho - R) / r
            expect = np.sort([1 / r, cphi / (R + r * cphi)])
            errs.append(np.max(np.abs(j["k"] - expect)))
            signs += int(cphi < -0.1)
    rows["torus_principal_max_abs_error_per_mm"] = float(max(errs))
    rows["torus_points"] = len(errs)
    rows["torus_points_with_a_concave_principal_curvature"] = signs
    # Torus sections by Z = h: two circles, outer convex (+1/(R + w)), inner concave (-1/(R - w)), w = sqrt(r^2 - h^2).
    sec_rows = []
    for h in (0.0, 1.0, 1.9):
        w = math.sqrt(r * r - h * h)
        O, n = np.array([0, 0, h]), np.array([0.0, 0, 1])
        sec = BRepAlgoAPI_Section(fr.face, gp_Pln(gp_Pnt(0, 0, h), gp_Dir(0, 0, 1)), False)
        sec.ComputePCurveOn1(True)
        sec.Approximation(True)
        sec.Build()
        edges, _ = A.shapes_of(sec.Shape(), TopAbs_EDGE)
        for e in edges:
            seg = A.Segment(fr, TopoDS.Edge_s(e), O, n)
            rows_ = seg.sample()
            rho = np.array([math.hypot(x["P"][0], x["P"][1]) for x in rows_])
            kap = np.array([x["kappa"] for x in rows_])
            outer = rho > R
            exp = np.where(outer, 1 / (R + w), -1 / (R - w))
            sec_rows.append({"z": h, "branch": "outer" if outer.mean() > 0.5 else "inner",
                             "max_abs_error_per_mm": float(np.nanmax(np.abs(kap - exp))),
                             "expected_per_mm": float(exp[0]), "samples": len(kap)})
    rows["torus_sections"] = sec_rows
    # Oblique plane through the sphere centre tilted 40 degrees: a great circle, radius 3 whatever the tilt.
    sphere_face = A.Face(0, TopoDS.Face_s(A.shapes_of(sphere, TopAbs_FACE)[0][0]))
    n = A.unit(np.array([math.sin(math.radians(40)), 0, math.cos(math.radians(40))]))
    sec = BRepAlgoAPI_Section(sphere_face.face, gp_Pln(gp_Pnt(0, 0, 0.5), gp_Dir(*n)), False)
    sec.ComputePCurveOn1(True)
    sec.Approximation(True)
    sec.Build()
    d = 0.5 * n[2]
    errs = []
    for e in A.shapes_of(sec.Shape(), TopAbs_EDGE)[0]:
        seg = A.Segment(sphere_face, TopoDS.Edge_s(e), np.array([0, 0, 0.5]), n)
        kap = np.array([x["kappa"] for x in seg.sample()])
        errs.append(np.nanmax(np.abs(kap - 1 / math.sqrt(9 - d * d))))
    rows["sphere_oblique_section"] = {"plane_offset_mm": d, "expected_radius_mm": math.sqrt(9 - d * d),
                                      "max_abs_error_per_mm": float(max(errs))}
    return rows


# ---------------------------------------------------------------------------
# 2. Revolved-cap whole models
# ---------------------------------------------------------------------------

def edge_mid(edge):
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    c = BRepAdaptor_Curve(edge)
    return np.array(c.Value(0.5 * (c.FirstParameter() + c.LastParameter())).Coord())


def build_revolved(variant):
    case = A.Case("flat", variant)
    shape = A.read_step(A.PARENTS / "flat" / variant / "parent.step")
    faces = [TopoDS.Face_s(f) for f in A.shapes_of(shape, TopAbs_FACE)[0]]
    keep, caps = [], []
    for f in faces:
        fr = A.Face(0, f)
        u0, u1, v0, v1 = fr.uv
        c = fr.point(0.5 * (u0 + u1), 0.5 * (v0 + v1))
        if fr.planar and abs(abs(c[1]) - 350) < 1e-6:
            edges = [TopoDS.Edge_s(e) for e in A.shapes_of(f, TopAbs_EDGE)[0]]
            upper = max(edges, key=lambda e: edge_mid(e)[2])
            te_seg = [e for e in edges if abs(edge_mid(e)[0] - 120) < 1e-3 and abs(edge_mid(e)[2]) < 1e-6]
            side = 1 if c[1] > 0 else -1
            axis = gp_Ax1(gp_Pnt(0, side * 350.0, 0), gp_Dir(-side, 0, 0))   # sweeps the upper curve outboard
            rev = BRepPrimAPI_MakeRevol(upper, axis, math.pi)
            cap = rev.Shape()
            caps.append(cap)
            if te_seg:
                ends = [TopoDS.Edge_s(e) for e in A.shapes_of(cap, TopAbs_EDGE)[0]]
                semi = [e for e in ends if abs(edge_mid(e)[0] - 120) < 1e-3 and abs(edge_mid(e)[2]) < 1e-6
                        and abs(abs(edge_mid(e)[1]) - 350) > 1e-3]
                wire = BRepBuilderAPI_MakeWire(semi[0], te_seg[0]).Wire()
                caps.append(BRepBuilderAPI_MakeFace(wire, True).Face())
        else:
            keep.append(f)
    sew = BRepBuilderAPI_Sewing(1e-6)
    for f in keep + caps:
        sew.Add(f)
    sew.Perform()
    shells = A.shapes_of(sew.SewedShape(), TopAbs_SHELL)[0]
    solid = BRepBuilderAPI_MakeSolid(TopoDS.Shell_s(shells[0])).Solid()
    fix = ShapeFix_Solid(solid)
    fix.Perform()
    solid = BRepBuilderAPI_NurbsConvert(fix.Solid(), True).Shape()
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"flat-{variant}-revolved.step"
    Interface_Static.SetCVal_s("write.step.schema", "AP214")
    w = STEPControl_Writer()
    w.Transfer(solid, STEPControl_AsIs)
    if w.Write(str(path)) != IFSelect_RetDone:
        raise IOError(path)
    return case, path, bool(BRepCheck_Analyzer(solid).IsValid())


def profile(case):
    """Exact attachment profile of the flat parent: X(s), h(s) and derivatives (upper skin, Y = 350)."""
    m = case.parent
    def at(s):
        j = m.skin_jet(np.array([s]), np.array([350.0]), +1)
        return (j["point"][0, 0], j["point"][0, 2], j["d_chord"][0, 0], j["d_chord"][0, 2],
                j["d_chord_chord"][0, 0], j["d_chord_chord"][0, 2])
    return at


def expected(case, at, s):
    X, h, X1, h1, X2, h2 = at(s)
    rx = h1 / X1
    rxx = (h2 * X1 - h1 * X2) / X1**3
    k_par = 1 / (h * math.sqrt(1 + rx * rx))
    k_mer = -rxx / (1 + rx * rx) ** 1.5
    kmax = max(abs(k_par), abs(k_mer))
    t = 2 * h
    req = 0.4 if t >= 0.8 else 0.45 * t
    return dict(X=X, h=h, k_par=k_par, k_mer=k_mer, R=1 / kmax, score=(1 / kmax) / req, t=t)


def s_of_X(at, X, s_end):
    return brentq(lambda s: at(s)[0] - X, 1e-15, s_end, xtol=1e-15)


def model_checks(variant):
    case, path, ok = build_revolved(variant)
    report, text = A.audit("flat", variant, path)
    text += "\n" + A.section_table(report)
    (OUT / f"flat-{variant}-revolved.json").write_text(json.dumps(report, indent=1) + "\n")
    (OUT / f"flat-{variant}-revolved.txt").write_text(text + "\n")
    at = profile(case)
    s_end = 1.0
    X_end, h_end = at(s_end)[0], at(s_end)[1]
    out = {"model": str(path.relative_to(A.REPO)), "built_valid": ok}
    # Sections: every cut is a half circle of radius h(X).
    sec = []
    for tip, cuts in report["sections"].items():
        for c in cuts:
            if "error" in c:
                sec.append({"tip": tip, "fraction": c["fraction"], "error": c["error"]})
                continue
            X = c["fraction"] * 120.0
            h = at(s_of_X(at, X, s_end))[1] if X < X_end else h_end
            sec.append({"tip": tip, "fraction": c["fraction"], "expected_radius_mm": h,
                        "min_radius_mm": c["min_radius"]["radius_mm"], "lip_radius_mm": c["lip"]["radius_mm"],
                        "relative_error": c["min_radius"]["radius_mm"] / h - 1, "ratio_or_margin_expected": 1.0,
                        "ratio": c.get("ratio_to_half_thickness"), "attachment_thickness_vs_2h": c["attachment_thickness_mm"] / (2 * h) - 1,
                        "reversals": len(c["reversals"])})
    good = [r for r in sec if "error" not in r]
    out["sections"] = {"cuts": len(sec), "traced": len(good),
                       "max_abs_relative_radius_error": max(abs(r["relative_error"]) for r in good),
                       "max_abs_thickness_relative_error": max(abs(r["attachment_thickness_vs_2h"]) for r in good),
                       "reversals": sum(r["reversals"] for r in good), "rows": sec}
    # Principal curvatures over the revolved faces, point by point.
    model = A.Model(path, case)
    errs, n = [], 0
    for fr in model.cap_faces():
        if fr.planar:
            continue
        u0, u1, v0, v1 = fr.uv
        for u in np.linspace(u0, u1, 41)[1:-1]:
            for v in np.linspace(v0, v1, 41)[1:-1]:
                j = fr.jet(u, v)
                if j["N"] is None or not fr.inside(u, v):
                    continue
                X = j["P"][0]
                if not (1e-3 < X < X_end - 1e-3):
                    continue
                e = expected(case, at, s_of_X(at, X, s_end))
                exp = np.sort([e["k_par"], e["k_mer"]])
                errs.append(np.max(np.abs(j["k"] - exp) / np.maximum(np.abs(exp), 1e-3)))
                n += 1
    out["principal_curvature"] = {"points": n, "max_relative_error": float(max(errs))}
    # Map worst against the closed form at the reported point and over the profile.
    cm = report["cap_map"]
    w = cm["worst"]
    e = expected(case, at, s_of_X(at, w["xyz"][0], s_end))
    ss = np.concatenate([np.linspace(1e-4, 0.02, 400), np.linspace(0.02, 0.98, 2000), 1 - np.geomspace(0.02, 1e-6, 600)])
    ss = ss[ss < s_end]
    grid = [expected(case, at, s) for s in ss]
    excl = cm["collapse_exclusion_mm"]
    def far(g):
        d_le = math.hypot(g["X"], g["h"])
        d_te = math.hypot(X_end - g["X"], g["h"]) if variant == "sharp" else np.inf
        return min(d_le, d_te) >= excl
    best = min((g for g in grid if far(g)), key=lambda g: g["score"])
    from scipy.optimize import minimize_scalar
    k = int(np.argmin([g["score"] if far(g) else np.inf for g in grid]))
    lo, hi = ss[max(k - 1, 0)], ss[min(k + 1, len(ss) - 1)]
    rr = minimize_scalar(lambda s: expected(case, at, s)["score"], bounds=(lo, hi), method="bounded", options={"xatol": 1e-14})
    best = min(best, expected(case, at, float(rr.x)), key=lambda g: g["score"])
    rmin = min(g["R"] for g in grid if far(g))
    k = int(np.argmin([g["R"] if far(g) else np.inf for g in grid]))
    lo, hi = ss[max(k - 1, 0)], ss[min(k + 1, len(ss) - 1)]
    rr = minimize_scalar(lambda s: expected(case, at, s)["R"], bounds=(lo, hi), method="bounded", options={"xatol": 1e-14})
    rmin = min(rmin, expected(case, at, float(rr.x))["R"])
    out["cap_map"] = {"tool_worst_radius_mm": w["radius_mm"], "closed_form_radius_at_that_X_mm": e["R"],
                      "relative_error_at_point": w["radius_mm"] / e["R"] - 1,
                      "tool_worst_score": w["score_radius_over_required"], "closed_form_min_score": best["score"],
                      "closed_form_min_score_at_X_mm": best["X"],
                      "tool_tightest_radius_mm": cm["tightest"]["radius_mm"],
                      "closed_form_min_radius_mm": rmin}
    # Joins: the mismatch is the round-the-circle curvature.
    j = report["joins"]["summary"]["parent skin to cap"]
    out["joins"] = {"gap_mm": j["gap_mm"]["max"], "angle_deg": j["angle_deg"]["max"],
                    "tool_max_curvature_mismatch_per_mm": j["curvature_mismatch_per_mm"]["max"],
                    "closed_form_k_par_at_that_point": expected(case, at, s_of_X(at, j["curvature_mismatch_per_mm"]["at_xyz"][0], s_end))["k_par"],
                    "closed_form_max_k_par_beyond_exclusion": max(g["k_par"] for g in grid if far(g))}
    # Outline: the lip sits at 350 + h(X); the locked outline at 350 + 200 p (1 - p) where X = 120 x(p).
    xh = lambda p: p**3 * (10 - 15 * p + 6 * p * p)
    fr_ = np.concatenate([[0.005, 0.01], np.linspace(0.02, 0.98, 97), [0.99, 0.995]])
    devs = []
    for f in fr_:
        X = 120 * f
        p = brentq(lambda p: xh(p) - f, 0, 1, xtol=1e-15)
        h = at(s_of_X(at, X, s_end))[1]
        devs.append(abs(h - 200 * p * (1 - p)))
    scan = report["outline"]["lip_scan"]
    out["outline"] = {"tool_max_deviation_mm": max(s["max_deviation_mm"] for s in scan.values()),
                      "closed_form_max_deviation_mm": max(devs)}
    # Volume and extent.
    vol_add = math.pi * quad(lambda s: at(s)[1] ** 2 * at(s)[2], 0, s_end, epsabs=1e-12, epsrel=1e-13, limit=400)[0]
    v = report["validity"]
    out["volume"] = {"tool_cap_volume_mm3": v["cap_volume_mm3"], "closed_form_cap_volume_mm3": vol_add,
                     "relative_error": v["cap_volume_mm3"] / vol_add - 1}
    from scipy.optimize import minimize_scalar
    s0 = max(np.linspace(0, s_end, 20001), key=lambda s: at(s)[1])
    hmax = -minimize_scalar(lambda s: -at(s)[1], bounds=(max(0, s0 - 1e-4), min(s_end, s0 + 1e-4)), method="bounded",
                            options={"xatol": 1e-15}).fun
    out["extent"] = {"tool_mm": v["extent"]["measured_mm"], "closed_form_mm": 700 + 2 * hmax}
    out["te_closure"] = report["te_closure"]["tips"][0] | {"approach_sections": [
        {k: a[k] for k in ("distance_to_te_mm", "attachment_thickness_mm", "min_radius_mm", "ratio_to_half_thickness")}
        for a in report["te_closure"]["tips"][0]["approach_sections"]]}
    return out, text


def main():
    result = {"primitives": primitive_checks()}
    for variant in ("sharp", "blunt"):
        result[f"flat-{variant}-revolved"], text = model_checks(variant)
        print(text)
    path = OUT / "synthetic_check.json"
    path.write_text(json.dumps(A.clean(result), indent=1) + "\n")
    print(json.dumps(A.clean({k: (v if k == "primitives" else {kk: vv for kk, vv in v.items() if kk not in ("sections",)})
                              for k, v in result.items()}), indent=1))


if __name__ == "__main__":
    main()
