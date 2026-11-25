"""Build wing.step: the swept55 / sharp parent with a rounded cap at each tip.

Run from the project root:
  .venv/bin/python runs/cap-radius/builds/swept55/sharp/build_cap.py

Construction (one B-spline face per tip, through the lip):
  u = the parent's chord parameter q, the parent's own u knots, degree raised
      5 -> 7 (exact) so the cross-boundary derivative D(u) * P_v(u) is exact.
  w in [0, 1]: w = 0 upper cut edge, w = 1/2 lip (locked outline), w = 1 lower
      cut edge; theta = theta_of_w(w, q) in cap_model (pi w near the edges and
      the lip, shifted in between so the start of the roll lies on a knot line).
  Fixed poles: rows 0 / N-1 = the parent's own (degree-raised) cut-edge poles
      (C0, exact); rows 1 / N-2 from the parent's span derivative times
      pi D(u) (exact tangent plane); rows 2 / N-3 make the second w
      derivative equal pi^2 D^2 P_vv (curvature match, exact for this STEP);
      columns 0 and last = the leading- and trailing-edge points (collapsed).
  Free poles: separable least squares to the analytic cap (cap_model.py,
      knob table design_table.json from design.py).
The left cap is the mirror image of the right cap in Y. The parent's skin
faces are used unchanged from parent.step; its two flat end faces are dropped.
"""
import hashlib
import json
import pathlib
import sys
import warnings

import numpy as np
from scipy.interpolate import BSpline

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cap_model as cm  # noqa: E402

from OCP.BRep import BRep_Tool  # noqa: E402
from OCP.BRepAdaptor import BRepAdaptor_Surface  # noqa: E402
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeSolid,  # noqa: E402
                                BRepBuilderAPI_Sewing)
from OCP.BRepCheck import BRepCheck_Analyzer  # noqa: E402
from OCP.BRepGProp import BRepGProp  # noqa: E402
from OCP.GProp import GProp_GProps  # noqa: E402
from OCP.Geom import Geom_BSplineSurface  # noqa: E402
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Plane  # noqa: E402
from OCP.IFSelect import IFSelect_RetDone  # noqa: E402
from OCP.Interface import Interface_Static  # noqa: E402
from OCP.STEPControl import STEPControl_AsIs, STEPControl_Reader, STEPControl_Writer  # noqa: E402
from OCP.ShapeFix import ShapeFix_Solid  # noqa: E402
from OCP.TColgp import TColgp_Array2OfPnt  # noqa: E402
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal  # noqa: E402
from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL  # noqa: E402
from OCP.TopExp import TopExp_Explorer  # noqa: E402
from OCP.TopoDS import TopoDS  # noqa: E402
from OCP.gp import gp_Pnt  # noqa: E402

warnings.simplefilter("ignore")
PARENT_STEP = cm.PARENT_DIR / "parent.step"
OUT = HERE / "wing.step"
DEG_U, DEG_W = 7, 5


# ----------------------------------------------------------------- parent
def read_step(path):
    r = STEPControl_Reader()
    if r.ReadFile(str(path)) != IFSelect_RetDone:
        raise RuntimeError(f"cannot read {path}")
    r.TransferRoots()
    return r.OneShape()


def faces_of(shape):
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    out = []
    while ex.More():
        out.append(TopoDS.Face_s(ex.Current()))
        ex.Next()
    return out


def bspline_data(bs):
    nu, nv = bs.NbUPoles(), bs.NbVPoles()
    P = np.array([[[bs.Pole(i, j).X(), bs.Pole(i, j).Y(), bs.Pole(i, j).Z()] for j in range(1, nv + 1)]
                  for i in range(1, nu + 1)])
    uk = [bs.UKnot(i) for i in range(1, bs.NbUKnots() + 1)]
    um = [bs.UMultiplicity(i) for i in range(1, bs.NbUKnots() + 1)]
    vk = [bs.VKnot(i) for i in range(1, bs.NbVKnots() + 1)]
    vm = [bs.VMultiplicity(i) for i in range(1, bs.NbVKnots() + 1)]
    return P, uk, um, vk, vm, bs.UDegree(), bs.VDegree()


def make_bspline(P, uk, um, vk, vm, du, dv):
    nu, nv = P.shape[:2]
    poles = TColgp_Array2OfPnt(1, nu, 1, nv)
    for i in range(nu):
        for j in range(nv):
            poles.SetValue(i + 1, j + 1, gp_Pnt(*map(float, P[i, j])))

    def arr(vals, cls):
        a = cls(1, len(vals))
        for k, x in enumerate(vals):
            a.SetValue(k + 1, x)
        return a
    return Geom_BSplineSurface(poles, arr(uk, TColStd_Array1OfReal), arr(vk, TColStd_Array1OfReal),
                               arr(um, TColStd_Array1OfInteger), arr(vm, TColStd_Array1OfInteger), du, dv)


def full_knots(k, m):
    return np.repeat(np.asarray(k, float), m)


def parent_skins():
    shape = read_step(PARENT_STEP)
    skins, planes = [], []
    for f in faces_of(shape):
        a = BRepAdaptor_Surface(f)
        if a.GetType() == GeomAbs_BSplineSurface:
            bs = a.BSpline()
            p = bs.Value(0.5, 390.0)
            skins.append({"face": f, "bs": bs, "Y": p.Y(), "Z": p.Z()})
        elif a.GetType() == GeomAbs_Plane:
            planes.append(f)
    return shape, skins, planes


# ----------------------------------------------------------------- cap fit
def w_knots():
    """Knots in w, mirrored about the lip: coarse near the cut edges, dense over the roll,
    a C2 knot (multiplicity 3) on the roll-start line w = W_ON, C3 knots at the ends of
    the parameter shift (W_A, W_C)."""
    half = np.unique(np.concatenate([np.linspace(0, cm.W_A, 14), np.linspace(cm.W_A, cm.W_ON, 19),
                                     np.linspace(cm.W_ON, 0.5, 32)]))
    k = np.unique(np.round(np.concatenate([half, 1 - half]), 12))
    mult = np.ones(len(k), int)
    for target, m in ((cm.W_A, 2), (cm.W_ON, 3), (cm.W_C, 2)):
        for w in (target, 1 - target):
            i = int(np.argmin(np.abs(k - w)))
            if abs(k[i] - w) > 1e-9:
                k = np.insert(k, np.searchsorted(k, w), w)
                mult = np.insert(mult, np.searchsorted(k, w) - 0, 1)
                i = int(np.argmin(np.abs(k - w)))
            mult[i] = m
    mult[0] = mult[-1] = DEG_W + 1
    return list(k), list(mult)


def boundary_rows(up, lo):
    """Rows 0,1,2 (upper) and the mirrored lower rows from the parent's degree-raised surfaces."""
    rows = {}
    for key, bs0 in (("up", up), ("lo", lo)):
        P, uk, um, vk, vm, du, dv = bspline_data(bs0)
        s = make_bspline(P, uk, um, vk, vm, du, dv)
        s.IncreaseDegree(DEG_U, dv)
        P7, uk7, um7, vk7, vm7, du7, dv7 = bspline_data(s)
        assert abs(vk7[0]) < 1e-12 and abs(vk7[-1] - 390.0) < 1e-12 and dv7 == 4 and len(vk7) == 2
        edge = P7[:, -1]
        Pv = 4 * (P7[:, -1] - P7[:, -2]) / 390.0                    # d/dv poles at v = 390
        Pvv = 12 * (P7[:, -1] - 2 * P7[:, -2] + P7[:, -3]) / 390.0**2
        rows[key] = (edge, Pv, Pvv)
        rows["uknots"], rows["umults"] = uk7, um7
    return rows


def spline_product_poles(tu, poles, factor_fn):
    """Poles (same degree-7 space) of factor(u) * spline(u); exact when the product lies in the space."""
    us = np.unique(np.concatenate([np.linspace(0, 1, 2001), tu]))
    Bm = BSpline.design_matrix(us, tu, DEG_U).toarray()
    vals = (Bm @ poles) * factor_fn(us)[:, None]
    sol, *_ = np.linalg.lstsq(Bm, vals, rcond=None)
    return sol, float(np.abs(Bm @ sol - vals).max())


def tensor_eval(A, P, Bm):
    """sum_ij A[a,i] P[i,j,c] B[b,j] with plain matrix products."""
    return np.stack([A @ P[..., c] @ Bm.T for c in range(P.shape[-1])], axis=-1)


def build_right_cap(rows, report):
    uk, um = rows["uknots"], rows["umults"]
    tu = full_knots(uk, um)
    nu = len(tu) - DEG_U - 1
    wk, wm = w_knots()
    tw = full_knots(wk, wm)
    nw = len(tw) - DEG_W - 1
    D = lambda u: 240.0 * u * (1 - u)
    P = np.zeros((nu, nw, 3))
    fixed = np.zeros((nu, nw), bool)
    k = DEG_W
    a1, a2 = tw[k + 1], tw[k + 2]                     # first interior knots from w = 0
    b1, b2 = 1 - tw[-k - 2], 1 - tw[-k - 3]           # from w = 1 (mirror)
    for key, (j0, j1, j2, sgn, h1, h2) in {"up": (0, 1, 2, +1, a1, a2), "lo": (nw - 1, nw - 2, nw - 3, -1, b1, b2)}.items():
        edge, Pv, Pvv = rows[key]
        Sw, e1 = spline_product_poles(tu, Pv, lambda u: np.pi * D(u))          # dS/dw (into the cap)
        Sww, e2 = spline_product_poles(tu, Pvv, lambda u: (np.pi * D(u)) ** 2)  # second derivative
        report[f"product_fit_max_{key}"] = [e1, e2]
        # w runs 0 -> 1; at w = 1 the derivative into the cap is -d/dw
        P[:, j0] = edge
        P[:, j1] = edge + h1 / k * Sw
        # S''(end) = (k-1)/h1 * (Q1 - Q0), Q0 = k (P1-P0)/h1, Q1 = k (P2-P1)/h2  ->  P2 = P1 + h2 (Q0 + h1/(k-1) S'') / k
        Q0 = Sw
        Q1 = Q0 + h1 / (k - 1) * Sww
        P[:, j2] = P[:, j1] + h2 / k * Q1
        fixed[:, [j0, j1, j2]] = True
    le = cm.parent.skin(0.0, 390.0, +1)
    te = cm.parent.skin(1.0, 390.0, +1)
    report["le_point_vs_step_edge"] = float(np.abs(P[0, 0] - le).max())
    report["te_point_vs_step_edge"] = float(np.abs(P[-1, 0] - te).max())
    P[0, :] = P[0, 0]
    P[-1, :] = P[-1, 0]
    fixed[0, :] = fixed[-1, :] = True
    # samples
    qs = np.unique(np.concatenate([np.linspace(0, 1, 1601), tu[(tu > 0) & (tu < 1)],
                                   np.geomspace(1e-6, 0.02, 200), 1 - np.geomspace(1e-6, 0.02, 200)]))
    ws = np.unique(np.concatenate([np.linspace(0, 1, 2401), tw]))
    Q, Wg = np.meshgrid(qs, ws, indexing="ij")
    Th = cm.theta_of_w(Wg, Q)
    F = np.empty(Q.shape + (3,))
    for i0 in range(0, len(qs), 100):
        F[i0:i0 + 100] = cm.point(Q[i0:i0 + 100], Th[i0:i0 + 100])
    Bu = BSpline.design_matrix(qs, tu, DEG_U).toarray()
    Bw = BSpline.design_matrix(ws, tw, DEG_W).toarray()
    fu = np.arange(1, nu - 1)
    fw = np.arange(3, nw - 3)
    P0 = np.where(fixed[..., None], P, 0.0)
    resid = F - tensor_eval(Bu, P0, Bw)
    # rows weighted by 1 / (local thickness + 0.01 mm): the trailing-edge rows carry
    # sub-micron features and need the same relative accuracy as the thick rows
    wt = 1.0 / (cm.thickness_at_attachment(qs) + 0.01)
    Au = np.linalg.pinv(wt[:, None] * Bu[:, fu]) * wt[None, :]
    Aw = np.linalg.pinv(Bw[:, fw])
    free = tensor_eval(Au, resid, Aw)
    P[np.ix_(fu, fw)] = free
    fit = tensor_eval(Bu, P, Bw)
    err = np.linalg.norm(fit - F, axis=-1)
    rel = err / (cm.thickness_at_attachment(qs)[:, None] + 1e-6)
    report["fit_max_error_over_local_thickness"] = float(rel.max())
    i, j = np.unravel_index(np.argmax(err), err.shape)
    report["fit_max_error_mm"] = float(err.max())
    report["fit_max_error_at_q_w"] = [float(qs[i]), float(ws[j])]
    report["fit_rms_error_mm"] = float(np.sqrt((err**2).mean()))
    report["cap_poles_u_w"] = [int(nu), int(nw)]
    return P, uk, um, wk, wm


def face_from(P, uk, um, wk, wm):
    s = make_bspline(P, uk, um, wk, wm, DEG_U, DEG_W)
    mf = BRepBuilderAPI_MakeFace(s, 1e-7)
    if not mf.IsDone():
        raise RuntimeError("cap face failed")
    return mf.Face(), s


def main():
    report = {}
    shape, skins, planes = parent_skins()
    right = [s for s in skins if s["Y"] > 0]
    up = next(s for s in right if s["Z"] > 0)["bs"]
    lo = next(s for s in right if s["Z"] < 0)["bs"]
    rows = boundary_rows(up, lo)
    report["parent_Pvv_max_mm_per_mm2"] = float(max(np.abs(rows["up"][2]).max(), np.abs(rows["lo"][2]).max()))
    P, uk, um, wk, wm = build_right_cap(rows, report)
    Pl = P.copy()
    Pl[..., 1] *= -1
    cap_r, _ = face_from(P, uk, um, wk, wm)
    cap_l, _ = face_from(Pl, uk, um, wk, wm)
    sew = BRepBuilderAPI_Sewing(1e-6)
    for s in skins:
        sew.Add(s["face"])
    sew.Add(cap_r)
    sew.Add(cap_l)
    sew.Perform()
    sewn = sew.SewedShape()
    report["sewing_free_edges"] = sew.NbFreeEdges()
    report["sewing_multiple_edges"] = sew.NbMultipleEdges()
    report["sewing_degenerated_shapes"] = sew.NbDegeneratedShapes()
    ex = TopExp_Explorer(sewn, TopAbs_SHELL)
    shell = TopoDS.Shell_s(ex.Current())
    ms = BRepBuilderAPI_MakeSolid(shell)
    fix = ShapeFix_Solid(ms.Solid())
    fix.Perform()
    solid = fix.Solid()
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(solid, props, 1e-9)
    report["volume_mm3"] = props.Mass()
    report["brepcheck_valid"] = BRepCheck_Analyzer(solid).IsValid()
    Interface_Static.SetCVal_s("write.step.schema", "AP214IS")
    Interface_Static.SetCVal_s("write.step.unit", "MM")
    w = STEPControl_Writer()
    w.Transfer(solid, STEPControl_AsIs)
    if w.Write(str(OUT)) != IFSelect_RetDone:
        raise RuntimeError("STEP write failed")
    data = OUT.read_bytes()
    report["step_bytes"] = len(data)
    report["step_sha256"] = hashlib.sha256(data).hexdigest()
    (HERE / "build_report.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
