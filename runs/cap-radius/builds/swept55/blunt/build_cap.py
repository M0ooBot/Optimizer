"""Build the capped swept55 blunt wing and export wing.step.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/swept55/blunt/build_cap.py

Steps
  1. sample the analytic cap (cap_model.surface) on a dense (U, v) grid of the right-tip upper patch;
  2. fit one quintic B-spline patch with exact boundary structure:
       v = 0 row  -> the parent's upper cut edge at Y = 390 (fitted, sewn to the parent edge),
       v = 1 row  -> the lip, every pole at Z = 0 and the next row straight above it
                     (vertical tangent plane -> upper and mirrored lower patch share it: G1),
       U = 0      -> collapsed onto the leading-edge point,
       U = 1      -> quarter circle of radius 0.2 in the trailing-edge plane (all poles in it);
  3. mirror it to the lower patch (Z -> -Z) and close the end of the flat TE face with the
     planar half-round face bounded by the two quarter circles and the TE segment;
  4. mirror the tip to Y < 0, drop the parent's two flat end faces, sew, make the solid, export.
"""
import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cap_model as cm  # noqa: E402
from bspline_fit import basis_matrix, clamped_knots, constrained_lstsq  # noqa: E402

from OCP.BRep import BRep_Tool  # noqa: E402
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace,  # noqa: E402
                                BRepBuilderAPI_MakeSolid, BRepBuilderAPI_MakeWire,
                                BRepBuilderAPI_Sewing, BRepBuilderAPI_Transform)
from OCP.BRepCheck import BRepCheck_Analyzer  # noqa: E402
from OCP.BRepGProp import BRepGProp  # noqa: E402
from OCP.GProp import GProp_GProps  # noqa: E402
from OCP.Geom import Geom_BSplineSurface, Geom_Plane  # noqa: E402
from OCP.GeomAbs import GeomAbs_Plane  # noqa: E402
from OCP.BRepAdaptor import BRepAdaptor_Surface  # noqa: E402
from OCP.IFSelect import IFSelect_RetDone  # noqa: E402
from OCP.Interface import Interface_Static  # noqa: E402
from OCP.STEPControl import STEPControl_AsIs, STEPControl_Reader, STEPControl_Writer  # noqa: E402
from OCP.ShapeFix import ShapeFix_Shape, ShapeFix_Solid  # noqa: E402
from OCP.TColgp import TColgp_Array2OfPnt  # noqa: E402
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal  # noqa: E402
from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL  # noqa: E402
from OCP.TopExp import TopExp_Explorer  # noqa: E402
from OCP.TopoDS import TopoDS  # noqa: E402
from OCP.gp import gp_Ax2, gp_Pln, gp_Pnt, gp_Dir, gp_Trsf, gp_Ax1  # noqa: E402

DEG_U, DEG_V = 5, 5
N_U, N_V = 200, 40            # poles
SEW_TOL = 1e-5                # mm

# ---------------------------------------------------------------------------- parametrization
P0_LE, Q0_TE = 0.003, 1e-3    # grading scales at the LE corner (p) and the TE apex (P_APEX - p)


def p_of_U(U):
    """Lip parameter p as a logistic function of the patch parameter U: graded towards both ends."""
    Pa = cm.P_APEX
    L0 = np.log(P0_LE / (Pa + Q0_TE))
    L1 = np.log((Pa + P0_LE) / Q0_TE)
    L = L0 + np.asarray(U, float) * (L1 - L0)
    e = np.exp(L)
    return np.clip((e * (Pa + Q0_TE) - P0_LE) / (1 + e), 0.0, Pa)


def surface_U(U, v):
    return cm.surface(p_of_U(U) / cm.P_APEX, v)


def u_grade(s, alpha=0.85):
    """Map uniform s in [0, 1] to U, concentrating at both ends (spacing x (1 - alpha) there)."""
    s = np.asarray(s, float)
    return s - alpha * np.sin(2 * np.pi * s) / (2 * np.pi)


def v_breaks(n_int):
    """Interior v knots, denser towards the lip."""
    s = np.linspace(0, 1, n_int + 2)[1:-1]
    return 1 - (1 - s) ** 1.6


def start_derivative_rows(k, deg, n):
    """Rows r1, r2 with S'(0) = r1 . P and S''(0) = r2 . P for a clamped B-spline."""
    a = deg / (k[deg + 1] - k[1])
    b = (deg - 1) / (k[deg + 2] - k[2])
    r1 = np.zeros(n); r1[0], r1[1] = -a, a
    # S''(0) = (deg - 1)/(k[deg+1]-k[2]) * (Q1 - Q0), Q_i = first-derivative poles
    a1 = deg / (k[deg + 2] - k[2])
    c = (deg - 1) / (k[deg + 1] - k[2])
    r2 = np.zeros(n)
    r2[0] += c * a; r2[1] += -c * a - c * a1; r2[2] += c * a1
    return r1, r2


def end_derivative_rows(k, deg, n):
    """Rows with S'(1) = r1 . P and S''(1) = r2 . P (mirror of start_derivative_rows)."""
    kr = 1 - k[::-1]
    r1, r2 = start_derivative_rows(kr, deg, n)
    return -r1[::-1], r2[::-1]


def fit_patch(nu_s=1500, nv_s=220):
    Us = u_grade(np.linspace(0, 1, nu_s), 0.5)
    vs = 1 - (1 - np.linspace(0, 1, nv_s)) ** 1.5
    UU, VV = np.meshgrid(Us, vs, indexing="ij")
    with np.errstate(all="ignore"):
        D = surface_U(UU, VV)                      # (nu_s, nv_s, 3)
    D[0] = [cm.SWEEP * cm.ATT, cm.ATT, 0.0]          # collapsed LE corner exactly
    bad = ~np.isfinite(D).all(-1)
    if bad.any():
        raise RuntimeError(f"non-finite samples at U = {np.unique(UU[bad])[:10]}")
    D[:, -1, 2] = 0.0
    # cross derivatives at the attachment (v = 0): the cap there is the parent's own surface,
    # so match its first and second v-derivatives exactly (G1 and G2 to the parent skin)
    ku = clamped_knots(N_U, DEG_U)
    kv = clamped_knots(N_V, DEG_V, v_breaks(N_V - DEG_V - 1))
    Av = basis_matrix(vs, kv, DEG_V)
    Au = basis_matrix(Us, ku, DEG_U)
    Pv = v_fit(Us, vs, D, Av, kv)
    # U = 1 end (the half-round guide): match the analytic first and second U-derivatives there,
    # taken from v-fits of five rows next to it (backward 5-point stencils on the v-poles)
    hu = 2e-4
    Ue = 1 - hu * np.arange(5)
    with np.errstate(all="ignore"):
        De = surface_U(*np.meshgrid(Ue, vs, indexing="ij"))
    De[:, -1, 2] = 0.0
    Pe = v_fit(Ue, vs, De, Av, kv)
    dU1 = (25 * Pe[0] - 48 * Pe[1] + 36 * Pe[2] - 16 * Pe[3] + 3 * Pe[4]) / (12 * hu)
    dU2 = (35 * Pe[0] - 104 * Pe[1] + 114 * Pe[2] - 56 * Pe[3] + 11 * Pe[4]) / (12 * hu * hu)
    r1, r2 = end_derivative_rows(ku, DEG_U, N_U)
    # (tried: also constraining dU1, dU2 there -> the apex corner made the lip ring; not used)
    Cu = np.zeros((2, N_U)); Cu[0, 0] = 1; Cu[1, -1] = 1
    P = np.zeros((N_U, N_V, 3))
    for c in range(3):
        Y = Pv[:, :, c]
        dd = np.vstack([Y[0], Pe[0][:, c]])
        P[:, :, c] = constrained_lstsq(Au, Y, Cu, dd)
    return P, ku, kv, (Us, vs, D)


def v_fit(Us, vs, D, Av, kv):
    """Per U sample: v-direction fit with first pole = attachment point, G1 + G2 to the parent
    (analytic cross derivatives), lip row at z = 0 with the next row straight above it."""
    nu_s = len(Us)
    hv = 1e-4
    with np.errstate(all="ignore"):
        S5 = np.stack([surface_U(Us, np.full_like(Us, j * hv)) for j in range(5)], 0)
    S5[:, Us == 0] = D[0, 0]
    Dv0 = (-25 * S5[0] + 48 * S5[1] - 36 * S5[2] + 16 * S5[3] - 3 * S5[4]) / (12 * hv)
    Dvv0 = (35 * S5[0] - 104 * S5[1] + 114 * S5[2] - 56 * S5[3] + 11 * S5[4]) / (12 * hv * hv)
    nv = N_V
    d1, d2 = start_derivative_rows(kv, DEG_V, nv)
    Cxy = np.zeros((4, nv)); Cxy[0, 0] = 1; Cxy[1] = d1; Cxy[2] = d2; Cxy[3, -1] = 1; Cxy[3, -2] = -1
    Cz = np.zeros((4, nv)); Cz[0, 0] = 1; Cz[1] = d1; Cz[2] = d2; Cz[3, -1] = 1
    Pv = np.zeros((nu_s, nv, 3))
    for c in range(3):
        Y = D[:, :, c].T
        C = Cz if c == 2 else Cxy
        dd = np.vstack([D[:, 0, c], Dv0[:, c], Dvv0[:, c], np.zeros(nu_s)])
        Pv[:, :, c] = constrained_lstsq(Av, Y, C, dd).T
    return Pv


def eval_patch(P, ku, kv, U, v):
    Au = basis_matrix(np.atleast_1d(U), ku, DEG_U)
    Av = basis_matrix(np.atleast_1d(v), kv, DEG_V)
    return np.einsum("ia,abc,jb->ijc", Au, P, Av)


def to_occ_surface(P, ku, kv):
    nu, nv, _ = P.shape
    poles = TColgp_Array2OfPnt(1, nu, 1, nv)
    for i in range(nu):
        for j in range(nv):
            poles.SetValue(i + 1, j + 1, gp_Pnt(*map(float, P[i, j])))

    def knots_mults(k):
        vals, mults = np.unique(np.round(k, 15), return_counts=True)
        K = TColStd_Array1OfReal(1, len(vals))
        M = TColStd_Array1OfInteger(1, len(vals))
        for i, (a, m) in enumerate(zip(vals, mults)):
            K.SetValue(i + 1, float(a)); M.SetValue(i + 1, int(m))
        return K, M

    Ku, Mu = knots_mults(ku)
    Kv, Mv = knots_mults(kv)
    return Geom_BSplineSurface(poles, Ku, Kv, Mu, Mv, DEG_U, DEG_V)


def read_step(path):
    r = STEPControl_Reader()
    assert r.ReadFile(path) == IFSelect_RetDone
    r.TransferRoots()
    return r.OneShape()


def faces_of(shape):
    out = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        out.append(TopoDS.Face_s(ex.Current()))
        ex.Next()
    return out


def mirrored(shape, plane_normal):
    t = gp_Trsf()
    t.SetMirror(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(*plane_normal)))
    return BRepBuilderAPI_Transform(shape, t, True).Shape()


def build():
    P, ku, kv, (Us, vs, D) = fit_patch()
    # fit error on a validation grid (between the samples)
    Uc = np.r_[np.linspace(0.0005, 0.9995, 700), 1 - np.geomspace(1e-6, 1e-3, 60)]
    vc = np.r_[np.linspace(0.0005, 0.9995, 211), 1 - np.geomspace(1e-6, 1e-3, 20)]
    with np.errstate(all="ignore"):
        Dc = surface_U(*np.meshgrid(Uc, vc, indexing="ij"))
    Fc = eval_patch(P, ku, kv, Uc, vc)
    err = np.linalg.norm(Fc - Dc, axis=-1)
    i = np.unravel_index(np.nanargmax(err), err.shape)
    fit = {"max_mm": float(np.nanmax(err)), "rms_mm": float(np.sqrt(np.nanmean(err**2))),
           "at_U": float(Uc[i[0]]), "at_v": float(vc[i[1]]), "at_xyz": Fc[i].round(5).tolist()}
    print("fit", fit)

    srf = to_occ_surface(P, ku, kv)
    up = BRepBuilderAPI_MakeFace(srf, 1e-7).Face()
    lo = mirrored(up, (0, 0, 1))
    # half-round face in the TE plane: the TE segment plus the two quarter circles (U = 1 isolines)
    from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf  # noqa: F401
    c_up = srf.UIso(1.0)
    e_up = BRepBuilderAPI_MakeEdge(c_up).Edge()
    e_lo = mirrored(e_up, (0, 0, 1))
    pA = gp_Pnt(*map(float, P[-1, 0]))
    pB = gp_Pnt(pA.X(), pA.Y(), -pA.Z())
    e_seg = BRepBuilderAPI_MakeEdge(pB, pA).Edge()
    w = BRepBuilderAPI_MakeWire()
    for e in (e_up, TopoDS.Edge_s(e_lo), e_seg):
        w.Add(e)
    n_te = np.r_[cm.E_B, 0.0]
    pln = gp_Pln(pA, gp_Dir(*map(float, n_te)))
    disc = BRepBuilderAPI_MakeFace(pln, w.Wire(), True).Face()
    tip = [up, TopoDS.Face_s(lo), disc]
    tip_l = [TopoDS.Face_s(mirrored(f, (0, 1, 0))) for f in tip]

    parent = read_step(os.path.join(cm.PARENT_DIR, "parent.step"))
    keep = []
    for f in faces_of(parent):
        if BRepAdaptor_Surface(f).GetType() == GeomAbs_Plane:
            continue                                   # the two flat end faces at the attachments
        keep.append(f)
    sew = BRepBuilderAPI_Sewing(SEW_TOL)
    for f in keep + tip + tip_l:
        sew.Add(f)
    sew.Perform()
    sewn = sew.SewedShape()
    ex = TopExp_Explorer(sewn, TopAbs_SHELL)
    shell = TopoDS.Shell_s(ex.Current())
    ms = BRepBuilderAPI_MakeSolid(shell)
    solid = ms.Solid()
    sf = ShapeFix_Solid(solid)
    sf.Perform()
    solid = sf.Solid()
    ok = BRepCheck_Analyzer(solid).IsValid()
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(solid, props, 1e-9)
    info = {"fit": fit, "valid": bool(ok), "volume_mm3": props.Mass(), "faces": len(faces_of(solid)),
            "free_edges": sew.NbFreeEdges(), "multiple_edges": sew.NbMultipleEdges(),
            "degenerated_shapes": sew.NbDegeneratedShapes(),
            "poles": [N_U, N_V], "degree": [DEG_U, DEG_V], "sew_tol_mm": SEW_TOL}
    out = os.path.join(HERE, "wing.step")
    Interface_Static.SetCVal_s("write.step.schema", "AP214IS")
    Interface_Static.SetCVal_s("write.step.unit", "MM")
    wr = STEPControl_Writer()
    wr.Transfer(solid, STEPControl_AsIs)
    assert wr.Write(out) == IFSelect_RetDone
    info["step_bytes"] = os.path.getsize(out)
    info["step_sha256"] = hashlib.sha256(open(out, "rb").read()).hexdigest()
    np.savez_compressed(os.path.join(HERE, "cap_patch.npz"), poles=P, ku=ku, kv=kv)
    info["design"] = {k: (list(v) if isinstance(v, tuple) else v) for k, v in cm.DESIGN.items()}
    json.dump(info, open(os.path.join(HERE, "build_report.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in info.items() if k != "design"}, indent=1))
    return solid


if __name__ == "__main__":
    build()
