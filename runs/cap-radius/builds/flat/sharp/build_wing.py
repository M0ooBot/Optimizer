"""Build runs/cap-radius/builds/flat/sharp/wing.step: the parent wing with a rounded cap at each tip.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/flat/sharp/build_wing.py

Construction (see README.md):
  * The parent skins are taken unchanged from parent.step (two B-spline faces,
    u = native chord parameter, degree 5 with 29 knots, 33 poles).
  * Each cap is ONE B-spline face. Its u direction is the parent's u basis, so
    its first and last pole rows are exactly the parent's cut-edge poles at
    Y = +/-350 (upper edge, then lower edge). Its v direction wraps upper edge
    -> lip -> lower edge on the 33-pole degree-5 basis of wrap_design.py.
  * Cap poles: P[j, i] = (X_j, +/-(350 + r_j eta_i), h_j zeta_i), where X_j,
    h_j are the parent edge poles (X and upper Z), r_j are the B-spline
    coefficients of the locked outline reach r(u) = 200 u (1 - u) on the same u
    basis, and eta, zeta are the wrap design poles (eta^5 + zeta^2 = 1). Every wrap guide (iso-u) lies in
    the plane X = X(u), the lip iso-curve (v = 0.5) is exactly the locked outline,
    and the first three wrap rows differ from the edge only in Y (exact G1 + G2
    to the parent skin, whose spanwise curvature is zero).
  * The four faces are sewn into one shell and closed into one solid.
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
PARENT_DIR = ROOT / "runs/cap-radius/parents/flat/sharp"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PARENT_DIR))

import parent  # noqa: E402
import wrap_design  # noqa: E402

from OCP.BRep import BRep_Tool  # noqa: E402
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeSolid, BRepBuilderAPI_Sewing  # noqa: E402
from OCP.BRepCheck import BRepCheck_Analyzer  # noqa: E402
from OCP.BRepGProp import BRepGProp  # noqa: E402
from OCP.GProp import GProp_GProps  # noqa: E402
from OCP.Geom import Geom_BSplineSurface  # noqa: E402
from OCP.GeomAdaptor import GeomAdaptor_Surface  # noqa: E402
from OCP.GeomAbs import GeomAbs_BSplineSurface  # noqa: E402
from OCP.IFSelect import IFSelect_RetDone  # noqa: E402
from OCP.Interface import Interface_Static  # noqa: E402
from OCP.STEPControl import STEPControl_AsIs, STEPControl_Reader, STEPControl_Writer  # noqa: E402
from OCP.ShapeFix import ShapeFix_Solid  # noqa: E402
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal  # noqa: E402
from OCP.TColgp import TColgp_Array2OfPnt  # noqa: E402
from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL  # noqa: E402
from OCP.TopExp import TopExp_Explorer  # noqa: E402
from OCP.TopoDS import TopoDS  # noqa: E402
from OCP.gp import gp_Pnt  # noqa: E402

SEW_TOL = 1e-6


def read_step(path):
    reader = STEPControl_Reader()
    if reader.ReadFile(str(path)) != IFSelect_RetDone:
        raise RuntimeError(f"cannot read {path}")
    reader.TransferRoots()
    return reader.OneShape()


def faces(shape):
    out = []
    exp = TopExp_Explorer(shape, TopAbs_FACE)
    while exp.More():
        out.append(TopoDS.Face_s(exp.Current()))
        exp.Next()
    return out


def parent_skins(shape):
    """The parent's upper and lower B-spline skin faces and their surfaces."""
    skins = {}
    for f in faces(shape):
        srf = BRep_Tool.Surface_s(f)
        if GeomAdaptor_Surface(srf).GetType() != GeomAbs_BSplineSurface:
            continue
        bs = GeomAdaptor_Surface(srf).BSpline()
        z = bs.Value(0.5, 0.0).Z()
        skins["upper" if z > 0 else "lower"] = (f, bs)
    return skins


def knots_of(bs, direction):
    if direction == "u":
        n, k, m, deg = bs.NbUKnots(), bs.UKnot, bs.UMultiplicity, bs.UDegree()
    else:
        n, k, m, deg = bs.NbVKnots(), bs.VKnot, bs.VMultiplicity, bs.VDegree()
    kn = [k(i) for i in range(1, n + 1)]
    mu = [m(i) for i in range(1, n + 1)]
    return kn, mu, deg


def flat_knots(kn, mu):
    return np.repeat(kn, mu)


def edge_poles(bs, side_index):
    """Poles of the iso-V curve at the first (0) or last (-1) V pole column (clamped: exact)."""
    nv = bs.NbVPoles()
    col = 1 if side_index == 0 else nv
    return np.array([[bs.Pole(i, col).X(), bs.Pole(i, col).Y(), bs.Pole(i, col).Z()]
                     for i in range(1, bs.NbUPoles() + 1)])


def coefficients_on_basis(func, t, deg):
    """Exact B-spline coefficients of a polynomial on knot vector t (interpolate at Greville points)."""
    n = len(t) - deg - 1
    grev = np.array([t[i + 1:i + deg + 1].mean() for i in range(n)])
    N = BSpline.design_matrix(grev, t, deg).toarray()
    c = np.linalg.solve(N, func(grev))
    check = np.linspace(0, 1, 2001)
    err = np.abs(BSpline(t, c, deg)(check) - func(check)).max()
    return c, err


def make_surface(poles, ukn, umu, udeg, vkn, vmu, vdeg):
    nu, nv = poles.shape[:2]
    arr = TColgp_Array2OfPnt(1, nu, 1, nv)
    for j in range(nu):
        for i in range(nv):
            arr.SetValue(j + 1, i + 1, gp_Pnt(*map(float, poles[j, i])))

    def real_arr(vals):
        a = TColStd_Array1OfReal(1, len(vals))
        for i, v in enumerate(vals):
            a.SetValue(i + 1, float(v))
        return a

    def int_arr(vals):
        a = TColStd_Array1OfInteger(1, len(vals))
        for i, v in enumerate(vals):
            a.SetValue(i + 1, int(v))
        return a

    return Geom_BSplineSurface(arr, real_arr(ukn), real_arr(vkn), int_arr(umu), int_arr(vmu), udeg, vdeg)


def cap_poles(upper_edge, lower_edge, r_coef, eta, zeta, side):
    """Pole grid (nu x nv x 3) for the cap at Y = side * 350."""
    X = upper_edge[:, 0]
    h = upper_edge[:, 2]
    assert np.allclose(lower_edge[:, 0], X, atol=1e-12) and np.allclose(lower_edge[:, 2], -h, atol=1e-12)
    assert np.allclose(upper_edge[:, 1], side * parent.ATTACHMENT, atol=1e-12)
    P = np.empty((len(X), len(eta), 3))
    P[:, :, 0] = X[:, None]
    P[:, :, 1] = side * (parent.ATTACHMENT + r_coef[:, None] * eta[None, :])
    P[:, :, 2] = h[:, None] * zeta[None, :]
    # first and last rows are the parent's own edge poles, copied bit for bit
    P[:, 0, :] = upper_edge
    P[:, -1, :] = lower_edge
    return P


def volume(shape):
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props, 1e-9)
    return props.Mass()


def build():
    shape = read_step(PARENT_DIR / "parent.step")
    skins = parent_skins(shape)
    up_face, up = skins["upper"]
    lo_face, lo = skins["lower"]
    ukn, umu, udeg = knots_of(up, "u")
    assert (ukn, umu, udeg) == knots_of(lo, "u")
    t_u = flat_knots(ukn, umu)

    # r(u) on the parent's u basis (u is the native chord parameter of parent.py)
    r_coef, r_err = coefficients_on_basis(lambda u: 200 * u * (1 - u), t_u, udeg)
    eta, zeta, dev = wrap_design.design()
    vkn_flat = wrap_design.v_knots()
    vkn, vmu = np.unique(vkn_flat, return_counts=True)

    info = {"r_coefficient_fit_max_err_mm": float(r_err), "wrap": wrap_design.report(eta, zeta, dev)}
    cap_faces = {}
    for side, idx in ((+1, -1), (-1, 0)):
        ue, le = edge_poles(up, idx), edge_poles(lo, idx)
        P = cap_poles(ue, le, r_coef, eta, zeta, side)
        srf = make_surface(P, ukn, umu, udeg, vkn, vmu, wrap_design.DEG)
        mf = BRepBuilderAPI_MakeFace(srf, 1e-9)
        if not mf.IsDone():
            raise RuntimeError("cap face failed")
        cap_faces[side] = mf.Face()
        info[f"cap_poles_{'right' if side > 0 else 'left'}"] = list(P.shape[:2])

    sew = BRepBuilderAPI_Sewing(SEW_TOL)
    for f in (up_face, lo_face, cap_faces[+1], cap_faces[-1]):
        sew.Add(f)
    sew.Perform()
    sewn = sew.SewedShape()
    info["sewing"] = {"free_edges": sew.NbFreeEdges(), "multiple_edges": sew.NbMultipleEdges(),
                      "degenerated_shapes": sew.NbDegeneratedShapes()}
    exp = TopExp_Explorer(sewn, TopAbs_SHELL)
    shell = TopoDS.Shell_s(exp.Current())
    exp.Next()
    assert not exp.More(), "more than one shell"
    ms = BRepBuilderAPI_MakeSolid(shell)
    fix = ShapeFix_Solid(ms.Solid())
    fix.Perform()
    solid = fix.Solid()
    info["brepcheck_valid"] = bool(BRepCheck_Analyzer(solid).IsValid())
    info["volume_mm3"] = volume(solid)
    if info["volume_mm3"] < 0:
        raise RuntimeError("solid is inside out")

    Interface_Static.SetCVal_s("write.step.schema", "AP214IS")
    Interface_Static.SetCVal_s("write.step.unit", "MM")
    Interface_Static.SetCVal_s("write.step.product.name", "wing flat sharp rounded caps")
    writer = STEPControl_Writer()
    writer.Transfer(solid, STEPControl_AsIs)
    out = HERE / "wing.step"
    if writer.Write(str(out)) != IFSelect_RetDone:
        raise RuntimeError("STEP write failed")
    data = out.read_bytes()
    info["step_bytes"] = len(data)
    info["step_sha256"] = hashlib.sha256(data).hexdigest()
    (HERE / "build.json").write_text(json.dumps(info, indent=1))
    return info


if __name__ == "__main__":
    info = build()
    print(json.dumps({k: v for k, v in info.items() if k != "wrap"}, indent=1))
