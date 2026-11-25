"""Build wing.step for fin5 / sharp: the parent fin closed by a rounded tip cap.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/fin5/sharp/build_cap.py

The parent's two skin faces and its root face are taken unchanged from
parent.step; its flat attachment face (Z = 240) is replaced by two cap patches
(+Y and -Y) built here. See cap_design.py and README.md for the construction.
"""
import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cap_design as cd  # noqa: E402

from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCP.IFSelect import IFSelect_RetDone
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_SHELL
from OCP.TopoDS import TopoDS
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Plane
from OCP.Geom import Geom_BSplineSurface
from OCP.TColgp import TColgp_Array2OfPnt
from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
from OCP.gp import gp_Pnt
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepLib import BRepLib
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.Interface import Interface_Static

PARENT_STEP = cd.PARENT_DIR / "parent.step"
OUT_STEP = HERE / "wing.step"
SEW_TOL = 1e-7


def read_step(path):
    r = STEPControl_Reader()
    if r.ReadFile(str(path)) != IFSelect_RetDone:
        raise RuntimeError(f"cannot read {path}")
    r.TransferRoots()
    return r.OneShape()


def faces(shape):
    out = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        out.append(TopoDS.Face_s(ex.Current()))
        ex.Next()
    return out


def classify_parent(shape):
    found = {}
    for f in faces(shape):
        a = BRepAdaptor_Surface(f)
        if a.GetType() == GeomAbs_BSplineSurface:
            s = a.BSpline()
            side = "upper" if s.Value(0.5, cd.ATT).Y() > 0 else "lower"
            found[side] = (f, s)
        elif a.GetType() == GeomAbs_Plane:
            z = a.Plane().Location().Z()
            found["root" if abs(z) < 1e-6 else "top"] = (f, None)
    assert set(found) == {"upper", "lower", "root", "top"}, found.keys()
    return found


def poles(s):
    return np.array([[[s.Pole(i, j).X(), s.Pole(i, j).Y(), s.Pole(i, j).Z()]
                      for j in range(1, s.NbVPoles() + 1)] for i in range(1, s.NbUPoles() + 1)])


def u_knots(s):
    return ([s.UKnot(i) for i in range(1, s.NbUKnots() + 1)],
            [s.UMultiplicity(i) for i in range(1, s.NbUKnots() + 1)])


def make_surface(net, uk, um, udeg, vt, vdeg):
    nu, nv = net.shape[:2]
    arr = TColgp_Array2OfPnt(1, nu, 1, nv)
    for i in range(nu):
        for j in range(nv):
            arr.SetValue(i + 1, j + 1, gp_Pnt(*map(float, net[i, j])))
    vk, vm = np.unique(vt, return_counts=True)
    a_uk, a_um = TColStd_Array1OfReal(1, len(uk)), TColStd_Array1OfInteger(1, len(uk))
    for i, (k, m) in enumerate(zip(uk, um)):
        a_uk.SetValue(i + 1, float(k))
        a_um.SetValue(i + 1, int(m))
    a_vk, a_vm = TColStd_Array1OfReal(1, len(vk)), TColStd_Array1OfInteger(1, len(vk))
    for i, (k, m) in enumerate(zip(vk, vm)):
        a_vk.SetValue(i + 1, float(k))
        a_vm.SetValue(i + 1, int(m))
    return Geom_BSplineSurface(arr, a_uk, a_vk, a_um, a_vm, udeg, vdeg)


def precise_volume(shape):
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props, 1e-9)
    return props.Mass()


def build():
    parent = read_step(PARENT_STEP)
    pf = classify_parent(parent)
    up_s, lo_s = pf["upper"][1], pf["lower"][1]
    Pu, Pl = poles(up_s), poles(lo_s)
    uk, um = u_knots(up_s)
    assert u_knots(lo_s) == (uk, um)
    assert up_s.UDegree() == lo_s.UDegree()
    vlast = Pu.shape[1] - 1
    top_u, top_l = Pu[:, vlast], Pl[:, vlast]
    assert np.allclose(top_u[:, 2], cd.ATT, atol=0) and np.allclose(top_l[:, 2], cd.ATT, atol=0)
    mirror_err = float(np.max(np.abs(top_l * [1, -1, 1] - top_u)))

    # chord functions of the parent section in its own u parameter:
    # top row = (m0 + a0 cx_i, 0.12 a0 sy_i, 240)
    a0, m0 = cd.a_of(cd.ATT), cd.m_of(cd.ATT)
    cx = (top_u[:, 0] - m0) / a0
    sy = top_u[:, 1] / (cd.TR * a0)
    # check the whole parent net has this structure (rows at Z_j = 0, 60, ..., 240)
    struct_err = 0.0
    for j in range(Pu.shape[1]):
        z = Pu[0, j, 2]
        pred = np.stack([cd.m_of(z) + cd.a_of(z) * cx, cd.TR * cd.a_of(z) * sy, np.full_like(cx, z)], -1)
        struct_err = max(struct_err, float(np.max(np.abs(pred - Pu[:, j]))))

    vt, coef, fit_err = cd.fit_spanwise()
    Mx, Ax, Ay, Zc = coef["Mx"], coef["Ax"], coef["Ay"], coef["Z"]
    nv = len(Mx)
    net_u = np.empty((len(cx), nv, 3))
    for j in range(nv):
        net_u[:, j] = np.stack([Mx[j] + Ax[j] * cx, Ay[j] * sy, np.full_like(cx, Zc[j])], -1)
    row0_formula_err = float(np.max(np.abs(net_u[:, 0] - top_u)))
    net_u[:, 0] = top_u                      # attachment row = the parent's own edge poles
    net_l = net_u * [1, -1, 1]
    net_l[:, 0] = top_l
    net_u[:, -1] = net_u[0, -1]              # collapsed crown row, exactly one point
    net_l[:, -1] = net_u[0, -1]

    udeg = up_s.UDegree()
    cap_u = make_surface(net_u, uk, um, udeg, vt, cd.V_DEGREE)
    cap_l = make_surface(net_l, uk, um, udeg, vt, cd.V_DEGREE)

    fu = BRepBuilderAPI_MakeFace(cap_u, 1e-9).Face()
    fl = BRepBuilderAPI_MakeFace(cap_l, 1e-9).Face()

    sew = BRepBuilderAPI_Sewing(SEW_TOL)
    for f in (pf["upper"][0], pf["lower"][0], pf["root"][0], fu, fl):
        sew.Add(f)
    sew.Perform()
    sewn = sew.SewedShape()
    ex = TopExp_Explorer(sewn, TopAbs_SHELL)
    shell = TopoDS.Shell_s(ex.Current())
    ex.Next()
    assert not ex.More(), "more than one shell after sewing"
    solid = BRepBuilderAPI_MakeSolid(shell).Solid()
    BRepLib.OrientClosedSolid_s(solid)
    valid = BRepCheck_Analyzer(solid).IsValid()
    vol = precise_volume(solid)
    assert valid and vol > 0, (valid, vol)

    Interface_Static.SetCVal_s("write.step.schema", "AP214IS")
    Interface_Static.SetCVal_s("write.step.unit", "MM")
    w = STEPControl_Writer()
    w.Transfer(solid, STEPControl_AsIs)
    assert w.Write(str(OUT_STEP)) == IFSelect_RetDone

    info = {
        "parent_step": str(PARENT_STEP.relative_to(cd.ROOT)),
        "output": str(OUT_STEP.relative_to(cd.ROOT)),
        "construction": "elliptic horizontal sections: X = m(Z) - a(Z) sqrt(1-h^3) cos(theta), Y = +/-Yp(h) sin(theta), h = (Z-240)/5",
        "wrap_profile": {
            "ramp_fraction_of_arc": cd.RAMP,
            "peak_curvature_per_mm": cd.PROFILE["K"],
            "tightest_profile_radius_mm": 1 / cd.PROFILE["K"],
            "crown_profile_radius_mm": 1 / (cd.PROFILE["K"] * cd.PROFILE["rho"]),
            "arc_length_mm": cd.PROFILE["S"],
            "attachment_lean_rad": float(cd.LEAN),
        },
        "patch": {"u": "parent skin u (chord parameter s = theta/pi), parent knots and degree",
                  "u_degree": udeg, "u_poles": len(cx),
                  "v": "h = 1 - (1 - v)^2", "v_degree": cd.V_DEGREE, "v_poles": nv,
                  "v_knots": [float(k) for k in vt]},
        "spanwise_fit_max_error_mm": fit_err,
        "parent_net_structure_max_error_mm": struct_err,
        "parent_lower_vs_mirrored_upper_top_row_mm": mirror_err,
        "cap_row0_formula_vs_parent_edge_poles_mm": row0_formula_err,
        "sewing_tolerance_mm": SEW_TOL,
        "brepcheck_valid": bool(valid),
        "volume_mm3": vol,
    }
    (HERE / "build.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    build()
