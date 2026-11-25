"""Build the Wing B (sharp) rounded tip cap and write wing.step.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/wingb/sharp/build_cap.py

Construction (see README.md):
  * profiles: upper cut edge -> 3D perimeter (locked outline at the continued
    camber height) -> lower cut edge;
  * one spanwise wrap guide per outline parameter p, running straight in the
    continued planform coordinates (S, xi) from the attachment station
    xi_a(p) to the perimeter point O(p), carrying the parent's camber and a
    squared half-thickness Q = Qb - Q_O (1 - v/a)^m that closes with a
    designed lip radius;
  * the guide net is interpolated by one bicubic B-spline surface
    (degenerate at the leading-edge and trailing-edge points) and sewn to the
    parent's skin and root faces.
"""
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import capgeom as g  # noqa: E402

NP = 241        # guides (outline parameter nodes, including the two collapsed ends)
NHALF = 100     # wrap intervals per side (attachment -> lip)
EPS_FRACTION = 0.5  # wrap-parameter lip scale as a fraction of the guide's design lip radius


P_COS = 1.0  # share of cosine spacing in the guide distribution


def p_of_t(t):
    """Outline parameter of the surface u-parameter t: half cosine, half uniform spacing, so guide
    rows are somewhat denser at the collapsed ends while staying linear in t there."""
    t = np.asarray(t, float)
    return P_COS * 0.5 * (1.0 - np.cos(np.pi * t)) + (1 - P_COS) * t


def guide_grid(np_=NP, nhalf=NHALF):
    t = np.linspace(0.0, 1.0, np_)
    p = p_of_t(t)
    d = g.design_guides(p[1:-1])
    sig = np.linspace(0.0, 2.0, 2 * nhalf + 1)
    grid = np.zeros((np_, 2 * nhalf + 1, 3))
    ders = np.zeros((np_, 2, 2, 3))   # [row, end (attachment upper / lower), order-1, xyz]
    tab = g.attachment_tables()
    # exact attachment points
    xa = g.xi_attach(p)
    jet = g.attachment_jet(xa)
    Xa, Za = g.plan_xz(np.full_like(xa, g.S_ATT), xa)
    up_att = np.stack([Xa, jet[1][0], Za], -1)
    lo_att = np.stack([Xa, jet[-1][0], Za], -1)
    meta = []
    for k in range(1, np_ - 1):
        i = k - 1
        G = g.guide_dense(d, i, K=4000, q=4.0)
        upper = np.stack([G["X"], G["Cm"] + G["H"], G["Z"]], -1)
        lower = np.stack([G["X"], G["Cm"] - G["H"], G["Z"]], -1)
        loop = np.r_[upper, lower[::-1][1:]]          # upper attachment -> lip -> lower attachment
        i_lip = len(upper) - 1
        keep = np.r_[True, np.linalg.norm(np.diff(loop, axis=0), axis=1) > 1e-13]
        i_lip = int(np.sum(keep[:i_lip + 1]) - 1)
        loop = loop[keep]
        eps = EPS_FRACTION * max(d["R_lip_design"][i], 1e-6)
        tau = g.loop_parameter(loop, i_lip, eps)
        grid[k] = np.stack([np.interp(sig, tau, loop[:, c]) for c in range(3)], -1)
        # exact first/second sigma-derivatives at both attachment ends (base surface = parent's G2 continuation)
        seg = np.linalg.norm(np.diff(loop, axis=0), axis=1)
        arc = np.r_[0, np.cumsum(seg)]
        s_up, s_lo = arc[i_lip], arc[-1] - arc[i_lip]
        for end, (s_end, sgn) in enumerate(((s_up, +1.0), (s_lo, -1.0))):
            T, kN = g.guide_end_jet(d, i, +1 if end == 0 else -1)
            lin = 0.5
            phi = np.arcsinh(s_end / eps) + lin
            f1 = 1 / np.sqrt(s_end**2 + eps**2) + lin / s_end
            f2 = -s_end / (s_end**2 + eps**2) ** 1.5
            da = sgn * phi / f1
            dda = phi**2 * f2 / f1**3
            ders[k, end, 0] = T * da
            ders[k, end, 1] = kN * da * da + T * dda
        grid[k, 0] = up_att[k]
        grid[k, -1] = lo_att[k]
    grid[0, :] = up_att[0]          # leading-edge point (upper = lower at the nose)
    grid[-1, :] = up_att[-1]        # trailing-edge point
    return t, sig, grid, d, ders


VDEG = 5  # wrap direction: quintic, clamped to the base surface's first and second derivatives


def bspline_from_grid(u, v, grid, ders):
    """Interpolating B-spline: quintic in the wrap direction with the exact attachment derivatives
    (so the join is tangent- and curvature-continuous up to interpolation in u), cubic across guides."""
    from scipy.interpolate import make_interp_spline
    from OCP.Geom import Geom_BSplineSurface
    from OCP.TColgp import TColgp_Array2OfPnt
    from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
    from OCP.gp import gp_Pnt
    rows = []
    kv = None
    for k in range(grid.shape[0]):
        bc = ([(1, ders[k, 0, 0]), (2, ders[k, 0, 1])], [(1, ders[k, 1, 0]), (2, ders[k, 1, 1])])
        sp = make_interp_spline(v, grid[k], k=VDEG, bc_type=bc)
        rows.append(sp.c)
        kv = sp.t
    rows = np.array(rows)                       # (nu, nV, 3)
    spu = make_interp_spline(u, rows, k=3, axis=0)
    ku = spu.t
    C = spu.c                                   # (nU, nV, 3)
    nU, nV = C.shape[:2]
    poles = TColgp_Array2OfPnt(1, nU, 1, nV)
    for i in range(nU):
        for j in range(nV):
            poles.SetValue(i + 1, j + 1, gp_Pnt(*map(float, C[i, j])))

    def km(kn):
        uniq, mult = np.unique(np.round(kn, 14), return_counts=True)
        a = TColStd_Array1OfReal(1, len(uniq))
        b = TColStd_Array1OfInteger(1, len(uniq))
        for i, (x, m) in enumerate(zip(uniq, mult)):
            a.SetValue(i + 1, float(x))
            b.SetValue(i + 1, int(m))
        return a, b
    uk, um = km(ku)
    vk, vm = km(kv)
    return Geom_BSplineSurface(poles, uk, vk, um, vm, 3, VDEG), C


def build(out_step):
    from OCP.BRep import BRep_Tool
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeSolid,
                                    BRepBuilderAPI_Sewing)
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.GeomAbs import GeomAbs_Plane
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.Interface import Interface_Static
    from OCP.STEPControl import STEPControl_AsIs, STEPControl_Reader, STEPControl_Writer
    from OCP.ShapeFix import ShapeFix_Solid, ShapeFix_Shape
    from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    t0 = time.time()
    u, v, grid, d, ders = guide_grid()
    surf, C = bspline_from_grid(u, v, grid, ders)
    cap_face = BRepBuilderAPI_MakeFace(surf, 1e-7).Face()

    reader = STEPControl_Reader()
    assert reader.ReadFile(str(g.PARENT_DIR / "parent.step")) == IFSelect_RetDone
    reader.TransferRoots()
    parent = reader.OneShape()
    keep = []
    exp = TopExp_Explorer(parent, TopAbs_FACE)
    while exp.More():
        f = TopoDS.Face_s(exp.Current())
        s = BRepAdaptor_Surface(f)
        if s.GetType() == GeomAbs_Plane and abs(s.Plane().Location().X() + 540.0) < 1e-6 \
                and abs(abs(s.Plane().Axis().Direction().X()) - 1) < 1e-9:
            pass  # the flat attachment face is replaced by the cap
        else:
            keep.append(f)
        exp.Next()
    assert len(keep) == 2, len(keep)
    sew = BRepBuilderAPI_Sewing(1e-4)
    for f in keep:
        sew.Add(f)
    sew.Add(cap_face)
    sew.Perform()
    sewn = sew.SewedShape()
    exp = TopExp_Explorer(sewn, TopAbs_SHELL)
    shell = TopoDS.Shell_s(exp.Current())
    solid = BRepBuilderAPI_MakeSolid(shell).Solid()
    fix = ShapeFix_Solid(solid)
    fix.Perform()
    solid = fix.Solid()
    ok = BRepCheck_Analyzer(solid).IsValid()
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(solid, props, 1e-9)
    w = STEPControl_Writer()
    Interface_Static.SetCVal_s("write.step.unit", "MM")
    w.Transfer(solid, STEPControl_AsIs)
    assert w.Write(str(out_step)) == IFSelect_RetDone
    info = dict(valid=bool(ok), volume_mm3=props.Mass(), free_edges=sew.NbFreeEdges(),
                multiple_edges=sew.NbMultipleEdges(), degenerated=sew.NbDegeneratedShapes(),
                poles=list(C.shape[:2]), seconds=time.time() - t0)
    return info, (u, v, grid, d, surf)


if __name__ == "__main__":
    out = g.BUILD_DIR / "wing.step"
    info, _ = build(out)
    data = out.read_bytes()
    info["step_bytes"] = len(data)
    info["step_sha256"] = hashlib.sha256(data).hexdigest()
    info["design"] = dict(g.DESIGN, taper_level_off_mm=g.TAPER_RETURN_MM, lip_bead=g.LIP_BEAD,
                          xi_attach_p0=g.XI_ATTACH_P0, guides=NP, guide_spacing="cosine in p",
                          wrap_intervals_per_side=NHALF, wrap_eps_fraction_of_lip_radius=EPS_FRACTION,
                          degrees=dict(across_guides=3, wrap=VDEG))
    (g.BUILD_DIR / "build.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))
