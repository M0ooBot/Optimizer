"""Build wing.step for fin40/blunt: the parent solid plus the rounded tip cap.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/fin40/blunt/build_cap.py

Construction (see README.md):
1. The cap skin is two NURBS patches (+Y and -Y halves).  In the chord
   direction each is an exact rational-quadratic half ellipse (the parent's
   12 % section shape); in the span direction a quintic B-spline carries the
   section centre Xc(v), semi-chord A(v) and height z(v) from cap_def.py.
   Row v = 0 is the parent's (uncut) attachment ellipse, row v = 1 collapses
   to the crown point (280 mm).
2. The two patches and a flat face at Z = 240 are sewn into a cap solid.
3. The cap solid is trimmed by the parent's trailing-edge plane
   x = 160 - 10 z / 240 (a large box on the aft side is subtracted).
4. The trimmed cap is fused with parent.step and coplanar faces are unified,
   so the flat trailing-edge face continues into the cap and closes there.
"""
import hashlib
import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cap_def as cd  # noqa: E402

from OCP.gp import gp_Pnt, gp_Pln, gp_Dir, gp_Ax3, gp_Vec, gp_Trsf, gp_Ax1  # noqa: E402
from OCP.TColgp import TColgp_Array2OfPnt  # noqa: E402
from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger, TColStd_Array2OfReal  # noqa: E402
from OCP.Geom import Geom_BSplineSurface  # noqa: E402
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire,  # noqa: E402
                                BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid, BRepBuilderAPI_Transform)
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox  # noqa: E402
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse  # noqa: E402
from OCP.TopTools import TopTools_ListOfShape  # noqa: E402
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain  # noqa: E402
from OCP.ShapeFix import ShapeFix_Solid, ShapeFix_Shape  # noqa: E402
from OCP.BRepCheck import BRepCheck_Analyzer  # noqa: E402
from OCP.TopExp import TopExp_Explorer  # noqa: E402
from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_SOLID, TopAbs_SHELL  # noqa: E402
from OCP.TopoDS import TopoDS  # noqa: E402
from OCP.BRep import BRep_Tool  # noqa: E402
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface as BRepAdaptor_Surface_  # noqa: E402
from OCP.GeomAbs import GeomAbs_Plane as GeomAbs_Plane_  # noqa: E402
from OCP.Bnd import Bnd_Box  # noqa: E402
from OCP.BRepBndLib import BRepBndLib  # noqa: E402


def _bbox(shape):
    b = Bnd_Box()
    BRepBndLib.AddOptimal_s(shape, b, False, False)
    return list(b.Get())
from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs  # noqa: E402
from OCP.Interface import Interface_Static  # noqa: E402
from OCP.IFSelect import IFSelect_RetDone  # noqa: E402
from OCP.GProp import GProp_GProps  # noqa: E402
from OCP.BRepGProp import BRepGProp  # noqa: E402

PARENT_STEP = cd.PARENT_DIR / "parent.step"
OUT_STEP = HERE / "wing.step"

VDEG = 5
TDEG = 5


def v_of_w(w):
    """Inverse of w = 40 v (2 - v)."""
    return 1.0 - np.sqrt(1.0 - np.asarray(w, float) / cd.REACH)


def v_knots():
    """Span knots: dense through the trailing-edge face closure (w < WE), graded above; the
    knots at w = WE and w = W2 repeat VDEG - 2 times because the blend laws are only C2 there."""
    vwe, vw2 = float(v_of_w(cd.WE)), float(v_of_w(cd.W2))
    k1 = list(np.linspace(0, vwe, 17)[1:-1])
    w2 = np.concatenate([np.geomspace(cd.WE, 12.0, 14)[1:], np.linspace(12.0, 38.0, 27)[1:]])
    k2 = [float(k) for k in v_of_w(w2) if abs(k - vw2) > 1e-3]
    k3 = list(1.0 - np.geomspace(1.0 - v_of_w(38.0), 0.004, 10)[1:])
    interior = sorted(k1 + k2 + k3 + [vwe] * (VDEG - 2) + [vw2] * (VDEG - 2))
    return np.array([0.0] * (VDEG + 1) + interior + [1.0] * (VDEG + 1))


def t_knots():
    """Chordwise knots on t in [0, pi]: uniform, with extra density toward the trailing end."""
    a = list(np.linspace(0, np.pi, 41)[1:-1])
    b = list(np.pi - np.geomspace(0.02, 0.9, 14))
    interior = sorted(set(np.round(a + b, 12)))
    return np.array([0.0] * (TDEG + 1) + interior + [np.pi] * (TDEG + 1))


def skin_grid(t, v):
    """Analytic +Y skin points on the grid t (n) x v (m): arrays (n, m, 3)."""
    v = np.asarray(v, float)
    vj = cd.Jet.constant(v, 1)
    h = vj * (2.0 - vj)
    z = cd.Z0 + cd.REACH * h
    sq = (1.0 - vj) * cd.jsqrt(1.0 + h + h * h)
    e = cd.ends(z, sq)
    xL, XE, xT = e["xL"].v, e["XE"].v, e["xT"].v
    b0, b1, b2 = (c.v for c in e["C"])
    T = np.asarray(t, float)[:, None]
    A = 0.5 * (XE - xL)
    x = xL + A * (1.0 - np.cos(T))
    r = xT - x
    q = (r**3 + b2 * r * r + b1 * r + b0) / (r + cd.ELL) ** 3
    Y = cd.parent.THICKNESS_RATIO * A * np.sin(T) * np.sqrt(q)
    return np.stack([x, Y, np.broadcast_to(z.v, x.shape)], -1)


def fit_skin():
    """Tensor-product least-squares B-spline of the +Y skin; returns knots and poles (nt, nv, 3)."""
    from scipy.interpolate import make_lsq_spline
    tk, vk = t_knots(), v_knots()
    vwe, vw2 = float(v_of_w(cd.WE)), float(v_of_w(cd.W2))
    ts = np.unique(np.concatenate([np.linspace(0, np.pi, 300), (1 - np.cos(np.linspace(0, np.pi, 500))) * np.pi / 2,
                                   np.pi - np.geomspace(1e-4, 1.0, 200)]))
    vs = np.unique(np.concatenate([np.linspace(0, 1, 1500), np.linspace(0, vwe, 600),
                                   np.linspace(0, vw2, 600), 1 - np.geomspace(1e-6, 0.05, 300)]))
    S = skin_grid(ts, vs)                                     # (nt, nv, 3)
    nt, nv = len(ts), len(vs)
    # along v for every t sample, then along t for every v coefficient
    sv = make_lsq_spline(vs, S.transpose(1, 0, 2).reshape(nv, -1), vk, k=VDEG)
    Cv = sv.c.reshape(-1, nt, 3)                              # (ncv, nt, 3)
    st = make_lsq_spline(ts, Cv.transpose(1, 0, 2).reshape(nt, -1), tk, k=TDEG)
    C = st.c.reshape(-1, Cv.shape[0], 3)                      # (nct, ncv, 3)
    # exact structure: z depends on v only; Y = 0 on both lip columns; zero chordwise slope of x
    # at the lip columns (tangent along Y there, so the mirrored patch meets it G1); crown collapse
    zfit = make_lsq_spline(vs, S[0, :, 2], vk, k=VDEG).c
    zfit[0], zfit[-1], zfit[-2] = cd.Z0, cd.Z0 + cd.REACH, cd.Z0 + cd.REACH
    C[:, :, 2] = zfit[None, :]
    C[0, :, 1] = 0.0
    C[-1, :, 1] = 0.0
    # x is separable: x = xL(v) + (XE - xL)(v) c(t), c = (1 - cos t)/2.  Fit c(t) with its end
    # slopes held at zero (c0 = c1 = 0, c_{n-2} = c_{n-1} = 1) and xL, XE - xL along v.
    from scipy.interpolate import BSpline
    Bts = BSpline.design_matrix(ts, tk, TDEG).toarray()
    ct = (1 - np.cos(ts)) / 2
    fixed = np.zeros(Bts.shape[1])
    fixed[-2:] = 1.0
    inner = np.arange(2, Bts.shape[1] - 2)
    sol, *_ = np.linalg.lstsq(Bts[:, inner], ct - Bts @ fixed, rcond=None)
    cco = fixed.copy()
    cco[inner] = sol
    xl = make_lsq_spline(vs, S[0, :, 0], vk, k=VDEG).c
    xe = make_lsq_spline(vs, S[-1, :, 0], vk, k=VDEG).c
    C[:, :, 0] = xl[None, :] + (xe - xl)[None, :] * cco[:, None]
    crown = skin_grid(np.array([0.0]), np.array([1.0]))[0, 0]
    C[:, -1, :] = crown
    # check against the analytic skin on a fresh grid (same parameterisation)
    from scipy.interpolate import BSpline
    tt = np.unique(np.concatenate([np.linspace(0, np.pi, 701), np.pi - np.geomspace(1e-5, 0.5, 100)]))
    vv = np.unique(np.concatenate([np.linspace(0, 1, 2001), np.linspace(0, vwe, 400), 1 - np.geomspace(1e-7, 0.05, 100)]))
    Bt = BSpline.design_matrix(tt, tk, TDEG).toarray()
    Bv = BSpline.design_matrix(vv, vk, VDEG).toarray()
    F = np.einsum("ti,ijc,vj->tvc", Bt, C, Bv)
    err = np.linalg.norm(F - skin_grid(tt, vv), axis=-1)
    i, j = np.unravel_index(np.argmax(err), err.shape)
    info = {"fit_max_error_mm": float(err.max()), "fit_max_error_at_t_v": [float(tt[i]), float(vv[j])],
            "poles_t_v": [int(C.shape[0]), int(C.shape[1])], "degree_t_v": [TDEG, VDEG],
            "attachment_row_max_error_mm": float(err[:, 0].max())}
    return tk, vk, C, info


def make_patch(tk, vk, C, side):
    nt, nv = C.shape[0], C.shape[1]
    poles = TColgp_Array2OfPnt(1, nt, 1, nv)
    for i in range(nt):
        for j in range(nv):
            poles.SetValue(i + 1, j + 1, gp_Pnt(C[i, j, 0], side * C[i, j, 1], C[i, j, 2]))

    def knots(t):
        vals, counts = np.unique(np.round(t, 13), return_counts=True)
        k = TColStd_Array1OfReal(1, len(vals))
        m = TColStd_Array1OfInteger(1, len(vals))
        for i, (a, c) in enumerate(zip(vals, counts)):
            k.SetValue(i + 1, float(a))
            m.SetValue(i + 1, int(c))
        return k, m

    uk, um = knots(tk)
    vk_, vm = knots(vk)
    return Geom_BSplineSurface(poles, uk, vk_, um, vm, TDEG, VDEG, False, False)


def faces_of(shape):
    out = []
    exp = TopExp_Explorer(shape, TopAbs_FACE)
    while exp.More():
        out.append(TopoDS.Face_s(exp.Current()))
        exp.Next()
    return out


def count(shape, kind):
    n = 0
    exp = TopExp_Explorer(shape, kind)
    while exp.More():
        n += 1
        exp.Next()
    return n


def bottom_edge(face):
    """The boundary edge of a cap patch lying in the plane Z = 240."""
    exp = TopExp_Explorer(face, TopAbs_EDGE)
    while exp.More():
        e = TopoDS.Edge_s(exp.Current())
        if not BRep_Tool.Degenerated_s(e):
            c = BRepAdaptor_Curve(e)
            zs = [c.Value(c.FirstParameter() + f * (c.LastParameter() - c.FirstParameter())).Z()
                  for f in (0.0, 0.3, 0.5, 0.7, 1.0)]
            if max(abs(z - cd.Z0) for z in zs) < 1e-9:
                return e
        exp.Next()
    raise RuntimeError("no bottom edge")


def read_step(path):
    r = STEPControl_Reader()
    if r.ReadFile(str(path)) != IFSelect_RetDone:
        raise RuntimeError(f"cannot read {path}")
    r.TransferRoots()
    return r.OneShape()


def volume(shape):
    p = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, p, 1e-9)
    return p.Mass()


def _scipy_to_occ(spl, dim):
    """Convert a clamped scipy BSpline with vector coefficients into a Geom(2d)_BSplineCurve."""
    from OCP.TColgp import TColgp_Array1OfPnt, TColgp_Array1OfPnt2d
    from OCP.Geom import Geom_BSplineCurve
    from OCP.Geom2d import Geom2d_BSplineCurve
    from OCP.gp import gp_Pnt2d
    k = spl.k
    c = np.asarray(spl.c)
    vals, counts = np.unique(np.round(spl.t, 14), return_counts=True)
    kn = TColStd_Array1OfReal(1, len(vals))
    mu = TColStd_Array1OfInteger(1, len(vals))
    for i, (a, m) in enumerate(zip(vals, counts)):
        kn.SetValue(i + 1, float(a))
        mu.SetValue(i + 1, int(m))
    n = len(spl.t) - k - 1
    if dim == 2:
        poles = TColgp_Array1OfPnt2d(1, n)
        for i in range(n):
            poles.SetValue(i + 1, gp_Pnt2d(float(c[i, 0]), float(c[i, 1])))
        return Geom2d_BSplineCurve(poles, kn, mu, k)
    poles = TColgp_Array1OfPnt(1, n)
    for i in range(n):
        poles.SetValue(i + 1, gp_Pnt(float(c[i, 0]), float(c[i, 1]), float(c[i, 2])))
    return Geom_BSplineCurve(poles, kn, mu, k)


def _fit2d(uv, params):
    from scipy.interpolate import make_interp_spline
    return _scipy_to_occ(make_interp_spline(np.asarray(params, float), np.asarray(uv, float), k=5), 2)


def _fit3d(xyz, params):
    from scipy.interpolate import make_interp_spline
    return _scipy_to_occ(make_interp_spline(np.asarray(params, float), np.asarray(xyz, float), k=5), 3)


def fix_edges(wing):
    """Give the attachment and trailing-edge-face edges exact, same-parameter representations.

    The boolean leaves (a) the cap-side pcurve of each attachment edge parameterised
    differently from the parent's 3D edge curve and (b) an approximated intersection
    curve for the edge between each cap patch and the flat trailing-edge face.  Both
    pass the standard check but not the exact one.  Here (a) gets a new pcurve fitted to
    the 3D curve's own parameter and (b) is rebuilt from points solved exactly on the
    cap surface in the trailing-edge plane, with matching pcurves on both faces.
    """
    from scipy.optimize import brentq
    from OCP.TopExp import TopExp
    from OCP.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
    from OCP.BRep import BRep_Builder
    from OCP.GeomProjLib import GeomProjLib
    from OCP.Geom import Geom_Plane
    from OCP.ShapeAnalysis import ShapeAnalysis_Surface
    from OCP.BRepTools import BRepTools

    B = BRep_Builder()
    named = {}
    for f in faces_of(wing):
        a = BRepAdaptor_Surface_(f)
        bb = _bbox(f)
        if a.GetType() == GeomAbs_Plane_:
            named["root" if bb[5] < 1e-3 else "te"] = f
        elif bb[5] < cd.Z0 + 1e-3:
            named["parent_up" if bb[1] + bb[4] > 0 else "parent_lo"] = f
        else:
            named["cap_up" if bb[1] + bb[4] > 0 else "cap_lo"] = f
    emap = TopTools_IndexedDataMapOfShapeListOfShape()
    TopExp.MapShapesAndAncestors_s(wing, TopAbs_EDGE, TopAbs_FACE, emap)

    def owners(e):
        lst = emap.FindFromKey(e)
        out = []
        for f in lst:
            for k, g in named.items():
                if g.IsSame(f):
                    out.append(k)
        return sorted(set(out))

    report = {}
    te_plane = BRep_Tool.Surface_s(named["te"])
    te_pln = Geom_Plane.DownCast(te_plane) if hasattr(Geom_Plane, "DownCast") else te_plane
    for i in range(1, emap.Extent() + 1):
        e = TopoDS.Edge_s(emap.FindKey(i))
        if BRep_Tool.Degenerated_s(e):
            continue
        own = owners(e)
        c = BRepAdaptor_Curve(e)
        t0, t1 = c.FirstParameter(), c.LastParameter()
        for side in ("up", "lo"):
            cap = "cap_" + side
            if own == sorted(["parent_" + side, cap]):
                srf = BRep_Tool.Surface_s(named[cap])
                sas = ShapeAnalysis_Surface(srf)
                from scipy.interpolate import make_lsq_spline
                uproj = lambda t: sas.ValueOfUV(c.Value(t), 1e-12).X()
                ts = np.linspace(t0, t1, 3001)
                us = np.array([uproj(t) for t in ts])
                kint = np.linspace(t0, t1, 161)[1:-1]
                knots = np.concatenate([[t0] * 6, kint, [t1] * 6])
                sp = make_lsq_spline(ts, np.stack([us, np.zeros_like(us)], 1), knots, k=5)
                cc = sp.c.copy()
                cc[0, 0], cc[-1, 0] = us[0], us[-1]
                sp.c = cc
                pc = _scipy_to_occ(sp, 2)
                dev = 0.0
                for t in np.linspace(t0, t1, 4001):
                    q2 = pc.Value(t)
                    dev = max(dev, c.Value(t).Distance(srf.Value(q2.X(), q2.Y())))
                if dev > 1e-6:
                    import sys as _s
                    print("DEBUG attachment", side, t0, t1, c.Value(t0).Coord(), c.Value(t1).Coord(), uv[:3], uv[-3:], file=_s.stderr)
                    worst = max(range(4001), key=lambda k: 0)
                    for t in np.linspace(t0, t1, 11):
                        q2 = pc.Value(t); p = c.Value(t); q = sas.ValueOfUV(p, 1e-12)
                        print("  t", t, p.Coord(), "fit u", q2.X(), "proj u", q.X(), q.Y(), "gap", sas.Gap(), file=_s.stderr)
                    raise RuntimeError(f"attachment pcurve deviation {dev}")
                tol = max(1e-7, 1.5 * dev)
                B.UpdateEdge(e, pc, named[cap], tol)
                B.SameParameter(e, True)
                report[f"attachment_{side}_cap_pcurve_max_deviation_mm"] = dev
            if own == sorted([cap, "te"]):
                srf = BRep_Tool.Surface_s(named[cap])
                xz = lambda u, v: srf.Value(u, v)
                te_x = lambda z: 160.0 - 10.0 * z / 240.0
                umax = srf.Bounds()[1]
                gv = lambda v: xz(umax, v).X() - te_x(xz(umax, v).Z())
                vc = brentq(gv, 1e-9, 0.1, xtol=1e-15)
                taus = np.linspace(0.0, 1.0, 601)
                uv, xyz = [], []
                for tau in taus:
                    v = vc * (1.0 - (1.0 - tau) ** 2)
                    fu = lambda u: xz(u, v).X() - te_x(xz(u, v).Z())
                    u = umax if tau == 1.0 else brentq(fu, 0.5 * umax, umax, xtol=1e-15)
                    uv.append((u, v))
                    p = xz(u, v)
                    xyz.append((p.X(), p.Y(), p.Z()))
                # direction of the existing edge: does its curve start at the corner (z = 240)?
                starts_at_corner = abs(c.Value(t0).Z() - cd.Z0) < abs(c.Value(t1).Z() - cd.Z0)
                params = taus if starts_at_corner else 1.0 - taus
                c3 = _fit3d(xyz, params)
                pc_cap = _fit2d(uv, params)
                pc_te = GeomProjLib.Curve2d_s(c3, 0.0, 1.0, te_plane)
                dev = 0.0
                for t in np.linspace(0, 1, 4001):
                    a = c3.Value(t)
                    q = pc_cap.Value(t)
                    dev = max(dev, a.Distance(srf.Value(q.X(), q.Y())))
                    q = pc_te.Value(t)
                    dev = max(dev, a.Distance(te_plane.Value(q.X(), q.Y())))
                tol = max(1e-7, 1.5 * dev)
                B.UpdateEdge(e, c3, tol)
                B.UpdateEdge(e, pc_cap, named[cap], tol)
                B.UpdateEdge(e, pc_te, named["te"], tol)
                B.Range(e, 0.0, 1.0)
                B.SameRange(e, True)
                B.SameParameter(e, True)
                # move the tip vertex onto the exact tip
                tip = gp_Pnt(*xyz[-1])
                from OCP.TopExp import TopExp as TE_
                v1, v2 = TE_.FirstVertex_s(e), TE_.LastVertex_s(e)
                for vx in (v1, v2):
                    if abs(BRep_Tool.Pnt_s(vx).Z() - cd.Z0) > 1e-3:
                        B.UpdateVertex(vx, tip, 1e-7)
                report[f"te_arch_{side}_max_deviation_mm"] = dev
                report["tip_v_on_cap"] = vc
                report["tip_xyz"] = list(xyz[-1])
    # the lip edge on the trailing-edge branch starts at the tip: set its range to start at v_c
    vc = report["tip_v_on_cap"]
    for i in range(1, emap.Extent() + 1):
        e = TopoDS.Edge_s(emap.FindKey(i))
        if BRep_Tool.Degenerated_s(e) or owners(e) != ["cap_lo", "cap_up"]:
            continue
        c = BRepAdaptor_Curve(e)
        if abs(c.Value(c.FirstParameter()).X() - 150.0) < 1.0 or abs(c.Value(c.LastParameter()).X() - 150.0) < 1.0:
            t0, t1 = c.FirstParameter(), c.LastParameter()
            if t0 > 1e-3 and t0 < 0.1:
                B.Range(e, vc, t1)
            report["lip_te_edge_range"] = [vc, t1]
    return report


def main():
    tk, vk, C, fit_info = fit_skin()
    up = make_patch(tk, vk, C, +1)
    lo = make_patch(tk, vk, C, -1)
    f_up = BRepBuilderAPI_MakeFace(up, 1e-9).Face()
    f_lo = BRepBuilderAPI_MakeFace(lo, 1e-9).Face()
    wire = BRepBuilderAPI_MakeWire(bottom_edge(f_up), bottom_edge(f_lo)).Wire()
    f_bot = BRepBuilderAPI_MakeFace(wire, True).Face()

    sew = BRepBuilderAPI_Sewing(1e-6)
    for f in (f_up, f_lo, f_bot):
        sew.Add(f)
    sew.Perform()
    shell = TopoDS.Shell_s(sew.SewedShape())
    cap = BRepBuilderAPI_MakeSolid(shell).Solid()
    fix = ShapeFix_Solid(cap)
    fix.Perform()
    cap = fix.Solid()
    cap_volume_full = volume(cap)
    if cap_volume_full < 0:
        raise RuntimeError("cap solid is inside out")

    # trailing-edge plane x = 160 - 10 z/240: subtract a box on its aft side
    box = BRepPrimAPI_MakeBox(gp_Pnt(0.0, -50.0, 200.0), 100.0, 100.0, 120.0).Shape()
    # place the box face x' = 0 on the plane: rotate about Y then translate
    ang = np.arctan2(10.0, 240.0)  # plane normal (1, 0, 10/240) -> tilt about +Y by -ang
    tr = gp_Trsf()
    tr.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)), -ang)
    t2 = gp_Trsf()
    t2.SetTranslation(gp_Vec(160.0, 0.0, 0.0))
    box = BRepBuilderAPI_Transform(box, t2.Multiplied(tr), True).Shape()
    # sanity: box near face must lie on the plane x = 160 - 10 z / 240
    cut = BRepAlgoAPI_Cut(cap, box)
    cut.SetFuzzyValue(1e-7)
    cut.Build()
    if not cut.IsDone():
        raise RuntimeError("TE trim failed")
    cap_cut = cut.Shape()

    parent = read_step(PARENT_STEP)
    fuse = BRepAlgoAPI_Fuse(parent, cap_cut)
    fuse.SetFuzzyValue(1e-6)
    fuse.Build()
    if not fuse.IsDone():
        raise RuntimeError("fuse failed")
    unify = ShapeUpgrade_UnifySameDomain(fuse.Shape(), False, True, False)  # faces only, keep skins apart
    unify.Build()
    wing = unify.Shape()
    sf = ShapeFix_Shape(wing)
    sf.Perform()
    wing = sf.Shape()
    edge_fix = fix_edges(wing)

    ana = BRepCheck_Analyzer(wing)
    ana_exact = BRepCheck_Analyzer(wing, True, False, True)
    info = {
        "valid": bool(ana.IsValid()),
        "exact_valid_before_export": bool(ana_exact.IsValid()),
        "edge_fix": edge_fix,
        "solids": count(wing, TopAbs_SOLID),
        "shells": count(wing, TopAbs_SHELL),
        "faces": count(wing, TopAbs_FACE),
        "volume_mm3": volume(wing),
        "parent_volume_mm3": volume(parent),
        "cap_volume_untrimmed_mm3": cap_volume_full,
        "cap_volume_trimmed_mm3": volume(cap_cut),
        "skin_fit": fit_info,
        "design_constants_mm": {"WE": cd.WE, "W2": cd.W2, "ELL": cd.ELL},
        "te_face_closure_height_above_attachment_mm": cd.te_closure_height(),
    }
    Interface_Static.SetCVal_s("write.step.schema", "AP214IS")
    Interface_Static.SetCVal_s("write.step.unit", "MM")
    w = STEPControl_Writer()
    w.Transfer(wing, STEPControl_AsIs)
    if w.Write(str(OUT_STEP)) != IFSelect_RetDone:
        raise RuntimeError("STEP write failed")
    info["step_bytes"] = OUT_STEP.stat().st_size
    info["step_sha256"] = hashlib.sha256(OUT_STEP.read_bytes()).hexdigest()
    (HERE / "build.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
