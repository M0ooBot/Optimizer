"""Build the rounded tip caps for the twisted/blunt parent and write wing.step.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/twisted/blunt/build_cap.py

Steps: per-guide design (cap_model.design) -> dense samples of three cap faces
per tip (upper flat region, lip band, lower flat region) -> tensor-product
least-squares B-spline fits with shared boundary curves -> faces sewn to the
parent's skins (parent.step with its two planar end faces removed) -> one solid
-> STEP.  The left cap is the mirror image of the right cap in the XZ plane.
"""
import json
import os
import sys
import time
import pickle
import hashlib
import numpy as np
from scipy.interpolate import BSpline

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cap_model as cm  # noqa: E402

from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCP.IFSelect import IFSelect_RetDone
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID
from OCP.TopoDS import TopoDS
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_Plane
from OCP.TColgp import TColgp_Array2OfPnt
from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
from OCP.gp import gp_Pnt, gp_Trsf, gp_Ax2, gp_Dir
from OCP.Geom import Geom_BSplineSurface
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid,
                                BRepBuilderAPI_Transform)
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.ShapeFix import ShapeFix_Solid, ShapeFix_Shape
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp
from OCP.Interface import Interface_Static
from OCP.Geom2d import Geom2d_BSplineCurve
from OCP.TColgp import TColgp_Array1OfPnt2d
from OCP.gp import gp_Pnt2d
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
from OCP.BRepLib import BRepLib
from OCP.ShapeFix import ShapeFix_Face
from OCP.GeomAdaptor import GeomAdaptor_Surface
from OCP.GeomAbs import GeomAbs_BSplineSurface

OUT_STEP = os.path.join(HERE, "wing.step")
OUT_FIT = os.path.join(HERE, "build_fit.json")
CACHE = os.path.join(HERE, "_design_cache.pkl")

NG, NG_CP = 1201, 401          # samples / control points along the rim parameter g
NL, NL_CP = 161, 61           # lip: samples / control points across the wrap
NF, NF_CP = 61, 21            # flat faces: samples / control points from band edge to attachment
DEG = 3
SEW_TOL = 2e-4
Y_TOP = 501.0                 # the continued-skin surfaces of the flat faces span Y = 440 .. Y_TOP
Y_KNOTS = cm.CUT + np.array([0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.23, 0.27, 0.32, 0.4, 0.5,
                             0.65, 0.85, 1.1, 1.4, 1.8, 2.3, 2.9, 3.6, 4.5, 6.0, 8.0, 11.0, 15.0, 21.0, 29.0, 40.0])


def clamped_knots(n_cp, deg):
    inner = np.linspace(0, 1, n_cp - deg + 1)
    return np.concatenate([[0.0] * deg, inner, [1.0] * deg])


def basis(x, t, deg):
    return BSpline.design_matrix(np.clip(x, 0, 1), t, deg).toarray()


def lsq_curve(x, P, t, deg, fix_ends=True):
    """1-D least squares with the end control points fixed to the end data."""
    Bm = basis(x, t, deg)
    n = Bm.shape[1]
    C = np.zeros((n, P.shape[-1]))
    if fix_ends:
        C[0], C[-1] = P[0], P[-1]
        rhs = P - np.outer(Bm[:, 0], C[0]) - np.outer(Bm[:, -1], C[-1])
        C[1:-1] = np.linalg.lstsq(Bm[:, 1:-1], rhs, rcond=None)[0]
    else:
        C = np.linalg.lstsq(Bm, P, rcond=None)[0]
    return C


def lsq_surface(xg, xs, X, tg, ts, deg, fixed):
    """Tensor least squares, X (Ng, Ns, 3); fixed = dict of boundary rows/cols of control points."""
    Bg, Bs = basis(xg, tg, deg), basis(xs, ts, deg)
    ng, ns = Bg.shape[1], Bs.shape[1]
    C = np.zeros((ng, ns, 3))
    C[0, :], C[-1, :] = fixed["g0"], fixed["g1"]
    C[:, 0], C[:, -1] = fixed["s0"], fixed["s1"]
    Cb = C.copy()
    Cb[1:-1, 1:-1] = 0.0
    Xr = X - np.einsum("ai,ijk,bj->abk", Bg, Cb, Bs)
    Pg = np.linalg.pinv(Bg[:, 1:-1]); Ps = np.linalg.pinv(Bs[:, 1:-1])
    C[1:-1, 1:-1] = np.einsum("ia,abk,jb->ijk", Pg, Xr, Ps)
    fitX = np.einsum("ai,ijk,bj->abk", Bg, C, Bs)
    return C, np.linalg.norm(fitX - X, axis=-1)


def occ_surface(C, tg, ts, deg):
    ng, ns = C.shape[:2]
    poles = TColgp_Array2OfPnt(1, ng, 1, ns)
    for i in range(ng):
        for j in range(ns):
            poles.SetValue(i + 1, j + 1, gp_Pnt(*map(float, C[i, j])))

    def kn(t):
        u, c = np.unique(t, return_counts=True)
        K = TColStd_Array1OfReal(1, len(u)); M = TColStd_Array1OfInteger(1, len(u))
        for k, (a, b) in enumerate(zip(u, c)):
            K.SetValue(k + 1, float(a)); M.SetValue(k + 1, int(b))
        return K, M

    Ku, Mu = kn(tg); Kv, Mv = kn(ts)
    return Geom_BSplineSurface(poles, Ku, Kv, Mu, Mv, deg, deg)


def occ_curve2d(C, t, deg):
    n = len(C)
    poles = TColgp_Array1OfPnt2d(1, n)
    for i in range(n):
        poles.SetValue(i + 1, gp_Pnt2d(float(C[i, 0]), float(C[i, 1])))
    u, c = np.unique(t, return_counts=True)
    K = TColStd_Array1OfReal(1, len(u)); M = TColStd_Array1OfInteger(1, len(u))
    for k, (a, b) in enumerate(zip(u, c)):
        K.SetValue(k + 1, float(a)); M.SetValue(k + 1, int(b))
    return Geom2d_BSplineCurve(poles, K, M, deg)


def parent_skin_rows(parent_shape):
    """s-knots (degree 5) and the poles of the upper/lower right skins along the attachment (|Y| = 440)."""
    rows = {}
    for f in faces_of(parent_shape):
        ad = GeomAdaptor_Surface(BRep_Tool.Surface_s(f))
        if ad.GetType() != GeomAbs_BSplineSurface:
            continue
        b = ad.BSpline()
        if b.NbUPoles() < 100:
            continue
        P = np.array([[(b.Pole(i, j).X(), b.Pole(i, j).Y(), b.Pole(i, j).Z()) for j in range(1, b.NbVPoles() + 1)]
                      for i in range(1, b.NbUPoles() + 1)])
        if P[..., 1].mean() < 0:
            continue
        t = []
        for k in range(1, b.NbUKnots() + 1):
            t += [b.UKnot(k)] * b.UMultiplicity(k)
        if abs(b.VKnot(b.NbVKnots()) - cm.CUT) > 1e-9:
            raise RuntimeError("unexpected parent parametrisation")
        rows[float(P[:, -1, 2].mean())] = (np.array(t), b.UDegree(), P[:, -1])
    zs = sorted(rows)
    out = {1: rows[zs[-1]], -1: rows[zs[0]]}
    for side in (1, -1):
        t, d, row = out[side]
        chk = BSpline(t, row, d)(np.array([0.0, 0.3, 1.0]))
        if np.abs(chk - cm.pa.skin(np.array([0.0, 0.3, 1.0]), cm.CUT, side)).max() > 1e-7:
            raise RuntimeError("parent attachment row does not match parent.py")
    return out


def flat_surface(side, srows):
    """Continued (and TE-compensated) parent skin over s in [0, 1], Y in [440, Y_TOP]; row Y = 440 = parent poles."""
    ts, ds, row440 = srows[side]
    dy = 5
    ty = np.concatenate([[cm.CUT] * (dy + 1), Y_KNOTS, [Y_TOP] * (dy + 1)])
    us = np.unique(ts)
    sx = np.concatenate([np.linspace(us[i], us[i + 1], 5)[:-1] for i in range(len(us) - 1)] + [[1.0]])
    uy = np.unique(ty)
    yx = np.concatenate([np.linspace(uy[i], uy[i + 1], 7)[:-1] for i in range(len(uy) - 1)] + [[Y_TOP]])
    X = np.stack([cm.cap_point(sx, np.full_like(sx, y), np.full_like(sx, float(side))) for y in yx], axis=1)
    Bs = BSpline.design_matrix(sx, ts, ds).toarray()
    By = BSpline.design_matrix(yx, ty, dy).toarray()
    ns, ny = Bs.shape[1], By.shape[1]
    C = np.zeros((ns, ny, 3))
    C[:, 0] = row440
    Xr = X - np.einsum("ai,ik,b->abk", Bs, C[:, 0], By[:, 0])
    C[:, 1:] = np.einsum("ia,abk,jb->ijk", np.linalg.pinv(Bs), Xr, np.linalg.pinv(By[:, 1:]))
    fitX = np.einsum("ai,ijk,bj->abk", Bs, C, By)
    err = np.linalg.norm(fitX - X, axis=-1)
    return C, ts, ty, ds, dy, float(err.max()), float(np.unravel_index(np.argmax(err), err.shape)[1])


def occ_surface_k(C, tu, tv, du, dv):
    ng, ns = C.shape[:2]
    poles = TColgp_Array2OfPnt(1, ng, 1, ns)
    for i in range(ng):
        for j in range(ns):
            poles.SetValue(i + 1, j + 1, gp_Pnt(*map(float, C[i, j])))

    def kn(t):
        u, c = np.unique(t, return_counts=True)
        K = TColStd_Array1OfReal(1, len(u)); M = TColStd_Array1OfInteger(1, len(u))
        for k, (a, b) in enumerate(zip(u, c)):
            K.SetValue(k + 1, float(a)); M.SetValue(k + 1, int(b))
        return K, M

    Ku, Mu = kn(tu); Kv, Mv = kn(tv)
    return Geom_BSplineSurface(poles, Ku, Kv, Mu, Mv, du, dv)


def trimmed_flat_face(surf, band_sY, g, tg):
    """Face on the continued-skin surface bounded by the attachment iso-line and the band-edge pcurve."""
    # band edge pcurve (s, Y) along g, reversed so the wire runs (0,440)->(1,440)->(back over the tip)
    Cb = lsq_curve(g, band_sY, tg, DEG)
    Cb[0] = (0.0, cm.CUT); Cb[-1] = (1.0, cm.CUT)
    c2 = occ_curve2d(Cb[::-1].copy(), (1 - tg[::-1]).copy(), DEG)
    line = occ_curve2d(np.array([[0.0, cm.CUT], [1.0, cm.CUT]]), np.array([0.0, 0.0, 1.0, 1.0]), 1)
    e1 = BRepBuilderAPI_MakeEdge(line, surf).Edge()
    e2 = BRepBuilderAPI_MakeEdge(c2, surf).Edge()
    BRepLib.BuildCurves3d_s(e1); BRepLib.BuildCurves3d_s(e2)
    mw = BRepBuilderAPI_MakeWire(e1, e2)
    face = BRepBuilderAPI_MakeFace(surf, mw.Wire(), True).Face()
    fx = ShapeFix_Face(face); fx.Perform()
    return fx.Face(), Cb


def load_parent():
    r = STEPControl_Reader()
    if r.ReadFile(os.path.join(cm.PARENT_DIR, "parent.step")) != IFSelect_RetDone:
        raise RuntimeError("cannot read parent.step")
    r.TransferRoots()
    return r.OneShape()


def faces_of(shape):
    out = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        out.append(TopoDS.Face_s(ex.Current()))
        ex.Next()
    return out


def volume(shape):
    pr = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, pr, 1e-9)
    return pr.Mass()


def main():
    t0 = time.time()
    g = np.linspace(0, 1, NG)
    p = cm.p_of_g(g)
    model_src = open(os.path.join(HERE, "cap_model.py"), "rb").read()
    model_src = model_src.split(b"# ---- sampling of the three faces")[0]      # the design part only
    key = hashlib.sha256(model_src + p.tobytes()).hexdigest()
    dz = None
    if os.path.exists(CACHE):
        with open(CACHE, "rb") as fh:
            k, d = pickle.load(fh)
        if k == key:
            dz = d
    if dz is None:
        dz = cm.design(p)
        with open(CACHE, "wb") as fh:
            pickle.dump((key, dz), fh)
    print("design %.1f s" % (time.time() - t0)); sys.stdout.flush()
    sl = np.linspace(-1, 1, NL)
    tf = np.linspace(0, 1, NF)
    XL, infoL = cm.lip_points(dz, sl)
    XU, _ = cm.flat_points(dz, tf, +1)
    XD, _ = cm.flat_points(dz, tf, -1)
    # the leading-edge pole (p = 0): every sample is the leading-edge point
    LE = cm.pa.skin(0.0, cm.CUT, +1)
    for X in (XL, XU, XD):
        X[g <= 0] = LE
    print("samples %.1f s" % (time.time() - t0)); sys.stdout.flush()
    tg = clamped_knots(NG_CP, DEG)
    tl = clamped_knots(NL_CP, DEG)
    tfk = clamped_knots(NF_CP, DEG)
    # shared boundary curves along g
    edgeU = lsq_curve(g, XL[:, -1], tg, DEG)      # lip sigma=+1 == upper flat t=0
    edgeD = lsq_curve(g, XL[:, 0], tg, DEG)       # lip sigma=-1 == lower flat t=0
    attU = lsq_curve(g, XU[:, -1], tg, DEG)       # upper attachment edge
    attD = lsq_curve(g, XD[:, -1], tg, DEG)       # lower attachment edge
    # rows at g = 0 (pole) and g = 1 (trailing edge)
    def row(Xrow, x, t):
        return lsq_curve(x, Xrow, t, DEG)
    lipg0 = np.repeat(LE[None], NL_CP, 0)
    lipg1 = row(XL[-1], (sl + 1) / 2, tl)
    lipg1[0], lipg1[-1] = edgeD[-1], edgeU[-1]
    fU0 = np.repeat(LE[None], NF_CP, 0); fD0 = fU0.copy()
    fU1 = np.repeat(XU[-1, :1], NF_CP, 0); fD1 = np.repeat(XD[-1, :1], NF_CP, 0)
    fU1[:] = edgeU[-1]; fD1[:] = edgeD[-1]
    CL, eL = lsq_surface(g, (sl + 1) / 2, XL, tg, tl, DEG, dict(g0=lipg0, g1=lipg1, s0=edgeD, s1=edgeU))
    parent = load_parent()
    srows = parent_skin_rows(parent)
    FU = flat_surface(+1, srows); FD = flat_surface(-1, srows)
    print("flat surfaces fit max error (mm): upper %.2e lower %.2e" % (FU[5], FD[5])); sys.stdout.flush()
    # band-edge chart points (s, Y) of each side
    bandU = np.stack([cm.s_of_a(dz["aE_up"]), dz["YE_up"]], -1); bandU[g <= 0] = (0.0, cm.CUT)
    bandD = np.stack([cm.s_of_a(dz["aE_lo"]), dz["YE_lo"]], -1); bandD[g <= 0] = (0.0, cm.CUT)
    bandU[:, 0] = np.clip(bandU[:, 0], 0, 1); bandD[:, 0] = np.clip(bandD[:, 0], 0, 1)
    eU = np.zeros((1, 1)); eD = np.zeros((1, 1))
    fit = {"lip_max_mm": float(eL.max()), "lip_worst_p": float(p[np.unravel_index(np.argmax(eL), eL.shape)[0]]),
           "upper_flat_surface_max_mm": FU[5], "lower_flat_surface_max_mm": FD[5],
           "poles": {"lip": list(CL.shape[:2]), "flat": list(FU[0].shape[:2])}}
    print("fit", json.dumps(fit)); sys.stdout.flush()
    with open(os.path.join(HERE, "_fit_poles.pkl"), "wb") as fh:
        pickle.dump(dict(CL=CL, tg=tg, tl=tl, g=g, p=p, sl=sl, XL=XL, wL=infoL["w"], bandU=bandU, bandD=bandD), fh)
    lip_face = BRepBuilderAPI_MakeFace(occ_surface(CL, tg, tl, DEG), 1e-7).Face()
    sU = occ_surface_k(FU[0], FU[1], FU[2], FU[3], FU[4]); sD = occ_surface_k(FD[0], FD[1], FD[2], FD[3], FD[4])
    fU, cbU = trimmed_flat_face(sU, bandU, g, tg)
    fD, cbD = trimmed_flat_face(sD, bandD, g, tg)
    right = [fU, lip_face, fD]
    tr = gp_Trsf(); tr.SetMirror(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)))
    left = [TopoDS.Face_s(BRepBuilderAPI_Transform(f, tr, True).Shape()) for f in right]
    pfaces = [f for f in faces_of(parent) if BRepAdaptor_Surface(f).GetType() != GeomAbs_Plane]
    sew = BRepBuilderAPI_Sewing(SEW_TOL)
    for f in pfaces + right + left:
        sew.Add(f)
    sew.Perform()
    sewed = sew.SewedShape()
    ex = TopExp_Explorer(sewed, TopAbs_SHELL)
    shell = TopoDS.Shell_s(ex.Current())
    ms = BRepBuilderAPI_MakeSolid(shell)
    solid = ms.Solid()
    sf = ShapeFix_Solid(solid); sf.Perform()
    solid = TopoDS.Solid_s(sf.Solid()) if sf.Solid().ShapeType() == TopAbs_SOLID else solid
    vol = volume(solid)
    if vol < 0:
        solid.Reverse(); vol = -vol
    ok = BRepCheck_Analyzer(solid).IsValid()
    nfree = sew.NbFreeEdges()
    print("solid valid", ok, "free edges", nfree, "faces", len(faces_of(solid)), "volume", vol); sys.stdout.flush()
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    for i in range(1, nfree + 1):
        e = sew.FreeEdge(i)
        ca = BRepAdaptor_Curve(e)
        u0, u1 = ca.FirstParameter(), ca.LastParameter()
        pts = [ca.Value(u0 + (u1 - u0) * t) for t in (0, 0.5, 1)]
        print("  free edge", [(round(q.X(), 4), round(q.Y(), 4), round(q.Z(), 4)) for q in pts])
    for nm, e in (("lip", eL),):
        em = e.max(axis=1)
        idx = np.argsort(em)[::-1][:6]
        print("  worst", nm, [(round(float(p[i]), 6), float(em[i]), int(np.argmax(e[i]))) for i in idx])
    for nm, e in (("lip", eL),):
        out = []
        for k in range(10):
            blk = e[int(k * NG / 10):int((k + 1) * NG / 10)]
            i, j = np.unravel_index(np.argmax(blk), blk.shape)
            out.append((float(blk[i, j]), int(j)))
        print("  fit by g decile", nm, out)
    print("  p at g deciles", [float(p[int(k * (NG - 1) / 10)]) for k in range(11)])
    Interface_Static.SetCVal_s("write.step.schema", "AP214IS")
    w = STEPControl_Writer()
    w.Transfer(solid, STEPControl_AsIs)
    if w.Write(OUT_STEP) != IFSelect_RetDone:
        raise RuntimeError("STEP write failed")
    fit.update({"brepcheck_valid_before_export": bool(ok), "sewing_free_edges": int(nfree),
                "volume_mm3": vol, "sew_tolerance_mm": SEW_TOL,
                "samples": {"g": NG, "lip": NL, "flat": NF}, "control_points": {"g": NG_CP, "lip": NL_CP, "flat": NF_CP},
                "degree": DEG, "step_bytes": os.path.getsize(OUT_STEP), "seconds": time.time() - t0})
    with open(OUT_FIT, "w") as fh:
        json.dump(fit, fh, indent=1)
    print("wrote", OUT_STEP, os.path.getsize(OUT_STEP), "bytes; %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
