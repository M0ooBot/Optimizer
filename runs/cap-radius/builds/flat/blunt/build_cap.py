"""Build wing.step: parent (from parent.step) plus a rounded cap at each tip.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/flat/blunt/build_cap.py

The cap surface is defined analytically in design.py. It is written as an exact
tensor-product B-spline: in the wrap direction w the functions 1, E(w), S(w) are
polynomials (degree 5 Bezier), in the chord direction s the parent's own edge
curve (iso-curve of the parent skin at |Y| = 350, with extra knots inserted) gives
X(s) and h(s) exactly, and the reach R(s) is interpolated at the Greville points
of the same knot vector. The cap's cut-edge rows of poles are the parent edge's
poles, so cap and parent share the same edge geometry.
"""
import json, os, sys, hashlib
import numpy as np
from scipy.interpolate import BSpline

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import design as D  # noqa: E402

from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCP.Interface import Interface_Static
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL
from OCP.TopoDS import TopoDS
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Plane
from OCP.Geom import Geom_BSplineSurface, Geom_BSplineCurve, Geom_Plane
from OCP.TColgp import TColgp_Array2OfPnt, TColgp_Array1OfPnt
from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
from OCP.gp import gp_Pnt, gp_Dir, gp_Ax2, gp_Trsf, gp_Pln, gp_Ax3
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire,
                                BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid, BRepBuilderAPI_Transform)
from OCP.ShapeFix import ShapeFix_Solid, ShapeFix_Shape
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp

PARENT_STEP = os.path.join(D.PARENT_DIR, 'parent.step')
OUT_STEP = os.path.join(HERE, 'wing.step')
OUT_JSON = os.path.join(HERE, 'build.json')
DEG_S = 5


def read_step(path):
    r = STEPControl_Reader()
    if r.ReadFile(path) != 1:
        raise RuntimeError(path)
    r.TransferRoots()
    return r.OneShape()


def faces(shape):
    out = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        out.append(TopoDS.Face_s(ex.Current()))
        ex.Next()
    return out


def bspline_knots_full(knots, mults):
    return np.repeat(np.asarray(knots, float), np.asarray(mults, int))


def bezier_coeffs(func, deg=5):
    """Bernstein coefficients on w in [-1, 1] of a polynomial of degree <= deg."""
    w = np.linspace(-1, 1, deg + 1)
    t = (w + 1) / 2
    from math import comb
    Bm = np.array([[comb(deg, k) * tt ** k * (1 - tt) ** (deg - k) for k in range(deg + 1)] for tt in t])
    return np.linalg.solve(Bm, func(w))


def main():
    parent = read_step(PARENT_STEP)
    skin_faces, te_face, end_faces = {}, None, []
    for f in faces(parent):
        ad = BRepAdaptor_Surface(f)
        if ad.GetType() == GeomAbs_Plane:
            end_faces.append(f)
            continue
        bs = ad.BSpline()
        if bs.UDegree() == 5:
            z = bs.Value(0.5, 0.0).Z()
            skin_faces['upper' if z > 0 else 'lower'] = (f, bs)
        else:
            te_face = f
    assert set(skin_faces) == {'upper', 'lower'} and te_face is not None and len(end_faces) == 2

    # ---- parent upper cut edge at Y = +350 in the parent chord parameter s ----
    bs_up = skin_faces['upper'][1]
    edge = bs_up.VIso(D.ATT)
    edge = Geom_BSplineCurve.DownCast(edge) if not isinstance(edge, Geom_BSplineCurve) else edge
    # extra knots: geometric towards the trailing edge (the reach R(s) varies on the scale of
    # the 2e-3 mm aft shift there) and a few towards the leading edge
    pk = [edge.Knot(i) for i in range(1, edge.NbKnots() + 1)]
    extra = list(1 - 4e-4 * 1.15 ** np.arange(0, 45))
    extra += list(np.geomspace(2e-3, 0.02, 4))
    extra = [k for k in extra if 0 < k < 1 and min(abs(k - q) for q in pk) > 2e-5]
    for k in sorted(extra):
        edge.InsertKnot(float(k), 1, 1e-12)
    nk = edge.NbKnots()
    knots = np.array([edge.Knot(i) for i in range(1, nk + 1)])
    mults = np.array([edge.Multiplicity(i) for i in range(1, nk + 1)])
    npl = edge.NbPoles()
    EP = np.array([[edge.Pole(i).X(), edge.Pole(i).Y(), edge.Pole(i).Z()] for i in range(1, npl + 1)])
    tfull = bspline_knots_full(knots, mults)

    # ---- reach R(s) interpolated at the Greville points ----
    grev = np.array([tfull[i + 1:i + 1 + DEG_S].mean() for i in range(npl)])
    R_at = D.reach(grev)[0]
    colloc = BSpline.design_matrix(grev, tfull, DEG_S).toarray()
    Rpol = np.linalg.solve(colloc, R_at)
    Rspl = BSpline(tfull, Rpol, DEG_S)
    sc = np.unique(np.concatenate([np.linspace(0, 1, 20001), 1 - np.geomspace(1e-9, 0.05, 4000), np.geomspace(1e-9, 0.05, 2000)]))
    fit_err = np.abs(Rspl(sc) - D.reach(sc)[0])
    # X and h of the edge B-spline against parent.py
    Xspl = BSpline(tfull, EP[:, 0], DEG_S); Hspl = BSpline(tfull, EP[:, 2], DEG_S)
    ref = D.PP.skin(sc, D.ATT, +1)
    edge_err = max(np.abs(Xspl(sc) - ref[:, 0]).max(), np.abs(Hspl(sc) - ref[:, 2]).max())

    # ---- wrap-direction Bezier coefficients ----
    one = bezier_coeffs(lambda w: np.ones_like(w))
    e = bezier_coeffs(lambda w: (1 - w * w) * (1 + D.A * w * w))
    sco = bezier_coeffs(lambda w: (15 * w - 10 * w ** 3 + 3 * w ** 5) / 8)
    assert abs(e[0]) < 1e-14 and abs(e[-1]) < 1e-14 and abs(sco[0] + 1) < 1e-14 and abs(sco[-1] - 1) < 1e-14
    e[0] = e[-1] = 0.0; sco[0] = -1.0; sco[-1] = 1.0

    def cap_surface(side):
        poles = TColgp_Array2OfPnt(1, npl, 1, 6)
        for i in range(npl):
            for j in range(6):
                X = EP[i, 0]
                Y = side * (D.ATT + Rpol[i] * e[j])
                Z = EP[i, 2] * sco[j]
                if i == 0:  # collapsed leading-edge row
                    X, Y, Z = EP[0, 0], side * D.ATT, 0.0
                poles.SetValue(i + 1, j + 1, gp_Pnt(float(X), float(Y), float(Z)))
        uk = TColStd_Array1OfReal(1, nk); um = TColStd_Array1OfInteger(1, nk)
        for i in range(nk):
            uk.SetValue(i + 1, float(knots[i])); um.SetValue(i + 1, int(mults[i]))
        vk = TColStd_Array1OfReal(1, 2); vm = TColStd_Array1OfInteger(1, 2)
        vk.SetValue(1, -1.0); vk.SetValue(2, 1.0); vm.SetValue(1, 6); vm.SetValue(2, 6)
        return Geom_BSplineSurface(poles, uk, vk, um, vm, DEG_S, 5)

    built = []
    for side in (+1, -1):
        srf = cap_surface(side)
        cap = BRepBuilderAPI_MakeFace(srf, 0.0, 1.0, -1.0, 1.0, 1e-7).Face()
        # closure of the flat trailing-edge face: plane X = 120 bounded by the TE face's end
        # segment and the cap's s = 1 iso-curve
        iso = srf.UIso(1.0)
        e1 = BRepBuilderAPI_MakeEdge(iso).Edge()
        pa = gp_Pnt(D.CHORD, side * D.ATT, -D.H_TE); pb = gp_Pnt(D.CHORD, side * D.ATT, D.H_TE)
        e2 = BRepBuilderAPI_MakeEdge(pb, pa).Edge()
        wire = BRepBuilderAPI_MakeWire(e1, e2).Wire()
        pln = gp_Pln(gp_Pnt(D.CHORD, 0, 0), gp_Dir(1, 0, 0))
        clo = BRepBuilderAPI_MakeFace(pln, wire, True).Face()
        built += [cap, clo]

    sew = BRepBuilderAPI_Sewing(1e-6)
    for f in [skin_faces['upper'][0], skin_faces['lower'][0], te_face] + built:
        sew.Add(f)
    sew.Perform()
    sewn = sew.SewedShape()
    ex = TopExp_Explorer(sewn, TopAbs_SHELL)
    shells = []
    while ex.More():
        shells.append(TopoDS.Shell_s(ex.Current())); ex.Next()
    assert len(shells) == 1, len(shells)
    ms = BRepBuilderAPI_MakeSolid(shells[0])
    fixer = ShapeFix_Solid(ms.Solid()); fixer.Perform()
    solid = fixer.Solid()
    ok = BRepCheck_Analyzer(solid).IsValid()
    props = GProp_GProps(); BRepGProp.VolumeProperties_s(solid, props, 1e-9)

    Interface_Static.SetCVal_s('write.step.schema', 'AP214IS')
    Interface_Static.SetCVal_s('write.step.unit', 'MM')
    w = STEPControl_Writer(); w.Transfer(solid, STEPControl_AsIs); w.Write(OUT_STEP)
    sha = hashlib.sha256(open(OUT_STEP, 'rb').read()).hexdigest()
    info = {
        'design': {'A': D.A, 'K': D.K, 'law': D.LAW, 'eps0_mm': D.EPS0, 'eps1_mm': D.EPS1,
                   'E(w)': '(1-w^2)(1+A w^2)', 'S(w)': '(15w-10w^3+3w^5)/8'},
        'cap_poles_s_by_w': [npl, 6], 'cap_knots_s': knots.tolist(), 'cap_mults_s': mults.tolist(),
        'reach_fit_max_error_mm': float(fit_err.max()),
        'reach_fit_max_error_at_s': float(sc[np.argmax(fit_err)]),
        'edge_vs_parent_py_max_mm': float(edge_err),
        'sewing_tolerance_mm': 1e-6,
        'faces': len(faces(solid)), 'valid': bool(ok), 'volume_mm3': props.Mass(),
        'step_sha256': sha, 'step_bytes': os.path.getsize(OUT_STEP),
    }
    json.dump(info, open(OUT_JSON, 'w'), indent=1)
    print(json.dumps({k: v for k, v in info.items() if k not in ('cap_knots_s', 'cap_mults_s')}, indent=1))


if __name__ == '__main__':
    main()
