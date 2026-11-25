"""Build wing.step: the locked parent (parent.step skins, unchanged) plus a rounded cap at each tip.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/twisted/sharp/build_wing.py

Each cap is one B-spline face fitted to the analytic cap of cap_model.py:
  u direction  the parent's own chord-parameter basis (degree 5, the parent skin's knot vector), so
               the cap's first and last pole rows are copied from the parent skin's attachment pole
               rows and the cap shares the parent's attachment edges exactly;
  w direction  degree 5, upper attachment (w = 0) -> rim (w = 1) -> lower attachment (w = 2).
Fitting: interpolation in u at the Greville abscissae of the parent basis; constrained least squares
in w with value, first and second derivative imposed at both attachments (G2 to the parent).
"""
import json
import pathlib
import sys
import warnings

import numpy as np
from scipy.interpolate import BSpline

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cap_model as C  # noqa: E402

from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCP.Interface import Interface_Static
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL
from OCP.TopoDS import TopoDS
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_BSplineSurface
from OCP.Geom import Geom_BSplineSurface
from OCP.TColgp import TColgp_Array2OfPnt
from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
from OCP.gp import gp_Pnt
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid
from OCP.ShapeFix import ShapeFix_Solid
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps

warnings.filterwarnings("ignore", category=RuntimeWarning)
PARENT_STEP = C.PARENT_DIR / "parent.step"
OUT_STEP = HERE / "wing.step"
DEG_W = 5
# w knots (upper half; mirrored for the lower half): spacing W_COARSE up to W_RAMP[0], W_FINE beyond
# W_RAMP[1], cosine ramp between.  The roll's onset moves from w = 0 (forward) to w = 0.71 (aft), and
# near the trailing edge the guides run almost parallel to the u lines, which amplifies w-fit error
# in the chordwise curvature, so the rim half is knotted densely.
W_COARSE, W_FINE, W_RAMP = 0.04, 0.007, (0.35, 0.65)


def read_step(path):
    r = STEPControl_Reader()
    if r.ReadFile(str(path)) != 1:
        raise RuntimeError(f"cannot read {path}")
    r.TransferRoots()
    return r.OneShape()


def faces_of(shape):
    out = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        out.append(TopoDS.Face_s(ex.Current()))
        ex.Next()
    return out


def classify_parent(shape):
    """Return {('right'|'left', 'upper'|'lower'): (face, Geom_BSplineSurface)} and the planar faces."""
    skins, planes = {}, []
    for f in faces_of(shape):
        a = BRepAdaptor_Surface(f)
        if a.GetType() != GeomAbs_BSplineSurface:
            planes.append(f)
            continue
        s = a.BSpline()
        tip = s.Value(0.5, s.VKnot(s.NbVKnots()))
        mid_up = s.Value(0.5, 0.0)
        side = "right" if tip.Y() > 0 else "left"
        skin = "upper" if mid_up.Z() > 3.0 else "lower"
        skins[(side, skin)] = (f, s)
    assert len(skins) == 4 and len(planes) == 2, (len(skins), len(planes))
    return skins, planes


def u_basis(s):
    knots = np.array([s.UKnot(i) for i in range(1, s.NbUKnots() + 1)])
    mults = np.array([s.UMultiplicity(i) for i in range(1, s.NbUKnots() + 1)])
    flat = np.repeat(knots, mults)
    return s.UDegree(), knots, mults, flat


def poles_of(s):
    nu, nv = s.NbUPoles(), s.NbVPoles()
    P = np.empty((nu, nv, 3))
    for i in range(nu):
        for j in range(nv):
            p = s.Pole(i + 1, j + 1)
            P[i, j] = (p.X(), p.Y(), p.Z())
    return P


def w_knots():
    ks, w = [], 0.0
    while True:
        t = min(max((w - W_RAMP[0]) / (W_RAMP[1] - W_RAMP[0]), 0.0), 1.0)
        w += W_COARSE + (W_FINE - W_COARSE) * (0.5 - 0.5 * np.cos(np.pi * t))
        if w >= 1 - 1e-9:
            break
        ks.append(w)
    half = np.array(ks) / w
    inner = np.concatenate([half, [1.0], 2 - half[::-1]])
    knots = np.concatenate([[0.0], inner, [2.0]])
    mults = np.array([DEG_W + 1] + [1] * len(inner) + [DEG_W + 1])
    return knots, mults, np.repeat(knots, mults)


def w_samples():
    s = np.linspace(0, 1, 1601)
    upper = np.unique(np.concatenate([s, 1 - (1 - s) ** 2, 1 - (1 - s) ** 3]))
    return np.unique(np.concatenate([upper, 2 - upper]))


def fit_cap():
    """Poles (nu, nw, 3) of the right-tip cap in the parent's u basis."""
    shape = read_step(PARENT_STEP)
    skins, _ = classify_parent(shape)
    su = skins[("right", "upper")][1]
    sl = skins[("right", "lower")][1]
    deg_u, uk, um, uflat = u_basis(su)
    assert np.allclose(uflat, u_basis(sl)[3]) and su.VKnot(2) == 440.0
    nu = su.NbUPoles()
    greville = np.array([uflat[i + 1:i + deg_u + 1].mean() for i in range(nu)])
    greville[0], greville[-1] = 0.0, 1.0

    wk, wm, wflat = w_knots()
    nw = len(wflat) - DEG_W - 1
    ws = w_samples()
    Aw = BSpline.design_matrix(ws, wflat, DEG_W).toarray()
    eye = np.eye(nw)
    cons = []
    for w0 in (0.0, 2.0):
        for der in (0, 1, 2):
            cons.append([BSpline(wflat, eye[k], DEG_W).derivative(der)(w0) if der else
                         BSpline(wflat, eye[k], DEG_W)(w0) for k in range(nw)])
    Cm = np.array(cons)
    K = np.block([[Aw.T @ Aw, Cm.T], [Cm, np.zeros((6, 6))]])
    Kinv = np.linalg.inv(K)

    # w-direction fit at each Greville abscissa
    Cw = np.empty((nu, nw, 3))
    for i, u in enumerate(greville):
        uu = np.full_like(ws, u)
        D = C.cap(uu, ws)
        rhs_c = []
        for w0 in (0.0, 2.0):
            j = C.cap_jet(np.array([u]), np.array([w0]))
            rhs_c += [j["point"][0], j["d_w"][0], j["d_w_w"][0]]
        rhs = np.vstack([Aw.T @ D, np.array(rhs_c)])
        Cw[i] = (Kinv @ rhs)[:nw]
    # u-direction interpolation in the parent's basis
    Au = BSpline.design_matrix(greville, uflat, deg_u).toarray()
    Pc = np.linalg.solve(Au, Cw.reshape(nu, -1)).reshape(nu, nw, 3)

    # exact shared boundary with the parent and exact collapsed ends
    Pu, Pl = poles_of(su), poles_of(sl)
    report = {
        "boundary_row_change_upper_mm": float(np.abs(Pc[:, 0] - Pu[:, -1]).max()),
        "boundary_row_change_lower_mm": float(np.abs(Pc[:, -1] - Pl[:, -1]).max()),
    }
    Pc[:, 0] = Pu[:, -1]
    Pc[:, -1] = Pl[:, -1]
    report["le_column_change_mm"] = float(np.abs(Pc[0] - Pu[0, -1]).max())
    report["te_column_change_mm"] = float(np.abs(Pc[-1] - Pu[-1, -1]).max())
    Pc[0, :] = Pu[0, -1]
    Pc[-1, :] = Pu[-1, -1]
    basis = dict(deg_u=deg_u, uk=uk, um=um, deg_w=DEG_W, wk=wk, wm=wm)
    return Pc, basis, report, skins, shape


def make_surface(Pc, b):
    nu, nw, _ = Pc.shape
    poles = TColgp_Array2OfPnt(1, nu, 1, nw)
    for i in range(nu):
        for j in range(nw):
            poles.SetValue(i + 1, j + 1, gp_Pnt(*Pc[i, j]))

    def arr(vals, cls):
        a = cls(1, len(vals))
        for k, v in enumerate(vals):
            a.SetValue(k + 1, type(a.Value(1))(v) if False else v)
        return a
    uk = TColStd_Array1OfReal(1, len(b["uk"]))
    um = TColStd_Array1OfInteger(1, len(b["um"]))
    for k, (kv, mv) in enumerate(zip(b["uk"], b["um"])):
        uk.SetValue(k + 1, float(kv)); um.SetValue(k + 1, int(mv))
    wk = TColStd_Array1OfReal(1, len(b["wk"]))
    wm = TColStd_Array1OfInteger(1, len(b["wm"]))
    for k, (kv, mv) in enumerate(zip(b["wk"], b["wm"])):
        wk.SetValue(k + 1, float(kv)); wm.SetValue(k + 1, int(mv))
    return Geom_BSplineSurface(poles, uk, wk, um, wm, b["deg_u"], b["deg_w"])


def fit_quality(surf):
    """Sampled distance and normal-curvature fidelity of the B-spline cap against the analytic cap."""
    from OCP.GeomLProp import GeomLProp_SLProps
    p = np.concatenate([np.linspace(0.002, 0.05, 25), np.linspace(0.05, 0.95, 181), np.linspace(0.95, 0.998, 25)])
    w = np.concatenate([np.linspace(0, 0.9, 91), np.linspace(0.9, 1.1, 201), np.linspace(1.1, 2, 91)])
    PP, WW = np.meshgrid(p, w, indexing="ij")
    ana = C.principal(C.cap_jet(PP, WW))["kmax"]
    pts = C.cap(PP, WW)
    kb = np.zeros_like(ana)
    dist = np.zeros_like(ana)
    for i in range(len(p)):
        for j in range(len(w)):
            props = GeomLProp_SLProps(surf, p[i], w[j], 2, 1e-9)
            q = props.Value()
            dist[i, j] = np.linalg.norm(np.array([q.X(), q.Y(), q.Z()]) - pts[i, j])
            kb[i, j] = max(abs(props.MaxCurvature()), abs(props.MinCurvature()))
    # error relative to the tightest analytic curvature of the same guide
    rel = np.abs(kb - ana) / ana.max(axis=1, keepdims=True)
    per_p = rel.max(axis=1)
    out = {"max_distance_to_analytic_cap_mm": float(dist.max()),
           "max_kmax_error_relative_to_guide_peak": float(per_p.max()),
           "at_p": float(p[per_p.argmax()]),
           "grid": f"{len(p)} x {len(w)} (p 0.002..0.998, w 0..2, dense through the lip)"}
    for lim in (0.9, 0.95, 0.98):
        out[f"max_kmax_error_relative_to_guide_peak_p_le_{lim}"] = float(per_p[p <= lim].max())
    rb, ra = 1 / kb.max(axis=1), 1 / ana.max(axis=1)
    out["max_relative_error_of_guide_min_radius"] = float(np.abs(rb / ra - 1).max())
    return out


def main():
    Pc, basis, report, skins, parent_shape = fit_cap()
    right = make_surface(Pc, basis)
    PcL = Pc.copy()
    PcL[..., 1] *= -1
    PlU, PlL = poles_of(skins[("left", "upper")][1]), poles_of(skins[("left", "lower")][1])
    report["left_boundary_mirror_mismatch_mm"] = float(max(np.abs(PcL[:, 0] - PlU[:, -1]).max(),
                                                         np.abs(PcL[:, -1] - PlL[:, -1]).max()))
    PcL[:, 0], PcL[:, -1] = PlU[:, -1], PlL[:, -1]
    left = make_surface(PcL, basis)
    report["fit_quality_right_cap"] = fit_quality(right)

    sew = BRepBuilderAPI_Sewing(1e-6)
    for key in [("right", "upper"), ("right", "lower"), ("left", "upper"), ("left", "lower")]:
        sew.Add(skins[key][0])
    for s in (right, left):
        mf = BRepBuilderAPI_MakeFace(s, 1e-7)
        sew.Add(mf.Face())
    sew.Perform()
    sewn = sew.SewedShape()
    ex = TopExp_Explorer(sewn, TopAbs_SHELL)
    shell = TopoDS.Shell_s(ex.Current())
    ms = BRepBuilderAPI_MakeSolid(shell)
    fix = ShapeFix_Solid(ms.Solid())
    fix.Perform()
    solid = fix.Solid()
    report["sewing"] = {"free_edges": sew.NbFreeEdges(), "multiple_edges": sew.NbMultipleEdges(),
                        "degenerated_shapes": sew.NbDegeneratedShapes()}
    report["brepcheck_valid"] = bool(BRepCheck_Analyzer(solid).IsValid())
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(solid, props, 1e-9)
    report["volume_mm3"] = props.Mass()

    Interface_Static.SetCVal_s("write.step.schema", "AP214")
    Interface_Static.SetCVal_s("write.step.unit", "MM")
    w = STEPControl_Writer()
    w.Transfer(solid, STEPControl_AsIs)
    if w.Write(str(OUT_STEP)) != 1:
        raise RuntimeError("STEP write failed")
    report["step_bytes"] = OUT_STEP.stat().st_size
    report["cap_poles"] = list(Pc.shape[:2])
    report["design"] = {"A": C.A_SHAPE, "B": C.B_SHAPE, "C": C.C_SHAPE, "K": C.K_FLAT, "E_TE": C.E_TE,
                        "P_A": C.P_A, "P_B": C.P_B, "w_knots": {"coarse": W_COARSE, "fine": W_FINE, "ramp": W_RAMP}}
    (HERE / "build.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
