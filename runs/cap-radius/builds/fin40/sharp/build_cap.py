"""Build the rounded tip cap for the fin40 / sharp parent and write wing.step.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/fin40/sharp/build_cap.py

Construction (see README.md):
  profiles  : lower cut edge -> 3D perimeter (locked outline, Y = 0) -> upper cut edge,
              the cut edges taken pole-for-pole from the parent skins at Z = 240
  guides    : one planar wrap guide per chord station u; in side view it runs on a
              straight line from the cut point A(u) to the perimeter point L(u); its
              half-thickness follows the cap thickness law
                  Y^2 = 0.0144 * [ a(Z)^2 (1 - h^3) - (X - m(Z))^2 ],  h = (Z - 240)/40
              which is the parent's own law (12 % ellipse, Y^2 = 0.0144 (X - X_le)(X_te - X))
              with the chord ends moved onto the locked outline. It differs from the
              parent by an O(h^3) term, so it matches position, tangent plane and
              curvature at the attachment, and it closes with a finite lip radius.
  surface   : one B-spline face, u along the chord (parent knots, collapsed at the LE
              and TE points), v from the lower cut edge (-1) through the lip (0) to the
              upper cut edge (+1); boundary rows are the parent's own cut-edge poles,
              interior poles are a least-squares fit of the analytic guide family.
"""
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import numpy as np
from scipy.interpolate import BSpline

from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCP.Interface import Interface_Static
from OCP.IFSelect import IFSelect_RetDone
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_SHELL
from OCP.TopoDS import TopoDS
from OCP.BRep import BRep_Tool
from OCP.GeomAdaptor import GeomAdaptor_Surface
from OCP.GeomAbs import GeomAbs_Plane, GeomAbs_BSplineSurface
from OCP.Geom import Geom_BSplineSurface
from OCP.TColgp import TColgp_Array2OfPnt
from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
from OCP.gp import gp_Pnt
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.ShapeFix import ShapeFix_Solid
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
PARENT_DIR = ROOT / "runs/cap-radius/parents/fin40/sharp"
sys.path.insert(0, str(PARENT_DIR))
import parent  # noqa: E402  (exact parent definition, locked outline)

K = parent.THICKNESS_RATIO          # 0.12: full thickness / chord
ZA = parent.ATTACHMENT              # 240
REACH = parent.REACH                # 40

# --------------------------------------------------------------------------
# Design parameters (the only free choices of this cap)
# --------------------------------------------------------------------------
ALPHA = -0.1        # perimeter mapping, see w_of()
EPS = 0.02          # small linear share of the perimeter mapping (keeps the corner parametrisation polar)
V_DEGREE = 5
V_SPANS = 24        # uniform B-spline spans in v over [-1, 1]
U_INSERT = 1        # knots inserted per parent u-span (exact refinement)
N_FIT_U, N_FIT_V = 721, 601


def m_of(z):
    return 80.0 + 27.5 * z / ZA


def a_of(z):
    return 80.0 - 37.5 * z / ZA


def xi_of(u):
    """Parent chord fraction at chord parameter u (theta = pi u)."""
    return 0.5 * (1.0 - np.cos(np.pi * u))


def w_of(u):
    """Perimeter parameter of the guide at chord station u: -1 = LE end, 0 = crown, +1 = TE end.

    w(u) = -(1 - EPS) cos(pi q) + EPS (2u - 1),  q = u + ALPHA sin(2 pi u).
    Away from the corners the perimeter point trails the cut point roughly as u^2
    (like the cut point itself), so each side-view guide line meets the perimeter
    almost at a right angle and the guide plane is close to the lip's normal plane.
    The small linear share EPS makes the map polar (perimeter speed > 0) inside the
    last few percent of chord, which keeps the B-spline well conditioned at the
    collapsed corners.
    """
    q = u + ALPHA * np.sin(2.0 * np.pi * u)
    return -(1.0 - EPS) * np.cos(np.pi * q) + EPS * (2.0 * u - 1.0)


def perimeter(w):
    """Locked outline point at perimeter parameter w (h = 1 - w^2). Same curve as
    parent.tip_outline, written in w so it is regular through the crown."""
    h = 1.0 - w**2
    z = ZA + REACH * h
    x = m_of(z) + a_of(z) * w * np.sqrt(1.0 + h + h * h)
    return x, z


def thickness_G(x, z):
    """Cap thickness law: half-thickness = K * sqrt(G)."""
    h = (z - ZA) / REACH
    return a_of(z) ** 2 * (1.0 - h**3) - (x - m_of(z)) ** 2


def cap_point(u, v):
    """Analytic cap S(u, v). u in [0,1] chord station, v in [-1,1] (lower cut, lip, upper cut)."""
    u, v = np.broadcast_arrays(np.asarray(u, float), np.asarray(v, float))
    xa = 65.0 + 85.0 * xi_of(u)
    xl, zl = perimeter(w_of(u))
    r = 1.0 - v**2                       # 0 on the cut edge, 1 on the perimeter
    x = xa + r * (xl - xa)
    z = ZA + r * (zl - ZA)
    g = np.clip(thickness_G(x, z), 0.0, None)
    y = np.sign(v) * K * np.sqrt(g)
    return np.stack([x, y, z], -1)


# --------------------------------------------------------------------------
# Parent
# --------------------------------------------------------------------------

def read_step(path):
    r = STEPControl_Reader()
    if r.ReadFile(str(path)) != IFSelect_RetDone:
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
    skins, planes = {}, {}
    for f in faces_of(shape):
        s = BRep_Tool.Surface_s(f)
        ad = GeomAdaptor_Surface(s)
        if ad.GetType() == GeomAbs_BSplineSurface:
            p = s.Value(0.5, 120.0)
            skins[+1 if p.Y() > 0 else -1] = (f, ad.BSpline())
        elif ad.GetType() == GeomAbs_Plane:
            z = ad.Plane().Location().Z()
            planes["root" if abs(z) < 1e-6 else "top"] = f
    return skins, planes


def cut_row(bs):
    """Poles of the skin's iso-curve at Z = 240 (last v pole row: the surface is Bezier in v)."""
    assert abs(bs.LastVKnotIndex() and bs.VKnot(bs.LastVKnotIndex()) - ZA) < 1e-12
    n = bs.NbUPoles()
    nv = bs.NbVPoles()
    return np.array([[bs.Pole(i, nv).X(), bs.Pole(i, nv).Y(), bs.Pole(i, nv).Z()] for i in range(1, n + 1)])


def u_knots(bs):
    k = [bs.UKnot(i) for i in range(1, bs.NbUKnots() + 1)]
    m = [bs.UMultiplicity(i) for i in range(1, bs.NbUKnots() + 1)]
    return np.array(k), np.array(m), bs.UDegree()


def refine_curve(poles, knots, mults, deg, inserts):
    """Exact knot refinement of a non-rational B-spline curve (Boehm insertion)."""
    t = np.repeat(knots, mults).astype(float)
    P = poles.copy()
    new = []
    for a, b in zip(knots[:-1], knots[1:]):
        new += list(a + (b - a) * np.arange(1, inserts + 1) / (inserts + 1))
    for x in new:
        k = np.searchsorted(t, x, side="right") - 1
        Q = np.empty((len(P) + 1, 3))
        for i in range(len(P) + 1):
            if i <= k - deg:
                Q[i] = P[i]
            elif i > k:
                Q[i] = P[i - 1]
            else:
                al = (x - t[i]) / (t[i + deg] - t[i])
                Q[i] = (1 - al) * P[i - 1] + al * P[i]
        P = Q
        t = np.insert(t, k + 1, x)
    return P, t


# --------------------------------------------------------------------------
# Fit
# --------------------------------------------------------------------------

def fit_cap(row_lo, row_up, tu, pu):
    pv = V_DEGREE
    inner = np.linspace(-1, 1, V_SPANS + 1)
    tv = np.concatenate([[-1.0] * pv, inner, [1.0] * pv])
    nu = len(tu) - pu - 1
    nv = len(tv) - pv - 1
    # sample: cosine-clustered in u (like the parent), uniform in v
    us = 0.5 * (1 - np.cos(np.pi * np.linspace(0, 1, N_FIT_U)))
    us = 0.5 * (us + np.linspace(0, 1, N_FIT_U))
    vs = np.linspace(-1, 1, N_FIT_V)
    D = cap_point(us[:, None], vs[None, :])
    Bu = BSpline.design_matrix(us, tu, pu).toarray()
    Bv = BSpline.design_matrix(vs, tv, pv).toarray()
    C = np.zeros((nu, nv, 3))
    C[:, 0] = row_lo
    C[:, -1] = row_up
    C[0, :] = row_lo[0]
    C[-1, :] = row_lo[-1]
    Bui, Bvi = Bu[:, 1:-1], Bv[:, 1:-1]
    Pu, Pv = np.linalg.pinv(Bui), np.linalg.pinv(Bvi)
    for c in range(3):
        Dp = D[..., c] - Bu @ C[..., c] @ Bv.T
        C[1:-1, 1:-1, c] = Pu @ Dp @ Pv.T
    F = np.einsum("ia,abc,jb->ijc", Bu, C, Bv)
    err = np.linalg.norm(F - D, axis=-1)
    return C, tv, pv, float(err.max())


def make_bspline(C, tu, pu, tv, pv):
    nu, nv = C.shape[:2]
    poles = TColgp_Array2OfPnt(1, nu, 1, nv)
    for i in range(nu):
        for j in range(nv):
            poles.SetValue(i + 1, j + 1, gp_Pnt(*C[i, j]))

    def km(t):
        k, m = np.unique(np.round(t, 15), return_counts=True)
        ka = TColStd_Array1OfReal(1, len(k))
        ma = TColStd_Array1OfInteger(1, len(k))
        for i, (a, b) in enumerate(zip(k, m)):
            ka.SetValue(i + 1, float(a))
            ma.SetValue(i + 1, int(b))
        return ka, ma
    uk, um = km(tu)
    vk, vm = km(tv)
    return Geom_BSplineSurface(poles, uk, vk, um, vm, pu, pv)


def main():
    shape = read_step(PARENT_DIR / "parent.step")
    skins, planes = classify_parent(shape)
    up_face, up_bs = skins[+1]
    lo_face, lo_bs = skins[-1]
    k, m, pu = u_knots(up_bs)
    k2, m2, _ = u_knots(lo_bs)
    assert np.allclose(k, k2) and np.all(m == m2)
    row_up, tu = refine_curve(cut_row(up_bs), k, m, pu, U_INSERT)
    row_lo, _ = refine_curve(cut_row(lo_bs), k, m, pu, U_INSERT)
    C, tv, pv, fit_err = fit_cap(row_lo, row_up, tu, pu)
    surf = make_bspline(C, tu, pu, tv, pv)
    cap = BRepBuilderAPI_MakeFace(surf, 1e-7).Face()

    sew = BRepBuilderAPI_Sewing(1e-6)
    for f in (up_face, lo_face, planes["root"], cap):
        sew.Add(f)
    sew.Perform()
    sewn = sew.SewedShape()
    ex = TopExp_Explorer(sewn, TopAbs_SHELL)
    shell = TopoDS.Shell_s(ex.Current())
    solid = BRepBuilderAPI_MakeSolid(shell).Solid()
    fix = ShapeFix_Solid(solid)
    fix.Perform()
    solid = fix.Solid()
    ok = BRepCheck_Analyzer(solid).IsValid()

    Interface_Static.SetCVal_s("write.step.schema", "AP214")
    Interface_Static.SetCVal_s("write.step.unit", "MM")
    w = STEPControl_Writer()
    w.Transfer(solid, STEPControl_AsIs)
    out = HERE / "wing.step"
    if w.Write(str(out)) != IFSelect_RetDone:
        raise RuntimeError("STEP write failed")

    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(solid, props, 1e-9)
    info = {
        "fit_max_error_mm": fit_err,
        "cap_poles": [int(C.shape[0]), int(C.shape[1])],
        "cap_degrees": [int(pu), int(pv)],
        "alpha": ALPHA, "eps": EPS, "v_spans": V_SPANS, "u_insert": U_INSERT,
        "sewing_free_edges": sew.NbFreeEdges(),
        "sewing_multiple_edges": sew.NbMultipleEdges(),
        "brepcheck_valid_before_export": bool(ok),
        "volume_mm3": props.Mass(),
        "step_bytes": out.stat().st_size,
    }
    (HERE / "build.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
