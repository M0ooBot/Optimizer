"""Build the fin5/blunt rounded tip cap and the whole-fin solid.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/fin5/blunt/build_cap.py

Construction (see README.md):
  1. Analytic cap surface (cap_analytic.py): Y^2 = G_par * psi, G2 to both skins,
     lip exactly on the locked outline.
  2. Wrap sections ("guides"): for each attachment chord parameter s the cap is cut by
     the plane through the Y axis direction and the side-view line from the attachment
     point A(s) to its perimeter point P(s) (equal arc-length fractions).
  3. One B-spline face, u = s (the parent's own chord parameter and knots, refined at
     the ends), t = normalised arc length around each wrap section
     (t = 0 on the -Y cut edge, 0.5 on the lip, 1 on the +Y cut edge).
     Boundary poles are fixed: the parent's top pole rows (exact shared edges), the LE
     point (collapsed u = 0 edge) and the trailing-edge segment (u = 1).
     Interior poles: one sparse least-squares solve against exact analytic points plus
     heavily weighted analytic first and second cross-boundary derivatives at both cut
     edges (tangent and curvature continuity with the skins).
  4. Replace the parent's flat attachment face by the cap face, sew, make the solid,
     export wing.step.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cap_analytic as ca  # noqa: E402

from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs  # noqa: E402
from OCP.TopExp import TopExp_Explorer  # noqa: E402
from OCP.TopAbs import TopAbs_FACE  # noqa: E402
from OCP.BRep import BRep_Tool  # noqa: E402
from OCP.TopoDS import TopoDS  # noqa: E402
from OCP.GeomAdaptor import GeomAdaptor_Surface  # noqa: E402
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Plane  # noqa: E402
from OCP.Geom import Geom_BSplineSurface, Geom_BSplineCurve  # noqa: E402
from OCP.TColgp import TColgp_Array2OfPnt, TColgp_Array1OfPnt  # noqa: E402
from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger  # noqa: E402
from OCP.gp import gp_Pnt  # noqa: E402

DEG = 5
PARENT_STEP = ca.PARENT_DIR / "parent.step"
OUT_STEP = HERE / "wing.step"


# --------------------------------------------------------------------------- parent
def read_step(path):
    r = STEPControl_Reader()
    if r.ReadFile(str(path)) != 1:
        raise RuntimeError(f"cannot read {path}")
    r.TransferRoots()
    return r.OneShape()


def faces_of(shape):
    out = []
    e = TopExp_Explorer(shape, TopAbs_FACE)
    while e.More():
        out.append(TopoDS.Face_s(e.Current()))
        e.Next()
    return out


def classify_parent(shape):
    """Return dict: skin+ / skin- (BSpline faces), te (ruled BSpline), root, top (planes)."""
    d = {}
    for f in faces_of(shape):
        s = BRep_Tool.Surface_s(f)
        a = GeomAdaptor_Surface(s)
        if a.GetType() == GeomAbs_Plane:
            z = a.Plane().Location().Z()
            d["top" if abs(z - ca.Z0) < 1e-6 else "root"] = f
        elif a.GetType() == GeomAbs_BSplineSurface:
            b = a.BSpline()
            if b.NbVPoles() == 2 or b.VDegree() == 1:
                d["te"] = f
            else:
                y = b.Value(0.5, 120.0).Y()
                d["skin+" if y > 0 else "skin-"] = f
    assert set(d) == {"skin+", "skin-", "te", "root", "top"}, d.keys()
    return d


def top_row_curve(face):
    """The skin's iso-curve at the attachment (last V pole row) as a Geom_BSplineCurve."""
    b = GeomAdaptor_Surface(BRep_Tool.Surface_s(face)).BSpline()
    assert abs(b.VKnot(b.NbVKnots()) - ca.Z0) < 1e-9
    nu = b.NbUPoles()
    poles = TColgp_Array1OfPnt(1, nu)
    for i in range(1, nu + 1):
        poles.SetValue(i, b.Pole(i, b.NbVPoles()))
    kn = TColStd_Array1OfReal(1, b.NbUKnots())
    mu = TColStd_Array1OfInteger(1, b.NbUKnots())
    for i in range(1, b.NbUKnots() + 1):
        kn.SetValue(i, b.UKnot(i))
        mu.SetValue(i, b.UMultiplicity(i))
    return Geom_BSplineCurve(poles, kn, mu, b.UDegree())


def curve_data(c):
    poles = np.array([[c.Pole(i).X(), c.Pole(i).Y(), c.Pole(i).Z()] for i in range(1, c.NbPoles() + 1)])
    flat = []
    for i in range(1, c.NbKnots() + 1):
        flat += [c.Knot(i)] * c.Multiplicity(i)
    return poles, np.array(flat)


# --------------------------------------------------------------------------- sampling
def section_points(s, t, want_length=False):
    """Exact cap points on the wrap section at chord parameter s, at normalised arc length t."""
    t = np.asarray(t, float)
    if s <= 0.0:
        out = np.tile(parent_le(), (t.size, 1))
        return (out, 0.0) if want_length else out
    if s >= 1.0:
        out = te_segment(t)
        return (out, 0.2) if want_length else out
    local = s > TE_LOCAL_FROM
    if local:
        tau, P, A, Pl, (q, z, qA, qP, zP) = ca.section_te(s)
        dX = -(10.0 / ca.Z0) * np.diff(z) - np.diff(q)
        L = np.concatenate([[0.0], np.cumsum(np.sqrt(dX**2 + np.diff(P[:, 1])**2 + np.diff(z)**2))])
    else:
        tau, P, A, Pl = ca.section(s)
        L = ca._cumlen(P)
    Lh = L[-1]
    out = np.empty((t.size, 3))
    for side, mask in ((+1, t >= 0.5), (-1, t < 0.5)):
        tt = t[mask] if side > 0 else 1.0 - t[mask]
        dist_from_A = Lh - (tt - 0.5) * 2.0 * Lh
        ta = np.interp(dist_from_A, L, tau)
        if local:  # here tau holds sig = 1 - tau
            qq = qP + ta * (qA - qP)
            zz = (1.0 - ta) * zP
            Y = np.sqrt(ca.G_local(qq, zz))
            Z = ca.Z0 + zz
            X = ca.X_te(Z) - qq
        else:
            X = A[0] + ta * (Pl[0] - A[0])
            Z = A[2] + ta * (Pl[2] - A[2])
            Y = np.sqrt(np.clip(ca.G_value(X, Z), 0, None))
        Y = np.where(ta == (1.0 if local else 0.0), A[1], Y)
        out[mask] = np.stack([X, side * Y, Z], -1)
    return (out, Lh) if want_length else out


def edge_derivative_targets(s, Lh):
    """Analytic first and second t-derivatives of the cap at the +Y cut edge (t = 1).

    The wrap section leaves A(s) along the side-view line A -> P(s) with Y = sqrt(G);
    G and its derivatives equal the parent's at the attachment (psi = 1 + O(h^3)),
    so these targets carry the parent skin's tangent plane and curvature.
    """
    A = ca.parent.skin(np.array([s]), ca.Z0, +1)[0]
    P = ca.lip_point(np.array([s]))[0]
    dX, dZ = P[0] - A[0], P[2] - A[2]
    g = ca.G(np.array([A[0]]), np.array([ca.Z0]))
    G = g["v"][0]; Y = np.sqrt(G)
    G1 = g["X"][0] * dX + g["Z"][0] * dZ
    G2 = g["XX"][0] * dX**2 + 2 * g["XZ"][0] * dX * dZ + g["ZZ"][0] * dZ**2
    Y1 = G1 / (2 * Y)
    Y2 = G2 / (2 * Y) - G1**2 / (4 * Y**3)
    Pt = np.array([dX, Y1, dZ]); Ptt = np.array([0.0, Y2, 0.0])
    sp = np.linalg.norm(Pt); T = Pt / sp
    k = (Ptt - Ptt.dot(T) * T) / sp**2
    return -2 * Lh * T, (2 * Lh) ** 2 * k


TE_LOCAL_FROM = 0.9  # use the well-conditioned local TE formulation beyond this chord parameter


def parent_le():
    return np.array([65.0, 0.0, ca.Z0])


def te_segment(t):
    lo, hi = TE_CORNERS
    return lo[None, :] + np.asarray(t, float)[:, None] * (hi - lo)[None, :]


TE_CORNERS = None


# --------------------------------------------------------------------------- fit
def clamped(interior, deg=DEG):
    return np.concatenate([[0.0] * (deg + 1), np.sort(interior), [1.0] * (deg + 1)])


def greville(knots, deg=DEG):
    n = len(knots) - deg - 1
    return np.array([knots[i + 1:i + deg + 1].mean() for i in range(n)])


def fit_joint(u_knots, t_knots, row_minus, row_plus, u_samp, t_samp, log=print,
              w1=1e3, w2=1e3, cache=None):
    """Least squares for all free poles at once: exact analytic points on the sample grid
    plus weighted analytic first and second cross-boundary derivatives at both cut edges
    (normalised by the section length). Boundary poles fixed as in fit()."""
    import scipy.sparse as sp
    from scipy.sparse.linalg import splu
    nu = len(u_knots) - DEG - 1
    nt = len(t_knots) - DEG - 1
    if cache is not None:
        S, Lh = cache
    else:
        S = np.empty((u_samp.size, t_samp.size, 3))
        Lh = np.empty(u_samp.size)
        for i, u in enumerate(u_samp):
            S[i], Lh[i] = section_points(u, t_samp, want_length=True)
    globals()["_LAST_SAMPLES"] = (S, Lh)
    Bu = BSpline.design_matrix(u_samp, u_knots, DEG)
    Bt = BSpline.design_matrix(t_samp, t_knots, DEG)
    Gu = (Bu.T @ Bu).tocsr(); Gt = (Bt.T @ Bt).tocsr()
    N = sp.kron(Gu, Gt, format="csr")
    rhs = np.stack([(Bu.T @ (Bt.T @ S[:, :, k].T).T) for k in range(3)], -1).reshape(nu * nt, 3)
    # derivative targets
    m = (u_samp > 0) & (u_samp < DERIV_U_MAX)
    ud = u_samp[m]; L2 = 2 * Lh[m]
    D1 = np.empty((ud.size, 3)); D2 = np.empty((ud.size, 3))
    for k, (u, L) in enumerate(zip(ud, Lh[m])):
        D1[k], D2[k] = edge_derivative_targets(u, L)
    eye = np.eye(nt)
    d1_1 = BSpline(t_knots, eye, DEG).derivative(1)(1.0); d2_1 = BSpline(t_knots, eye, DEG).derivative(2)(1.0)
    d1_0 = BSpline(t_knots, eye, DEG).derivative(1)(0.0); d2_0 = BSpline(t_knots, eye, DEG).derivative(2)(0.0)
    Bd = BSpline.design_matrix(ud, u_knots, DEG)
    mirror = np.array([1.0, -1.0, 1.0])
    s1 = np.full(ud.size, float(w1)); s2 = np.full(ud.size, float(w2))
    for dvec, target, wscale in ((d1_1, D1, s1), (d2_1, D2, s2), (d1_0, -D1 * mirror, s1), (d2_0, D2 * mirror, s2)):
        c = np.abs(dvec).max()          # express each equation in pole units (coefficients ~ 1)
        dvec = dvec / c
        target = target / c
        Wd = sp.diags(wscale) @ Bd
        Op = sp.kron(Wd, sp.csr_matrix(dvec[None, :]), format="csr")   # rows: u samples
        N = N + (Op.T @ Op)
        rhs = rhs + Op.T @ (wscale[:, None] * target)
    N = N.tocsr()
    # fixed poles
    P = np.zeros((nu, nt, 3))
    fixed = np.zeros((nu, nt), bool)
    P[:, 0] = row_minus; P[:, -1] = row_plus; fixed[:, [0, -1]] = True
    P[0, :] = parent_le(); fixed[0, :] = True
    P[-1, :] = te_segment(greville(t_knots)); fixed[-1, :] = True
    f = fixed.ravel(); free = ~f
    Pv = P.reshape(nu * nt, 3)
    Nff = (N[free][:, free] + 1e-12 * sp.identity(int(free.sum()))).tocsc()
    Nfx = N[free][:, f]
    b = rhs[free] - Nfx @ Pv[f]
    lu = splu(Nff)
    Pv[free] = np.stack([lu.solve(b[:, k]) for k in range(3)], -1)
    P = Pv.reshape(nu, nt, 3)
    fitted = np.stack([(Bu @ (Bt @ P[:, :, k].T).T) for k in range(3)], -1)
    err = np.linalg.norm(fitted - S, axis=-1)
    i = np.unravel_index(np.argmax(err), err.shape)
    log(f"joint fit: poles {nu}x{nt}, samples {u_samp.size}x{t_samp.size}, max err {err.max():.3e} mm at u={u_samp[i[0]]:.6f} t={t_samp[i[1]]:.4f}, rms {np.sqrt((err**2).mean()):.3e}")
    return P, err, S


def make_surface(P, u_knots, t_knots):
    nu, nt, _ = P.shape
    arr = TColgp_Array2OfPnt(1, nu, 1, nt)
    for i in range(nu):
        for j in range(nt):
            arr.SetValue(i + 1, j + 1, gp_Pnt(*P[i, j]))

    def km(k):
        vals, mults = np.unique(np.round(k, 14), return_counts=True)
        a = TColStd_Array1OfReal(1, len(vals)); b = TColStd_Array1OfInteger(1, len(vals))
        for i, (v, m) in enumerate(zip(vals, mults)):
            a.SetValue(i + 1, float(v)); b.SetValue(i + 1, int(m))
        return a, b
    uk, um = km(u_knots); tk, tm = km(t_knots)
    return Geom_BSplineSurface(arr, uk, tk, um, tm, DEG, DEG)


# --------------------------------------------------------------------------- knots / samples
U_EXTRA = [0.001, 0.002, 0.003, 0.004, 0.0055, 0.007, 0.009, 0.012, 0.015,
           0.98, 0.985, 0.9875, 0.99, 0.9925, 0.994, 0.9955, 0.997, 0.998, 0.9985, 0.999,
           0.99925, 0.9995, 0.9997, 0.9998, 0.9999, 0.99995, 0.99998, 0.99999]
NT = 80
T_EXTRA = [0.002, 0.0045, 0.0075, 0.0105]
DERIV_U_MAX = 1 - 1e-4  # derivative targets are not imposed in the singular TE corner
S_MAX = 1 - 1e-10


def knots_and_rows(parent_faces):
    cp = top_row_curve(parent_faces["skin+"])
    cm = top_row_curve(parent_faces["skin-"])
    k0 = np.array([cp.Knot(i) for i in range(1, cp.NbKnots() + 1)])
    ins = sorted(set(list((k0[1:] + k0[:-1]) / 2) + U_EXTRA))
    for c in (cp, cm):
        for k in ins:
            c.InsertKnot(float(k), 1)
    pp, ku = curve_data(cp)
    pm, ku2 = curve_data(cm)
    assert np.allclose(ku, ku2)
    T_EXTRA_sym = sorted(set(T_EXTRA + [1 - x for x in T_EXTRA]))
    tk = clamped(np.unique(np.round(np.concatenate([np.linspace(0, 1, NT + 1)[1:-1], T_EXTRA_sym]), 12)))
    return ku, tk, pm, pp


def samples():
    us = np.unique(np.concatenate([np.linspace(0, S_MAX, 900), np.geomspace(1e-5, 0.03, 120),
                                   1 - np.geomspace(1e-10, 0.03, 300), [0.0, 1.0]]))
    us = us[(us <= S_MAX) | (us == 1.0)]
    ts = np.unique(np.concatenate([np.linspace(0, 1, 801), np.geomspace(1e-6, 0.03, 120),
                                   1 - np.geomspace(1e-6, 0.03, 120)]))
    return us, ts


# --------------------------------------------------------------------------- solid
def build_solid(parent, cap_surface, tol=1e-6):
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid
    from OCP.ShapeFix import ShapeFix_Solid
    from OCP.TopAbs import TopAbs_SHELL
    cap_face = BRepBuilderAPI_MakeFace(cap_surface, 1e-7).Face()
    sew = BRepBuilderAPI_Sewing(tol)
    for k in ("skin+", "skin-", "te", "root"):
        sew.Add(parent[k])
    sew.Add(cap_face)
    sew.Perform()
    sewed = sew.SewedShape()
    e = TopExp_Explorer(sewed, TopAbs_SHELL)
    shell = TopoDS.Shell_s(e.Current())
    solid = BRepBuilderAPI_MakeSolid(shell).Solid()
    fix = ShapeFix_Solid(solid)
    fix.Perform()
    return fix.Solid(), sew


def write_step(shape, path):
    from OCP.Interface import Interface_Static
    Interface_Static.SetCVal_s("write.step.schema", "AP214IS")
    Interface_Static.SetCVal_s("write.step.unit", "MM")
    w = STEPControl_Writer()
    w.Transfer(shape, STEPControl_AsIs)
    if w.Write(str(path)) != 1:
        raise RuntimeError("STEP write failed")


def main():
    global TE_CORNERS
    parent_shape = read_step(PARENT_STEP)
    parent = classify_parent(parent_shape)
    ku, tk, pm, pp = knots_and_rows(parent)
    TE_CORNERS = (pm[-1], pp[-1])
    us, ts = samples()
    P, err, S = fit_joint(ku, tk, pm, pp, us, ts)
    surf = make_surface(P, ku, tk)
    solid, sew = build_solid(parent, surf)
    from OCP.BRepCheck import BRepCheck_Analyzer
    ok = BRepCheck_Analyzer(solid).IsValid()
    write_step(solid, OUT_STEP)
    info = {
        "cap_poles": list(P.shape[:2]),
        "degree": [DEG, DEG],
        "u_knots_distinct": int(len(np.unique(ku))),
        "t_knots_distinct": int(len(np.unique(tk))),
        "fit_samples": [int(us.size), int(ts.size)],
        "fit_max_error_mm": float(err.max()),
        "fit_rms_error_mm": float(np.sqrt((err ** 2).mean())),
        "fit_max_error_mm_u_le_0p99": float(err[us <= 0.99].max()),
        "eps": ca.EPS,
        "brepcheck_valid_before_export": bool(ok),
        "step_bytes": OUT_STEP.stat().st_size,
    }
    (HERE / "build.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
