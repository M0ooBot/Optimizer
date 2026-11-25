"""Build the Wing B (blunt) rounded tip cap and write wing.step.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/wingb/blunt/build_cap.py
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.interpolate import BSpline

sys.path.insert(0, str(Path(__file__).resolve().parent))
import capgeom as G  # noqa: E402

OUT = Path(__file__).resolve().parent
DEG = 3
T_INTERIOR = 30
SPLIT_PASSES = 0
SPLIT_RATIO = 4.0


def clamp_knots(interior, deg=DEG):
    return np.concatenate([[0.0] * (deg + 1), np.asarray(interior, float), [1.0] * (deg + 1)])


def _solve_s(f, target, lo=0.0, hi=1.0, it=60):
    """Bisection for an increasing function f(s) = target."""
    lo, hi = float(lo), float(hi)
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if f(np.array([mid]))[0] < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def lip_split_s():
    """Station where the trailing-edge strip starts (outline Delta = L_TE)."""
    f = lambda s: -(-G.outline_frame(G.station_p(s))[0][:, 0] - G.ATT)
    return _solve_s(f, -G.L_TE, 0.9, 1.0)


def s_knots(s_split):
    T = G._taylor()
    uk = np.unique(T["tu"])
    un = G.u_nose0()
    parent = np.concatenate([(uk[(uk > un) & (uk < 1)] - un) / (1 - un), (un - uk[(uk > 0) & (uk < un)]) / un])
    parent = np.unique(parent)
    extra = [s_split] + [_solve_s(G.station_p, pk) for pk in (0.18, 0.82)]
    extra += list(1 - np.geomspace(3e-2, 3e-5, 24) * (1 - s_split) / 2e-2)
    extra = np.array([e for e in extra if 0 < e < 1 and np.min(np.abs(parent - e)) > 2e-6])
    k = np.unique(np.concatenate([parent, extra]))
    # split spans much longer than a neighbour
    for _ in range(SPLIT_PASSES):
        kk = np.concatenate([[0.0], k, [1.0]])
        span = np.diff(kk)
        left = np.concatenate([[np.inf], span[:-1]])
        right = np.concatenate([span[1:], [np.inf]])
        long_ = span > SPLIT_RATIO * np.minimum(left, right)
        if not long_.any():
            break
        k = np.unique(np.concatenate([k, 0.5 * (kk[:-1] + kk[1:])[long_]]))
    return k[(k > 0) & (k < 1)]


def sample_points(knots, per_span):
    g, _ = np.polynomial.legendre.leggauss(per_span)
    k = np.unique(knots)
    a, b = k[:-1], k[1:]
    pts = ((0.5 * (b - a))[:, None] * (g[None, :] + 1) + a[:, None]).ravel()
    return np.unique(np.concatenate([pts, k]))


def basis(knots, x, deg=DEG):
    return BSpline.design_matrix(np.clip(x, 0, 1), knots, deg).tocsr()


def lsq_1d(knots, x, y, zero=()):
    A = basis(knots, x)
    keep = [i for i in range(A.shape[1]) if i not in set(zero)]
    Ak = A[:, keep]
    c = spla.spsolve((Ak.T @ Ak).tocsc(), Ak.T @ y)
    out = np.zeros(A.shape[1])
    out[keep] = c
    return out


def constrained_lsq(A, B, data, fixed, ties=()):
    """min || A P B^T - data || over the poles P (ns x nt);
    fixed: {(i, j): value}; ties: groups of (i, j) sharing one unknown."""
    ns, nt = A.shape[1], B.shape[1]
    N = ns * nt
    owner = np.arange(N)
    for grp in ties:
        idx = [i * nt + j for i, j in grp]
        owner[idx] = idx[0]
    fixed_idx = {i * nt + j: v for (i, j), v in fixed.items()}
    free = [k for k in range(N) if owner[k] == k and k not in fixed_idx]
    col = {k: c for c, k in enumerate(free)}
    rows, cols, g0 = [], [], np.zeros(N)
    for k in range(N):
        o = owner[k]
        if o in fixed_idx:
            g0[k] = fixed_idx[o]
        else:
            rows.append(k)
            cols.append(col[o])
    Gm = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(N, len(free)))
    K = sp.kron(A, B, format="csr")
    KG = K @ Gm
    rhs = data.ravel() - K @ g0
    q = spla.spsolve((KG.T @ KG).tocsc(), KG.T @ rhs)
    return (Gm @ q + g0).reshape(ns, nt)


def fit(log=print):
    t0 = time.time()
    s_split = lip_split_s()
    sk = clamp_knots(s_knots(s_split))
    tu_ = np.linspace(0, 1, T_INTERIOR + 2)[1:-1]
    tk = clamp_knots(np.sort(np.concatenate([tu_, tu_[0] * np.array([0.25, 0.5, 0.75])])))
    ns, nt = len(sk) - DEG - 1, len(tk) - DEG - 1
    s = sample_points(sk, 3)
    t = sample_points(tk, 6)
    A, Bt = basis(sk, s), basis(tk, t)
    log(f"s poles {ns}, t poles {nt}; samples {len(s)} x {len(t)}; s_split {s_split:.6f}")

    # shared lip row: outline position, mean height and (trailing-edge strip) half-height
    sd = np.unique(np.concatenate([np.linspace(0, 1, 20001), sample_points(sk, 6)]))
    o = G.outline_frame(G.station_p(sd))[0]
    m_o, h_o, *_ = G.cap_fields(o[:, 0], o[:, 1])
    zero = [i for i in range(ns) if sk[i] < s_split - 1e-12]
    lip = {"X": lsq_1d(sk, sd, o[:, 0]), "Z": lsq_1d(sk, sd, o[:, 1]),
           "m": lsq_1d(sk, sd, m_o), "h": lsq_1d(sk, sd, np.where(sd < s_split, 0.0, h_o), zero)}
    for k_, v in (("X", o[:, 0]), ("Z", o[:, 1]), ("m", m_o)):
        log(f"  lip row {k_}: max error {np.abs(basis(sk, sd) @ lip[k_] - v).max():.2e}")
    log(f"  lip row h: max error {np.abs(basis(sk, sd) @ lip['h'] - np.where(sd < s_split, 0.0, h_o)).max():.2e}")
    patches = {}
    for side in (+1, -1):
        name = "upper" if side > 0 else "lower"
        zc, yc, _ = G.cut(sd, side)
        cutrow = {"X": np.full(ns, -G.ATT), "Z": lsq_1d(sk, sd, zc), "Y": lsq_1d(sk, sd, yc)}
        log(f"  {name} cut row: max error Z {np.abs(basis(sk, sd) @ cutrow['Z'] - zc).max():.2e}, "
            f"Y {np.abs(basis(sk, sd) @ cutrow['Y'] - yc).max():.2e}")
        liprow = {"X": lip["X"], "Z": lip["Z"], "Y": lip["m"] + side * lip["h"]}
        # rows 1 and 2 from the exact first and second t-derivatives at the cut, so the
        # cap leaves the parent with its tangent plane and curvature (G2 by construction)
        hstep = 2e-3
        tt = np.arange(5) * hstep
        qd, _ = G.plan_map(sd[:, None], tt[None, :], side)
        md, hd, *_ = G.cap_fields(qd[..., 0], qd[..., 1])
        V = np.vander(tt, 5, increasing=True)          # quartic through the 5 samples
        k1, k2 = tk[DEG + 1], tk[DEG + 2]
        rows12 = {}
        for f, d in (("X", qd[..., 0]), ("Y", md + side * hd), ("Z", qd[..., 1])):
            c = np.linalg.solve(V, d.T)                 # coefficients (5, n)
            d1, d2 = c[1], 2 * c[2]
            P0 = basis(sk, sd) @ cutrow[f]
            P1 = P0 + k1 / 3.0 * d1
            P2 = P1 + k2 / 3.0 * (d1 + k1 / 2.0 * d2)
            rows12[f] = (lsq_1d(sk, sd, P1), lsq_1d(sk, sd, P2))
        q, p = G.plan_map(s[:, None], t[None, :], side)
        m, h, *_ = G.cap_fields(q[..., 0], q[..., 1])
        data = {"X": q[..., 0], "Z": q[..., 1], "Y": m + side * h}
        log(f"  {name} sampled ({time.time() - t0:.0f}s)")
        poles = {}
        for f, d in data.items():
            fixed = {}
            for j in range(nt):
                fixed[(0, j)] = cutrow[f][0]
                fixed[(ns - 1, j)] = cutrow[f][-1]
            for i in range(1, ns - 1):
                fixed[(i, 0)] = cutrow[f][i]
                fixed[(i, 1)] = rows12[f][0][i]
                fixed[(i, 2)] = rows12[f][1][i]
                fixed[(i, nt - 1)] = liprow[f][i]
            ties = [[(i, nt - 2), (i, nt - 1)] for i in range(1, ns - 1)] if f in ("X", "Z") else ()
            poles[f] = constrained_lsq(A, Bt, d, fixed, ties)
            res = A @ poles[f] @ Bt.T - d
            i, j = np.unravel_index(np.argmax(np.abs(res)), res.shape)
            log(f"    {f}: max residual {np.abs(res).max():.2e} at s={s[i]:.5f} t={t[j]:.4f}")
        patches[name] = np.stack([poles["X"], poles["Y"], poles["Z"]], -1)
    log(f"fitted in {time.time() - t0:.0f}s")
    return {"sk": sk, "tk": tk, "s_split": s_split, "upper": patches["upper"], "lower": patches["lower"]}


# ----------------------------------------------------------------------------
# Open CASCADE: surfaces, faces, sewing with the parent, STEP export
# ----------------------------------------------------------------------------
def occ_bspline_surface(poles, uk, vk, du=DEG, dv=DEG):
    from OCP.Geom import Geom_BSplineSurface
    from OCP.TColgp import TColgp_Array2OfPnt
    from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
    from OCP.gp import gp_Pnt
    nu, nv, _ = poles.shape
    arr = TColgp_Array2OfPnt(1, nu, 1, nv)
    for i in range(nu):
        for j in range(nv):
            arr.SetValue(i + 1, j + 1, gp_Pnt(*map(float, poles[i, j])))

    def km(k):
        u, m = np.unique(np.round(k, 15), return_counts=True)
        a, b = TColStd_Array1OfReal(1, len(u)), TColStd_Array1OfInteger(1, len(u))
        for i, (x, c) in enumerate(zip(u, m)):
            a.SetValue(i + 1, float(x))
            b.SetValue(i + 1, int(c))
        return a, b

    ua, um = km(uk)
    va, vm = km(vk)
    return Geom_BSplineSurface(arr, ua, va, um, vm, du, dv)


def parent_faces():
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.GeomAbs import GeomAbs_Plane
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_Reader
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    r = STEPControl_Reader()
    assert r.ReadFile(str(G.P.STEP_PATH)) == IFSelect_RetDone
    r.TransferRoots()
    keep, dropped = [], []
    ex = TopExp_Explorer(r.OneShape(), TopAbs_FACE)
    while ex.More():
        f = TopoDS.Face_s(ex.Current())
        ad = BRepAdaptor_Surface(f)
        if ad.GetType() == GeomAbs_Plane:
            pr = GProp_GProps()
            BRepGProp.SurfaceProperties_s(f, pr)
            c = pr.CentreOfMass()
            if abs(c.X() + G.ATT) < 1e-6:
                dropped.append(f)
                ex.Next()
                continue
        keep.append(f)
        ex.Next()
    assert len(dropped) == 1 and len(keep) == 3
    return keep


def build_solid(F, log=print):
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.ShapeFix import ShapeFix_Solid
    from OCP.TopAbs import TopAbs_SHELL
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    sk, tk = F["sk"], F["tk"]
    up = occ_bspline_surface(F["upper"], sk, tk)
    lo = occ_bspline_surface(F["lower"], sk, tk)
    strip_poles = np.stack([F["upper"][:, -1, :], F["lower"][:, -1, :]], 1)
    st = occ_bspline_surface(strip_poles, sk, np.array([0.0, 0.0, 1.0, 1.0]), DEG, 1)
    st.Segment(float(F["s_split"]), 1.0, 0.0, 1.0)
    faces = [BRepBuilderAPI_MakeFace(s, 1e-7).Face() for s in (up, lo, st)]
    sew = BRepBuilderAPI_Sewing(SEW_TOL)
    for f in parent_faces() + faces:
        sew.Add(f)
    sew.Perform()
    log(f"sewing: free edges {sew.NbFreeEdges()}, multiple edges {sew.NbMultipleEdges()}, "
        f"degenerated {sew.NbDegeneratedShapes()}")
    sh = sew.SewedShape()
    ex = TopExp_Explorer(sh, TopAbs_SHELL)
    shells = []
    while ex.More():
        shells.append(TopoDS.Shell_s(ex.Current()))
        ex.Next()
    assert len(shells) == 1, f"{len(shells)} shells"
    ms = BRepBuilderAPI_MakeSolid(shells[0])
    fix = ShapeFix_Solid(ms.Solid())
    fix.Perform()
    solid = fix.Solid()
    log(f"BRepCheck valid: {BRepCheck_Analyzer(solid).IsValid()}")
    return solid, {"upper": up, "lower": lo, "strip": st}


def write_step(shape, path):
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.Interface import Interface_Static
    from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer
    w = STEPControl_Writer()
    Interface_Static.SetCVal_s("write.step.unit", "MM")
    Interface_Static.SetCVal_s("write.step.schema", "AP214IS")
    w.Transfer(shape, STEPControl_AsIs)
    assert w.Write(str(path)) == IFSelect_RetDone


SEW_TOL = 1e-6

if __name__ == "__main__":
    F = fit()
    solid, _ = build_solid(F)
    write_step(solid, OUT / "wing.step")
    print("wrote", OUT / "wing.step", (OUT / "wing.step").stat().st_size, "bytes")
