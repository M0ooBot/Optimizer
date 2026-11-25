"""Design of the rounded tip cap for fin5 / sharp (shared by build and verification).

Construction (see README.md): every horizontal section of the cap is an ellipse
centred on the parent's mid-chord line X = m(Z), with
    semi-chord      a_s(Z) = a(Z) * sqrt(1 - h^3)       (exactly the locked outline)
    half thickness  b_s(Z) = Yp(h)                      (the designed wrap profile)
    h = (Z - 240) / 5,  m = 80 + 27.5 Z/240,  a = 80 - 37.5 Z/240,
and the ellipse is parametrised by the parent's own chord angle theta, so
    X = m(Z) - a_s(Z) cos(theta),  Y = +/- b_s(Z) sin(theta).
At Z = 240 this is the parent section itself. The cap's meridians (theta = const)
are the spanwise wrap guides; theta = 0 and theta = pi are the 3D perimeter (the
locked outline, lying in Y = 0); the crown (h = 1) is a collapsed pole.

Yp(h) is the mid-chord wrap profile: a plane curve that leaves the parent skin
with its slope and zero curvature (G2), ramps its curvature up with a smoothstep
over the first S1 fraction of its arc length, then relaxes smoothly to a
slightly flatter crown, arriving horizontal at (Y = 0, Z = 245).

Spanwise parameter of the CAD patches: v in [0, 1] with h = 1 - (1 - v)^2, which
makes every coordinate a smooth function of v up to the crown pole.
"""
import importlib.util
import pathlib

import numpy as np
from scipy.interpolate import CubicSpline, BSpline
from scipy.optimize import fsolve, brentq
from scipy.integrate import cumulative_simpson

ROOT = pathlib.Path(__file__).resolve().parents[5]
PARENT_DIR = ROOT / "runs/cap-radius/parents/fin5/sharp"


def load_parent():
    spec = importlib.util.spec_from_file_location("fin5_sharp_parent", PARENT_DIR / "parent.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


par = load_parent()
ATT = par.ATTACHMENT          # 240
REACH = par.REACH             # 5
TOP = ATT + REACH             # 245
DM = 27.5 / ATT               # dm/dZ
DA = -37.5 / ATT              # da/dZ
TR = par.THICKNESS_RATIO      # 0.12 (half thickness = 0.06 c = 0.12 a)


def m_of(z):
    return 80.0 + DM * np.asarray(z, float)


def a_of(z):
    return 80.0 + DA * np.asarray(z, float)


# ---------------------------------------------------------------- wrap profile
RAMP = 0.12          # fraction of the profile arc length over which curvature ramps up
W0 = TR * a_of(ATT)  # 5.1 mm, attachment half thickness at mid-chord
LEAN = np.arctan(-TR * DA)   # parent skin leans inward by atan(0.01875)
N_PROFILE = 200001


def _smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3 - 2 * x)


def _psi(sig, rho):
    up = _smoothstep(sig / RAMP)
    down = 1 - (1 - rho) * _smoothstep((sig - RAMP) / (1 - RAMP))
    return np.where(sig < RAMP, up, down)


def _integrate(K, S, rho, n=N_PROFILE):
    sig = np.linspace(0.0, 1.0, n)
    k = K * _psi(sig, rho)
    phi = LEAN + cumulative_simpson(k * S, x=sig, initial=0.0)
    y = W0 - cumulative_simpson(np.sin(phi) * S, x=sig, initial=0.0)
    z = cumulative_simpson(np.cos(phi) * S, x=sig, initial=0.0)
    return sig, k, phi, y, z


def solve_profile():
    def F(p):
        K, S, rho = p
        _, _, phi, y, z = _integrate(K, S, rho, 20001)
        return [phi[-1] - np.pi / 2, y[-1], z[-1] - REACH]
    sol = fsolve(F, [0.23, 8.0, 0.75], xtol=1e-13)

    def F_fine(p):
        _, _, phi, y, z = _integrate(*p)
        return [phi[-1] - np.pi / 2, y[-1], z[-1] - REACH]
    K, S, rho = fsolve(F_fine, sol, xtol=1e-14)
    sig, k, phi, y, z = _integrate(K, S, rho)
    return {"K": K, "S": S, "rho": rho, "sig": sig, "k": k, "phi": phi, "y": y, "z": z}


PROFILE = solve_profile()
_zs = CubicSpline(PROFILE["sig"], PROFILE["z"])
_ys = CubicSpline(PROFILE["sig"], PROFILE["y"])


def profile_y_of_height(zeta):
    """Mid-chord half thickness Yp at height zeta = Z - 240 (0..5), from the designed curve."""
    zeta = np.atleast_1d(np.asarray(zeta, float))
    out = np.empty_like(zeta)
    for i, zz in enumerate(zeta):
        if zz <= 0:
            out[i] = W0
        elif zz >= REACH:
            out[i] = 0.0
        else:
            s = brentq(lambda q: _zs(q) - zz, 0.0, 1.0, xtol=1e-15, rtol=1e-15, maxiter=500)
            out[i] = _ys(s)
    return out


# ------------------------------------------------------------ spanwise functions
def h_of_v(v):
    v = np.asarray(v, float)
    return 1 - (1 - v) ** 2


def z_of_v(v):
    return ATT + REACH * h_of_v(v)


def s1_of_v(v):
    """sqrt(1 - h^3) written smoothly: (1 - v) sqrt(1 + h + h^2)."""
    h = h_of_v(v)
    return (1 - np.asarray(v, float)) * np.sqrt(1 + h + h * h)


def spanwise_functions(v):
    """Mx, Ax, Ay, Z at spanwise parameter v (exact design values)."""
    v = np.asarray(v, float)
    z = z_of_v(v)
    return {"Mx": m_of(z), "Ax": a_of(z) * s1_of_v(v), "Ay": profile_y_of_height(z - ATT), "Z": z}


# derivatives at v = 0 (G2 to the ruled parent) : Z' = 10, Z'' = -10, S1 and Yp'' vanish
Z1, Z2 = 2 * REACH, -2 * REACH
END0 = {
    "Mx": (m_of(ATT), DM * Z1, DM * Z2),
    "Ax": (a_of(ATT), DA * Z1, DA * Z2),
    "Ay": (W0, TR * DA * Z1, TR * DA * Z2),
    "Z": (ATT, Z1, Z2),
}

# ------------------------------------------------------------ B-spline fitting
V_DEGREE = 5
V_INTERIOR = 32


def ramp_end_v():
    """Spanwise parameter where the profile's curvature ramp ends (a C3 joint of the profile)."""
    i = int(round(RAMP * (N_PROFILE - 1)))
    h = PROFILE["z"][i] / REACH
    return float(1 - np.sqrt(1 - h))


def v_knots():
    """Uniform interior knots plus a double knot at the ramp end, so the C3 joint is represented exactly."""
    inner = np.linspace(0, 1, V_INTERIOR + 2)[1:-1]
    vr = ramp_end_v()
    inner = inner[np.abs(inner - vr) > 0.5 / (V_INTERIOR + 1)]
    inner = np.sort(np.concatenate([inner, [vr, vr]]))
    return np.concatenate([np.zeros(V_DEGREE + 1), inner, np.ones(V_DEGREE + 1)])


def _basis_rows(t, p, x, nu):
    n = len(t) - p - 1
    rows = np.zeros((len(x), n))
    for j in range(n):
        c = np.zeros(n)
        c[j] = 1.0
        b = BSpline(t, c, p)
        rows[:, j] = b.derivative(nu)(x) if nu else b(x)
    return rows


def fit_spanwise(n_samples=1201):
    """Constrained least squares of Mx, Ax, Ay, Z on a common degree-5 basis in v.

    Exact at v = 0 for value, first and second derivative (G2 to the parent) and
    exact at v = 1 for the value (Ax = Ay = 0: collapsed crown, Z = 245).
    """
    t = v_knots()
    p = V_DEGREE
    n = len(t) - p - 1
    xs = 0.5 - 0.5 * np.cos(np.linspace(0, np.pi, n_samples))
    vals = spanwise_functions(xs)
    Bm = _basis_rows(t, p, xs, 0)
    C0 = np.vstack([_basis_rows(t, p, np.array([0.0]), d) for d in (0, 1, 2)])
    C1 = _basis_rows(t, p, np.array([1.0]), 0)
    coefs, err = {}, {}
    for key in ("Mx", "Ax", "Ay", "Z"):
        end1 = {"Mx": m_of(TOP), "Ax": 0.0, "Ay": 0.0, "Z": TOP}[key]
        Cm = np.vstack([C0, C1])
        d = np.array(list(END0[key]) + [end1])
        kkt = np.block([[Bm.T @ Bm, Cm.T], [Cm, np.zeros((len(d), len(d)))]])
        rhs = np.concatenate([Bm.T @ vals[key], d])
        sol = np.linalg.solve(kkt, rhs)
        c = sol[:n]
        coefs[key] = c
        check = np.linspace(0, 1, 4001)
        err[key] = float(np.max(np.abs(BSpline(t, c, p)(check) - spanwise_functions(check)[key])))
    return t, coefs, err
