"""Analytic definition of the rounded tip cap for the flat/blunt parent.

The cap is one surface per tip, parameterised by the parent's own chord
parameter s (0 = leading edge, 1 = trailing-edge corner) and a wrap parameter
w in [-1, 1] (w = +1 on the upper cut edge, w = 0 on the lip, w = -1 on the
lower cut edge). Every wrap lies in a constant-X plane (a chord station):

    X  = X_parent(s)                       (parent skin chordwise position)
    |Y| = 350 + R(s) * E(w)                (spanwise reach to the lip)
    Z  = h(s) * S(w)                       (h = parent half thickness at s)

E(w) = (1 - w^2) (1 + A w^2) and S(w) = (15 w - 10 w^3 + 3 w^5) / 8.
S'(1) = S''(1) = 0 and E(+-1) = 0, so the cap meets the parent skin with
matching tangent plane and matching curvature (G2) along the cut edges.
E'(0) = 0, S'(0) != 0, so each wrap crosses the lip vertically and the
upper and lower halves are one smooth surface.

R(s) is the locked outline's reach at the same X, R_lock, smoothly maxed with
K * h(s) so that at the trailing edge (where the locked outline has zero
reach) the wrap becomes a rounded closure of the flat trailing-edge face with
apex radius h_TE = 0.2 mm (half the trailing-edge thickness). The smooth max
moves the lip by less than 1e-5 mm in plan (measured by the verification).

Needs NumPy and runs/cap-radius/parents/flat/blunt/parent.py.
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.normpath(os.path.join(HERE, '..', '..', '..', 'parents', 'flat', 'blunt'))
sys.path.insert(0, PARENT_DIR)
import parent as PP  # noqa: E402
from parent import Jet  # noqa: E402

CHORD = PP.CHORD
ATT = PP.ATTACHMENT
H_TE = PP.TE_THICKNESS / 2.0

# ---- design parameters (see README / skill-notes for how they were chosen) ----
A = 0.0          # E(w) = (1 - w^2)(1 + A w^2)
NMAX = 12        # exponent of the smooth max between R_lock and K h
LIP_TE = H_TE    # apex radius of the closure of the trailing-edge face (face plane)
# slice lip radius of the wrap = (225/128) h^2 / ((1 - A) R)  ->  K = R_TE / h_TE
K = (225.0 / 128.0) * H_TE / ((1.0 - A) * LIP_TE)


def _xpoly(q):
    return q * q * q * (10 - 15 * q + 6 * q * q)


def _xinv(t):
    """Inverse of the smoothstep x on [0, 1] (vectorised bisection + Newton)."""
    t = np.asarray(t, float)
    lo = np.zeros_like(t); hi = np.ones_like(t)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        m = _xpoly(mid) < t
        lo = np.where(m, mid, lo); hi = np.where(m, hi, mid)
    q = 0.5 * (lo + hi)
    return q


def _pow_frac(j, a):
    v = j.v
    return j.chain(v ** a, a * v ** (a - 1), a * (a - 1) * v ** (a - 2))


def r_lock_jet(eps, X=None):
    """Locked outline reach at chord station eps = 120 - X (jet).

    Uses q = 1 - p from eps on the aft half and p from X on the forward half
    (each is the accurate one near its own end).
    """
    if X is not None:
        fwd = X.v < CHORD / 2
        if np.any(fwd):
            T = X / CHORD
            p = np.clip(_xinv(np.clip(T.v, 0, 1)), 1e-300, 1.0)
            d1 = 30 * p * p * (1 - p) ** 2
            d2 = 60 * p * (1 - p) * (1 - 2 * p)
            with np.errstate(all='ignore'):
                pj = T.chain(p, 1.0 / d1, -d2 / d1 ** 3)
            rf = 200.0 * pj * (1 - pj)
            ra = r_lock_jet(eps)
            return PP.jwhere(fwd, rf, ra)
    T = eps / CHORD
    tv = np.clip(T.v, 0.0, 1.0)
    q = _xinv(tv)
    q = np.clip(q, 1e-300, 1.0)
    d1 = 30 * q * q * (1 - q) ** 2
    d2 = 60 * q * (1 - q) * (1 - 2 * q)
    with np.errstate(all='ignore'):
        qj = T.chain(q, 1.0 / d1, -d2 / d1 ** 3)
    return 200.0 * qj * (1 - qj)


LAW = 'local'
EPS0 = 2e-3      # aft shift of the locked outline at the TE (mm), see r_mod_jet
EPS1 = 6e-3      # decay length of that shift and offset (mm)


def r_mod_jet(eps, h, X=None):
    """Reach of the lip at chord station eps = 120 - X."""
    if LAW == 'max':
        rl = r_lock_jet(eps)
        small = eps.v < 1e-13
        rl_n = rl ** NMAX
        rl_n = PP.jwhere(small, Jet.constant(np.zeros(eps.v.shape), eps.n), rl_n)
        rc_n = (K * h) ** NMAX
        return _pow_frac(rl_n + rc_n, 1.0 / NMAX)
    if LAW == 'shift':
        # locked outline shifted aft by EPS0, then offset so the reach at the TE is K h_TE
        r0 = r_lock_jet(Jet.constant(np.full(eps.v.shape, EPS0), eps.n))
        return r_lock_jet(eps + EPS0) - r0 + K * H_TE
    if LAW == 'local':
        # near the TE only: shift the locked outline aft by EPS0 * exp(-eps/EPS1) and add a
        # decaying offset so that the reach at the TE is exactly K h_TE (closure apex radius h_TE)
        arg = eps * (-1.0 / EPS1)
        decay = arg.chain(np.exp(arg.v), np.exp(arg.v), np.exp(arg.v))
        r0 = r_lock_jet(Jet.constant(np.full(eps.v.shape, EPS0), eps.n))
        return r_lock_jet(eps + EPS0 * decay, X) + (K * H_TE - r0) * decay
    raise ValueError(LAW)


def surface_jet(s, w, side=+1):
    """Cap point with exact first/second derivatives in (s, w). side +1 = right tip (+Y)."""
    s = np.asarray(s, float); w = np.asarray(w, float)
    s, w = np.broadcast_arrays(s, w)
    sj = Jet.variable(s, 0, 2); wj = Jet.variable(w, 1, 2)
    span = Jet.constant(np.full(s.shape, ATT), 2)
    Xp, _, Zp = PP._skin(sj, span, +1)
    h = Zp
    eps = CHORD - Xp
    Rm = r_mod_jet(eps, h, Xp)
    E = (1 - wj * wj) * (1 + A * wj * wj)
    S = (15 * wj - 10 * wj ** 3 + 3 * wj ** 5) / 8.0
    X = Xp + 0 * wj
    Y = side * (ATT + Rm * E)
    Z = h * S
    return PP._pack([X, Y, Z], 's', 'w')


def reach(s):
    """R(s) values (mm) and the locked reach at the same X."""
    j = surface_jet(s, 0.0)
    X = j['point'][..., 0]
    return j['point'][..., 1] - ATT, X


def half_thickness(s):
    return PP.skin(np.asarray(s, float), ATT, +1)[..., 2]


def principal(jet):
    """Principal curvatures (k1 >= k2 by |k|), normal, and the metric for a packed jet."""
    Ps, Pw = jet['d_s'], jet['d_w']
    Pss, Psw, Pww = jet['d_s_s'], jet['d_s_w'], jet['d_w_w']
    n = np.cross(Ps, Pw)
    nn = np.linalg.norm(n, axis=-1, keepdims=True)
    n = n / nn
    E1 = (Ps * Ps).sum(-1); F1 = (Ps * Pw).sum(-1); G1 = (Pw * Pw).sum(-1)
    L = (Pss * n).sum(-1); M = (Psw * n).sum(-1); N = (Pww * n).sum(-1)
    det = E1 * G1 - F1 * F1
    Hm = (E1 * N - 2 * F1 * M + G1 * L) / (2 * det)
    Kg = (L * N - M * M) / det
    disc = np.sqrt(np.maximum(Hm * Hm - Kg, 0))
    k1 = Hm + disc; k2 = Hm - disc
    kmax = np.maximum(np.abs(k1), np.abs(k2))
    return {'k1': k1, 'k2': k2, 'kmax': kmax, 'n': n, 'H': Hm, 'K': Kg}
