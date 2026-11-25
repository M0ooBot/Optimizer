"""Analytic rounded-tip cap for the swept55 / sharp parent (right tip, +Y).

The cap is written on the parent's own extended planform. A cap point is
indexed by the parent's chord parameter q (0 = LE, 1 = TE; chord fraction
U(q) = 3q^2 - 2q^3) and a wrap angle theta in [0, pi] (0 = upper cut edge,
pi/2 = lip, pi = lower cut edge):

    reach      D(q) = 240 q (1 - q)            (locked outline: r_lip = 390 + D)
    station    r = 390 + D(q) sin(theta)
    plan       X = tan55 r + c(r) U(q),  Y = r,  c(r) = 210 - (115/390) r
    height     Z = t(q) C(q, x) zeta(theta)
    thickness  C = c(r) + k_u(q) * TAPER * D * (x - h(x)),  h' = (1 - (x/x0)^2)^3 up to x0, 0 after
               (k_u = 0: the parent's taper continues; k_u = 1: the taper is eased
               out by the start of the roll so the lip keeps more of the attachment thickness)
    profile    x = sin(theta); zeta = sign(cos theta) * 1                     for x <= x0
                                zeta = sign(cos theta) * sqrt(1 - xi^b),     xi = (x - x0)/(1 - x0)

Each wrap guide is a straight plan line at constant chord fraction (the line
the parent's own constant-chord-fraction lines continue on), so at theta = 0
the cap is the parent's extended skin plus a vertical offset that is O(x^3) (the
taper easing) and zero up to x0 for the roll, so the cap meets the parent
with matched position, tangent plane and curvature (G2). At theta = pi/2 it reaches the
locked outline (p = q) with a vertical tangent.

Knobs along the chord: x0(q) (design table). x0 = 0 is a cubic superellipse
spread over the whole reach; a larger x0 keeps the parent's extended skin
further out and rounds only the outer (1 - x0) of the reach, which opens the
lip radius. b(q) = 3 (smooth onset) blends to 2.25 (rounder nose for the same
thickness) and k_u(q) blends 0 -> 1 only toward the thin trailing edge, where
the lip radius is limited by the local thickness.
"""

import importlib.util
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[5]
PARENT_DIR = ROOT / "runs/cap-radius/parents/swept55/sharp"

_spec = importlib.util.spec_from_file_location("parent_swept55_sharp", PARENT_DIR / "parent.py")
parent = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(parent)
Jet = parent.Jet

ATT = parent.ATTACHMENT            # 390
TAPER = parent.TAPER
SWEEP = parent.SWEEP
ROOT_CHORD = parent.ROOT_CHORD


def jexp(x):
    e = np.exp(x.v)
    return x.chain(e, e, e)


def jlog(x):
    return x.chain(np.log(x.v), 1.0 / x.v, -1.0 / x.v**2)


def jlog1p(x):
    return x.chain(np.log1p(x.v), 1.0 / (1 + x.v), -1.0 / (1 + x.v) ** 2)


def jexpm1(x):
    e = np.exp(x.v)
    return x.chain(np.expm1(x.v), e, e)


# ---------------------------------------------------------------- design knob
B_EXP = 3.0           # superellipse exponent of the roll (> 2 keeps the join G2)
B_TE = 2.5            # exponent toward the trailing edge
B_BLEND = (0.70, 0.85)
KU_BLEND = (0.60, 0.80)


def _smoothstep_jet(qj, lo, hi):
    u = (qj - lo) / (hi - lo)
    uc = np.clip(u.v, 0, 1)
    u = Jet(uc, np.where(((u.v > 0) & (u.v < 1))[..., None], u.d, 0 * u.d), np.where(((u.v > 0) & (u.v < 1))[..., None, None], u.h, 0 * u.h))
    return u * u * u * (10 - 15 * u + 6 * u * u)


def b_jet(qj):
    return B_EXP + (B_TE - B_EXP) * _smoothstep_jet(qj, *B_BLEND)


def ku_jet(qj):
    return _smoothstep_jet(qj, *KU_BLEND)


def b_of_q(q):
    return b_jet(Jet.constant(np.asarray(q, float), 1)).v


def ku_of_q(q):
    return ku_jet(Jet.constant(np.asarray(q, float), 1)).v
R_FLOOR = 0.5         # design lip radius in the thick part (rule minimum 0.4 + 25 %)


def thickness_at_attachment(q):
    """Parent thickness (mm) at the attachment at chord parameter q."""
    qj = Jet.constant(np.asarray(q, float), 1)
    return 2 * 95.0 * parent._t(qj).v


def _load_table():
    import json
    f = pathlib.Path(__file__).with_name("design_table.json")
    if f.exists():
        d = json.loads(f.read_text())
        return np.array(d["q"]), np.array(d["x0"])
    return np.array([0.0, 1.0]), np.array([0.0, 0.0])


DESIGN_TABLE = _load_table()


def x0_of_q(q, table=None):
    """x0(q): the same natural cubic spline the jets use."""
    return _x0_jet(Jet.constant(np.asarray(q, float), 1), DESIGN_TABLE if table is None else table).v


# ------------------------------------------------- CAD parameterisation of the wrap
# The B-spline face uses w in [0, 1] (0 upper cut edge, 1/2 lip, 1 lower cut edge).
# theta = pi w near the cut edges (so the edge rows of poles stay exact) and around
# the lip; in between, theta(w, q) is shifted so that the start of the roll
# (x = x0(q)) falls on the fixed knot line w = W_ON wherever the roll starts
# late (q > ~0.6). The geometry is unchanged; only its parameterisation moves.
W_A, W_ON, W_C = 0.20, 0.344, 0.45
GAMMA_BLEND = (0.55, 0.68)


def _s5(s):
    s = np.clip(s, 0.0, 1.0)
    return s**3 * (10 - 15 * s + 6 * s * s)


def _bump(w):
    w = np.asarray(w, float)
    up = _s5((w - W_A) / (W_ON - W_A))
    down = _s5((W_C - w) / (W_C - W_ON))
    return np.where(w <= W_ON, up, down)


def theta_of_w(w, q):
    w, q = np.broadcast_arrays(np.asarray(w, float), np.asarray(q, float))
    x0 = np.clip(x0_of_q(q), 0.0, 0.999)
    gamma = _s5((q - GAMMA_BLEND[0]) / (GAMMA_BLEND[1] - GAMMA_BLEND[0]))
    delta = gamma * (np.arcsin(x0) - np.pi * W_ON)
    half = np.minimum(w, 1 - w)
    th = np.pi * half + delta * _bump(half)
    return np.where(w <= 0.5, th, np.pi - th)


def _x0_jet(qj, table):
    from scipy.interpolate import CubicSpline
    qn, xn = table
    cs = CubicSpline(qn, xn, bc_type="natural")
    v = qj.v
    return qj.chain(cs(v), cs(v, 1), cs(v, 2))


def _g(xi, b):
    """(1 - xi^b) / (1 - xi) for jets xi in [0, 1] and b; finite at xi = 1 (value b)."""
    v = np.clip(xi.v, 1e-6, 1.0)
    xi = Jet(v, xi.d, xi.h)
    e = 1 - xi
    near = e.v < 1e-4
    # far: (1 - exp(b log xi)) / (1 - xi)
    es = Jet(np.where(near, 1.0, e.v), e.d, e.h)
    xs = Jet(np.where(near, 0.5, xi.v), xi.d, xi.h)
    far = (1 - jexp(b * jlog(xs))) / es
    # near: series in e = 1 - xi
    ser = b - b * (b - 1) * 0.5 * e + b * (b - 1) * (b - 2) * (1 / 6) * e * e - b * (b - 1) * (b - 2) * (b - 3) * (1 / 24) * e * e * e
    return parent.jwhere(near, ser, far)


def _h(x, x1):
    """h(x) = int_0^x (1 - (s/x1)^2)^3 ds up to x1, constant (16/35 x1) after."""
    u = x / x1
    inside = x1 * (u - u * u * u + 0.6 * u**5 - (1 / 7) * u**7)
    after = x1 * (16.0 / 35.0) + 0 * x
    return parent.jwhere(u.v < 1, inside, after)


def surface_jet(q, theta, table=None, x0_const=None, b_const=None, ku_const=None):
    """Point and exact first/second derivatives in (q, theta). Returns a dict like parent.skin_jet."""
    if table is None:
        table = DESIGN_TABLE
    qj, th = parent._seed(q, theta)
    U = parent._U(qj)
    t = parent._t(qj)
    D = 240.0 * qj * (1 - qj)
    s = parent.jsin(th)
    co = parent.jcos(th)
    r = ATT + D * s
    c = ROOT_CHORD - TAPER * r
    if x0_const is not None:
        x0 = Jet.constant(np.broadcast_to(np.asarray(x0_const, float), qj.v.shape).copy(), 2)
    else:
        x0 = _x0_jet(qj, table)
    b = b_jet(qj) if b_const is None else _const(b_const, qj)
    ku = ku_jet(qj) if ku_const is None else _const(ku_const, qj)
    xi = (s - x0) / (1 - x0)
    # outside the flat run: zeta = cos * sqrt((1 + xi + xi^2) / ((1 + sin)(1 - x0)))  (no cancellation at the lip)
    # 1 - xi^b = (1 - xi) g(xi), 1 - xi = cos^2 / ((1 + sin)(1 - x0)): no cancellation at the lip
    round_ = co * parent.jsqrt(_g(xi, b) / ((1 + s) * (1 - x0)))
    sign = np.where(co.v >= 0, 1.0, -1.0)
    flat = Jet.constant(sign, 2)
    zeta = parent.jwhere(s.v > x0.v, round_, flat)
    x1 = Jet(np.maximum(x0.v, 0.05), x0.d, x0.h)
    C = c + ku * TAPER * D * (s - _h(s, x1))
    Z = C * t * zeta
    X = SWEEP * r + c * U
    return parent._pack([X, r + 0 * qj, Z], "q", "theta")


def _const(value, like):
    return Jet.constant(np.broadcast_to(np.asarray(value, float), like.v.shape).copy(), like.n)


def point(q, theta, **kw):
    return surface_jet(q, theta, **kw)["point"]


def principal(jet):
    """Principal curvatures from a jet dict: k1, k2, kmax = max |k|, mean H, Gauss K, unit normal n."""
    Su, Sv = jet["d_q"], jet["d_theta"]
    Suu, Suv, Svv = jet["d_q_q"], jet["d_q_theta"], jet["d_theta_theta"]
    n = np.cross(Su, Sv)
    n = n / np.linalg.norm(n, axis=-1, keepdims=True)
    E = (Su * Su).sum(-1); F = (Su * Sv).sum(-1); G = (Sv * Sv).sum(-1)
    L = (Suu * n).sum(-1); M = (Suv * n).sum(-1); N = (Svv * n).sum(-1)
    den = E * G - F * F
    H = (E * N - 2 * F * M + G * L) / (2 * den)
    K = (L * N - M * M) / den
    disc = np.sqrt(np.maximum(H * H - K, 0))
    k1, k2 = H + disc, H - disc
    kmax = np.maximum(np.abs(k1), np.abs(k2))
    return {"k1": k1, "k2": k2, "kmax": kmax, "H": H, "K": K, "n": n}
