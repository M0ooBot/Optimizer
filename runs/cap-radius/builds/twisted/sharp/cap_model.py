"""Analytic rounded tip cap for the twisted/sharp parent (right tip; the left tip mirrors it in XZ).

Construction (see README.md):
  profiles   upper attachment cut edge (w = 0) -> locked 3D rim (w = 1) -> lower cut edge (w = 2)
  guides     one spanwise wrap guide per chord parameter p (the parent's native chord parameter s = u = p)
  surface    S(p, w) = place(Y, airfoil(p, sigma)), i.e. the parent's own section formula with the
             unit thickness side (+1 upper, -1 lower) replaced by a signed thickness fraction sigma,
             placed at the station Y with the parent's taper, twist, dihedral and sweep continued.
             Y = 440 + L(p) eta(w), L(p) = 60 sin(pi p) (the locked outline reach),
             eta = 1 - v^2, v = |1 - w|, sigma = sign(1 - w) F(v, p) (lip profile below).
  At w = 0 / 2 the surface is the parent skin (sigma = +-1, Y = 440); at w = 1 it is the locked rim
  (sigma = 0 is the camber line, Y = 440 + 60 sin(pi p), native u = p).  F = 1 + O(eta^3), so the
  cap continues the parent skin to second order (G2) at the attachment.  The roll share E(p) is the
  design variable for the printable radius rule.

Lengths in mm.  Needs NumPy and the parent definition (parent.py).
"""
import importlib.util
import pathlib
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[5]
PARENT_DIR = ROOT / "runs/cap-radius/parents/twisted/sharp"
_spec = importlib.util.spec_from_file_location("parent_twisted_sharp", PARENT_DIR / "parent.py")
P = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(P)
Jet = P.Jet

CUT = P.CUT
REACH = P.REACH

# ---- lip profile ------------------------------------------------------------------------------
# Along each guide, e = 1 - eta = v^2 is the fraction of the reach still to go to the rim.
# The thickness fraction is F = sqrt(G(x)), x = e / E(p), F = 1 for x >= 1, with
#   G(x) = 1 - (1 - x)^K (1 + A x + B x^2 + C x^3)
# G(1) = 1 with G'(1) = G''(1) = 0, so the cap meets the slab (or, for E = 1, the parent skin) G2.
# Near the rim F^2 ~ (3 - A) x - (3 - 3A) x^2: A_SHAPE = sqrt(12) - 3 makes the lip osculate a circle
# to second order; a little more makes the lip the tightest point of the turn (progressive).
# E(p) is the share of the reach used by the roll: 1 (spread over the whole reach) forward,
# falling smoothly to E_TE aft, where the section is too thin for a spread roll to keep 0.4 mm.
A_SHAPE = 1.4631
B_SHAPE = 1.3569
C_SHAPE = 1.9858
K_FLAT = 4      # G = 1 - (1 - x)^K_FLAT q(x): the roll meets the slab C^(K_FLAT - 1)
E_TE = 0.0837
P_A, P_B = 0.32, 0.88


def _smooth(t):
    t = P.jwhere(t.v <= 0, P._const(0.0, t), P.jwhere(t.v >= 1, P._const(1.0, t), t))
    return t**3 * (10 - 15 * t + 6 * t * t)


def E_of_p(pj):
    return E_TE + (1.0 - E_TE) * (1 - _smooth((pj - P_A) * (1.0 / (P_B - P_A))))


def E_value(p):
    return E_of_p(Jet.constant(np.asarray(p, float), 0)).v


def g_coefficients():
    """Coefficients (ascending) of g(x) = G(x) / x, G = 1 - (1 - x)^3 (1 + A x + B x^2)."""
    q = np.polynomial.Polynomial([1.0, A_SHAPE, B_SHAPE, C_SHAPE])
    G = 1 - np.polynomial.Polynomial([1.0, -1.0]) ** K_FLAT * q
    return G.coef[1:]


def _g(x):
    c = g_coefficients()
    out = P._const(c[-1], x)
    for k in range(len(c) - 2, -1, -1):
        out = out * x + c[k]
    return out


def profile_F(v, E):
    """Thickness fraction F(v) = v sqrt(g(v^2/E)/E) inside the roll, 1 outside (jets)."""
    x = v * v / E
    inside = x.v < 1.0
    xi = P.jwhere(inside, x, P._const(0.5, x))
    Fin = v * P.jsqrt(_g(xi) / E)
    return P.jwhere(inside, Fin, P._const(1.0, x))


def _cap_coords(pj, wj):
    upper = wj.v <= 1.0
    v = P.jwhere(upper, 1 - wj, wj - 1)
    eta = 1 - v * v
    F = profile_F(v, E_of_p(pj))
    sigma = P.jwhere(upper, F, -F)
    L = REACH * P.jsin(pj * np.pi)
    Y = CUT + L * eta
    xn, zn = P._airfoil(pj, sigma)
    return P._place(Y, xn, zn)


def cap_jet(p, w):
    """Right-tip cap point and exact derivatives in (p, w)."""
    pj, wj = P._seed(p, w)
    return P._pack(_cap_coords(pj, wj), "p", "w")


def cap(p, w):
    return cap_jet(p, w)["point"]


def principal(jet):
    """Principal curvatures (k_max >= k_min by absolute value) and unit normal from a jet dict."""
    Su, Sv = jet["d_p"], jet["d_w"]
    Suu, Suv, Svv = jet["d_p_p"], jet["d_p_w"], jet["d_w_w"]
    N = np.cross(Su, Sv)
    nn = np.linalg.norm(N, axis=-1, keepdims=True)
    N = N / nn
    E = (Su * Su).sum(-1); Fm = (Su * Sv).sum(-1); G = (Sv * Sv).sum(-1)
    L = (Suu * N).sum(-1); M = (Suv * N).sum(-1); Nn = (Svv * N).sum(-1)
    det = E * G - Fm * Fm
    K = (L * Nn - M * M) / det
    H = (E * Nn - 2 * Fm * M + G * L) / (2 * det)
    disc = np.sqrt(np.maximum(H * H - K, 0))
    k1, k2 = H + disc, H - disc
    kmax = np.maximum(np.abs(k1), np.abs(k2))
    return {"kmax": kmax, "k1": k1, "k2": k2, "normal": N}


def attachment_thickness(xi):
    """Parent attachment thickness (mm) at chord fraction xi (parent.section_thickness at Y = 440)."""
    xi = np.clip(np.asarray(xi, float), 0.0, 1.0)
    return np.nan_to_num(P.section_thickness(xi, CUT), nan=0.0)


def chord_fraction(p):
    """Chord fraction of guide p: the camber-line x of native parameter u = p (the rim point's fraction)."""
    return np.sin(np.pi * np.asarray(p, float) / 2) ** 2


def required_radius(p):
    """Printable rule: 0.4 mm where the attachment section is at least 0.8 mm thick, else half the thickness."""
    t = attachment_thickness(chord_fraction(p))
    return np.where(t >= 0.8, 0.4, t / 2), t
