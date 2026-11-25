"""Analytic definition of the fin5/blunt rounded tip cap.

The fin is symmetric about Y = 0, so each skin is a graph Y = +/- sqrt(G(X, Z)).
For the parent, G_par = 0.0144 a^2 (1 + xi)(xi_f - xi), with
    a(Z)  = half chord = 80 - 37.5 Z/240,  m(Z) = mid chord = 80 + 27.5 Z/240,
    xi    = (X - m)/a  (-1 at the leading edge, +1 at the trailing edge),
    xi_f  = 2*scale(Z) - 1  (the virtual end of the keep-chord scaled ellipse).
The cap (240 <= Z <= 245) is
    G_cap = G_par(X, Z) * psi,   psi = 1 - h^3/(1 - xi^2),   h = (Z - 240)/5.
psi = 1 at the attachment to second order in h (so the cap is curvature
continuous with both skins) and psi = 0 exactly on the locked outline
xi^2 + h^3 = 1, which is parent.tip_outline written in the same coordinates.
The lip (Y = 0) is therefore the locked outline.  The implicit surface
F = Y^2 - G_cap = 0 is smooth across the lip, so the two sides join with a
common tangent plane and matching curvature.
"""
import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[5]
PARENT_DIR = ROOT / "runs/cap-radius/parents/fin5/blunt"
_spec = importlib.util.spec_from_file_location("fin5_blunt_parent", PARENT_DIR / "parent.py")
parent = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(parent)
Jet = parent.Jet

Z0 = parent.ATTACHMENT
REACH = parent.REACH
K = parent.THICKNESS_RATIO ** 2  # (B/A)^2 = 0.0144


EPS = 0.25  # wrap-profile ramp parameter (see psi below); chosen in the design study (README)


def _cbrt(x):
    r = np.cbrt(x.v)
    return x.chain(r, 1.0 / (3.0 * r * r), -2.0 / (9.0 * r**5))


def _G_jets(X, Z, eps=None):
    """Jets (in X, Z) of G_par and G_cap."""
    eps = EPS if eps is None else eps
    Xj, Zj = parent._seed(X, Z)
    a = 80.0 - (37.5 / Z0) * Zj
    m = 80.0 + (27.5 / Z0) * Zj
    xi = (Xj - m) / a
    _, scale = parent._cut(Zj)
    xi_f = scale * 2.0 - 1.0
    gpar = a * a * (xi + 1.0) * (xi_f - xi) * K
    h = (Zj - Z0) / REACH
    x3 = h ** 3 / (1.0 - xi * xi)          # (h/h_o)^3, h_o = outline height at this xi
    if eps == float("inf"):
        psi = 1.0 - x3
    else:
        x = h * _cbrt(1.0 / (1.0 - xi * xi))  # h/h_o, smooth in h
        psi = 1.0 - x3 * (1.0 + eps) / (x + eps)
    return gpar, gpar * psi, xi, h


def G(X, Z, eps=None):
    """G_cap and its first/second derivatives: dict with v, X, Z, XX, XZ, ZZ."""
    _, g, _, _ = _G_jets(X, Z, eps)
    return {"v": g.v, "X": g.d[..., 0], "Z": g.d[..., 1], "XX": g.h[..., 0, 0],
            "XZ": g.h[..., 0, 1], "ZZ": g.h[..., 1, 1]}


def xi_h(X, Z):
    X, Z = np.broadcast_arrays(np.asarray(X, float), np.asarray(Z, float))
    a = 80.0 - (37.5 / Z0) * Z
    m = 80.0 + (27.5 / Z0) * Z
    return (X - m) / a, (Z - Z0) / REACH


def principal_curvatures(X, Z, eps=None):
    """Principal curvatures (k1 >= k2 in magnitude order kmax=|k|max) of F = Y^2 - G = 0 at Y = +sqrt(G)."""
    g = G(X, Z, eps)
    Y = np.sqrt(np.clip(g["v"], 0, None))
    grad = np.stack([-g["X"], 2 * Y, -g["Z"]], -1)
    H = np.zeros(grad.shape + (3,))
    H[..., 0, 0] = -g["XX"]; H[..., 0, 2] = H[..., 2, 0] = -g["XZ"]; H[..., 2, 2] = -g["ZZ"]; H[..., 1, 1] = 2.0
    n2 = np.einsum("...i,...i", grad, grad)
    n = np.sqrt(n2)
    # project Hessian onto tangent plane; shape operator = P H P / |grad|
    unit = grad / n[..., None]
    P = np.eye(3) - unit[..., :, None] * unit[..., None, :]
    S = P @ H @ P / n[..., None, None]
    ev = np.linalg.eigvalsh(S)  # one eigenvalue ~0 (normal direction)
    # drop the eigenvalue closest to zero along the normal: use the two with largest |.| ... safer: sort by |ev|
    idx = np.argsort(np.abs(ev), -1)
    ev = np.take_along_axis(ev, idx, -1)[..., 1:]
    return ev, Y, unit


def attachment_thickness(X):
    """Parent attachment-section thickness at the chord fraction of X (nan outside)."""
    xi = (np.asarray(X, float) - 65.0) / 85.0
    return parent.section_thickness(np.clip(xi, 0, 1), Z0)


def required_radius(X):
    t = attachment_thickness(X)
    return np.where(t >= 0.8, 0.4, t / 2)


def G_value(X, Z, eps=None):
    """Values only (fast) of G_cap at points (X, Z) inside the cap domain."""
    eps = EPS if eps is None else eps
    X, Z = np.broadcast_arrays(np.asarray(X, float), np.asarray(Z, float))
    a = 80.0 - (37.5 / Z0) * Z
    m = 80.0 + (27.5 / Z0) * Z
    xi = (X - m) / a
    _, scale = parent._cut(Jet.constant(Z, 0))
    xi_f = 2.0 * scale.v - 1.0
    gpar = K * a * a * (1.0 + xi) * (xi_f - xi)
    h = (Z - Z0) / REACH
    om = np.clip(1.0 - xi * xi, 1e-300, None)
    x3 = h**3 / om
    if eps == float("inf"):
        psi = 1.0 - x3
    else:
        x = h / np.cbrt(om)
        psi = 1.0 - x3 * (1.0 + eps) / (x + eps)
    return gpar * psi


# ---------------------------------------------------------------------------
# Perimeter and wrap-section construction (the "guides")
# ---------------------------------------------------------------------------

def outline_point(xi_p):
    """Locked outline point at outline coordinate xi_p in [-1, 1] (-1 = LE end, +1 = TE end)."""
    xi_p = np.asarray(xi_p, float)
    h = np.cbrt(np.clip(1.0 - xi_p**2, 0, None))
    Z = Z0 + REACH * h
    a = 80.0 - (37.5 / Z0) * Z
    m = 80.0 + (27.5 / Z0) * Z
    return np.stack([m + a * xi_p, np.zeros_like(Z), Z], -1)


def _cumlen(P):
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=-1))])


class Perimeter:
    """Arc-length parametrisation of the locked outline (LE end -> crown -> TE end)."""

    def __init__(self, n=400001):
        phi = np.linspace(0.0, np.pi, n)
        xi = -np.cos(phi)
        # cluster towards both ends, where h ~ (1 - xi^2)^(1/3) is steep in xi
        t = np.linspace(0, 1, n)
        xi = -np.cos(np.pi * (3 * t**2 - 2 * t**3))
        self.xi = xi
        self.L = _cumlen(outline_point(xi))
        self.total = self.L[-1]

    def xi_at(self, frac):
        return np.interp(np.asarray(frac) * self.total, self.L, self.xi)


class CutEdge:
    """Arc-length fraction of the +Y attachment edge as a function of the chord parameter s."""

    def __init__(self, n=200001):
        self.s = np.linspace(0, 1, n)
        self.L = _cumlen(parent.skin(self.s, Z0, +1))
        self.total = self.L[-1]

    def frac(self, s):
        return np.interp(s, self.s, self.L) / self.total


_PER = None
_CUT = None


def _tables():
    global _PER, _CUT
    if _PER is None:
        _PER, _CUT = Perimeter(), CutEdge()
    return _PER, _CUT


def lip_point(s):
    """Perimeter point paired with attachment chord parameter s (equal arc-length fractions)."""
    per, cut = _tables()
    return outline_point(per.xi_at(cut.frac(s)))


def section(s, n_tau=6000):
    """Exact points of the +Y half of the wrap section at chord parameter s.

    The section lies in the plane that contains the Y axis and the side-view
    line from the attachment point A(s) to its perimeter point P(s).
    Returns tau (0 at A, 1 at P) and XYZ points (Y >= 0).
    """
    A = parent.skin(np.array([s]), Z0, +1)[0]
    P = lip_point(np.array([s]))[0]
    tau = np.unique(np.concatenate([np.linspace(0, 1, n_tau), 1 - np.geomspace(1e-15, 0.5, n_tau)]))
    X = A[0] + tau * (P[0] - A[0])
    Z = A[2] + tau * (P[2] - A[2])
    g = G_value(X, Z)
    Y = np.sqrt(np.clip(g, 0, None))
    Y[0] = A[1]
    Y[-1] = 0.0
    return tau, np.stack([X, Y, Z], -1), A, P


# ---------------------------------------------------------------------------
# Trailing-edge end in local coordinates (well conditioned near the TE segment)
#   p = X_TE(Z) - X  (distance forward of the parent's trailing-edge line)
#   z = Z - 240
# ---------------------------------------------------------------------------
def _scale_minus_one(Z):
    theta_end, _ = parent._cut(Jet.constant(np.asarray(Z, float), 0))
    te = theta_end.v
    xi_end = 0.5 * (1.0 - np.cos(te))
    return 0.5 * (1.0 + np.cos(te)) / xi_end


def G_local(p, z, eps=None):
    """G_cap from local TE coordinates; identical to G_value but without cancellation."""
    eps = EPS if eps is None else eps
    p, z = np.broadcast_arrays(np.asarray(p, float), np.asarray(z, float))
    Z = Z0 + z
    a = 80.0 - (37.5 / Z0) * Z
    w = p / a                       # 1 - xi
    om = w * (2.0 - w)              # 1 - xi^2
    gpar = K * a * a * (2.0 - w) * (2.0 * _scale_minus_one(Z) + w)
    h = z / REACH
    om = np.clip(om, 1e-300, None)
    x3 = h**3 / om
    x = h / np.cbrt(om)
    psi = 1.0 - x3 * (1.0 + eps) / (x + eps)
    return gpar * np.clip(psi, 0.0, None)


def X_te(Z):
    return 160.0 - 10.0 * np.asarray(Z, float) / Z0


class TEtables:
    """Arc lengths measured from the trailing-edge end (cut edge and outline)."""

    def __init__(self):
        # cut edge: s from 1 down to 0.5, clustered at 1
        d = np.unique(np.concatenate([[0.0], np.geomspace(1e-13, 0.5, 40001)]))
        self.ds = d
        s = 1.0 - d
        te = parent.cut(Z0)["theta_end"]
        c = float(parent.chord(Z0)); sc = float(parent.cut(Z0)["scale"])
        th = s * te
        q = c * sc * np.sin(0.5 * (te + th)) * np.sin(0.5 * (te - th))   # X_TE - X_A (stable)
        y = 0.5 * parent.THICKNESS_RATIO * c * sc * np.sin(th)
        self.Lc = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(q), np.diff(y)))])
        # outline TE branch: h from 0 up to 0.95, clustered at 0
        h = np.unique(np.concatenate([[0.0], np.geomspace(1e-13, 0.95, 60001)]))
        self.h = h
        Z = Z0 + REACH * h
        a = 80.0 - (37.5 / Z0) * Z
        p = a * h**3 / (1.0 + np.sqrt(1.0 - h**3))      # X_TE - X_outline (stable)
        dX = -(10.0 / Z0) * np.diff(Z) - np.diff(p)
        self.Lo = np.concatenate([[0.0], np.cumsum(np.hypot(dX, np.diff(Z)))])

    def q_y_A(self, s):
        te = parent.cut(Z0)["theta_end"]
        c = float(parent.chord(Z0)); sc = float(parent.cut(Z0)["scale"])
        th = s * te
        return (c * sc * np.sin(0.5 * (te + th)) * np.sin(0.5 * (te - th)),
                0.5 * parent.THICKNESS_RATIO * c * sc * np.sin(th))

    def h_P(self, s, cut_total, out_total):
        Lc = np.interp(1.0 - s, self.ds, self.Lc)
        return np.interp(Lc / cut_total * out_total, self.Lo, self.h)


_TET = None


def section_te(s, n_tau=6000):
    """Like section(), in local TE coordinates; for s close to 1.

    Returns sig (1 at A, 0 at P), XYZ, A, P and the local data (q, z, qA, qP, zP)."""
    global _TET
    if _TET is None:
        _TET = TEtables()
    per, cut = _tables()
    qA, yA = _TET.q_y_A(s)
    hP = float(_TET.h_P(s, cut.total, per.total))
    zP = REACH * hP
    aP = 80.0 - (37.5 / Z0) * (Z0 + zP)
    qP = aP * hP**3 / (1.0 + np.sqrt(1.0 - hP**3))
    # sig = 1 - tau, the fraction of the side-view line still to go to P; sampled down to
    # 1e-34 so that the sliver next to the outline is resolved even next to the TE segment
    sig = np.unique(np.concatenate([np.linspace(0, 1, n_tau), np.geomspace(1e-34, 0.5, 2 * n_tau)]))[::-1]
    tau = sig  # returned in place of tau (decreasing from 1 at A to 0 at P)
    q = qP + sig * (qA - qP)
    z = (1.0 - sig) * zP
    Y = np.sqrt(G_local(q, z))
    Y[0] = yA
    Y[-1] = 0.0
    Z = Z0 + z
    X = X_te(Z) - q
    A = np.array([X_te(Z0) - qA, yA, Z0])
    P = np.array([X_te(Z0 + zP) - qP, 0.0, Z0 + zP])
    return tau, np.stack([X, Y, Z], -1), A, P, (q, z, qA, qP, zP)
