"""Analytic definition of the rounded tip cap for fin40/blunt.

Notation (lengths in mm, z = height, w = z - 240, h = w / 40):
  outline   xL(z) = m - a sqrt(1 - h^3)  (LE branch),  xT(z) = m + a sqrt(1 - h^3)  (TE branch)
            m = 80 + 27.5 z/240,  a = 80 - 37.5 z/240              (parent.tip_outline)
  parent    Y^2 = K2 (x - LE)(E - x),  K2 = 0.12^2, LE = 65 z/240, E = LE + c(z) scale(z)
            (12 % ellipse, scaled by 'scale' and cut by the trailing-edge plane x = TE(z))
            E - TE = dE(z) ~ 0.0327 mm at the attachment.

The cap surface is Y^2 = G(x, z) with

  G = K2 [ (x - xL)(xT - x) + d(w) (x - xL) (1 - b(w) (l / (r + l))^3) ],   r = xT - x,

  d(w) = dE(z) mu(w): the parent's 'scale' excess, faded out slowly (mu = 1 - S(w/W2));
  b(w) = S(w/WE): closes the flat trailing-edge face; the weight (l/(r + l))^3 confines it
         to a few l of the TE branch (l = 3 mm);
  S(u) = 10u^3 - 15u^4 + 6u^5 (quintic smoothstep: S = O(u^3) at 0, C2 at 1).

At w = 0 the bracket equals (x - LE)(E - x) + O(w^3) for every x (outline ends are O(h^3)
from the parent's linear LE/TE lines; mu and b are 1 - O(w^3) and O(w^3)), so the cap
matches the parent skin in position, tangent plane and curvature (G2).  For w >= WE the
bracket vanishes exactly on both outline branches, so the lip is the locked outline.
Below that the trailing end of each section lies aft of the trailing-edge plane and is cut
by it: the flat trailing-edge face continues a little way into the cap and closes where
the section end crosses the plane (about 1.29 mm above the attachment).

Each horizontal section factorises as G = K2 (x - xL)(XE - x) q with XE the section's
trailing zero and q = C(r)/(r + l)^3 > 0 smooth, so with x = xL + (XE - xL)(1 - cos t)/2 the skin is
Y = 0.12 A sin t sqrt(q), A = (XE - xL)/2: a smooth parameterisation for the B-spline fit.
"""
import importlib.util
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[5]
PARENT_DIR = ROOT / "runs" / "cap-radius" / "parents" / "fin40" / "blunt"

_spec = importlib.util.spec_from_file_location("fin40_blunt_parent", PARENT_DIR / "parent.py")
parent = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(parent)
Jet, jsin, jcos, jsqrt, jwhere = parent.Jet, parent.jsin, parent.jcos, parent.jsqrt, parent.jwhere

Z0 = parent.ATTACHMENT          # 240
REACH = parent.REACH            # 40
K2 = parent.THICKNESS_RATIO**2  # 0.0144
WE = 1.5                        # height over which the flat trailing-edge face is closed (mm)
W2 = 20.0                       # height over which the parent's scale excess fades out (mm)
ELL = 3.0                       # chordwise width of the trailing-edge closure term (mm)


def _smoothstep(w, length):
    """S(w/length) clamped to [0, 1] as a jet."""
    u = w * (1.0 / length)
    u = jwhere(u.v < 1.0, u, u * 0.0 + 1.0)
    u = jwhere(u.v > 0.0, u, u * 0.0)
    return 10 * u**3 - 15 * u**4 + 6 * u**5


def te_plane_x(z):
    return 160.0 - 10.0 * z / 240.0


def _excess(z):
    """dE(z) = E - TE from the parent's keep-chord cut, as a jet."""
    theta_end, scale = parent._cut(z)
    c = 160.0 - 75.0 * z / 240.0
    return 65.0 * z / 240.0 + c * scale - te_plane_x(z)


def _closure_poly(d, b):
    """Coefficients of P(r) = (r + d)(r + l)^3 - d b l^3 (monic quartic), lowest first."""
    l = ELL
    return [d * (l**3) * (1.0 - b), l**3 + 3 * l * l * d, 3 * l * l + 3 * l * d, 3 * l + d]


def ends(z, sq=None):
    """Jets: xL, xT, XE (trailing zero of the section), d, b, r1, cubic cofactor at height jet z.

    With r = xT - x the section is G = K2 (x - xL) P(r) / (r + l)^3; r1 is the root of P in
    [-d, 0] (XE = xT - r1) and P(r) = (r - r1) C(r) with C > 0 on the section.
    """
    h = (z - Z0) * (1.0 / REACH)
    m = 80.0 + 27.5 * z / 240.0
    a = 80.0 - 37.5 * z / 240.0
    if sq is None:
        sq = jsqrt(1.0 - h**3)
    xL = m - a * sq
    xT = m + a * sq
    w = z - Z0

    def dbz(zz):
        ww = zz - Z0
        return _excess(zz) * (1.0 - _smoothstep(ww, W2)), _smoothstep(ww, WE)

    d, b = dbz(z)
    # real root r1 in [-d, 0] by bisection + Newton, then exact derivatives (implicit function)
    dv, bv = np.asarray(d.v, float), np.asarray(b.v, float)
    poly = lambda r: (r + dv) * (r + ELL) ** 3 - dv * bv * ELL**3
    r1v = parent._bisect(poly, -dv - 1e-15, np.zeros_like(dv) + 1e-15)

    def residual(unknowns, outer):
        dd, bb = dbz(outer)
        r = unknowns[0]
        return [(r + dd) * (r + ELL) ** 3 - dd * bb * ELL**3]

    if np.ndim(z.v) == 0:
        r1v = np.asarray(r1v)
    r1 = parent._implicit(residual, [r1v], z)[0] if z.n > 0 else Jet.constant(r1v, 0)
    c0, c1, c2, c3 = _closure_poly(d, b)
    b2 = c3 + r1
    b1 = c2 + r1 * b2
    b0 = c1 + r1 * b1
    return {"xL": xL, "xT": xT, "XE": xT - r1, "d": d, "b": b, "r1": r1, "C": (b0, b1, b2),
            "TE": te_plane_x(z)}


def implicit_G(x, z):
    """G(x, z) as a jet (surface: Y^2 = G).  Smooth everywhere including the crown."""
    h = (z - Z0) * (1.0 / REACH)
    m = 80.0 + 27.5 * z / 240.0
    a = 80.0 - 37.5 * z / 240.0
    base = a * a * (1.0 - h**3) - (x - m) ** 2          # (x - xL)(xT - x)
    w = z - Z0
    if np.any(np.asarray(w.v) < W2):
        hh = jwhere(h.v < 0.9, h, h * 0.0 + 0.9)       # the extra term is zero above W2 anyway
        sq = jsqrt(1.0 - hh**3)
        xL = m - a * sq
        xT = m + a * sq
        d = _excess(z) * (1.0 - _smoothstep(w, W2))
        b = _smoothstep(w, WE)
        r = xT - x
        base = base + d * (x - xL) * (1.0 - b * (ELL / (r + ELL)) ** 3)
    return K2 * base


def surface_jet(t, v):
    """Cap skin (+Y side) point and exact derivatives in (t, v).

    t in [0, pi]: 0 = LE-branch lip, pi = the section's trailing zero XE.
    v in [0, 1]: z = 240 + 40 v (2 - v); v = 1 is the crown (collapsed).
    The trailing-edge plane trim (x <= TE(z)) is not applied here.
    """
    tj, vj = parent._seed(t, v)
    h = vj * (2.0 - vj)
    z = Z0 + REACH * h
    sq = (1.0 - vj) * jsqrt(1.0 + h + h * h)
    e = ends(z, sq)
    A = (e["XE"] - e["xL"]) * 0.5
    x = e["xL"] + A * (1.0 - jcos(tj))
    r = e["xT"] - x
    b0, b1, b2 = e["C"]
    q = (r**3 + b2 * r * r + b1 * r + b0) / (r + ELL) ** 3
    Y = parent.THICKNESS_RATIO * A * jsin(tj) * jsqrt(q)
    return parent._pack([x, Y, z], "t", "v")


def surface_points(t, v):
    """Vectorised points (no derivatives) of the +Y skin on a (t, v) grid."""
    t, v = np.broadcast_arrays(np.asarray(t, float), np.asarray(v, float))
    vj = Jet.constant(v, 1)
    h = vj * (2.0 - vj)
    z = Z0 + REACH * h
    sq = (1.0 - vj) * jsqrt(1.0 + h + h * h)
    e = ends(z, sq)
    xL, XE, xT = e["xL"].v, e["XE"].v, e["xT"].v
    b0, b1, b2 = (c.v for c in e["C"])
    A = 0.5 * (XE - xL)
    x = xL + A * (1.0 - np.cos(t))
    r = xT - x
    q = (r**3 + b2 * r * r + b1 * r + b0) / (r + ELL) ** 3
    return np.stack([x, parent.THICKNESS_RATIO * A * np.sin(t) * np.sqrt(q), z.v], -1)


def ends_values(z):
    e = ends(Jet.constant(np.asarray(z, float), 1))
    return e["xL"].v, e["XE"].v, e["TE"].v


def te_closure_height():
    """Height above the attachment where the section's trailing zero XE crosses the TE plane."""
    f = lambda w: ends_values(Z0 + w)[1] - ends_values(Z0 + w)[2]
    lo, hi = 1e-6, WE
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def te_arch_halfwidth(z):
    """Half-width of the flat trailing-edge face at height z (NaN above its closure)."""
    zj = Jet.constant(np.asarray(z, float), 1)
    xj = Jet.constant(te_plane_x(np.asarray(z, float)), 1)
    g = implicit_G(xj, zj).v
    return np.where(g >= 0, np.sqrt(np.abs(g)), np.nan)
