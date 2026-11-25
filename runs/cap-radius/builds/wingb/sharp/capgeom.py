"""Shared geometry for the Wing B (sharp) rounded tip cap.

Coordinates: STEP axes (X = -span, Y up, Z forward). The cap is described in
the extended planform coordinates (S, xi): S = projected span (X = -S), xi =
chord fraction of the parent's analytic section frame continued past the
attachment (leading-edge reference Z = -(35 f + 45 f^2), chord 200 - 125 f,
f = S/600). The locked outline uses exactly this continuation for its
leading- and trailing-edge ends.
"""
import importlib.util
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[5]
PARENT_DIR = ROOT / "runs/cap-radius/parents/wingb/sharp"
BUILD_DIR = Path(__file__).resolve().parent


@lru_cache(maxsize=1)
def parent():
    spec = importlib.util.spec_from_file_location("wingb_sharp_parent", PARENT_DIR / "parent.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S_ATT = 540.0


def frame(S):
    """Leading-edge Z, chord, and their first/second S-derivatives of the (continued) section frame."""
    S = np.asarray(S, float)
    f = S / 600.0
    zle = -(35.0 * f + 45.0 * f * f)
    zle1 = -(35.0 + 90.0 * f) / 600.0
    zle2 = np.full_like(S, -90.0 / 600.0**2)
    c = 200.0 - 125.0 * f
    c1 = np.full_like(S, -125.0 / 600.0)
    return zle, zle1, zle2, c, c1


def plan_xz(S, xi):
    zle, _, _, c, _ = frame(S)
    return -np.asarray(S, float), zle - np.asarray(xi, float) * c


def xi_of(S, Z):
    zle, _, _, c, _ = frame(S)
    return (zle - Z) / c


def _u_at_z(z_target, v, side, nose):
    """Native u on one side of the section at native v whose Z equals z_target (bisection)."""
    par = parent()
    skin = par.surfaces()[0]
    end = par.SURFACE_DATA["u_end_upper"] if side > 0 else par.SURFACE_DATA["u_end_lower"]
    z_target = np.atleast_1d(z_target)
    lo = np.full(z_target.shape, nose)
    hi = np.full(z_target.shape, end)
    for _ in range(64):
        mid = 0.5 * (lo + hi)
        z = np.array([skin.Value(float(m), v).Z() for m in mid])
        ahead = z > z_target
        lo, hi = np.where(ahead, mid, lo), np.where(ahead, hi, mid)
    return 0.5 * (lo + hi)


def attachment_jet(xi, S=S_ATT):
    """Upper/lower Y and their first and second S-derivatives at fixed frame chord fraction xi.

    Exact chain rule through the parent's B-spline derivatives (parent.skin_jet).
    Returns dict side -> (Y, Y_S, Y_SS) arrays, plus u per side.
    """
    par = parent()
    xi = np.atleast_1d(np.asarray(xi, float))
    zle, zle1, zle2, c, c1 = frame(S)
    zt = zle - xi * c
    zt1 = zle1 - xi * c1
    zt2 = zle2  # c'' = 0
    v = float(par.v_at_span(S))
    nose = float(par.nose_u(S))
    out = {}
    for side in (+1, -1):
        u = _u_at_z(zt, v, side, nose)
        j = par.skin_jet(u, np.full(u.shape, float(S)))
        P, Pu, PS = j["point"], j["d_chord"], j["d_span"]
        Puu, PuS, PSS = j["d_chord_chord"], j["d_chord_span"], j["d_span_span"]
        Zu, ZS = Pu[:, 2], PS[:, 2]
        u1 = (zt1 - ZS) / Zu
        u2 = (zt2 - Puu[:, 2] * u1**2 - 2 * PuS[:, 2] * u1 - PSS[:, 2]) / Zu
        Y = P[:, 1]
        Y1 = Pu[:, 1] * u1 + PS[:, 1]
        Y2 = Puu[:, 1] * u1**2 + 2 * PuS[:, 1] * u1 + PSS[:, 1] + Pu[:, 1] * u2
        out[side] = (Y, Y1, Y2)
        out["u%+d" % side] = u
        out["z%+d" % side] = P[:, 2]
    return out


# ----------------------------------------------------------------------------
# Attachment data as smooth functions of xi
# ----------------------------------------------------------------------------

XI_NOSE = None  # set by attachment_tables()


@lru_cache(maxsize=1)
def attachment_tables(n=2401):
    """Cubic splines in xi of camber Cm, squared half-thickness Q and their S-derivatives at S = 540.

    Cm = (Yu + Yl)/2, Q = ((Yu - Yl)/2)^2. Q is smooth through the round nose
    (Q ~ 2 R_LE c (xi - xi_nose)) and the sharp trailing edge, unlike the
    half-thickness itself.
    """
    from scipy.interpolate import CubicSpline
    global XI_NOSE
    par = parent()
    nose = float(par.nose_u(S_ATT))
    skin = par.surfaces()[0]
    v = float(par.v_at_span(S_ATT))
    xi_n = float(xi_of(S_ATT, skin.Value(nose, v).Z()))
    XI_NOSE = xi_n
    t = np.linspace(0.0, 1.0, n)
    xi = xi_n + (1.0 - xi_n) * (1 - np.cos(np.pi * t)) / 2
    xi[0] = xi_n + 1e-12
    xi[-1] = 1.0
    a = attachment_jet(xi)
    (Yu, Yu1, Yu2), (Yl, Yl1, Yl2) = a[1], a[-1]
    Cm, Cm1, Cm2 = (Yu + Yl) / 2, (Yu1 + Yl1) / 2, (Yu2 + Yl2) / 2
    h, h1, h2 = (Yu - Yl) / 2, (Yu1 - Yl1) / 2, (Yu2 - Yl2) / 2
    Q, Q1, Q2 = h * h, 2 * h * h1, 2 * (h1 * h1 + h * h2)
    # exact zeros at the ends
    Q[0] = 0.0; Q[-1] = 0.0
    tab = {k: CubicSpline(xi, val) for k, val in
           dict(Cm=Cm, Cm1=Cm1, Cm2=Cm2, Q=Q, Q1=Q1, Q2=Q2).items()}
    tab["xi_nose"] = xi_n
    tab["xi_samples"] = xi
    tab["raw"] = dict(Yu=Yu, Yl=Yl, Cm=Cm, Q=Q, Q1=Q1, Q2=Q2, Cm1=Cm1, Cm2=Cm2)
    return tab


# ----------------------------------------------------------------------------
# Base surface beyond the attachment
# ----------------------------------------------------------------------------

TAPER_RETURN_MM = 6.0  # D: span scale over which the parent's thickness taper is brought to a stop
LIP_BEAD = dict(h_gain=0.12, t_full=0.8, t_zero=1.3, extra_gain=0.08, extra_t=(0.1, 0.5), width=3.0)
# Thin trailing-edge band only: within a short distance of the lip (width x the closure length) the
# guide's half-thickness swells by up to 12 % (fading out between 0.8 and 1.3 mm attachment
# thickness; a further 8 % below 0.5 mm, complete at 0.1 mm) so that a curvature-continuous lip
# can still be as round as half the attachment thickness (see README).


def lip_bead_gain(t):
    """Gain on Q (squared half-thickness) of the lip bead for attachment thickness t at the lip's chord station."""
    t = np.asarray(t, float)
    tf, tz = LIP_BEAD["t_full"], LIP_BEAD["t_zero"]
    s = np.clip((tz - t) / (tz - tf), 0.0, 1.0)
    e0, e1 = LIP_BEAD["extra_t"]
    se_ = np.clip((e1 - t) / (e1 - e0), 0.0, 1.0)
    hg = LIP_BEAD["h_gain"] * s**3 * (10 - 15 * s + 6 * s * s) + LIP_BEAD["extra_gain"] * se_**3 * (10 - 15 * se_ + 6 * se_ * se_)
    return (1 + hg) ** 2 - 1


def base_fields(S, xi):
    """Camber Cm, squared half-thickness Qb and their (S, xi) partial derivatives.

    Cm is the exact second-order continuation of the parent (rise and camber).
    Qb = Q0 + Q1 Dr tanh(d/Dr) + Q2 (Dr tanh(d/Dr))^2 / 2, d = S - 540:
    identical to the parent's second-order continuation up to O(d^3) (so
    curvature-continuous at the attachment); the planform taper of the
    thickness then levels off monotonically (net -1.4 % half-thickness)
    instead of continuing to the tip.
    """
    tab = attachment_tables()
    S = np.asarray(S, float)
    xi = np.clip(np.asarray(xi, float), tab["xi_nose"], 1.0)
    d = S - S_ATT
    D = TAPER_RETURN_MM
    th = np.tanh(d / D)
    g1, g1d = D * th, 1 - th * th                       # d - d^3/(3D^2) + ...
    g2, g2d = (D * th) ** 2, 2 * D * th * (1 - th * th)  # d^2 - (2/3) d^4/D^2 + ...
    Q0, Q1, Q2 = tab["Q"](xi), tab["Q1"](xi), tab["Q2"](xi)
    Q0x, Q1x, Q2x = tab["Q"](xi, 1), tab["Q1"](xi, 1), tab["Q2"](xi, 1)
    Qb = Q0 + Q1 * g1 + 0.5 * Q2 * g2
    Qb_S = Q1 * g1d + 0.5 * Q2 * g2d
    Qb_x = Q0x + Q1x * g1 + 0.5 * Q2x * g2
    C0, C1, C2 = tab["Cm"](xi), tab["Cm1"](xi), tab["Cm2"](xi)
    Cm = C0 + C1 * d + 0.5 * C2 * d * d
    Cm_S = C1 + C2 * d
    Cm_x = tab["Cm"](xi, 1) + tab["Cm1"](xi, 1) * d + 0.5 * tab["Cm2"](xi, 1) * d * d
    return dict(Cm=Cm, Cm_S=Cm_S, Cm_x=Cm_x, Qb=Qb, Qb_S=Qb_S, Qb_x=Qb_x)


def plan_gradient(S, xi, dS, dx):
    """Convert (d/dS, d/dxi) partials to the plan gradient (d/dX, d/dZ)."""
    zle, zle1, _, c, c1 = frame(S)
    # X = -S, Z = zle - xi c  ->  dX = -dS, dZ = (zle1 - xi c1) dS - c dxi
    dSdX = -1.0
    dxidZ = -1.0 / c
    dxidX = (zle1 - xi * c1) / c  # from dZ = 0 when moving in X: dxi = (zle1 - xi c1) dS / c, dS = -dX
    dxidX = -dxidX
    gX = dS * dSdX + dx * dxidX
    gZ = dx * dxidZ
    return np.stack([gX, gZ], -1)


# ----------------------------------------------------------------------------
# Outline and guides
# ----------------------------------------------------------------------------

def outline(p):
    """Locked outline: plan (X, Z), (S, xi), unit tangent (d/dp direction) in plan."""
    par = parent()
    p = np.atleast_1d(np.asarray(p, float))
    O = par.tip_outline(p)
    h = 1e-6
    dO = (par.tip_outline(np.clip(p + h, 0, 1)) - par.tip_outline(np.clip(p - h, 0, 1)))
    T = dO / np.linalg.norm(dO, axis=-1, keepdims=True)
    S = -O[:, 0]
    xi = xi_of(S, O[:, 1])
    return O, S, xi, T


DESIGN = dict(
    R_cap=0.6,       # lip radius (mm) wherever the section is thick enough
    m_min=2.1,       # smallest closure exponent (>2 keeps the closure curvature-continuous)
    alpha=0.95,      # fraction of the guide's reach over which the roll is spread
    softmin_power=12.0,  # smooth minimum between R_cap and the achievable radius (p-norm)
    m_slope=3.0,     # closure exponent growth with section thickness (keeps the lip circle-like, k <= 1)
)


def _softmin(a, b, eps):
    return 0.5 * (a + b - np.sqrt((a - b) ** 2 + eps * eps))


def _softmax(a, b, eps):
    return 0.5 * (a + b + np.sqrt((a - b) ** 2 + eps * eps))


XI_ATTACH_P0 = 0.02


def xi_attach(p):
    """Attachment chord fraction of guide p: about p, but quadratic at the leading-edge end
    (xi - xi_nose ~ p^2 / p0) so the guides next to the nose leave the cut spanwise instead of
    running around the nose along the cut edge."""
    tab = attachment_tables()
    p = np.asarray(p, float)
    F = p * p * (1 + XI_ATTACH_P0) / (p + XI_ATTACH_P0)
    return tab["xi_nose"] + (1.0 - tab["xi_nose"]) * F


def _guide_geometry(p):
    """Geometric inputs of guide p: outline point, inboard normal, path direction and obliquity, plan length."""
    p = np.atleast_1d(np.asarray(p, float))
    O, S_O, xi_O, T = outline(p)
    xa = xi_attach(p)
    eps = 1e-6
    Sm, xm = S_O + (S_ATT - S_O) * eps, xi_O + (xa - xi_O) * eps
    Xm, Zm = plan_xz(Sm, xm)
    e = np.stack([Xm - O[:, 0], Zm - O[:, 1]], -1)
    e /= np.linalg.norm(e, axis=-1, keepdims=True)
    N = np.stack([-T[:, 1], T[:, 0]], -1)
    N *= np.sign(np.sum(N * e, -1))[:, None]
    sinb = np.sum(N * e, -1)
    tt = np.linspace(0, 1, 401)[:, None]
    Sp = S_ATT + tt * (S_O - S_ATT)
    xp = xa + tt * (xi_O - xa)
    Xp, Zp = plan_xz(Sp, xp)
    L = np.sum(np.hypot(np.diff(Xp, axis=0), np.diff(Zp, axis=0)), 0)
    return dict(p=p, O=O, S_O=S_O, xi_O=xi_O, xi_a=xa, T=T, N=N, e=e, sinb=sinb, L=L)


def design_guides(p, design=DESIGN):
    """Guide design laws for outline parameters p (0 < p < 1).

    Each guide runs straight in (S, xi) from the attachment station
    (540, xi_a(p)) to the outline point O(p); v is the plan distance from
    the lip. Squared half-thickness along it:
        Q(v) = Qb(path(v)) - Q_O * (1 - v/a)_+^m
    which vanishes at the lip with plan-normal gradient G, i.e. a lip radius
    G/2 before the camber slope shears the section. The laws are smooth in p
    (smooth minima, Gaussian-averaged obliquity) so neighbouring guides change
    shape gradually.
    """
    geo = _guide_geometry(p)
    p = geo["p"]
    # obliquity averaged over neighbouring guides (the outline is only C2 where its fairing ends)
    w = np.exp(-0.5 * np.arange(-4, 5) ** 2 / 2.0**2)
    w /= w.sum()
    sb = np.zeros_like(p)
    for k, wk in zip(range(-4, 5), w):
        sb += wk * _guide_geometry(np.clip(p + 0.006 * k, 1e-6, 1 - 1e-6))["sinb"]
    O, S_O, xi_O, N, L = geo["O"], geo["S_O"], geo["xi_O"], geo["N"], geo["L"]
    bf = base_fields(S_O, xi_O)
    t_O = thickness_att(xi_O)
    bead = np.where(xi_O > 0.5, lip_bead_gain(t_O), 0.0)
    Q_O = np.maximum(bf["Qb"], 1e-30) * (1 + bead)       # lip-region squared half-thickness incl. the bead
    gb = plan_gradient(S_O, xi_O, bf["Qb_S"], bf["Qb_x"])
    G_b = np.sum(gb * N, -1) * (1 + bead)
    gc = plan_gradient(S_O, xi_O, bf["Cm_S"], bf["Cm_x"])
    slope = np.sum(gc * N, -1)          # camber rise per mm moving inboard across the lip
    slope_d = np.clip(slope, -0.5, 0.5)  # (the micro-nose of the parent gives meaningless values)
    mmin, Rcap, alpha, kq = design["m_min"], design["R_cap"], design["alpha"], design["softmin_power"]
    Gam_max = np.sqrt(2 * Q_O * mmin / (mmin - 1))
    Gam_lo = mmin * Q_O / (alpha * L * sb)

    def softplus(x, wdt):
        return wdt * np.logaddexp(0.0, x / wdt)

    eta_min = 1 - 1 / mmin

    def laws(Gam):
        Gam = (Gam**4 + Gam_lo**4) ** 0.25
        m_spread = Gam * alpha * L * sb / Q_O    # closure spread over alpha * L
        eta = 2 * Q_O / Gam**2                   # k <= 1 (circle-like lip) needs m <= 1/(1 - eta) when eta < 1
        m_lin = mmin + design["m_slope"] * softplus(eta - eta_min, 0.02)
        x1 = softplus(m_spread - mmin, 0.3)
        x2 = m_lin - mmin
        m = mmin + (x1**-2 + np.maximum(x2, 1e-12) ** -2) ** -0.5
        return Gam, m, m * Q_O / Gam

    Gm, mm_, ap = laws(Gam_max)
    R_ach = lip_section_min_radius(Q_O, G_b, slope_d, Gm, mm_, ap)
    Gn, mn_, apn = laws(np.zeros_like(Q_O))
    R_nat = lip_section_min_radius(Q_O, G_b, slope_d, Gn, mn_, apn)
    R_t = (Rcap ** -kq + (0.995 * R_ach) ** -kq) ** (-1.0 / kq)
    # where the base surface alone already gives a rounder lip, blend smoothly into it
    R_t = np.minimum((R_t**6 + (0.98 * R_nat) ** 6) ** (1 / 6), 0.995 * R_ach)
    lo = np.full_like(Q_O, 0.0)
    hi = Gam_max.copy()
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        Gx, mx, ax = laws(mid)
        R = lip_section_min_radius(Q_O, G_b, slope_d, Gx, mx, ax)
        big = R > R_t
        lo, hi = np.where(big, lo, mid), np.where(big, mid, hi)
    Gam, m, a_p = laws(0.5 * (lo + hi))
    a = a_p / sb
    b_bead = LIP_BEAD["width"] * a                     # along the guide
    return dict(geo, sinb_smoothed=sb, Q_O=Q_O, bead=bead, b_bead=b_bead, Q_O_base=Q_O / (1 + bead), G_b=G_b, slope=slope, R_t=R_t, R_ach=R_ach, Gam=Gam, m=m, a=a,
                a_p=a_p, R_lip_design=0.5 * (G_b + m * Q_O / a_p),
                R_phys_design=lip_section_min_radius(Q_O, G_b, slope_d, Gam, m, a_p),
                t_att_O=t_O)


def lip_section_min_radius(Q_O, G_b, slope, Gam, m, a_p, n=1500):
    """Smallest radius of the model lip section normal to the outline (both sides).

    Model: Q(d) = Q_O + G_b d - Q_O (1 - d/a_p)^m (d = plan distance inboard
    of the lip), Y = slope d +- sqrt(Q). The camber slope shears the section,
    which tightens one side of the lip (the vertical-tangent radius G/2 is
    not the smallest radius once the section is inclined).
    """
    Q_O, G_b, slope, Gam, m, a_p = [np.asarray(x, float)[..., None] for x in (Q_O, G_b, slope, Gam, m, a_p)]
    x = np.linspace(0.0, 1.0, n) ** 2
    dd = a_p * x
    K = (1 - x) ** m
    Q = Q_O + G_b * dd - Q_O * K
    Q1 = G_b + Q_O * m / a_p * (1 - x) ** (m - 1)
    Q2 = -Q_O * m * (m - 1) / a_p**2 * (1 - x) ** np.maximum(m - 2, 0)
    w = np.sqrt(np.maximum(Q, 0))
    best = np.full(Q_O.shape[:-1], np.inf)
    for s in (+1, -1):
        ww = s * w
        Fy = 2 * ww
        Fd = -2 * ww * slope - Q1
        Fdd = 2 * slope**2 - Q2
        Fdy = -2 * slope
        k = np.abs(Fdd * Fy**2 - 2 * Fdy * Fd * Fy + 2 * Fd**2) / (Fd**2 + Fy**2) ** 1.5
        best = np.minimum(best, 1 / np.max(k, -1))
    return best


def guide_dense(d, i, K=3000, q=4.0):
    """Dense samples of guide i of a design table d, ordered from the attachment (v = L) to the lip (v = 0).

    Returns plan X, Z, camber Cm, half-thickness H, the plan distance v from
    the lip, and (S, xi) of every sample.
    """
    L = d["L"][i]
    r = np.linspace(0.0, 1.0, K + 1)
    v = L * r**q                      # clustered at the lip
    v = v[::-1]                       # attachment -> lip
    # plan distance -> path parameter t (straight in (S, xi)); invert the plan arclength
    tt = np.linspace(0, 1, 4001)
    S_O, xi_O, xa = d["S_O"][i], d["xi_O"][i], d["xi_a"][i]
    Sp = S_ATT + tt * (S_O - S_ATT)
    xp = xa + tt * (xi_O - xa)
    Xp, Zp = plan_xz(Sp, xp)
    seg = np.hypot(np.diff(Xp), np.diff(Zp))
    arc = np.r_[0, np.cumsum(seg)]
    dist_from_lip = arc[-1] - arc
    t = np.interp(v, dist_from_lip[::-1], tt[::-1])
    S = S_ATT + t * (S_O - S_ATT)
    xi = xa + t * (xi_O - xa)
    X, Z = plan_xz(S, xi)
    bf = base_fields(S, xi)
    vv = v * (L / arc[-1])
    a, m, Q_O = d["a"][i], d["m"][i], d["Q_O"][i]
    K_ = np.where(vv < a, np.clip(1 - vv / a, 0, 1) ** m, 0.0)
    bead, bb = d["bead"][i], d["b_bead"][i]
    x = np.clip(vv / max(bb, 1e-300), 0, 1)
    Sb = 1 - x**3 * (10 - 15 * x + 6 * x * x)          # 1 at the lip, 0 beyond the bead width
    Q = bf["Qb"] * (1 + bead * Sb) - Q_O * K_
    H = np.sqrt(np.maximum(Q, 0.0))
    H[-1] = 0.0
    return dict(X=X, Z=Z, Cm=bf["Cm"], H=H, Q=Q, v=vv, S=S, xi=xi)


def guide_end_jet(d, i, side, h=2e-4):
    """Unit tangent (towards the lip) and curvature vector of guide i at its attachment end.

    Near the attachment the guide lies on the base surface (closure term is
    exactly zero for v > a), so this is evaluated from the base fields along
    the straight (S, xi) path with a one-sided polynomial fit.
    """
    S_O, xi_O, xa = d["S_O"][i], d["xi_O"][i], d["xi_a"][i]
    t = h * np.arange(9)
    S = S_ATT + t * (S_O - S_ATT)
    xi = xa + t * (xi_O - xa)
    X, Z = plan_xz(S, xi)
    bf = base_fields(S, xi)
    Y = bf["Cm"] + side * np.sqrt(np.maximum(bf["Qb"], 0))
    P = np.stack([X, Y, Z], -1)
    V = np.vander(t / h, 7, increasing=True)
    coef, *_ = np.linalg.lstsq(V[:, :7], P, rcond=None)
    d1 = coef[1] / h
    d2 = 2 * coef[2] / h**2
    sp = np.linalg.norm(d1)
    T = d1 / sp
    kN = (d2 - np.dot(d2, T) * T) / sp**2
    return T, kN


def loop_parameter(loop, i_lip, eps, lin=0.5):
    """Wrap parameter in [0, 2] along a closed-over guide loop (upper attachment -> lip -> lower attachment).

    sigma = 1 -+ phi(s)/phi(s_end), phi(s) = asinh(s/eps) + lin*s/s_end, with s
    the arc length from the lip: uniform spacing through the lip at the scale
    eps, geometric growth away from it, so neighbouring interpolation nodes
    never jump in spacing (which made the bicubic net ring).
    """
    seg = np.linalg.norm(np.diff(loop, axis=0), axis=1)
    arc = np.r_[0, np.cumsum(seg)]
    s = arc - arc[i_lip]
    s_up, s_lo = -s[0], s[-1]
    sig = np.empty_like(s)
    up = s <= 0
    phi = lambda x, send: np.arcsinh(x / eps) + lin * x / send
    sig[up] = 1.0 - phi(-s[up], s_up) / phi(s_up, s_up)
    sig[~up] = 1.0 + phi(s[~up], s_lo) / phi(s_lo, s_lo)
    return sig


def thickness_att(xi):
    """Parent attachment thickness (mm, along Y) at frame chord fraction xi (from the Q table)."""
    tab = attachment_tables()
    xi = np.clip(np.asarray(xi, float), tab["xi_nose"], 1.0)
    return 2.0 * np.sqrt(np.maximum(tab["Q"](xi), 0.0))
