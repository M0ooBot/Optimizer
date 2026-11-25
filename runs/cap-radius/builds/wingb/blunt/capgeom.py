"""Geometry definition of the Wing B (blunt) rounded tip cap, used by build_cap.py.

Axes as the parent STEP: X = -span (outboard is -X), Y = up, Z = forward.
d (Delta) = projected span beyond the attachment (S - 540 mm).

Design in one paragraph: the parent skin is continued past the attachment by its
exact second-order Taylor expansion in span (E, G2 with the parent). Over the plan
region inside the locked outline the cap is the graph Y = m +/- h of E's mean line m
and a half-thickness h that closes to zero on the outline with a G2 round lip
(lip_profile) whose apex radius and width are scheduled along the outline
(schedule). The thin trailing-edge region keeps the attachment thickness at the
same chord position (+TH_GAIN) so its full half-round lip meets the printable rule;
the flat 0.4 mm trailing-edge face fades into that half-round over L_TE of span.
"""
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline, CubicSpline
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
PARENT_DIR = ROOT / "runs" / "cap-radius" / "parents" / "wingb" / "blunt"
sys.path.insert(0, str(PARENT_DIR))
import parent as P  # noqa: E402

ATT = P.ATTACHMENT          # 540
NOZZLE = 0.4                 # printable minimum radius (mm)
THIN = 0.8                   # thickness threshold of the printable rule (mm)


# ----------------------------------------------------------------------------
# E: the parent skin continued past the attachment by its exact second-order
# Taylor expansion in projected span at fixed native u (G2 with the parent).
# ----------------------------------------------------------------------------
@lru_cache(maxsize=1)
def _taylor():
    sk, _ = P.surfaces()
    nu, nv = sk.NbUPoles(), sk.NbVPoles()
    pu, pv = sk.UDegree(), sk.VDegree()

    def knots(n, knot, mult):
        out = []
        for i in range(1, n + 1):
            out += [knot(i)] * mult(i)
        return np.array(out)

    tu = knots(sk.NbUKnots(), sk.UKnot, sk.UMultiplicity)
    tv = knots(sk.NbVKnots(), sk.VKnot, sk.VMultiplicity)
    poles = np.array([[sk.Pole(i, j).Coord() for j in range(1, nv + 1)] for i in range(1, nu + 1)])
    w = np.array([[sk.Weight(i, j) for j in range(1, nv + 1)] for i in range(1, nu + 1)])
    assert np.allclose(w, w[:, :1]), "weights must be constant along v"
    w = w[:, 0]
    va = tv[-1]
    assert abs(va - P.SURFACE_DATA["v_attachment"]) < 1e-9
    q = BSpline(tv, np.transpose(poles, (1, 0, 2)).reshape(nv, -1), pv, extrapolate=True)
    q0 = q(va).reshape(nu, 3)
    q1 = q.derivative(1)(va).reshape(nu, 3)
    q2 = q.derivative(2)(va).reshape(nu, 3)
    xv, xvv = q1[:, 0], q2[:, 0]
    assert np.ptp(xv) < 1e-9 * abs(xv).max() + 1e-12 and np.ptp(q0[:, 0]) < 1e-9
    v1 = -1.0 / xv[0]
    v2 = -xvv[0] * v1 ** 2 / xv[0]
    c0, c1, c2 = q0, q1 * v1, q2 * v1 ** 2 + q1 * v2
    num = [BSpline(tu, w[:, None] * c, pu) for c in (c0, c1, c2)]
    den = BSpline(tu, w, pu)
    return {"num": num, "den": den, "num_d": [n.derivative() for n in num], "den_d": den.derivative(),
            "tu": tu}


def E(u, d, der=0):
    """Continued skin point at native u and Delta d (arrays broadcast). der=1 also returns d/du."""
    T = _taylor()
    u, d = np.broadcast_arrays(np.asarray(u, float), np.asarray(d, float))
    dd = d[..., None]
    n = T["num"][0](u) + dd * T["num"][1](u) + 0.5 * dd ** 2 * T["num"][2](u)
    W = T["den"](u)[..., None]
    p = n / W
    if not der:
        return p
    nd = T["num_d"][0](u) + dd * T["num_d"][1](u) + 0.5 * dd ** 2 * T["num_d"][2](u)
    return p, (nd - p * T["den_d"](u)[..., None]) / W


def _bisect(f, lo, hi, it=60):
    flo = f(lo)
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        same = np.sign(fm) == np.sign(flo)
        lo, flo = np.where(same, mid, lo), np.where(same, fm, flo)
        hi = np.where(same, hi, mid)
    return 0.5 * (lo + hi)


def _nose_u_exact(d):
    d = np.asarray(d, float)
    lo, hi = np.full(d.shape, 0.4155), np.full(d.shape, 0.5755)
    return _bisect(lambda u: E(u, d, 1)[1][..., 2], lo, hi)


@lru_cache(maxsize=1)
def _nose_spline():
    dg = np.linspace(-2.0, 80.0, 821)
    return CubicSpline(dg, _nose_u_exact(dg))


def nose_u(d):
    return _nose_spline()(np.asarray(d, float))


def section_ends(d):
    """(Z_nose, Y_nose, Z_te, Y_te_mid) of the continued section at Delta d."""
    d = np.asarray(d, float)
    un = nose_u(d)
    pn = E(un, d)
    p0, p1 = E(np.zeros(d.shape), d), E(np.ones(d.shape), d)
    return pn[..., 2], pn[..., 1], 0.5 * (p0[..., 2] + p1[..., 2]), 0.5 * (p0[..., 1] + p1[..., 1])


def side_y(d, z, side):
    """Y of the continued skin at Delta d and plan Z on one side (+1 upper, -1 lower)."""
    d, z = np.broadcast_arrays(np.asarray(d, float), np.asarray(z, float))
    un = nose_u(d)
    end = np.ones(d.shape) if side > 0 else np.zeros(d.shape)
    zn, zt = E(un, d)[..., 2], E(end, d)[..., 2]
    z = np.clip(z, zt, zn)
    # plain bisection: Newton is ill-conditioned next to the nose (dZ/du -> 0)
    u = _bisect(lambda u: E(u, d)[..., 2] - z, un, end, 46)
    return E(u, d)[..., 1], u


def smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x ** 3 * (10 - 15 * x + 6 * x ** 2)


def step7(x):
    """C3 step 0 -> 1 on [0, 1]."""
    x = np.clip(x, 0.0, 1.0)
    return x ** 4 * (35 - 84 * x + 70 * x ** 2 - 20 * x ** 3)


def flat_rise(x):
    """C-infinity rise from 0 (with zero 1st and 2nd derivatives at x = 0) towards 1."""
    x = np.maximum(np.asarray(x, float), 0.0)
    return 1.0 - np.exp(-x ** 3)


def soft_step(x):
    """C-infinity step 0 -> 1 around x = 0.5 (logistic)."""
    return 0.5 * (1.0 + np.tanh(6.0 * (np.asarray(x, float) - 0.5)))


@lru_cache(maxsize=1)
def attachment_ends():
    zn, yn, zt, yt = section_ends(np.array(0.0))
    return float(zn), float(zt)


def xi_hat(d, z):
    """Chord position nose (0) to trailing edge (1) of the continued section at Delta d."""
    zn, _, zt, _ = section_ends(d)
    return (zn - z) / (zn - zt)


def _att_half_thickness_exact(xi):
    zn, zt = attachment_ends()
    z = zn - np.asarray(xi, float) * (zn - zt)
    yu, _ = side_y(np.zeros_like(z), z, +1)
    yl, _ = side_y(np.zeros_like(z), z, -1)
    return 0.5 * (yu - yl)


@lru_cache(maxsize=1)
def _att_spline():
    eta = np.linspace(0.0, 1.0, 8001)
    return CubicSpline(eta, _att_half_thickness_exact(eta ** 2))


def att_half_thickness(xi):
    return _att_spline()(np.sqrt(np.clip(np.asarray(xi, float), 0.0, 1.0)))


# ----------------------------------------------------------------------------
# Base field: mean line and thickness of E; in the thin trailing-edge region
# (chord position > ~0.9) the thickness is blended (C2-flat at the attachment)
# into the attachment section's own thickness at the same chord position plus
# TH_GAIN, so the thin lip can meet the printable rule judged against it.
# ----------------------------------------------------------------------------
L_BLEND = 6.0
TH_GAIN = 0.08     # the thin trailing-edge region (chord position > ~0.9) is 8 % thicker away from the attachment


def base_fields(x, z):
    """(mean Y, half-thickness, chord position xi) of the uncapped base at plan points."""
    x, z = np.broadcast_arrays(np.asarray(x, float), np.asarray(z, float))
    d = -x - ATT
    yu, _ = side_y(d, z, +1)
    yl, _ = side_y(d, z, -1)
    xi = xi_hat(d, z)
    he = 0.5 * (yu - yl)
    # thin trailing-edge region only: replace the parent's natural taper by the attachment
    # thickness at the same chord position, thickened by TH_GAIN (elsewhere E's own taper)
    w = soft_step((xi - 0.85) / 0.1)
    ha = att_half_thickness(np.clip(xi, 0.0, 1.0)) * (1.0 + TH_GAIN)
    b = flat_rise(d / L_BLEND)
    return 0.5 * (yu + yl), he + (ha - he) * b * w, xi


def rule_radius(xi):
    """Printable rule: minimum radius at chord position xi (attachment thickness)."""
    h = att_half_thickness(np.clip(xi, 0.0, 1.0))
    return np.where(2 * h >= THIN, NOZZLE, h)


# ----------------------------------------------------------------------------
# Locked outline (plan X, Z) as a dense C2 spline in p, plus projection.
# ----------------------------------------------------------------------------
@lru_cache(maxsize=1)
def outline():
    p = np.linspace(0.0, 1.0, 40001)
    xz = P.tip_outline(p)
    cs = CubicSpline(p, xz, bc_type="not-a-knot")
    seg = np.linalg.norm(np.diff(xz, axis=0), axis=1)
    arc = np.concatenate([[0], np.cumsum(seg)])
    return {"p": p, "xz": xz, "cs": cs, "arc": arc / arc[-1], "length": arc[-1],
            "tree": cKDTree(xz), "p_of_arc": CubicSpline(arc / arc[-1], p)}


def outline_frame(p):
    """Point, unit tangent (increasing p) and unit inward normal in plan (X, Z)."""
    o = outline()
    pt, d1, d2 = o["cs"](p), o["cs"](p, 1), o["cs"](p, 2)
    sp = np.linalg.norm(d1, axis=-1, keepdims=True)
    t = d1 / sp
    # the outline runs LE -> tip -> TE; with X outboard negative the interior is on the right
    n = np.stack([-t[..., 1], t[..., 0]], -1)
    k = (d1[..., 0] * d2[..., 1] - d1[..., 1] * d2[..., 0]) / sp[..., 0] ** 3
    return pt, t, n, k


def p_of_arc(a):
    return np.clip(outline()["p_of_arc"](np.clip(a, 0.0, 1.0)), 0.0, 1.0)


def project(x, z):
    """Foot parameter p and signed inward distance of plan points to the outline."""
    o = outline()
    q = np.stack(np.broadcast_arrays(np.asarray(x, float), np.asarray(z, float)), -1)
    _, idx = o["tree"].query(q.reshape(-1, 2))
    p = o["p"][idx].reshape(q.shape[:-1])
    for _ in range(8):
        c, c1, c2 = o["cs"](p), o["cs"](p, 1), o["cs"](p, 2)
        r = c - q
        g = np.sum(r * c1, -1)
        h = np.sum(c1 * c1, -1) + np.sum(r * c2, -1)
        p = np.clip(p - g / h, 0.0, 1.0)
    pt, t, n, _ = outline_frame(p)
    return p, np.sum((q - pt) * n, -1)


# ----------------------------------------------------------------------------
# Lip schedule along the outline.
#   closed profile  h_c^2 = H^2 - H_out(p)^2 * F(delta/ell(p))   (F: lip_profile)
#   (delta: plan distance inside the outline, p: foot point on the outline)
#   apex radius     r = g + k0 H_out^2 / (2 ell),   g = (1/2) dH^2/ddelta at the outline.
# Near the trailing-edge corner the flat TE face fades into the half-round over
# L_TE of span: h = H - B (H - h_c), B = 0 at the attachment, 1 beyond L_TE.
# ----------------------------------------------------------------------------
L_TE = 2.0          # span over which the flat trailing-edge face fades into the half-round
SPREAD = 0.6        # fraction of the medial distance the roll may use
SOFT_K = 10.0
MEAN_ZONE = 4.0      # the mean surface is levelled within MEAN_ZONE * ell of the lip
NOSE_ZONE = 0.05     # mm; leading-edge apex regularisation zone
L_NOSE = 1.0         # mm of span over which that regularisation fades in


THETA_RAMP = 0.5     # turning angle (rad) over which the lip curvature ramps up from the flat


@lru_cache(maxsize=1)
def lip_profile():
    """Normalised lip profile: a G2 round whose curvature ramps from 0 (at the flat,
    as theta^(3/4)) to a constant over the first THETA_RAMP of turning and then stays
    constant to the apex. Returns F(X) = 1 - Phi(X)^2 (X: distance from the apex over
    the profile run, Phi: height over the rise), the apex slope k0 = -F'(0), the
    aspect a0 = run / rise and the constant-curvature radius over the rise."""
    v = np.linspace(0.0, (np.pi / 2) ** 0.25, 400001)
    th = v ** 4
    u = np.minimum(th / THETA_RAMP, 1.0 - 1e-16)
    ramp = np.where(th < THETA_RAMP, (-np.expm1(4 * np.log1p(-u))) ** 0.75, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        fx = np.where(v > 0, np.cos(th) * 4 * v ** 3 / ramp, 0.0)
        fy = np.where(v > 0, np.sin(th) * 4 * v ** 3 / ramp, 0.0)
    fx[0], fy[0] = fx[1], 0.0
    from scipy.integrate import cumulative_trapezoid
    x = cumulative_trapezoid(fx, v, initial=0.0)      # from the flat junction
    y = cumulative_trapezoid(fy, v, initial=0.0)
    run, rise = x[-1], y[-1]
    X = (run - x) / run                                # 0 at the apex, 1 at the flat
    Phi = (rise - y) / rise
    F = 1.0 - Phi ** 2
    order = np.argsort(X)
    Xs, Fs = X[order], F[order]
    keep = np.concatenate([[True], np.diff(Xs) > 1e-9])
    spl = CubicSpline(Xs[keep], Fs[keep])
    k0 = 2.0 * run / rise ** 2
    return {"F": spl, "k0": k0, "a0": run / rise, "r_over_rise": 1.0 / rise}


def lip_F(x):
    x = np.clip(np.asarray(x, float), 0.0, 1.0)
    return np.clip(lip_profile()["F"](x), 0.0, 1.0) * (x < 1.0)


def _soft_max(a, b, k=SOFT_K):
    return (a ** k + b ** k) ** (1.0 / k)


def _soft_min(a, b, k=SOFT_K):
    return (a ** -k + b ** -k) ** (-1.0 / k)


def medial_distance(p, smax=90.0, ns=361):
    """Distance along the inward normal from O(p) until the nearest outline point
    stops being O(p) or the attachment plane is reached."""
    pt, t, n, _ = outline_frame(p)
    s = np.linspace(0.0, smax, ns)[1:]
    q = pt[:, None, :] + s[None, :, None] * n[:, None, :]
    pf, _ = project(q[..., 0], q[..., 1])
    bad = (np.abs(pf - p[:, None]) > 2e-3) | (q[..., 0] > -ATT)
    first = np.where(bad.any(1), bad.argmax(1), ns - 2)
    return s[first]


@lru_cache(maxsize=1)
def schedule():
    from scipy.ndimage import gaussian_filter1d, minimum_filter1d
    p = np.linspace(0.0, 1.0, 2001)
    pt, t, n, k = outline_frame(p)
    m, H, xi = base_fields(pt[:, 0], pt[:, 1])
    H = np.maximum(H, 0.0)
    eps = 0.01
    m2, H2, _ = base_fields(pt[:, 0] + eps * n[:, 0], pt[:, 1] + eps * n[:, 1])
    g = (H2 ** 2 - H ** 2) / (2 * eps)
    mu = (m2 - m) / eps                               # slope of the mean surface across the lip
    d = -pt[:, 0] - ATT
    zn, _, _, _ = section_ends(d)
    _, hh, _ = base_fields(pt[:, 0], zn - 0.002)
    rho = hh ** 2 / (2 * 0.002)                      # base nose radius at this span
    pc = np.linspace(0.0, 1.0, 401)
    M = np.interp(p, pc, medial_distance(pc))
    Ms = gaussian_filter1d(minimum_filter1d(M, 81, mode="nearest"), 25, mode="nearest")
    ell_max = SPREAD * Ms
    # leading-edge side and tip: spread the roll as far as the medial distance allows
    ell_le = ell_max
    # trailing-edge side: thin, flat-thickness section; aim the apex radius at the base
    # nose radius but never above the half-thickness (a full half-round where thin)
    k0 = lip_profile()["k0"]
    r_round = lip_profile()["r_over_rise"] * H          # roundest G2 lip the thickness allows
    r_spread = k0 * H ** 2 / (2 * ell_max)
    r = _soft_min(_soft_max(rho, r_spread), np.maximum(r_round, 1e-9))
    ell_te = _soft_min(k0 * H ** 2 / (2 * r) + 1e-12, ell_max)
    w_te = soft_step((p - 0.55) / 0.2)
    ell = np.exp((1 - w_te) * np.log(ell_le) + w_te * np.log(ell_te))
    r_eff = g + k0 * H ** 2 / (2 * ell)
    return {"p": p, "pt": pt, "H": H, "m": m, "xi": xi, "g": g, "rho": rho, "M": M,
            "ell_max": ell_max, "r_target": r, "r": r_eff, "ell": ell, "rule": rule_radius(xi), "d": d,
            "mu": mu, "ell_cs": CubicSpline(p, ell), "H_cs": CubicSpline(p, H), "mu_cs": CubicSpline(p, mu),
            "m_cs": CubicSpline(p, m), "r_cs": CubicSpline(p, r_eff)}


def lip_ell(p):
    return np.maximum(schedule()["ell_cs"](p), 1e-9)


def te_blend(d, p):
    """B = 1 away from the trailing-edge corner; 0 at the attachment near the TE."""
    b = smoothstep(d / L_TE)
    w = smoothstep((p - 0.6) / 0.15)
    return 1.0 - (1.0 - b) * w


def cap_fields(x, z):
    """Mean Y and half-thickness h of the finished cap at plan points inside the outline."""
    m, H, xi = base_fields(x, z)
    H = np.maximum(H, 0.0)
    p, delta = project(x, z)
    ell = lip_ell(p)
    Hout = np.maximum(schedule()["H_cs"](p), 0.0)
    xx = np.clip(np.maximum(delta, 0.0) / ell, 0.0, 1.0)
    hc = np.sqrt(np.maximum(H ** 2 - Hout ** 2 * lip_F(xx), 0.0))
    d = -np.asarray(x) - ATT
    B = te_blend(d, p)
    h = H - B * (H - hc)
    S = schedule()
    dl = np.maximum(delta, 0.0)
    mu = S["mu_cs"](p)
    # leading-edge side: within NOSE_ZONE of the lip replace the base's micro-structure
    # (the parent's nose B-spline has a ~0.07 mm radius spot within a few micrometres
    # of its nose) by the smooth apex parabola and a linear mean; faded in from the
    # attachment over L_NOSE so the parent is still matched there
    wn = (1.0 - step7(dl / NOSE_ZONE)) * flat_rise(d / L_NOSE) * (1.0 - soft_step((p - 0.35) / 0.2))
    h = np.sqrt(np.maximum(h ** 2 + (2.0 * S["r_cs"](p) * dl - h ** 2) * wn, 0.0))
    m = m + (S["m_cs"](p) + mu * dl - m) * wn
    # trailing-edge side: level the mean surface across the lip so the thin half-round
    # is symmetric (otherwise one side turns 90 deg + the camber slope within H); the
    # slope is removed exactly within the lip zone and restored beyond it
    lm = np.minimum(MEAN_ZONE * ell, np.interp(p, S["p"], S["ell_max"]))
    b = step7((dl / lm - 0.3) / 0.7)
    # only where the lip is thin (half-thickness below ~0.7 mm), i.e. on the trailing-edge side
    w = soft_step((1.0 - np.maximum(S["H_cs"](p), 0.0)) / 0.6) * (p > 0.5)
    m = m - B * w * mu * dl * (1.0 - b)
    return m, h, xi, p, delta


HERMITE = 1.0
LIP_SHARE = 4.0


# ----------------------------------------------------------------------------
# Patch parametrisation (shared by the upper and lower patches):
#   s in [0, 1] along the attachment cut, s = upper-skin native u mapped linearly
#   (nose 0, trailing edge 1); t in [0, 1] from the cut (0) to the lip (1).
#   The plan path from the cut point C(s) to the outline point O(s) is a cubic
#   Hermite curve leaving C spanwise and reaching O along the outline normal;
#   omega = 1 - (1 - t)^2 so that the half-thickness (~ sqrt of the distance to
#   the lip) is linear in t at the lip.
# ----------------------------------------------------------------------------


@lru_cache(maxsize=1)
def u_nose0():
    return float(_nose_u_exact(np.array(0.0)))


def cut(s, side=+1):
    """Cut station s on one skin (+1 upper: u = un + s (1 - un); -1 lower: u = un - s un).
    Returns plan Z, Y, chord position xi."""
    s = np.asarray(s, float)
    un = u_nose0()
    u = un + s * (1.0 - un) if side > 0 else un - s * un
    pt = E(u, np.zeros(s.shape))
    zn, zt = attachment_ends()
    return pt[..., 2], pt[..., 1], np.clip((zn - pt[..., 2]) / (zn - zt), 0.0, 1.0)


def station_p(s):
    """Outline parameter reached by the guides from station s: equal arc fractions
    (chord position of the upper cut point = arc fraction of the outline)."""
    return p_of_arc(cut(s, +1)[2])


def plan_map(s, t, side=+1):
    s, t = np.broadcast_arrays(np.asarray(s, float), np.asarray(t, float))
    zc = cut(s, side)[0]
    c = np.stack([np.full(zc.shape, -ATT), zc], -1)
    p = station_p(s)
    o, tan, n, _ = outline_frame(p)
    dist = np.linalg.norm(o - c, axis=-1)
    L = dist[..., None] * HERMITE
    # the lip zone (delta < ell) takes a fixed share of t
    # (trailing-edge side only, where the lip zone is narrow; plain sigma^2 elsewhere)
    lam = _soft_min(np.ones(dist.shape), np.maximum(LIP_SHARE * lip_ell(p) / np.maximum(dist, 1e-12), 1e-12))
    lam = 1.0 - soft_step((p - 0.55) / 0.2) * (1.0 - lam)
    sig = 1.0 - t
    w = 1.0 - sig ** 2 * (lam + (1.0 - lam) * sig)
    h00 = 2 * w ** 3 - 3 * w ** 2 + 1
    h10 = w ** 3 - 2 * w ** 2 + w
    h01 = -2 * w ** 3 + 3 * w ** 2
    h11 = w ** 3 - w ** 2
    tc = np.array([-1.0, 0.0])
    q = h00[..., None] * c + h10[..., None] * L * tc + h01[..., None] * o + h11[..., None] * L * (-n)
    return q, p
