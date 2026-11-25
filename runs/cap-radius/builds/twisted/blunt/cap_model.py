"""Analytic model of the rounded tip cap for the twisted/blunt parent (right tip, +Y).

The cap is written as   point = M(s, Y) + w * lam(s, Y) * D(s, Y)
where M = (upper + lower)/2 and D = (upper - lower)/2 are taken from the
parent's own skin definition (parent.py skin_jet, continued past the
attachment), s is the parent's chord parameter, Y the span station and w the
thickness fraction (+1 upper skin, -1 lower skin, 0 the mid surface).

Outside a band along the lip w = +-1, so the cap there IS the parent skin
continued analytically (curvature continuous with the parent by
construction).  Inside the band (width l(p) measured inward from the lip
apex, perpendicular to the rim) w follows a rounded profile that reaches 0 at
the apex.  See README.md for the design rules.

Charts: a = s^2 (2 - s) is used instead of s for in-plane interpolation (it is
proportional to chordwise distance at both the leading and trailing edge).
"""
import os
import sys
import numpy as np
from scipy.interpolate import CubicSpline

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, *[os.pardir] * 5))
PARENT_DIR = os.path.join(ROOT, "runs", "cap-radius", "parents", "twisted", "blunt")
if PARENT_DIR not in sys.path:
    sys.path.insert(0, PARENT_DIR)
import parent as pa  # noqa: E402

CUT = float(pa.CUT)
H_TE = float(pa.TE_THICKNESS) / 2.0
RULE_MIN_R = 0.4          # printable minimum radius (nozzle width)
RULE_T = 0.8              # thickness threshold of the rule

# ---- design constants (see README) -------------------------------------------
ETA_OMEGA = 3.0           # span over which the TE-region thickness compensation fades in (mm)
THIN_MARGIN = 0.035       # extra thickness (fraction) in the TE thin zone (s > 0.90), faded in over the TE corner
                          # zone, so the inscribed-circle lip clears t/2 with a small margin
S_MARGIN = (0.90, 0.94)   # chord-parameter ramp of THIN_MARGIN
S_COMP = (0.80, 0.90)     # chord-parameter ramp of the thickness compensation
ETA_F = 0.20              # span of the trailing-edge corner transition (crease -> half-round), mm
REACH_FRAC = 0.35         # band width as a fraction of the apex's span reach (spreads the roll)
MARGIN = 1.05             # design margin on the printable radius
P_TE = (0.78, 0.925)      # rim-parameter ramp (TE side): superelliptic m=3 profile -> inscribed circle (m=2)
NOSE_FAC = 8.0            # near the LE the band is NOSE_FAC times the rim's offset from the continued LE
P_NOSE = (0.04, 0.12)     # rim-parameter ramp from the shifted-nose profile to the superelliptic one


def smooth5(t):
    t = np.clip(t, 0.0, 1.0)
    return t ** 3 * (10 - 15 * t + 6 * t * t)


def smooth7(t):
    t = np.clip(t, 0.0, 1.0)
    return t ** 4 * (35 - 84 * t + 70 * t * t - 20 * t ** 3)


def a_of_s(s):
    s = np.asarray(s, float)
    return s * s * (2 - s)


def s_of_a(a):
    """Inverse of a(s); a slightly above 1 (lip apex aft of the continued TE, < 0.02 mm) maps to s > 1."""
    a = np.clip(np.asarray(a, float), 0.0, 1.15)
    s = np.where(a < 0.5, np.sqrt(a / 2.0), a)
    for _ in range(60):
        f = s * s * (2 - s) - a
        fp = np.maximum(4 * s - 3 * s * s, 1e-12)
        s = np.clip(s - f / fp, 0.0, 1.25)
    return s


def dads(s):
    s = np.asarray(s, float)
    return 4 * s - 3 * s * s


# ---- parent evaluation ----------------------------------------------------------
def skins(s, Y):
    up = pa.skin_jet(s, Y, +1)
    lo = pa.skin_jet(s, Y, -1)
    return up, lo


def att_half(s):
    """|D| at the attachment for chord parameter s."""
    s = np.asarray(s, float)
    return 0.5 * np.linalg.norm(pa.skin(s, CUT, +1) - pa.skin(s, CUT, -1), axis=-1)


def lam_of(s, Y, Dn):
    """Thickness compensation factor (TE region only): blends |D| to its attachment value."""
    s = np.asarray(s, float)
    Y = np.asarray(Y, float)
    ks = smooth5((s - S_COMP[0]) / (S_COMP[1] - S_COMP[0]))
    w = smooth5((Y - CUT) / ETA_OMEGA) * ks
    marg = 1 + THIN_MARGIN * smooth5((s - S_MARGIN[0]) / (S_MARGIN[1] - S_MARGIN[0])) * smooth5((Y - CUT) / ETA_F)
    out = np.ones(np.broadcast(s, Y).shape) * marg
    sel = w > 0
    if np.any(sel):
        ss = np.broadcast_to(s, out.shape)[sel]
        ratio = att_half(ss) / np.maximum(Dn[sel], 1e-15)
        out[sel] = (1 + w[sel] * (ratio - 1)) * np.broadcast_to(marg, out.shape)[sel]
    return out


def MD(s, Y):
    """Mid surface M, compensated half-thickness vector lam*D, and the raw parts."""
    up, lo = skins(s, Y)
    M = 0.5 * (up["point"] + lo["point"])
    D = 0.5 * (up["point"] - lo["point"])
    Dn = np.linalg.norm(D, axis=-1)
    lam = lam_of(s, Y, Dn)
    return M, D * lam[..., None], dict(up=up, lo=lo, lam=lam)


def cap_point(s, Y, w):
    M, LD, _ = MD(s, Y)
    return M + np.asarray(w, float)[..., None] * LD


def section_frame_v(Y):
    Y = np.abs(np.asarray(Y, float))
    th = np.deg2rad(3.0) * Y / CUT
    origin = np.stack([90.0 * Y / CUT, np.asarray(Y, float), Y * np.tan(np.deg2rad(5.0))], -1)
    cd = np.stack([np.cos(th), 0 * th, np.sin(th)], -1)
    ud = np.stack([-np.sin(th), 0 * th, np.cos(th)], -1)
    return origin, cd, ud, 180.0 - 80.0 * Y / CUT


def chord_fraction(P, Y):
    """Chord fraction of a point measured in its own (continued) section frame."""
    o, cd, _, c = section_frame_v(Y)
    return np.sum((P - o) * cd, -1) / c


# ---- attachment thickness table and the printable rule ---------------------------
_tab = None


def _thickness_table():
    global _tab
    if _tab is None:
        g = np.linspace(0, 1, 4001)
        xi = np.unique(np.concatenate([0.5 - 0.5 * np.cos(np.pi * g),
                                       np.geomspace(1e-7, 1e-2, 400),
                                       1 - np.geomspace(1e-7, 3e-2, 400)]))
        xi = xi[(xi >= 0) & (xi <= 0.9999989206029309)]
        t = pa.section_thickness(xi, CUT)
        ok = np.isfinite(t)
        _tab = (xi[ok], t[ok])
    return _tab


def t_att(chi):
    xi, t = _thickness_table()
    chi = np.asarray(chi, float)
    return np.interp(np.clip(chi, 0, xi[-1]), xi, t)


def rule_radius(chi):
    t = t_att(chi)
    return np.where(t >= RULE_T, RULE_MIN_R, t / 2.0)


# ---- rim and its foot on the mid surface ------------------------------------------
def rim(p):
    j = pa.tip_outline_jet(np.asarray(p, float))
    return j["point"], j["d_p"]


def foot(p, iters=40):
    """Solve rim(p) = M(s, Y) + kappa * lamD(s, Y) for (s, Y, kappa)."""
    p = np.asarray(p, float)
    R, _ = rim(p)
    s = np.clip(p.copy(), 1e-6, 1.0)
    Y = R[..., 1].copy()
    k = np.zeros_like(p)
    h = 1e-7
    for it in range(iters):
        def F(s_, Y_, k_):
            M, LD, _ = MD(s_, Y_)
            return M + k_[..., None] * LD - R
        f0 = F(s, Y, k)
        Fs = (F(np.minimum(s + h, 1.0 + 0 * s), Y, k) - F(np.minimum(s + h, 1.0) - h, Y, k)) / h
        FY = (F(s, Y + h, k) - f0) / h
        _, LD, _ = MD(s, Y)
        J = np.stack([Fs, FY, LD], -1)
        step = np.linalg.solve(J, f0[..., None])[..., 0]
        s = np.clip(s - step[..., 0], 1e-9, 1.0)
        Y = np.maximum(Y - step[..., 1], CUT)
        k = k - step[..., 2]
        if np.max(np.abs(step)) < 1e-13:
            break
    res = np.linalg.norm(F(s, Y, k), axis=-1)
    return s, Y, k, res


# ---- per-guide design ------------------------------------------------------------
def _feature(pp):
    return (0.03 ** -4 + (0.25 * (pp + 2e-4)) ** -4 + (0.2 * (1 - pp) + 1e-4) ** -4) ** -0.25


_GSPL = None


def p_of_g(g):
    """Rim parameter p as a function of the surface parameter g in [0, 1].

    Samples are spread in proportion to 1/feature(p), feature(p) = smooth min of
    (0.03, 0.25 (p + 2e-4), 0.2 (1 - p) + 1e-4): log-like clustering into both
    corners (leading-edge pole, trailing-edge corner transition).  G(p) is a C2
    spline of the cumulative density, inverted by Newton steps so that p(g) is
    smooth (a piecewise-linear inverse leaves kinks that show up in the fit).
    """
    global _GSPL
    if _GSPL is None:
        pp = np.unique(np.concatenate([np.geomspace(1e-10, 0.5, 6000), 1 - np.geomspace(1e-10, 0.5, 6000), [0.0, 1.0]]))
        dens = 1.0 / _feature(pp)
        G = np.concatenate([[0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(pp))])
        _GSPL = (CubicSpline(pp, G / G[-1]), pp, G / G[-1])
    spl, pp, G = _GSPL
    g = np.clip(np.asarray(g, float), 0.0, 1.0)
    p = np.interp(g, G, pp)
    for _ in range(30):
        p = np.clip(p - (spl(p) - g) / np.maximum(spl(p, 1), 1e-300), 0.0, 1.0)
    return p


def phi_of_eta(eta):
    """Trailing-edge corner transition: 0 at the attachment (sharp corners of the blunt TE), 1 from ETA_F."""
    t = np.clip(np.asarray(eta, float) / ETA_F, 0.0, 1.0)
    return 1 - (1 - t * t) ** 3


def Phi(v, m):
    v = np.clip(v, 0.0, 1.0)
    return np.sqrt(np.maximum(1 - (1 - v * v) ** m, 0.0))


def Phi_inv_v2(w, m):
    """v^2 such that Phi(v, m) = w."""
    return 1 - (1 - np.clip(w, 0, 1) ** 2) ** (1.0 / m)


def M_jac(s, Y):
    up, lo = skins(s, Y)
    Ms = 0.5 * (up["d_chord"] + lo["d_chord"])
    MY = 0.5 * (up["d_span"] + lo["d_span"])
    return Ms, MY


def _radius_ratio(b, w, Hfun, bgrid, reqgrid, v):
    """Planar model of a band section: curve (b, w * H(b)); returns min over v of radius / rule."""
    H = Hfun(b)
    req = np.interp(b, bgrid, reqgrid)
    z = w * H
    db = np.gradient(b, v, axis=-1); dz = np.gradient(z, v, axis=-1)
    ddb = np.gradient(db, v, axis=-1); ddz = np.gradient(dz, v, axis=-1)
    k = np.abs(db * ddz - dz * ddb) / np.maximum(db * db + dz * dz, 1e-300) ** 1.5
    r = 1.0 / np.maximum(k[..., 2:-2], 1e-300)
    return np.min(r / req[..., 2:-2], axis=-1)


V_GRID = np.linspace(1e-3, 1.0, 401)


class Guide:
    """Profile data of one guide (one rim parameter p).

    Coordinates of the profile are physical, in the band plane: x is the distance
    from the lip apex inward along the band line, y the half-thickness measured
    perpendicular to it.  The parent's thickness offset D is not perpendicular to
    the band line (sh = D^.nu, up to 0.08 from dihedral and camber slope); the
    realisation in lip_points shifts each sample along the band so that the
    physical section is the designed curve (otherwise a circle comes out as a
    leaning ellipse whose smallest radius is ~10 % too small).
    """

    def __init__(self, bgrid, Hgrid, reqgrid, m, muN, muC, betaLE, sh=0.0, st=0.0):
        self.bgrid, self.Hgrid, self.reqgrid = bgrid, Hgrid, reqgrid
        self.m, self.muN, self.muC, self.betaLE = m, muN, muC, betaLE
        self.sh, self.st = sh, st          # components of the offset direction along the band line and the rim
        self.c = np.sqrt(1 - sh * sh - st * st)
        self.b_r = 0.0
        self._H2 = CubicSpline(bgrid, Hgrid ** 2)     # H^2 is smooth through the LE nose, H is not
        self.R_c, self.bt_c = None, None

    def __getstate__(self):
        d = dict(self.__dict__)
        d.pop("_H2", None)
        return d

    def __setstate__(self, d):
        self.__dict__.update(d)
        self._H2 = CubicSpline(self.bgrid, self.Hgrid ** 2)

    def Hfoot(self, beta):
        """|lam D| along the band line at chart distance beta from the rim foot."""
        return np.sqrt(np.maximum(self._H2(beta), 0.0))

    def H(self, x):
        """Design skin height (perpendicular to the band line) at distance x from the apex."""
        return self.c * self.Hfoot(np.asarray(x, float) - self.b_r)

    def dH(self, x):
        h = 1e-6
        return (self.H(x + h) - self.H(x - h)) / (2 * h)

    def zeta_skin(self, x, sgn, sh_eff):
        """Offset magnitude along D of the skin at physical band position x (side sgn)."""
        x = np.asarray(x, float)
        z = self.Hfoot(x - self.b_r)
        if np.all(np.asarray(sh_eff) == 0):
            return z
        for _ in range(30):
            z = 0.5 * z + 0.5 * self.Hfoot(x - sgn * z * sh_eff - self.b_r)
        return z

    def Z(self, x, W, sh_eff):
        """Realised signed offset along D for design sample (x, W), W in [-1, 1] (signed fraction).

        Z = W * zs_sym + W^2 * dz with zs_sym, dz the mean and half-difference of the upper and
        lower skin offsets at physical position x: exactly the skin at W = +-1, and smooth
        (no kink) through the apex W = 0 although the two skins differ under the shear.
        """
        zu = self.zeta_skin(x, 1.0, sh_eff); zl = self.zeta_skin(x, -1.0, sh_eff)
        W = np.asarray(W, float)
        return W * 0.5 * (zu + zl) + W * W * 0.5 * (zu - zl)

    def inscribed_circle(self):
        """Circle centred on the mid surface, touching the apex line and tangent to the skin."""
        H0 = float(self.H(0.0))
        R, bt = H0, H0
        for _ in range(100):
            Ht = float(self.H(bt)); dH = float(self.dH(bt))
            R_new = Ht * np.sqrt(1 + dH * dH)
            bt_new = R_new - Ht * dH
            done = abs(R_new - R) < 1e-14 and abs(bt_new - bt) < 1e-14
            R, bt = R_new, bt_new
            if done:
                break
        self.R_c, self.bt_c = R, bt

    def curve(self, ell, v):
        """(x, w) of the blended profile for band width ell at parameters v (w: fraction of the skin height)."""
        ell = np.asarray(ell, float)
        v = np.asarray(v, float)
        muN, muC, m = self.muN, self.muC, self.m
        bP = ell * v * v
        wP = Phi(v, m)
        w = wP * (1 - muN - muC)
        b = bP * (1 - muC)
        if muN > 0:
            Psi = 1 - smooth5(bP / np.maximum(ell, 1e-300))
            wN = self.H(bP - (self.betaLE - self.b_r) * Psi) / np.maximum(self.H(bP), 1e-300)
            w = w + muN * wN
        if muC > 0:
            bC = self.bt_c * v * v
            zC = self.R_c * np.sqrt(np.maximum(1 - (1 - bC / self.R_c) ** 2, 0.0))
            wC = zC / np.maximum(self.H(bC), 1e-300)
            w = w + muC * wC
            b = b + muC * bC
        return b, np.clip(w, 0.0, 1.0)

    def ratio(self, ell):
        ell = np.atleast_1d(np.asarray(ell, float))
        b, w = self.curve(ell[:, None], V_GRID[None, :])
        return _radius_ratio(b, w, self.H, self.bgrid + self.b_r, self.reqgrid, V_GRID)


def design(p):
    """Per-guide data for rim parameters p (0 <= p <= 1; p = 0 is the leading-edge pole)."""
    p = np.asarray(p, float)
    pin = np.clip(p, 1e-12, 1.0)
    s_f, Y_f, kap, res = foot(pin)
    if np.max(res) > 1e-9:
        raise RuntimeError("rim foot not converged: %g" % np.max(res))
    M_f, LD_f, _ = MD(s_f, Y_f)
    H_f = np.linalg.norm(LD_f, axis=-1)
    # foot-curve tangent by central differences in p
    h = 1e-6
    sa_, Ya_, _, _ = foot(np.clip(pin - h, 1e-12, 1)); sb_, Yb_, _, _ = foot(np.clip(pin + h, 0, 1))
    T = (cap_point(sb_, Yb_, 0 * sb_) - cap_point(sa_, Ya_, 0 * sa_))
    T /= np.linalg.norm(T, axis=-1, keepdims=True)
    Ms, MY = M_jac(s_f, Y_f)
    N = np.cross(Ms, MY); N /= np.linalg.norm(N, axis=-1, keepdims=True)
    # inward (into the cap): the rim runs LE -> TE with the cap on its right seen from above (N up)
    nu = np.cross(T, N); nu /= np.linalg.norm(nu, axis=-1, keepdims=True)
    Ma = Ms / np.maximum(dads(s_f), 1e-12)[..., None]
    J = np.stack([Ma, MY], -1)
    JTJ = np.einsum("...ki,...kj->...ij", J, J)
    JTn = np.einsum("...ki,...k->...i", J, nu)
    dch = np.linalg.solve(JTJ, JTn[..., None])[..., 0]
    a_f = a_of_s(s_f)
    eta_a = Y_f - CUT
    tTE = smooth7((p - P_TE[0]) / (P_TE[1] - P_TE[0]))
    m = np.where(p > 0.5, 3.0 - tTE, 2.0 + smooth5((H_f - 0.42) / (0.90 - 0.42)))
    # distance from the apex forward along the band line to the continued leading edge (a = 0)
    betaLE = np.where((p < 0.5) & (dch[:, 0] > 0), a_f / np.maximum(dch[:, 0], 1e-300), np.inf)
    ell_d = np.where(p < 0.5, (np.maximum(REACH_FRAC * eta_a, 1e-300) ** -4.0 + np.maximum(NOSE_FAC * betaLE, 1e-300) ** -4.0) ** -0.25,
                     REACH_FRAC * eta_a)
    ell_d = np.maximum(ell_d, 1e-3 * H_f + 1e-12)
    muN = np.where(p < 0.5, 1 - smooth5((p - P_NOSE[0]) / (P_NOSE[1] - P_NOSE[0])), 0.0)
    muC = np.where(p > 0.5, tTE, 0.0)
    # sample H and the rule along each band line (extended forward to the LE where needed)
    bmax = np.maximum(3.0 * np.maximum(ell_d, H_f), 1e-6)
    bmin = np.where(muN > 0, -np.minimum(betaLE, 1e3), -0.05 * bmax)
    nb = 200
    tb = np.linspace(0, 1, nb) ** 2
    nneg = 80
    tn = 0.5 + 0.5 * np.cos(np.pi * np.linspace(0, 1, nneg, endpoint=False))
    bgrid = np.concatenate([bmin[:, None] * tn[None, :], bmax[:, None] * tb[None, :]], axis=1)
    ach = a_f[:, None] + bgrid * dch[:, None, 0]
    Ych = Y_f[:, None] + bgrid * dch[:, None, 1]
    sch = s_of_a(ach)
    Mb, LDb, _ = MD(sch, np.maximum(Ych, CUT))
    Hgrid = np.linalg.norm(LDb, axis=-1)
    chib = chord_fraction(Mb, Ych)
    reqgrid = rule_radius(chib)
    # band width: smooth minimum of the design width and the width at which the superelliptic
    # profile's apex radius drops to MARGIN * 0.4 mm (irrelevant where the nose / circle profiles rule)
    ell_apex = m * H_f ** 2 / (2 * RULE_MIN_R * MARGIN) / np.maximum(1 - muN, 1e-12)
    ell_s = (ell_d ** -6.0 + ell_apex ** -6.0) ** (-1 / 6.0)
    Dhat_f = LD_f / np.maximum(H_f, 1e-300)[..., None]
    shear = np.sum(Dhat_f * nu, -1)
    stilt = np.sum(Dhat_f * T, -1)
    JTt = np.einsum("...ki,...k->...i", J, T)
    dTch = np.linalg.solve(JTJ, JTt[..., None])[..., 0]
    guides = []
    ell = np.empty_like(p)
    status = np.empty(p.shape, dtype=object)
    fach = np.empty_like(p)
    for i in range(len(p)):
        bg, Hg = bgrid[i], Hgrid[i]
        g = Guide(bg, Hg, reqgrid[i], m[i], muN[i], muC[i], betaLE[i], shear[i], stilt[i])
        g.inscribed_circle()              # inscribed circle tangent to the skin (trailing-edge thin zone)
        guides.append(g)
        if p[i] <= 0:
            ell[i] = 0.0; status[i] = "pole"; fach[i] = np.nan
            continue
        ell[i] = ell_s[i]
        status[i] = "apex-limited" if ell_apex[i] < ell_d[i] else "design"
        fach[i] = g.ratio(ell[i])[0]
    # trailing-edge corner transition (only at the trailing-edge end of the rim)
    phi = np.where(p > 0.5, phi_of_eta(eta_a), 1.0)
    end = p >= 1.0
    if np.any(end) and np.any(~end):
        j = np.argmax(np.where(~end, p, -1))
        ell[end] = ell[j]; fach[end] = fach[j]
    # outward shift of the apex (b_r: chart distance apex -> rim foot) so that the rim lies on the
    # lower side of the realised profile; iterated with the inscribed circle (both depend on b_r)
    ak = np.abs(kap)
    wall = 1 - phi
    b_r = np.zeros_like(p)
    for i in range(len(p)):
        g = guides[i]
        if p[i] <= 0:
            continue
        sh_eff = phi[i] * (1 - g.muN) * g.sh
        zt = ak[i] * H_f[i]                      # the rim's offset magnitude along D at its foot
        for it in range(6):
            if zt <= -float(g.Z(0.0, -wall[i], sh_eff)):
                x_r = 0.0
            else:
                # realised lower-side offset -Z(x(v), -w(v)) must equal |kappa| H_f
                lo, hi = 0.0, 1.0
                for _ in range(60):
                    vm = 0.5 * (lo + hi)
                    bb, ww = g.curve(ell[i], vm)
                    xx = phi[i] * float(bb)
                    zz = -float(g.Z(xx, -(wall[i] + phi[i] * float(ww)), sh_eff))
                    if zz < zt: lo = vm
                    else: hi = vm
                bb, ww = g.curve(ell[i], 0.5 * (lo + hi))
                x_r = phi[i] * float(bb)
            # lower side (sgn = -1): chart position b' = x + zeta * sh_eff with zeta = |kappa| H_f
            g.b_r = x_r + zt * sh_eff
            g.inscribed_circle()
        b_r[i] = g.b_r
        fach[i] = g.ratio(ell[i])[0]
    # band edge and the attachment pairing (flat-region path at constant chord fraction)
    bend = np.array([guides[i].curve(ell[i], 1.0)[0] for i in range(len(p))]) * phi
    sg = 0.5 - 0.5 * np.cos(np.pi * np.linspace(0, 1, 40001))
    chig = chord_fraction(cap_point(sg, np.full_like(sg, CUT), 0 * sg), np.full_like(sg, CUT))
    side_edge = {}
    for side in (+1, -1):
        she = phi * (1 - muN) * shear
        ste = phi * (1 - muN) * stilt
        zse = np.array([guides[i].zeta_skin(bend[i], side, she[i]) if p[i] > 0 else 0.0 for i in range(len(p))])
        bpe = bend - side * zse * she
        Eb = bpe - b_r
        tsh = -side * zse * ste
        aE = a_f + Eb * dch[:, 0] + tsh * dTch[:, 0]; YE = np.maximum(Y_f + Eb * dch[:, 1] + tsh * dTch[:, 1], CUT)
        sE = s_of_a(aE)
        ME = cap_point(sE, YE, 0 * sE)
        chiE = chord_fraction(ME, YE)
        sA = np.interp(np.clip(chiE, chig[0], chig[-1]), chig, sg)
        side_edge[side] = (aE, YE, sA)
    aE, YE, sA = side_edge[+1]
    out = dict(p=p, s_f=s_f, Y_f=Y_f, kappa=kap, a_f=a_f, dch=dch, nu=nu, H_f=H_f, eta_a=eta_a,
               m=m, ell=ell, ell_d=ell_d, status=status, fach=fach, stilt=stilt, dTch=dTch, phi=phi, b_r=b_r, aE=aE, YE=YE,
               sA=sA, aA=a_of_s(sA), muN=muN, muC=muC, betaLE=betaLE, guides=guides, bend=bend, shear=shear,
               aE_up=side_edge[1][0], YE_up=side_edge[1][1], aA_up=a_of_s(side_edge[1][2]),
               aE_lo=side_edge[-1][0], YE_lo=side_edge[-1][1], aA_lo=a_of_s(side_edge[-1][2]))
    pole = p <= 0
    if np.any(pole):
        for k_ in ("ell", "b_r", "eta_a", "H_f", "bend"):
            out[k_][pole] = 0.0
        out["phi"][pole] = 1.0
        out["sA"][pole] = 0.0; out["aA"][pole] = 0.0
        out["a_f"][pole] = 0.0; out["Y_f"][pole] = CUT; out["aE"][pole] = 0.0; out["YE"][pole] = CUT
        for k_ in ("aE_up", "aE_lo", "aA_up", "aA_lo"):
            out[k_][pole] = 0.0
        out["dTch"][pole] = 0.0
        out["YE_up"][pole] = CUT; out["YE_lo"][pole] = CUT
        out["dch"][pole] = 0.0
    return out


# ---- sampling of the three faces --------------------------------------------------
def lip_points(dz, sig):
    """Lip face: sig in [-1, 1] (+ upper side), apex at 0, band edges at +-1. Returns (Np, Ns, 3).

    Each designed sample (x, w) is realised as  M(b') + sgn * zeta * D^(b')  with
    zeta = w * (skin offset at physical position x) and b' = x - sgn * zeta * sh_eff,
    so the section in the band plane is the designed curve; at w = 1 the point is
    exactly the parent skin continued.
    """
    sig = np.asarray(sig, float)
    p = dz["p"]
    Np, Ns = len(p), len(sig)
    v = np.abs(sig)
    sg = np.where(sig >= 0, 1.0, -1.0)
    bch = np.zeros((Np, Ns)); fr = np.zeros((Np, Ns)); tch = np.zeros((Np, Ns))
    for i in range(Np):
        if p[i] <= 0:
            continue
        g = dz["guides"][i]
        phi = dz["phi"][i]; wall = 1 - phi
        sw = wall / max(wall + np.sqrt(g.m) * phi, 1e-300) if wall > 0 else 0.0
        inwall = v <= sw
        vp = np.where(inwall, 0.0, (v - sw) / max(1 - sw, 1e-300))
        bb, ww = g.curve(dz["ell"][i], vp)
        x = np.where(inwall, 0.0, phi * bb)
        w = np.where(inwall, wall * v / max(sw, 1e-300), wall + phi * ww)
        sh_eff = phi * (1 - g.muN) * g.sh       # no correction on the shifted parent nose (LE corner)
        Zs = g.Z(x, sg * w, sh_eff)             # signed realised offset along D
        bch[i] = x - Zs * sh_eff
        tch[i] = -Zs * phi * (1 - g.muN) * g.st  # shift along the rim so the section is perpendicular to it
        fr[i] = np.where(w >= 1.0, 1.0, np.abs(Zs) / np.maximum(g.Hfoot(bch[i] - g.b_r), 1e-300))
    beta = bch - dz["b_r"][:, None]
    a = dz["a_f"][:, None] + beta * dz["dch"][:, None, 0] + tch * dz["dTch"][:, None, 0]
    Y = dz["Y_f"][:, None] + beta * dz["dch"][:, None, 1] + tch * dz["dTch"][:, None, 1]
    s = s_of_a(a)
    # fr is zeta / |lam D| from the band-line spline; evaluate the exact |lam D| at the chart point
    M, LD, _ = MD(s, np.maximum(Y, CUT))
    Hex = np.linalg.norm(LD, axis=-1)
    Hsp = np.zeros_like(Hex)
    for i in range(Np):
        if p[i] > 0:
            Hsp[i] = dz["guides"][i].Hfoot(beta[i])
    wfac = np.where(fr >= 1.0 - 1e-12, 1.0, fr)
    P = M + (sg[None, :] * wfac)[..., None] * LD
    return P, dict(a=a, Y=Y, s=s, w=sg[None, :] * wfac, b=bch)


def flat_path(dz, t, side=+1):
    """Chart points of the flat-region path from the band edge (t=0) to the attachment (t=1): a straight chart line."""
    t = np.asarray(t, float)[None, :]
    k = "_up" if side > 0 else "_lo"
    a = (1 - t) * dz["aE" + k][:, None] + t * dz["aA" + k][:, None]
    Y = (1 - t) * dz["YE" + k][:, None] + t * CUT
    return a, Y


def flat_points(dz, t, side):
    a, Y = flat_path(dz, t, side)
    s = s_of_a(a)
    P = cap_point(s, np.maximum(Y, CUT), np.full(s.shape, float(side)))
    return P, dict(a=a, Y=Y, s=s)
