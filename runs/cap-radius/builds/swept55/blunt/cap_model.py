"""Analytic definition of the rounded tip cap (right tip, upper patch; lower = mirror in Z).

Swept55 blunt parent. Lengths in mm. Needs NumPy and the locked parent.py.

Construction (skill 'close-rounded-wingtip', as a boundary surface defined analytically):
  * first profile  = the parent's upper cut edge at Y = 390 (attachment chord fraction xi_a),
  * middle profile = the locked perimeter (lip) at Z = 0 (lip()),
  * guides         = one spanwise wrap curve per patch parameter u: in plan a straight line in
    (chord fraction, station) space from the attachment point to the lip point; sigma runs
    0 -> 1 along it.  Away from the trailing edge the guides follow the sweep (constant chord
    fraction of the continued planform).  Near the blunt TE they fan forward (guide_ends()),
    and the last guide (u = 1) runs along the TE line to the half-round apex.
  * height: a roll field tau (tau = 0 at the lip, >= 1 + h in the flat run next to the cut);
    tau = k (1 - sigma) on the guides, blended near the TE corner into te_field(), a
    parabolic-coordinate field whose tau = 1 - a/0.2 on the TE line gives the half-round
    of radius 0.2 that closes the end of the flat TE face.
    z = Zp(P) - H_lip * Delta(soft-max(1 - tau, 0)), Zp = the parent's own half thickness
    continued outboard (taper let out over D_taper), Delta = a stretched circular arc that is
    vertical at the lip and tangent to the parent's slope at the roll top: convex everywhere,
    exactly the parent (G2) in the flat run, lip radius R = H_lip^2 |grad tau| (1 + mu).
  * per guide, k (hence the roll width W = 1 / |grad tau|) is set from the lip-radius target
    min(R_design, kappa * H_lip) (profile_params()).
"""
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "parents", "swept55", "blunt"))


def _load_parent():
    spec = importlib.util.spec_from_file_location("parent_swept55_blunt", os.path.join(PARENT_DIR, "parent.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


parent = _load_parent()

ATT = parent.ATTACHMENT                 # 390
SWEEP = parent.SWEEP
TAPER = parent.TAPER
C_ATT = float(parent.chord(ATT))        # 95
TE_HALF = parent.TE_THICKNESS / 2       # 0.2
TE_SLOPE = SWEEP - TAPER                # dX/dr of the trailing-edge line
COS_TE = 1.0 / np.sqrt(1.0 + TE_SLOPE**2)

# ----------------------------------------------------------------------------- design
DESIGN = dict(
    R_design=0.6,      # target lip radius where the section is thick enough (mm)
    kappa=1.08,        # lip target where thinner: kappa * local half thickness (x (1 + mu) near the TE)
    kappa_te=1.0,      # TE-corner field: lip radius = kappa_te * local half thickness (the half-round)
    wedge_p=(0.75, 0.9),  # parent slope normal to the lip used by the arc (mu) for lip p beyond these
    m_max=0.5,         # cap on that slope
    D_taper=3.0,       # chord taper let-out length (mm)
    k_floor=1.3,       # minimum tail slope of the guide field (keeps a flat start -> G2 at the cut)
    h_max=0.3,         # half width (in tau) of the C2 flat-to-ellipse ramp
    F0=2.0,            # TE fan: attachment offset saturates at F0 mm forward of the TE
    fan_beta=0.06,     # TE fan keeps opening at this rate past the saturation (no stalled attachment)
    rho0=60.0,         # TE fan switched off between rho0 ...
    rho1=90.0,         # ... and rho1 (mm along the TE line from the corner; p = 0.5 is at 90.7)
    kappa0=0.36,       # TE corner field: parabolic coefficient at the half-round apex ...
    kappa_c=(0.0, 0.0),  # ... optional variation kappa0 exp(-c1 s + c2 s^2), s = (a - 0.2)/0.2 (unused)
    w_te=(2.0, 6.0),   # TE corner field used for lip distances < 2 mm, guide field beyond 6 mm
    h_te_d=0.5,        # ramp half width grows from 0 at the TE corner over this distance (mm)
)


def U(p):
    return p * p * (3 - 2 * p)


def Uinv(y):
    y = np.clip(y, 0.0, 1.0)
    return 0.5 - np.sin(np.arcsin(1 - 2 * y) / 3)


def smooth5(t):
    t = np.clip(t, 0.0, 1.0)
    return t**3 * (10 - 15 * t + 6 * t * t)


def chord(r):
    return parent.ROOT_CHORD - TAPER * r


def scale_at(r):
    """Keep-chord scale of the blunt cut at station r (parent definition, continued outboard)."""
    r = np.asarray(r, float)
    return np.asarray(parent.cut(r)["scale"], float)


def Zp(xi, r):
    """Parent half thickness (mm) at chord fraction xi and station r (continued past 390)."""
    xi, r = np.broadcast_arrays(np.asarray(xi, float), np.asarray(r, float))
    sc = scale_at(r)
    u = np.clip(xi / sc, 0.0, None)
    t = 0.45 * (0.2969 * np.sqrt(u) - 0.1260 * u - 0.3516 * u**2 + 0.2843 * u**3 - 0.1036 * u**4)
    return chord(r) * sc * t


def r_eff(r):
    D = DESIGN["D_taper"]
    return ATT + D * np.tanh((r - ATT) / D)


def plan_xy(xi, r):
    return SWEEP * r + chord(r) * xi, r


def _xir_eff(X, Y):
    xi, r = xy_to_xir(X, Y)
    return xi, r_eff(r)


def xy_to_xir(X, Y):
    return (X - SWEEP * Y) / chord(Y), Y


# ----------------------------------------------------------------------------- lip
R_APEX = ATT + TE_HALF * COS_TE                      # station of the half-round apex
P_APEX = float(0.5 + 0.5 * np.sqrt(1 - 4 * (R_APEX - ATT) / 240.0))
P_B = 1.0 - 2.0e-3                                   # lip nudged onto the TE plane after P_B


def lip(p):
    """Lip (perimeter) in (xi, r): the locked outline, nudged (< 1.2 um) onto the TE plane at the apex."""
    p = np.asarray(p, float)
    xi0 = U(p)
    s = smooth5((p - P_B) / (P_APEX - P_B))
    return xi0 + (1 - xi0) * s, ATT + 240 * p * (1 - p)


def lip_inward_normal(p, h=1e-7):
    """Unit plan normal of the lip pointing into the cap (right tip)."""
    p = np.asarray(p, float)
    pa, pb = np.clip(p - h, 0, P_APEX), np.clip(p + h, 0, P_APEX)
    Xa, Ya = plan_xy(*lip(pa))
    Xb, Yb = plan_xy(*lip(pb))
    t = np.stack([Xb - Xa, Yb - Ya], -1)
    t /= np.linalg.norm(t, axis=-1, keepdims=True)
    return np.stack([t[..., 1], -t[..., 0]], -1)


# ----------------------------------------------------------------------------- guides
E_A = np.array([TE_SLOPE, 1.0]) * COS_TE        # TE direction in plan (outboard)
E_B = np.array([-1.0, TE_SLOPE]) * COS_TE       # forward normal of the TE line in plan
A_DOT_EA, A_DOT_EB = -E_A[0], -E_B[0]           # attachment direction (-X) in that frame
P0_XY = np.array([SWEEP * ATT + C_ATT, ATT])    # TE corner (centre of the half-round) in plan

def guide_ends(u):
    """Attachment chord fraction and lip point of guide u (u = p / P_APEX)."""
    u = np.asarray(u, float)
    p = P_APEX * u
    xl, rl = lip(p)
    d = DESIGN
    rho_o = (rl - ATT) / COS_TE
    rho_cx = C_ATT * (1 - xl)
    lam = 4 * TE_HALF / (abs(A_DOT_EB) * d["kappa0"])               # makes the fan match the corner field
    e_o = np.maximum(rho_o - TE_HALF, 0)
    fan = d["F0"] * np.tanh(lam * e_o / d["F0"]) + d["fan_beta"] * e_o
    fan = fan * (1 - smooth5((rho_o - d["rho0"]) / (d["rho1"] - d["rho0"])))
    rho_a = (np.maximum(rho_cx, 0) ** 4 + fan**4) ** 0.25          # smooth max: monotone, = rho_cx far out
    xa = np.where(p > 0.5, 1 - rho_a / C_ATT, xl)
    return xa, xl, rl, p




def plan_map(u, sigma):
    xa, xl, rl, _ = guide_ends(u)
    xi = xa + sigma * (xl - xa)
    r = ATT + sigma * (rl - ATT)
    return xi, r


def grad_sigma_at_lip(u, h=1e-5):
    """|grad sigma| in plan at the lip end of guide u (1/mm): |dO/du| / |dO/du x dP/dsigma|."""
    u = np.asarray(u, float)
    p = P_APEX * u
    Xa, Ya = plan_xy(*lip(p - h))
    Xb, Yb = plan_xy(*lip(p + h))
    Pu = np.stack([Xb - Xa, Yb - Ya], -1) * (P_APEX / (2 * h))
    xa, xl, rl, _ = guide_ends(u)
    dr = rl - ATT
    dxi = xl - xa
    Ps = np.stack([SWEEP * dr - TAPER * dr * xl + chord(rl) * dxi, dr], -1)
    det = Pu[..., 0] * Ps[..., 1] - Pu[..., 1] * Ps[..., 0]
    return np.hypot(Pu[..., 0], Pu[..., 1]) / np.abs(det)


def profile_params(u):
    """Per-guide data: guide-field slope k (tau = k (1 - sigma)), roll width W, arc parameter mu,
    ramp half width h_g and the TE-corner weight w.

    Lip radius R = H_lip^2 k g (1 + mu), g = |grad sigma| at the lip, mu = m W / H_lip with m the
    parent's slope normal to the lip (only used near the TE, where the section thickens forward)."""
    u = np.asarray(u, float)
    d = DESIGN
    xa, xl, rl, p = guide_ends(u)
    zl = Zp(xl, r_eff(rl))
    g = grad_sigma_at_lip(u)
    # parent slope normal to the lip, into the cap (the TE wedge thickens forward); used near the TE
    OX, OY = plan_xy(xl, rl)
    n_in = lip_inward_normal(p)
    dn = 1e-3
    m = (Zp(*_xir_eff(OX + dn * n_in[..., 0], OY + dn * n_in[..., 1])) - zl) / dn
    chi = smooth5((p - d["wedge_p"][0]) / (d["wedge_p"][1] - d["wedge_p"][0]))
    m_eff = chi * np.clip(m, 0.0, d["m_max"])
    H = zl
    k = np.full(np.shape(u), 2.0)
    for _ in range(8):
        R_t = np.minimum(d["R_design"], d["kappa"] * H)
        W = 1 / (k * g)
        mu = np.minimum(m_eff * W / H, 0.9)
        k_t = R_t / (H**2 * g * (1 + mu))
        k = (np.maximum(k_t, 0) ** 4 + d["k_floor"] ** 4) ** 0.25        # smooth max(k_t, k_floor)
    W = 1 / (k * g)
    mu = np.minimum(m_eff * W / H, 0.9)
    sm = 1 - 1 / k
    h_g = np.minimum(d["h_max"], 0.8 * (k - 1))
    eps = np.maximum((rl - ATT) / COS_TE - TE_HALF, 0.0)
    w = np.where(p > 0.5, 1 - smooth5((eps - d["w_te"][0]) / (d["w_te"][1] - d["w_te"][0])), 0.0)
    return dict(k=k, sm=sm, W=W, mu=mu, m=m, h_g=h_g, w=w, R_target=R_t, R_lip_pred=k * H**2 * g * (1 + mu), zl=zl, H_top=H, g=g)


def phi_hat(x, h):
    """C2 'soft max(x, 0)': 0 for x <= -h, x for x >= h, smoothstep-integral ramp between."""
    hh = np.maximum(h, 1e-300)
    t = np.clip((x + hh) / (2 * hh), 0.0, 1.0)
    ramp = 2 * hh * (t**3 - 0.5 * t**4)
    return np.where(x >= hh, x, np.where(x <= -hh, 0.0, ramp))


def _lip_local(pp):
    xl, rl = lip(pp)
    X, Y = plan_xy(xl, rl)
    D = np.stack([X - P0_XY[0], Y - P0_XY[1]], -1)
    return D @ E_A, D @ E_B


def te_field(X, Y):
    """TE-corner tau field from plan coordinates (a along the TE from the corner, b' forward of the lip).

    tau = (sqrt(A^2 + kappa(a) b') - A) / (2 * 0.2),  A = a - 0.2:  tau = 1 - a/0.2 on the TE line
    inboard of the apex (the half-round), tau ~ kappa b' / (0.8 A) along the lip outboard of it.
    kappa grows with A so the lip radius settles at the target (half-round where thin)."""
    d = DESIGN
    D = np.stack([X - P0_XY[0], Y - P0_XY[1]], -1)
    a, b = D @ E_A, D @ E_B
    # lip offset at the same a (Newton on the lip parameter, p near 1)
    p = np.full(np.shape(a), P_APEX)
    for _ in range(30):
        la, _ = _lip_local(p)
        dp = 1e-7
        la2, _ = _lip_local(p - dp)
        dadp = (la - la2) / dp
        p = np.clip(p - (la - np.maximum(a, TE_HALF)) / dadp, 0.5, P_APEX)
    la, lb = _lip_local(p)
    bp = np.maximum(b - lb, 0.0)
    xl, rl = lip(p)
    zl = Zp(xl, r_eff(rl))
    R_t = np.minimum(d["R_design"], d["kappa_te"] * zl)
    A = a - TE_HALF
    s = A / TE_HALF
    k0 = d["kappa0"] * np.exp(-d["kappa_c"][0] * s + d["kappa_c"][1] * s * s)
    kap = (k0**4 + np.maximum(2 * TE_HALF * 2 * A * R_t / zl**2, 0) ** 4) ** 0.25
    tau = (np.sqrt(A * A + kap * bp) - A) / (2 * TE_HALF)
    dist0 = np.hypot(a, b)
    h = d["h_max"] * smooth5(dist0 / d["h_te_d"])
    return tau, h


def surface(u, v, params=None):
    """Upper cap point(s) at patch parameters u (chordwise, 0 = LE corner, 1 = TE half-round)
    and v (0 = attachment, 1 = lip).  sigma = 2v - v^2 keeps the patch regular at the lip."""
    u, v = np.broadcast_arrays(np.asarray(u, float), np.asarray(v, float))
    pp = profile_params(u) if params is None else params
    k, mu, zl, h_g, w = pp["k"], pp["mu"], pp["zl"], pp["h_g"], pp["w"]
    sigma = 2 * v - v * v
    xi, r = plan_map(u, sigma)
    X, Y = plan_xy(xi, r)
    tau = k * (1 - v) ** 2
    h = h_g
    if np.any(w > 0):
        tt, ht = te_field(X, Y)
        tau = np.where(w > 0, (1 - w) * tau + w * tt, tau)
        h = np.where(w > 0, (1 - w) * h_g + w * ht, h_g)
    # gap form: z = Zp(P) - H_lip * Delta(x), x = soft-max(1 - tau, 0).  Delta is the circular arc,
    # stretched, that is vertical at the lip and tangent to the parent's own slope m (mu = m W / H_lip)
    # at the roll top, so the turn is convex and may use the wedge of a thickening section.
    H = Zp(xi, r_eff(r))
    x = 1 - tau
    xe = phi_hat(x, h)
    c = 1 - mu
    N = np.sqrt(1 - mu * mu)
    te = 1 - xe                                        # = tau outside the ramp
    arc = np.sqrt(np.maximum(c * te * (2 - c * te), 0)) / N
    delta = 1 + mu * te - (1 + mu) * arc
    z = H - zl * delta
    z = np.where(v >= 1, 0.0, z)
    return np.stack([X, Y, z], -1)


def chord_fraction_of(X, Y):
    """Plan chord fraction of a point in the continued planform (|Y| used)."""
    r = np.abs(Y)
    return (X - SWEEP * r) / chord(r)
