"""1D design of the spanwise wrap guide (iteration 2: the p = 5 superellipse wrap).

At native chord parameter u the wrap guide lies in the plane X = X(u) and is

    Y - 350 = r(u) eta(v),    Z = h(u) zeta(v),    v in [0, 1]  (0 upper edge, 0.5 lip, 1 lower edge)

with r(u) = 200 u (1 - u) the locked outline reach and h(u) = u (1 - u) (34 - 20 u)
the parent's half thickness. The guide is the curve

    eta^P + zeta^2 = 1,   P = 5,

parametrised by t = pi v as eta = sin t, zeta = cos t * sqrt((1 + s + s^2 + s^3 + s^4) / (1 + s)),
s = sin t (analytic in t, including through the lip at t = pi / 2).

Properties, in units of (r, h):
  * at the attachment zeta = 1 - eta^5 / 2: the wrap leaves the cut edge tangent to the
    skin with zero spanwise curvature (G2; the first four derivatives of zeta vanish);
  * at the lip 1 - eta = zeta^2 / P...: a round lip with in-plane radius P h^2 / (2 r)
    and the in-plane curvature rising monotonically from the attachment into the lip;
  * the lip is at eta = 1 exactly: the perimeter is the locked outline;
  * near the leading- and trailing-edge corners, where the outline runs almost
    spanwise, the cap reproduces the parent's own LE / TE closure shifted along the
    outline (see README), so the surface's principal radii there stay at the
    parent's own values.

The functions eta(v), zeta(v) are represented by degree-5 B-spline functions on a
uniform knot vector with NV poles, mirrored about v = 0.5 (eta even, zeta odd).
Held exactly in the poles: eta_0 = 0, zeta_0 = zeta_1 = zeta_2 = 1 (G1 + G2 with the
parent skin, whose spanwise curvature is zero), eta(0.5) = 1, zeta_M = 0.
"""

import json
import numpy as np
from scipy.interpolate import BSpline

DEG = 5
M = 16                      # index of the middle (lip) pole
NV = 2 * M + 1              # poles across the full wrap
P_EXP = 5.0                 # superellipse exponent on eta


def v_knots():
    inner = np.linspace(0, 1, NV - DEG + 1)
    return np.concatenate([[0.0] * DEG, inner, [1.0] * DEG])


def basis(v):
    return BSpline.design_matrix(np.clip(np.asarray(v, float), 0, 1), v_knots(), DEG).toarray()


def exact(v):
    """eta(v), zeta(v) of the design curve for v in [0, 1]."""
    t = np.pi * np.asarray(v, float)
    s = np.abs(np.sin(t))
    q = np.sum([s**i for i in range(int(P_EXP))], axis=0)      # 1 + s + ... + s^(P-1)
    return s, np.cos(t) * np.sqrt(q / (1 + s))


def fit_poles(n=6001):
    v = np.linspace(0, 0.5, n)
    eta, zeta = exact(v)
    N = basis(v)
    Ev = np.zeros((n, M + 1)); Od = np.zeros((n, M + 1))
    for i in range(M + 1):
        Ev[:, i] = N[:, i] + (N[:, NV - 1 - i] if i < M else 0)
        Od[:, i] = N[:, i] - (N[:, NV - 1 - i] if i < M else 0)
    Nm = basis([0.5])[0]
    Evm = np.array([Nm[i] + (Nm[NV - 1 - i] if i < M else 0) for i in range(M + 1)])

    def constrained(A, y, fixed, eq_rows, eq_vals):
        free = [i for i in range(M + 1) if i not in fixed]
        fi, fv = list(fixed), np.array(list(fixed.values()), float)
        rhs = y - A[:, fi] @ fv
        Af = A[:, free]
        C = np.array([row[free] for row in eq_rows]).reshape(len(eq_rows), len(free))
        d = np.array([val - row[fi] @ fv for row, val in zip(eq_rows, eq_vals)])
        K = np.block([[Af.T @ Af, C.T], [C, np.zeros((len(C), len(C)))]])
        sol = np.linalg.solve(K, np.concatenate([Af.T @ rhs, d]))[:len(free)]
        out = np.zeros(M + 1); out[fi] = fv; out[free] = sol
        return out

    eh = constrained(Ev, eta, {0: 0.0}, [Evm], [1.0])
    zh = constrained(Od, zeta, {0: 1.0, 1: 1.0, 2: 1.0, M: 0.0}, [], [])
    full = lambda h, sgn: np.concatenate([h, sgn * h[-2::-1]])
    E, Z = full(eh, 1), full(zh, -1)
    vv = np.linspace(0, 1, 20001)
    e2, z2 = exact(vv)
    Nf = basis(vv)
    dev = {"eta": float(np.abs(Nf @ E - e2).max()), "zeta": float(np.abs(Nf @ Z - z2).max())}
    return E, Z, dev


def curve(v, k, E, Z, der=0):
    """Guide in units of h for aspect k = r / h: (k eta, zeta) and derivatives."""
    t = v_knots()
    se, sz = BSpline(t, k * E, DEG), BSpline(t, Z, DEG)
    return np.stack([se.derivative(der)(v) if der else se(v), sz.derivative(der)(v) if der else sz(v)], -1)


def analyse(k, E, Z, n=20001):
    v = np.linspace(0, 0.5, n)
    c, c1, c2 = (curve(v, k, E, Z, d) for d in range(3))
    sp = np.hypot(c1[:, 0], c1[:, 1])
    kap = -(c1[:, 0] * c2[:, 1] - c1[:, 1] * c2[:, 0]) / sp**3
    return {"v": v, "y": c[:, 0], "z": c[:, 1], "kappa": kap, "speed": sp}


def design():
    return fit_poles()


def report(E, Z, dev):
    rep = []
    for k in [5.88, 8.33, 10, 12, 14.29]:
        an = analyse(k, E, Z)
        kap = an["kappa"]
        dk = np.diff(kap[an["v"] > 0.02])
        rep.append({"k": k, "R_min_over_h": float(1 / kap.max()), "R_lip_over_h": float(1 / kap[-1]),
                    "R_lip_over_h_formula": P_EXP / (2 * k),
                    "min_kappa_h": float(kap.min()),
                    "largest_relative_kappa_drop_towards_lip": float(max(0.0, -(dk / kap.max()).min())),
                    "lip": [float(an["y"][-1]), float(an["z"][-1])],
                    "speed_ratio": float(an["speed"].max() / an["speed"].min())})
    return {"P": P_EXP, "NV": NV, "degree": DEG, "bspline_fit_max_dev": dev, "guides_units_of_h": rep}


if __name__ == "__main__":
    print(json.dumps(report(*design()), indent=1))
