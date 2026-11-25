"""Design-loop scan of the analytic cap: principal radii on a dense grid (finite differences).

Run from the project root:  .venv/bin/python runs/cap-radius/builds/swept55/blunt/design_scan.py
Prints the smallest radius against the printable rule.  (The release numbers come from
verify.py on the reimported STEP; this is only the design loop.)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cap_model as cm  # noqa: E402


def u_of_p(p):
    return p / cm.P_APEX


def u_of_rho(rho):
    dr = rho * cm.COS_TE
    p = 0.5 + 0.5 * np.sqrt(np.clip(1 - 4 * dr / 240.0, 0, None))
    return p / cm.P_APEX


def curvature(F, a, b, ha, hb):
    """Principal curvatures of F(a, b) -> xyz by central differences (arrays a, b)."""
    f = lambda da, db: F(a + da, b + db)
    P0 = f(0, 0)
    Pa = (f(ha, 0) - f(-ha, 0)) / (2 * ha)
    Pb = (f(0, hb) - f(0, -hb)) / (2 * hb)
    Paa = (f(ha, 0) - 2 * P0 + f(-ha, 0)) / ha**2
    Pbb = (f(0, hb) - 2 * P0 + f(0, -hb)) / hb**2
    Pab = (f(ha, hb) - f(ha, -hb) - f(-ha, hb) + f(-ha, -hb)) / (4 * ha * hb)
    n = np.cross(Pa, Pb)
    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    E, Fm, G = (Pa * Pa).sum(-1), (Pa * Pb).sum(-1), (Pb * Pb).sum(-1)
    L, M, N = (Paa * n).sum(-1), (Pab * n).sum(-1), (Pbb * n).sum(-1)
    den = E * G - Fm**2
    Hm = (E * N - 2 * Fm * M + G * L) / (2 * den)
    K = (L * N - M * M) / den
    disc = np.sqrt(np.maximum(Hm * Hm - K, 0))
    return P0, Hm + disc, Hm - disc, n


def scan(nv=200, na=300, verbose=True):
    v = 1 - (1 - np.linspace(0.0005, 0.9995, nv)) ** 1.0
    out = []
    # LE half: a = log p, p in [1e-4, 0.5]
    a = np.linspace(np.log(1e-4), np.log(0.5), na)
    A, V = np.meshgrid(a, v, indexing="ij")
    F = lambda aa, vv: cm.surface(u_of_p(np.exp(aa)), vv)
    out.append(("LE", A, V) + curvature(F, A, V, 3e-3, 1e-3))
    # TE half: a = log(rho - 0.2), rho in [0.2 + 1e-4, rho(p=0.5)]
    rho_mid = 60.0 / cm.COS_TE
    a = np.linspace(np.log(1e-4), np.log(rho_mid - 0.2), na)
    A, V = np.meshgrid(a, v, indexing="ij")
    F = lambda aa, vv: cm.surface(u_of_rho(0.2 + np.exp(aa)), vv)
    out.append(("TE", A, V) + curvature(F, A, V, 3e-3, 1e-3))
    res = {}
    for name, A, V, P, k1, k2, n in out:
        kmax = np.maximum(np.abs(k1), np.abs(k2))
        R = 1 / kmax
        xi = cm.chord_fraction_of(P[..., 0], P[..., 1])
        t_att = 2 * cm.Zp(np.clip(xi, 0, 1), 390.0)
        thick = t_att >= 0.8
        Rt = np.where(thick, R, np.nan)
        i = np.nanargmin(Rt)
        ratio = np.where(~thick, R / (t_att / 2), np.nan)
        j = np.nanargmin(ratio) if np.any(~thick) else None
        res[name] = dict(R_min_thick=float(Rt.flat[i]), at=P.reshape(-1, 3)[i].round(4).tolist(),
                         xi=float(xi.flat[i]), t_att=float(t_att.flat[i]), v=float(V.flat[i]),
                         a=float(A.flat[i]))
        if j is not None:
            res[name].update(thin_ratio_min=float(ratio.flat[j]), thin_at=P.reshape(-1, 3)[j].round(4).tolist(),
                             thin_xi=float(xi.flat[j]), thin_R=float(R.flat[j]), thin_v=float(V.flat[j]),
                             thin_a=float(A.flat[j]))
        res[name]["arrays"] = (A, V, P, k1, k2, R, xi, t_att)
    if verbose:
        for name in res:
            print(name, {k: v for k, v in res[name].items() if k != "arrays"})
    return res


if __name__ == "__main__":
    scan()
