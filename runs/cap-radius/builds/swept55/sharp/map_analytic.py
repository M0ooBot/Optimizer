"""Principal-radius map of the analytic cap (design check before fitting).

For each chord parameter q: smallest principal radius over the wrap and where,
lip radius, attachment thickness, rule verdict, bending reversals along the
guide (sign changes of the normal curvature in the wrap direction) and
saddle (K < 0) fraction. Run from the project root:
  .venv/bin/python runs/cap-radius/builds/swept55/sharp/map_analytic.py
"""
import pathlib
import sys
import warnings

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import cap_model as cm  # noqa: E402

warnings.simplefilter("ignore")


def main():
    th = np.concatenate([np.linspace(1e-4, np.pi / 2 - 1e-4, 4000), [np.pi / 2]])
    qs = np.concatenate([np.geomspace(1e-3, 0.02, 6), np.linspace(0.03, 0.97, 95), 1 - np.geomspace(0.025, 1e-3, 6)])
    worst = []
    for q in qs:
        J = cm.surface_jet(np.full_like(th, q), th)
        P = cm.principal(J)
        rad = 1 / P["kmax"]
        G = (J["d_theta"] ** 2).sum(-1)
        kn = (J["d_theta_theta"] * P["n"]).sum(-1) / G       # normal curvature along the guide
        sgn = np.sign(np.round(kn / np.abs(kn).max(), 9))
        rev = int(np.sum(sgn[1:] * sgn[:-1] < 0))
        t = float(cm.thickness_at_attachment(q))
        i = int(np.argmin(rad))
        need = 0.4 if t >= 0.8 else t / 2
        ok = rad[i] >= need * (1 - 1e-9)
        worst.append((q, t, rad[i], th[i], rad[-1], rev, float(np.mean(P["K"] < 0)), ok, rad[i] / need))
        print(f"q={q:.4f} t={t:.4f} minR={rad[i]:.4f} at x={np.sin(th[i]):.3f} lipR={rad[-1]:.4f} "
              f"need={need:.4f} ratio={rad[i]/need:.3f} rev={rev} saddle={np.mean(P['K'] < 0):.2f} {'OK' if ok else 'LOW'}")
    w = min(worst, key=lambda r: r[8])
    print("worst ratio", w)


if __name__ == "__main__":
    main()
