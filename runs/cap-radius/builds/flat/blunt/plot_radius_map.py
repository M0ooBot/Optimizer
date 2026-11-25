"""Plot the curvature-radius map of the right cap from the reimported STEP (radius / printable requirement).

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/flat/blunt/plot_radius_map.py
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import verify_cap as V  # noqa: E402
from OCP.TopAbs import TopAbs_FACE
from OCP.TopoDS import TopoDS
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_Plane


def main():
    shape = V.read(V.STEP)
    cap = None
    for f in V.explore(shape, TopAbs_FACE, TopoDS.Face_s):
        ad = BRepAdaptor_Surface(f)
        if ad.GetType() != GeomAbs_Plane:
            S = V.Surf(ad.BSpline(), f)
            if abs(S.vrange[0] + 1) < 1e-12 and S.P[..., 1].mean() > 0:
                cap = S
    s = np.unique(np.r_[np.linspace(1e-4, 1, 900), 1 - np.geomspace(1e-7, 0.05, 200)])
    w = np.linspace(0, 1, 400)  # upper half (symmetric)
    d = cap.grid(s, w)
    g = V.geometry(d)
    r = 1 / g['kmax']
    req, t = V.required_radius(d['P'][..., 0])
    ratio = r / req
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), dpi=110)
    X, Yp = d['P'][..., 0], d['P'][..., 1] - 350
    for ax, sel, ttl in ((axs[0], slice(None), 'whole cap (upper half), plan view'), (axs[1], None, 'aft end (s >= 0.92), X on a log scale; X = 120 drawn at 1e-9')):
        if sel is None:
            m = s >= 1 - 0.08
            Xs = np.log10(np.maximum(120 - X[m], 1e-9)); Ys = Yp[m]; R_ = ratio[m]
            pc = ax.pcolormesh(Ys, Xs, np.log10(R_), shading='auto', cmap='viridis', vmin=0, vmax=2)
            ax.set_ylabel('log10(120 - X)  [mm]'); ax.invert_yaxis()
        else:
            pc = ax.pcolormesh(Yp, X, np.log10(ratio), shading='auto', cmap='viridis', vmin=0, vmax=2)
            ax.set_ylabel('X [mm]'); ax.invert_yaxis()
        ax.set_xlabel('Y - 350 [mm]'); ax.set_title(ttl, fontsize=10)
    i = np.unravel_index(np.argmin(ratio), ratio.shape)
    axs[1].plot(Yp[i], np.log10(max(120 - X[i], 1e-9)), 'r+', ms=14)
    j = np.unravel_index(np.argmin(np.where(t >= 0.8, r, np.inf)), r.shape)
    axs[0].plot(Yp[j], X[j], 'r+', ms=14)
    axs[0].annotate(f'min radius where t >= 0.8 mm: {r[j]:.3f} mm', (Yp[j], X[j]), (2, 95), color='r', fontsize=9, arrowprops=dict(arrowstyle='-', color='r'))
    axs[1].annotate(f'min ratio {ratio[i]:.2f}: r = {r[i]:.3f} mm, t = {t[i]:.2f} mm', (Yp[i], np.log10(max(120 - X[i], 1e-9))), (5, -2.5), color='r', fontsize=9, arrowprops=dict(arrowstyle='-', color='r'))
    cb = fig.colorbar(pc, ax=axs, shrink=0.9); cb.set_label('log10( principal radius / required radius )   (0 = at the limit)')
    fig.suptitle('wing.step right cap: smallest principal radius over the printable requirement (0.4 mm, or t/2 where t < 0.8 mm)', fontsize=10)
    p = os.path.join(HERE, 'images', 'radius_ratio_map.png'); fig.savefig(p); print(p)


if __name__ == '__main__':
    main()
