"""Shaded close-ups of the trailing-edge corner of the right tip from the reimported wing.step.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/flat/blunt/render_te_corner.py

Structured triangles over each face's parameter window, CAD normals at the
vertices (exact from the STEP B-spline data), painter's order, Gouraud shading.
The flat trailing-edge face and the planar closure are drawn too.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import verify_cap as V  # noqa: E402
from OCP.TopAbs import TopAbs_FACE
from OCP.TopoDS import TopoDS
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_Plane


def grid_faces(s_lo, y_lo, n=260):
    shape = V.read(V.STEP)
    faces = []
    cap = None
    for f in V.explore(shape, TopAbs_FACE, TopoDS.Face_s):
        ad = BRepAdaptor_Surface(f)
        if ad.GetType() == GeomAbs_Plane:
            continue
        S = V.Surf(ad.BSpline(), f)
        if S.P[..., 1].mean() < 0:
            continue
        if S.ku == 5 and abs(S.vrange[0] + 1) < 1e-12:
            cap = S
            u = 1 - np.geomspace(1e-7, 1 - s_lo, n)[::-1]; u = np.r_[u, 1.0]
            e = np.linspace(0, 1, n); v = np.unique(np.r_[np.sqrt(1 - e), -np.sqrt(1 - e), np.linspace(-1, 1, n)])
        elif S.ku == 5:
            u = 1 - np.geomspace(1e-7, 1 - s_lo, n)[::-1]; u = np.r_[u, 1.0]; v = np.linspace(y_lo, 350, n)
        else:
            u = np.linspace(y_lo, 350, n); v = np.linspace(0, 1, 40)
        d = S.grid(u, v)
        faces.append((d['P'], np.cross(d['Pu'], d['Pv'])))
    # planar closure: (X, Y, Z) = (120, 350 + lam * R(1) E(w), h(1) S(w)) spans the region between the
    # trailing-edge face's end segment and the cap's s = 1 iso-curve (the closure boundary)
    w = np.linspace(-1, 1, n); lam = np.linspace(0, 1, 60)
    d = cap.grid(np.array([1.0]), w)
    Pc = d['P'][0]
    P = np.stack([np.full((60, n), 120.0), 350 + lam[:, None] * (Pc[None, :, 1] - 350), np.broadcast_to(Pc[None, :, 2], (60, n))], -1)
    faces.append((P, np.broadcast_to(np.array([1.0, 0, 0]), P.shape)))
    return faces


def draw(faces, view, up, center, width, path, title, light, px=(1100, 760), zebra_deg=None):
    view = np.asarray(view, float); view /= np.linalg.norm(view)
    up = np.asarray(up, float); right = np.cross(up, view); right /= np.linalg.norm(right); upv = np.cross(view, right)
    L = np.asarray(light, float); L /= np.linalg.norm(L)
    xs, ys, cs, tris, depth = [], [], [], [], []
    off = 0
    for P, N in faces:
        nu, nv = P.shape[:2]
        N = N / np.maximum(np.linalg.norm(N, axis=-1, keepdims=True), 1e-300)
        # orient toward the viewer-independent outside: +X for the TE region faces, else by Z sign / Y sign
        q = (P - center)
        x = q @ right; y = q @ upv; z = q @ view
        if zebra_deg is None:
            c = 0.18 + 0.78 * np.abs(N @ L)
        else:  # continuous reflection phase; striped per triangle below
            Nn = N * np.sign((N @ view))[..., None]
            r = -view - 2 * (Nn @ -view)[..., None] * Nn
            c = np.arcsin(np.clip(r @ L, -1, 1)) / np.radians(zebra_deg)
        idx = np.arange(nu * nv).reshape(nu, nv) + off
        a, b, cc, d = idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]
        t = np.concatenate([np.stack([a, b, cc], -1).reshape(-1, 3), np.stack([a, cc, d], -1).reshape(-1, 3)])
        xs.append(x.ravel()); ys.append(y.ravel()); cs.append(c.ravel()); depth.append(z.ravel()); tris.append(t)
        off += nu * nv
    x = np.concatenate(xs); y = np.concatenate(ys); c = np.concatenate(cs); z = np.concatenate(depth); t = np.concatenate(tris)
    h = width * px[1] / px[0]
    vis = (np.abs(x[t]).min(1) < width / 2 * 1.05) & (np.abs(y[t]).min(1) < h / 2 * 1.05)
    t = t[vis]
    t = t[np.argsort(z[t].mean(1))]
    fig = plt.figure(figsize=(px[0] / 100, px[1] / 100 + 0.35), dpi=100)
    ax = fig.add_axes([0, 0, 1, px[1] / (px[1] + 35)])
    ax.set_facecolor((0.55, 0.62, 0.7))
    if zebra_deg is None:
        ax.tripcolor(Triangulation(x, y, t), c, shading='gouraud', cmap='gray', vmin=0, vmax=1)
    else:
        fc = 0.12 + 0.8 * (np.floor(c[t].mean(1) + 0.25) % 2)
        ax.tripcolor(Triangulation(x, y, t), facecolors=fc, cmap='gray', vmin=0, vmax=1, edgecolors='none')
    ax.set_xlim(-width / 2, width / 2); ax.set_ylim(-h / 2, h / 2); ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    fig.text(0.01, 1 - 0.3 / (px[1] / 100 + 0.35), title, fontsize=10, va='center')
    fig.savefig(path); plt.close(fig)


def main():
    out = os.path.join(HERE, 'images'); os.makedirs(out, exist_ok=True)
    f = grid_faces(0.9, 330.0)
    draw(f, (0.8, 0.45, 0.4), (0, 0, 1), np.array([119.0, 352.0, 0]), 9.0, os.path.join(out, 'shaded_te_corner_0.png'),
         'right TE corner from aft-outboard-above, 9 mm wide (shaded, no edge lines)', light=(0.2, 0.6, 0.8))
    f = grid_faces(0.995, 349.0, n=420)
    draw(f, (1.0, 0.25, 0.15), (0, 0, 1), np.array([120.0, 350.2, 0]), 1.4, os.path.join(out, 'zebra_te_corner_1.png'),
         'right TE corner from aft, 1.4 mm wide, fine zebra (0.25 deg): flat TE face + planar half-round closure (left), cap (right); creases show as stripe jumps',
         light=(0.0, 0.6, 0.8), zebra_deg=0.25)
    print('ok')


if __name__ == '__main__':
    main()
