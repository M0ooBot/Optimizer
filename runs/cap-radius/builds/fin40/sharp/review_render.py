"""Render the reimported wing.step with CAD-normal zebra stripes and a radius map (PNG).

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/fin40/sharp/review_render.py

Each pixel's colour comes from the CAD surface normal of the STEP B-spline at that
point (dense parameter sampling, z-buffered), so there is no tessellation in the
shading and no edge lines. Camera, zoom and stripe lighting are fixed in VIEWS.
"""
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from surfkit import BSurf, read_step, faces_of, differential, principal, unit, zebra_phase  # noqa: E402
from OCP.BRep import BRep_Tool  # noqa: E402
from OCP.GeomAdaptor import GeomAdaptor_Surface  # noqa: E402
from OCP.GeomAbs import GeomAbs_BSplineSurface  # noqa: E402

OUT = HERE / "review"
WIDTH = np.radians(1.5)     # stripe width in reflected angle (finer than the numerical check's 5 deg)

# name: (view direction (eye -> scene), up, stripe direction angle deg, centre, half-height mm, size px, zmin of skins)
VIEWS = {
    "fin_side_+Y": ((0, -1, 0), (0, 0, 1), 0, (85, 0, 140), 150, (900, 1100), 0),
    "fin_oblique_front": ((1, -1.2, -0.35), (0, 0, 1), 20, (85, 0, 140), 150, (900, 1100), 0),
    "cap_oblique_front": ((1, -1.2, -0.35), (0, 0, 1), 20, (108, 0, 258), 34, (1200, 1000), 200),
    "cap_oblique_rear": ((-1, -1.0, -0.5), (0, 0, 1), 70, (110, 0, 258), 34, (1200, 1000), 200),
    "cap_top_down": ((0.15, -0.35, -1), (1, 0, 0), 45, (108, 0, 262), 30, (1200, 900), 200),
    "cap_te_corner": ((-0.6, -1, -0.25), (0, 0, 1), 30, (148, 0, 246), 7, (1000, 1000), 225),
    "cap_le_corner": ((1, -0.8, -0.3), (0, 0, 1), 30, (67, 0, 245), 7, (1000, 1000), 225),
}


def load():
    shape = read_step(HERE / "wing.step")
    surfs = []
    for f in faces_of(shape):
        if GeomAdaptor_Surface(BRep_Tool.Surface_s(f)).GetType() == GeomAbs_BSplineSurface:
            surfs.append(BSurf(f))
    return surfs


def sample(s, zmin, n_u, n_v, urange=(0.0, 1.0)):
    t = np.linspace(0, 1, n_u)
    t = 0.5 * (t + 0.5 * (1 - np.cos(np.pi * t)))
    us = urange[0] + (urange[1] - urange[0]) * t
    if s.nv == 5:     # parent skin: (s, Z)
        vs = np.linspace(zmin, s.vrange[1], n_v)
    else:             # cap: (u, v)
        us = np.clip(us, 1e-4, 1 - 1e-4)
        vs = np.linspace(-1, 1, n_v)
    out = []
    for i0 in range(0, n_u, 200):
        d = differential(s.jet(us[i0:i0 + 200], vs, grid=True), s.sign)
        out.append((d["P"].reshape(-1, 3), d["n"].reshape(-1, 3)))
    return np.concatenate([o[0] for o in out]), np.concatenate([o[1] for o in out])


def grid_samples(s, zmin, urange, n_u, n_v):
    t = np.linspace(0, 1, n_u)
    t = 0.5 * (t + 0.5 * (1 - np.cos(np.pi * t)))
    us = urange[0] + (urange[1] - urange[0]) * t
    if s.nv == 5:
        vs = np.linspace(zmin, s.vrange[1], n_v)
    else:
        us = np.clip(us, 1e-4, 1 - 1e-4)
        vs = np.linspace(-1, 1, n_v)
    d = differential(s.jet(us, vs, grid=True), s.sign)
    return d["P"], d["n"]


def render(name, surfs):
    """Gouraud-interpolated stripe phase over a dense CAD-normal grid; back faces culled; no edges."""
    from matplotlib.colors import ListedColormap, Normalize
    from matplotlib.tri import Triangulation
    V, up, ang, centre, half, (W, H), zmin = VIEWS[name][:7]
    urange = VIEWS[name][7] if len(VIEWS[name]) > 7 else (0.0, 1.0)
    V = unit(np.array(V, float))
    up = unit(np.array(up, float) - np.dot(up, V) * V)
    right = np.cross(V, up)
    dvec = np.cos(np.radians(ang)) * up + np.sin(np.radians(ang)) * right
    c = np.array(centre, float)
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-half * W / H, half * W / H)
    ax.set_ylim(-half, half)
    ax.set_axis_off()
    for s in surfs:
        n_u, n_v = (1800, 900) if s.nv != 5 else (1200, 600)
        P, n = grid_samples(s, zmin, urange, n_u, n_v)
        x = (P - c) @ right
        y = (P - c) @ up
        ph, _ = zebra_phase(n, V, dvec, WIDTH)
        front = (n @ V) < 0
        nu, nv = x.shape
        idx = np.arange(nu * nv).reshape(nu, nv)
        a, b, cc, dd = idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]
        tris = np.concatenate([np.stack([a, b, cc], -1).reshape(-1, 3), np.stack([a, cc, dd], -1).reshape(-1, 3)])
        fr = front.ravel()
        tris = tris[fr[tris].all(1)]
        xs, ys = x.ravel(), y.ravel()
        # drop degenerate triangles (collapsed cap ends)
        area = np.abs((xs[tris[:, 1]] - xs[tris[:, 0]]) * (ys[tris[:, 2]] - ys[tris[:, 0]])
                      - (xs[tris[:, 2]] - xs[tris[:, 0]]) * (ys[tris[:, 1]] - ys[tris[:, 0]]))
        tris = tris[area > 1e-14]
        if len(tris) == 0:
            continue
        lo, hi = np.floor(ph.min()), np.ceil(ph.max())
        nb = int(hi - lo)
        cols = [(0.93, 0.93, 0.93) if (int(lo) + k) % 2 == 0 else (0.10, 0.10, 0.10) for k in range(max(nb, 1))]
        ax.tripcolor(Triangulation(xs, ys, tris), ph.ravel(), cmap=ListedColormap(cols),
                     norm=Normalize(lo, lo + max(nb, 1)), shading="gouraud", rasterized=True)
    fig.savefig(OUT / f"zebra_{name}.png")
    plt.close(fig)


def radius_map(surfs):
    cap = [s for s in surfs if s.nv != 5][0]
    us = np.linspace(1e-3, 1 - 1e-3, 900)
    vs = np.linspace(0, 1, 500)          # +Y half of the cap (the -Y half is its mirror)
    d = differential(cap.jet(us, vs, grid=True), cap.sign)
    k1, k2, _, _ = principal(d)
    R = 1 / k1
    X, Z = d["P"][..., 0], d["P"][..., 2]
    fig, ax = plt.subplots(figsize=(9, 5.2), dpi=150)
    from matplotlib.colors import LogNorm
    pc = ax.pcolormesh(X, Z, np.clip(R, 0.4, 100.0), cmap="Blues_r", shading="gouraud", norm=LogNorm(0.4, 100.0))
    cs = ax.contour(X, Z, R, levels=[0.6, 0.8, 1.0, 2.0, 5.0], colors="#444444", linewidths=0.6)
    ax.clabel(cs, fmt="%.1f mm", fontsize=7)
    i = np.unravel_index(np.argmin(R), R.shape)
    ax.plot(X[i], Z[i], marker="o", ms=8, mfc="none", mec="#c0392b", mew=1.5)
    ax.annotate(f"smallest {R[i]:.3f} mm", (X[i], Z[i]), xytext=(-150, -25), textcoords="offset points",
                fontsize=8, color="#222222", arrowprops=dict(arrowstyle="-", color="#888888", lw=0.6))
    ax.set_aspect("equal")
    ax.set_xlabel("X (mm, aft)")
    ax.set_ylabel("Z (mm, up)")
    ax.set_title("Cap, +Y half: smaller principal radius 1/k1 (side view, log colour scale)", fontsize=9)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    cb = fig.colorbar(pc, ax=ax, shrink=0.8)
    cb.set_label("radius (mm)")
    fig.tight_layout()
    fig.savefig(OUT / "radius_map_side.png")


def main():
    OUT.mkdir(exist_ok=True)
    surfs = load()
    for name in VIEWS:
        render(name, surfs)
    radius_map(surfs)
    print("written", sorted(p.name for p in OUT.iterdir()))


if __name__ == "__main__":
    main()
