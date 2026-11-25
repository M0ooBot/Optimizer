"""Render zebra-stripe views of the reimported wing.step (tip region) from CAD normals.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/wingb/sharp/render_zebra.py

Orthographic cameras held fixed; stripes come from the reflected view ray
(same phase definition as verify_cap.zebra_join); no CAD edge lines drawn.
Writes zebra_*.png and radius_map.png next to this script.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import capgeom as g  # noqa: E402
import surfeval as se  # noqa: E402
import verify_cap as vc  # noqa: E402

NSTRIPES = 24


def sample(surface, u, v):
    U, V = np.meshgrid(u, v, indexing="ij")
    D = se.d2(surface, U, V)
    P, Pu, Pv = D[..., 0, :], D[..., 1, :], D[..., 2, :]
    n = np.cross(Pu, Pv)
    n /= np.maximum(np.linalg.norm(n, axis=-1, keepdims=True), 1e-300)
    return P, n


def main():
    _, faces, _ = vc.read_step(vc.STEP)
    cap_f, skin_f, _ = vc.classify(faces)
    par = g.parent()
    v_att = par.SURFACE_DATA["v_attachment"]
    v0 = float(par.v_at_span(470.0))
    skinP, skinN = sample(skin_f["surf"], np.linspace(0, 1, 1601), np.linspace(v0, v_att, 141))
    capP, capN = sample(cap_f["surf"], np.r_[np.linspace(0, 1, 801)], np.linspace(0, 2, 801))
    # orient normals outward: skin normal at the upper mid-chord should point +Y; cap lip normal at the tip -X
    if skinN[1200, -1, 1] < 0:
        skinN = -skinN
    if capN[400, 400, 0] > 0:
        capN = -capN
    views = {
        "top": (np.array([0.0, -1.0, 0.0]), np.array([1.0, 0.0, 0.0]), [1.0, 0.0, 0.3]),
        "oblique": (np.array([0.55, -0.6, -0.58]), np.array([0.0, 1.0, 0.0]), [1.0, 0.2, 0.0]),
        "bottom": (np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0]), [1.0, 0.0, -0.3]),
        "tip-on": (np.array([1.0, -0.25, 0.0]), np.array([0.0, 1.0, 0.0]), [0.0, 1.0, 0.2]),
    }
    for name, (d, up, axis) in views.items():
        d = d / np.linalg.norm(d)
        a = np.array(axis, float)
        a /= np.linalg.norm(a)
        e1 = np.cross(up, d)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(d, e1)
        fig, ax = plt.subplots(figsize=(12, 7), dpi=130)
        for P, n in ((skinP, skinN), (capP, capN)):
            dn = np.sum(n * d, -1)
            r = d - 2 * dn[..., None] * n
            phi = NSTRIPES * np.arccos(np.clip(np.sum(r * a, -1), -1, 1)) / np.pi
            col = np.floor(phi) % 2
            shade = 0.25 + 0.75 * col
            mask = dn > -0.02  # back-facing
            x, y = P @ e1, P @ e2
            ax.pcolormesh(x, y, np.ma.masked_where(mask, shade), cmap="gray", vmin=0, vmax=1, shading="gouraud",
                          rasterized=True)
        ax.set_aspect("equal")
        ax.set_title(f"zebra {name}: view {np.round(d, 2).tolist()}, stripe axis {np.round(a, 2).tolist()}, "
                     f"{NSTRIPES} stripes / 180 deg (CAD normals, reimported STEP)", fontsize=9)
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(g.BUILD_DIR / f"zebra_{name}.png")
        plt.close(fig)
    # plan-view map of the smallest principal radius against the printable rule
    U, V, D, C = vc.tabulate(cap_f["surf"], 401, 801)
    P = C["P"]
    xi = g.xi_of(-P[..., 0], P[..., 2])
    R = 1 / C["kmax"]
    ratio = R / vc.required_radius(xi)
    fig, axs = plt.subplots(1, 2, figsize=(13, 6), dpi=130)
    for ax, side, title in ((axs[0], C["n"][..., 1] >= 0, "upper half (sigma <= 1)"),
                            (axs[1], C["n"][..., 1] < 0, "lower half (sigma >= 1)")):
        half = V <= 1 if "upper" in title else V >= 1
        m = np.ma.masked_where(~half, np.log10(np.clip(ratio, 0.5, 50)))
        pc = ax.pcolormesh(-P[..., 0], -P[..., 2], m, cmap="viridis", vmin=np.log10(0.5), vmax=np.log10(50),
                           shading="gouraud", rasterized=True)
        ax.contour(-P[..., 0], -P[..., 2], np.ma.masked_where(~half, ratio), levels=[1.0], colors="r", linewidths=0.8)
        O = par.tip_outline(np.linspace(0, 1, 400))
        ax.plot(-O[:, 0], -O[:, 1], "k-", lw=0.6)
        ax.set_aspect("equal")
        ax.set_xlabel("projected span S (mm)")
        ax.set_ylabel("aft (mm)")
        ax.invert_yaxis()
        ax.set_title(title, fontsize=9)
    cb = fig.colorbar(pc, ax=axs, shrink=0.8)
    cb.set_label("log10(smallest principal radius / required radius); red: ratio = 1")
    fig.savefig(g.BUILD_DIR / "radius_map.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
