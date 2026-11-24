"""Render mirrored zebra lights using STEP samples and OCCT first derivatives."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from OCP.gp import gp_Pnt, gp_Vec

from check import HERE, load, unit


VIEW = unit(np.array([-0.38, 0., 0.93]))
HORIZONTAL = np.array([0., 1., 0.])
VERTICAL = np.cross(VIEW, HORIZONTAL)
STRIPE = np.cos(np.radians(50))*HORIZONTAL + np.sin(np.radians(50))*VERTICAL
PROJECTION = np.column_stack([HORIZONTAL, VERTICAL])


def sample(patch, u, v, orientation):
    point, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
    patch.D1(float(u), float(v), point, du, dv)
    normal = orientation * unit(np.cross(du.Coord(), dv.Coord()))
    return point.Coord(), normal


def collect(patches):
    polygons, depths, colors = [], [], []
    for (side, skin, cap), patch in patches.items():
        u = np.linspace(0, 1, 501)
        v = np.linspace(0, 1, 241 if cap else 101)
        p = np.array([[patch.Value(float(a), float(b)).Coord() for b in v] for a in u])
        cells = np.stack((p[:-1, :-1], p[1:, :-1], p[1:, 1:], p[:-1, 1:]), axis=2)
        centers, normals = map(np.array, zip(*(
            sample(patch, a, b, side*skin)
            for a in (u[:-1]+u[1:])*0.5 for b in (v[:-1]+v[1:])*0.5)))
        stripe = STRIPE * [1, side, 1]
        phase = 4*(2*(normals @ VIEW)*(normals @ stripe) - VIEW @ stripe)
        mix = 0.5 + 0.5*np.tanh(4*np.cos(2*np.pi*phase))
        dark = np.array([0.08, 0.15, 0.20])
        light = np.array([0.90, 0.94, 0.96])
        front = normals @ VIEW > 0
        polygons.append(cells.reshape(-1, 4, 3)[front] @ PROJECTION)
        depths.append(centers[front] @ VIEW)
        colors.append(dark + mix[front, None]*(light-dark))
    order = np.argsort(np.concatenate(depths))
    return np.concatenate(polygons)[order], np.concatenate(colors)[order]


def plan(ax, patches):
    samples = np.linspace(0, 1, 501)
    parent, cap = patches[1, 1, False], patches[1, 1, True]
    perimeter = [parent.Value(0., t).Coord() for t in samples]
    perimeter.extend(cap.Value(u, 1.).Coord() for u in samples)
    perimeter.extend(parent.Value(1., t).Coord() for t in samples[::-1])
    points = np.array(perimeter)
    for side in (-1, 1):
        ax.fill(side*points[:, 1], -points[:, 0], color="#285275")
        cut = np.array([parent.Value(u, 1.).Coord() for u in (0., 1.)])
        ax.plot(side*cut[:, 1], -cut[:, 0], color="white", lw=1.3, ls="--")
    ax.set_title("Plan view · dashed lines at cap starts", fontsize=12)
    ax.set_aspect("equal")
    ax.set_xlim(-475, 475)
    ax.set_ylim(-760, 35)
    ax.set_xlabel("Spanwise y (mm)")
    ax.set_ylabel("Aft x (mm)")
    ax.set_yticks([-700, -500, -300, -100, 0], [700, 500, 300, 100, 0])


def main():
    _, patches = load()
    polygons, colors = collect(patches)
    figure, axes = plt.subplots(1, 2, figsize=(15, 7), layout="constrained")
    plan(axes[0], patches)
    ax = axes[1]
    ax.add_collection(PolyCollection(polygons, facecolors=colors, edgecolors="none", antialiased=False))
    ax.set_aspect("equal")
    ax.set_xlim(-475, 475)
    ax.set_ylim(-710, 35)
    ax.set_axis_off()
    ax.set_title("Centered zebra view · mirrored lighting", fontsize=12)
    figure.suptitle("55° swept wing with rounded caps", fontsize=18)
    figure.savefig(HERE / "whole-wing-review.png", dpi=180)

    figure, axes = plt.subplots(1, 2, figsize=(12, 6), layout="constrained")
    for side, ax in zip((-1, 1), axes):
        selection = np.any((side * polygons[..., 0] > 365), axis=1)
        ax.add_collection(PolyCollection(polygons[selection], facecolors=colors[selection],
                                          edgecolors="none", antialiased=False))
        ax.set_aspect("equal")
        ax.set_xlim(*(np.array([-465, -368]) if side == -1 else [368, 465]))
        ax.set_ylim(-697, -480)
        ax.set_axis_off()
    figure.suptitle("Cap detail at the same camera and light angle", fontsize=16)
    figure.savefig(HERE / "zebra-detail.png", dpi=190)
    plt.close("all")


if __name__ == "__main__":
    main()
