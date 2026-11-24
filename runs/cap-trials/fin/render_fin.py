"""Render overview and fixed-view zebra evidence directly from the reimported STEP."""

from itertools import product
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.ticker import MultipleLocator
import numpy as np
from OCP.gp import gp_Pnt, gp_Vec

from build_fin import CUT, HEIGHT, HERE, OUTPUT
from check_fin import at, patches

INK = "#172232"
GRID = "#ced6de"
BLUE = "#165773"
GOLD = "#b85b27"


def grid(face, side, resolution):
    parameters = np.linspace(0, 1, resolution)
    points = np.empty((resolution, resolution, 3))
    normals = np.empty_like(points)
    for i, u in enumerate(parameters):
        for j, v in enumerate(parameters):
            p, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
            face.D1(u, v, p, du, dv)
            points[i, j] = p.Coord()
            n = np.cross(du.Coord(), dv.Coord())
            normals[i, j] = (side*n/np.linalg.norm(n) if np.linalg.norm(n) > 1e-12
                             else [0, 0, 1])
    vertices = [points[:-1, :-1], points[1:, :-1], points[1:, 1:], points[:-1, 1:]]
    normals = (normals[:-1, :-1] + normals[1:, :-1] + normals[1:, 1:] + normals[:-1, 1:])
    normals /= np.linalg.norm(normals, axis=-1, keepdims=True)
    return np.stack(vertices, axis=-2).reshape(-1, 4, 3), normals.reshape(-1, 3)


def pose(az, el):
    a, e = np.radians([az, el])
    eye = np.array([np.cos(e)*np.sin(a), np.cos(e)*np.cos(a), np.sin(e)])
    right = np.array([np.cos(a), -np.sin(a), 0])
    return eye, right, np.cross(right, eye)


def outline(ax, faces):
    for region, u in product(("parent", "cap"), (0, 1)):
        curve = np.array([at(faces[region, 1], u, v) for v in np.linspace(0, 1, 241)])
        ax.plot(curve[:, 0], curve[:, 2], lw=2.4, color=BLUE if region == "parent" else GOLD)
    ax.plot([0, 160], [0, 0], lw=2.4, color=BLUE)
    ax.plot([65, 150], [CUT, CUT], color=GOLD, lw=1, linestyle=(0, (4, 3)))
    ax.set_title("Fin outline and cap station", loc="left", color=INK)
    ax.set(xlim=(-10, 170), ylim=(-10, HEIGHT+10), aspect="equal")
    ax.xaxis.set_major_locator(MultipleLocator(40))
    ax.yaxis.set_major_locator(MultipleLocator(40))
    ax.set_xlabel("Aft (mm)")
    ax.set_ylabel("Height (mm)")
    ax.grid(color=GRID, linewidth=0.5)


def reflection(ax, mesh, eye_pose, angle, view_label):
    verts, normals = mesh
    eye, right, up = pose(*eye_pose)
    theta = np.radians(angle)
    axis = np.cos(theta)*right + np.sin(theta)*up
    on_screen = normals @ eye > 0
    verts, normals = verts[on_screen], normals[on_screen]
    near_to_far = np.argsort(verts.mean(axis=1) @ eye)
    verts, normals = verts[near_to_far], normals[near_to_far]
    coords = verts @ np.column_stack((right, up))
    phase = 3*(2*(normals @ eye)*(normals @ axis) - eye @ axis)
    gray = 0.18 + 0.74*np.clip(0.5 + 4*np.sin(2*np.pi*phase), 0, 1)
    colors = np.tile(gray[:, None], 3)
    collection = PolyCollection(coords, facecolors=colors, edgecolors="none", antialiased=False)
    ax.add_collection(collection)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.set_title(view_label, loc="left", color=INK)
    ax.axis("off")


def render_overview(faces):
    resolution = 321
    parts = [grid(face, side, resolution) for (_, side), face in faces.items()]
    mesh = tuple(np.concatenate(items) for items in zip(*parts))
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 6.1))
    fig.subplots_adjust(left=0.086, right=0.99, bottom=0.12, top=0.88, wspace=0.20)
    outline(axes[0], faces)
    reflection(axes[1], mesh, (15, 10), 51, "Zebra, side view")
    reflection(axes[2], mesh, (55, 25), 120, "Zebra, oblique view")
    fig.suptitle("Vertical fin with rounded cap", x=0.025, ha="left", fontsize=19, color=INK)
    fig.savefig(HERE / "fin-overview.png", dpi=180)
    plt.close(fig)


def render_sections():
    evidence = json.loads((HERE / "readback.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.65), layout="constrained")
    indices = [0, 2, 4, 7, 9]
    colors = plt.colormaps["viridis"](np.linspace(0.12, 0.83, len(indices)))
    for index, color in zip(indices, colors):
        section = evidence["sections"][index]
        z, y, curvature = np.array(section["trace_z_y_curvature"]).T
        label = f"{section['cut_chord_fraction_from_le']:.0%}"
        axes[0].plot(np.r_[y, -y[::-1]], np.r_[z, z[::-1]], color=color, label=label, lw=1.7)
        axes[1].plot(z, curvature, color=color, label=label, lw=1.7)
    axes[0].axhline(CUT, color=GOLD, lw=1, linestyle="--")
    axes[1].axvline(CUT, color=GOLD, lw=1, linestyle="--")
    axes[0].set(xlabel="Transverse position (mm)", ylabel="Height (mm)",
                title="Wrap sections through both skins", ylim=(220, HEIGHT+2))
    axes[1].set(xlabel="Height (mm)", ylabel="Section curvature (1/mm)",
                title="Transition into the lip", ylim=(-0.035, None))
    axes[1].legend(loc="upper left", title="Chord from LE", frameon=False)
    for ax in axes:
        ax.grid(color=GRID, linewidth=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(HERE / "fin-sections.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    render_overview(patches(OUTPUT)[0])
    render_sections()
