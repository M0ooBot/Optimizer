"""Offline 3D review of the Stage 1 winglet: a skin mesh built from the law's sections.

Run from the repository root (same Converter argument as build.py):

    uv run --locked --extra phase0 --extra cad python runs/winglet-b/stage1/viewer.py \
        --converter-src <dir>/src

Writes viewer.html next to this file: one self-contained page, no network needed.
The skin is the law's own sections joined by straight facets. It is a review mesh,
not the CAD surface (Stage 2 builds that).
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import plotly.graph_objects as go

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build

from planeopt.phase0_geometry import load_definition
from planeopt.phase0_winglet import BlendedWinglet, WingletDefinition

STRIPE_PERIOD_DEG = 7.0
LIGHT = np.array([0.35, -0.55, 0.76])  # display axes: outboard, aft, up
WING = "#7fa7c9"
WINGLET = "#e8a868"


def display(points):
    """Study axes (-X outboard, +Y up, +Z forward) to display axes (outboard, aft, up)."""
    p = np.asarray(points, float)
    return np.stack([-p[..., 0], -p[..., 2], p[..., 1]], axis=-1)


def grid(w):
    parent = [w.parent_frame(float(f)) for f in np.linspace(0.0, w.f0, 46)]
    bend = w.L * (1 - np.cos(np.linspace(0, math.pi, 121))) / 2
    leg = np.linspace(w.L, w.s_tip, 41)[1:]
    winglet = [w.frame(float(s)) for s in np.concatenate([bend, leg])]
    keep = np.r_[np.arange(0, 481, 2), 481]  # every other profile point, TE corners kept
    rings = [display(f.section()[keep]) for f in parent + winglet]
    return np.array(rings), len(parent)


def triangles(n_rings, n_ring, skip=()):
    """Quads between consecutive rings; rings in skip are not joined to the next one."""
    i, j, k = [], [], []
    for r in range(n_rings - 1):
        if r in skip:
            continue
        a, b = r * n_ring, (r + 1) * n_ring
        for p in range(n_ring - 1):
            i += [a + p, a + p + 1]
            j += [b + p, b + p]
            k += [a + p + 1, b + p + 1]
    return np.array(i), np.array(j), np.array(k)


def vertex_normals(verts, i, j, k):
    face = np.cross(verts[j] - verts[i], verts[k] - verts[i])
    normals = np.zeros_like(verts)
    for idx in (i, j, k):
        np.add.at(normals, idx, face)
    return normals / np.linalg.norm(normals, axis=1, keepdims=True)


def stripes(normals, light):
    light = light / np.linalg.norm(light)
    angle = np.degrees(np.arccos(np.clip(normals @ light, -1, 1)))
    return np.sin(2 * math.pi * angle / STRIPE_PERIOD_DEG)


def mesh(verts, i, j, k, name, color=None, intensity=None, visible=True, mirror=False):
    v = verts * (-1, 1, 1) if mirror else verts
    common = {"x": v[:, 0], "y": v[:, 1], "z": v[:, 2], "i": i, "j": k if mirror else j, "k": j if mirror else k,
                  "name": name, "hoverinfo": "skip", "showlegend": False, "visible": visible, "flatshading": False}
    if intensity is None:
        return go.Mesh3d(color=color, opacity=1.0,
                         lighting={"ambient": 0.62, "diffuse": 0.6, "specular": 0.25, "roughness": 0.5,
                                   "fresnel": 0.05},
                         lightposition={"x": 1500, "y": -2500, "z": 4000}, **common)
    return go.Mesh3d(intensity=intensity, intensitymode="vertex", showscale=False, cmin=-1, cmax=1,
                     colorscale=[[0, "#111"], [0.5, "#111"], [0.5, "#f4f4f4"], [1, "#f4f4f4"]],
                     lighting={"ambient": 1.0, "diffuse": 0.0, "specular": 0.0}, **common)


def line(points, name, color, width=5, visible=True, dash=None):
    p = np.asarray(points)
    return go.Scatter3d(x=p[:, 0], y=p[:, 1], z=p[:, 2], mode="lines", name=name, visible=visible,
                        line={"color": color, "width": width, "dash": dash}, hoverinfo="name",
                        showlegend=False)


def main():
    args = argparse.ArgumentParser()
    args.add_argument("--converter-src", required=True, type=Path)
    args.add_argument("--open-with", choices=("shaded", "stripes_a", "stripes_b"), default="shaded",
                      help="surface shown when the page opens (for screenshots)")
    args.add_argument("--output", type=Path, default=HERE / "viewer.html")
    opts = args.parse_args()
    sys.path.insert(0, str(opts.converter_src))
    from airfoil_converter import geometry as converter

    definition = WingletDefinition()
    parent_family, target_family = build.families(converter, definition)
    w = BlendedWinglet(load_definition(build.PARENT).wing, definition, parent_family, target_family)
    rings, n_parent = grid(w)
    n_rings, n_ring = rings.shape[:2]
    verts = rings.reshape(-1, 3)
    # The parent's last ring and the winglet's first are the same 540 mm section, kept twice so
    # each body has its own colour. Normals average both sides there, so shading stays seamless.
    i, j, k = triangles(n_rings, n_ring, skip={n_parent - 1})
    bridge = triangles(n_rings, n_ring)
    normals = vertex_normals(verts, *bridge)
    color = np.where(np.arange(n_rings) < n_parent, 0, 1).repeat(n_ring)
    zebra_a = stripes(normals, LIGHT)
    zebra_b = stripes(normals, np.array([-0.8, -0.2, 0.55]))

    traces, groups = [], {}

    def add(group, *items):
        groups.setdefault(group, [])
        for item in items:
            groups[group].append(len(traces))
            traces.append(item)

    shade = {"vertexcolor": np.where(color[:, None] == 0, [WING], [WINGLET]).ravel()}
    shaded = mesh(verts, i, j, k, "skin", color=WING)
    shaded.update(vertexcolor=shade["vertexcolor"], color=None)
    add("shaded", shaded)
    mirrored = mesh(verts, i, j, k, "mirrored skin", color=WING, mirror=True, visible=False)
    mirrored.update(vertexcolor=shade["vertexcolor"], color=None)
    add("mirror", mirrored)
    add("stripes_a", mesh(verts, i, j, k, "stripes A", intensity=zebra_a, visible=False))
    add("stripes_b", mesh(verts, i, j, k, "stripes B", intensity=zebra_b, visible=False))

    le = np.array([w.parent_frame(float(f)).le for f in np.linspace(0, w.f0, 60)]
                  + [w.frame(float(s)).le for s in np.linspace(0, w.s_tip, 160)])
    te = np.array([w.parent_frame(float(f)).te for f in np.linspace(0, w.f0, 60)]
                  + [w.frame(float(s)).te for s in np.linspace(0, w.s_tip, 160)])
    add("edges", line(display(le), "leading edge", "#c0392b"), line(display(te), "trailing edge", "#1e8449"))
    outlines = []
    for s in np.concatenate([np.linspace(0, w.L, 9), np.linspace(w.L, w.s_tip, 4)[1:]]):
        outlines.append(display(w.frame(float(s)).section()))
        outlines.append(np.full((1, 3), np.nan))
    add("sections", line(np.vstack(outlines), "winglet sections", "#1f4e79", width=3))
    attach = display(w.frame(0.0).section())
    add("sections", line(attach, "attachment, 540 mm", "#000", width=6))
    top = float(rings[..., 2].max()) + 15
    aft = rings[n_parent:, :, 1]
    plane = np.array([[600, aft.min() - 5, 30], [600, aft.max() + 5, 30], [600, aft.max() + 5, top],
                      [600, aft.min() - 5, top], [600, aft.min() - 5, 30]])
    add("span", line(plane, "600 mm projected half-span", "#7d3c98", width=4, dash="dash"))

    n = len(traces)

    def toggle(label, group):
        return {"label": label, "method": "restyle",
                    "args": [{"visible": True}, groups[group]], "args2": [{"visible": False}, groups[group]]}

    def surface(label, group):
        vis = [traces[t].visible for t in range(n)]
        for g in ("shaded", "stripes_a", "stripes_b"):
            for t in groups[g]:
                vis[t] = g == group
        return {"label": label, "method": "restyle", "args": [{"visible": [vis[t] for t in groups["shaded"] + groups["stripes_a"] + groups["stripes_b"]]},
                                                         groups["shaded"] + groups["stripes_a"] + groups["stripes_b"]]}

    ortho = {"type": "orthographic"}

    def camera(label, eye, up, center=(0, 0, 0)):
        return {"label": label, "method": "relayout", "args": [{
            "scene.camera": {"eye": dict(zip("xyz", eye)), "up": dict(zip("xyz", up)),
                             "center": dict(zip("xyz", center)), "projection": ortho}}]}

    extent = {
        "whole": {"x": [-10, 620], "y": [-10, 260], "z": [-40, top]},
        "tip": {"x": [500, 620], "y": [50, 190], "z": [30, 135]},
    }

    def zoom(label, key):
        r = extent[key]
        span = [r[a][1] - r[a][0] for a in "xyz"]
        largest = max(span)
        return {"label": label, "method": "relayout", "args": [{
            "scene.xaxis.range": r["x"], "scene.yaxis.range": r["y"], "scene.zaxis.range": r["z"],
            "scene.aspectmode": "manual",
            "scene.aspectratio": {a: s / largest * 2.2 for a, s in zip("xyz", span)}}]}

    menus = [
        {"type": "buttons", "direction": "right", "x": 0.0, "y": 1.11, "xanchor": "left", "buttons": [
            surface("Shaded", "shaded"), surface("Stripes A", "stripes_a"), surface("Stripes B", "stripes_b")]},
        {"type": "buttons", "direction": "right", "x": 0.0, "y": 1.05, "xanchor": "left", "buttons": [
            camera("Oblique", (1.4, -1.6, 0.9), (0, 0, 1)), camera("Front", (0, -2.5, 0), (0, 0, 1)),
            camera("Plan", (0, 0, 2.5), (0, -1, 0)), camera("From behind", (0, 2.5, 0.3), (0, 0, 1)),
            camera("Outboard", (2.5, 0, 0.2), (0, 0, 1)), camera("Below", (0.3, -0.6, -2.4), (0, -1, 0))]},
        {"type": "buttons", "direction": "right", "x": 0.62, "y": 1.05, "xanchor": "left", "buttons": [
            zoom("Whole half-wing", "whole"), zoom("Tip", "tip")]},
        {"type": "buttons", "direction": "right", "x": 0.62, "y": 1.11, "xanchor": "left", "buttons": [
            toggle("Edges", "edges"), toggle("Sections", "sections"), toggle("Span line", "span"),
            toggle("Mirror", "mirror")]},
    ]
    tip = extent["tip"]
    span = [tip[a][1] - tip[a][0] for a in "xyz"]
    largest = max(span)
    fig = go.Figure(traces)
    fig.update_layout(
        title={"text": "Wing B blended winglet · Stage 1 law (review mesh; winglet tip still open)<br>"
                       "<sup>Blue: Wing B · Amber: winglet from 540 mm · Stripes: isophotes, view-independent; "
                       "a kink or step in a stripe marks a break in the surface</sup>",
               "x": 0.01, "y": 0.975, "font": {"size": 16}},
        updatemenus=menus, margin={"l": 0, "r": 0, "t": 120, "b": 0}, paper_bgcolor="#eef2f5",
        scene={"xaxis": {"title": "outboard (mm)", "range": tip["x"]},
               "yaxis": {"title": "aft (mm)", "range": tip["y"]},
               "zaxis": {"title": "up (mm)", "range": tip["z"]},
               "aspectmode": "manual",
               "aspectratio": {a: s / largest * 2.2 for a, s in zip("xyz", span)},
               "camera": {"eye": {"x": 1.4, "y": -1.6, "z": 0.9}, "projection": ortho},
               "bgcolor": "#eef2f5"},
    )
    for group in ("shaded", "stripes_a", "stripes_b"):
        for t in groups[group]:
            traces[t].visible = group == opts.open_with
    fig.data = []
    fig.add_traces(traces)
    out = opts.output
    page = fig.to_html(include_plotlyjs=True, full_html=True, default_height="100vh", default_width="100%",
                       config={"displaylogo": False, "responsive": True})
    style = "<style>html, body { margin: 0; height: 100%; overflow: hidden; background: #eef2f5; }</style>"
    out.write_text(page.replace("<head>", "<head><meta charset=\"utf-8\"><title>Wing B winglet · Stage 1</title>"
                                + style, 1), encoding="utf-8")
    print(out, f"{out.stat().st_size / 1e6:.1f} MB", f"{len(verts)} vertices")


if __name__ == "__main__":
    main()
