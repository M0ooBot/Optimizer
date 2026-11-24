"""Offline 3D review of the Stage 1 winglet variants: a skin mesh built from each law's sections.

Run from the repository root (same Converter argument as build.py):

    uv run --locked --extra phase0 --extra cad python runs/winglet-b/stage1/viewer.py \
        --converter-src <dir>/src

Writes viewer.html next to this file: one self-contained page, no network needed. Each skin
is the law's own sections joined by flat facets. It is a review mesh, not the CAD surface
(Stage 2 builds that).
"""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import plotly.graph_objects as go

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build

STRIPE_PERIOD_DEG = 7.0
LIGHT = np.array([0.35, -0.55, 0.76])  # display axes: outboard, aft, up
WING = "#7fa7c9"
WINGLET = "#e8a868"
PROFILE_STEP = 3  # every third profile point; TE corners kept


def display(points):
    """Study axes (-X outboard, +Y up, +Z forward) to display axes (outboard, aft, up)."""
    p = np.asarray(points, float)
    return np.stack([-p[..., 0], -p[..., 2], p[..., 1]], axis=-1)


def rings(w):
    parent = [w.parent_frame(float(f)) for f in np.linspace(0.0, w.f0, 31)]
    winglet = [w.frame(float(s)) for s in w.stations(101)]
    n = len(parent[0].section())  # closed loop: upper TE ... lower TE, back to upper TE
    keep = np.unique(np.r_[np.arange(0, n - 1, PROFILE_STEP), n - 2, n - 1])
    return np.array([display(f.section()[keep]) for f in parent + winglet]), len(parent)


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
    return np.array(i, np.int32), np.array(j, np.int32), np.array(k, np.int32)


def vertex_normals(verts, i, j, k):
    face = np.cross(verts[j] - verts[i], verts[k] - verts[i])
    normals = np.zeros_like(verts)
    for idx in (i, j, k):
        np.add.at(normals, idx, face)
    return normals / np.linalg.norm(normals, axis=1, keepdims=True)


def variant_traces(name, w):
    ring_pts, n_parent = rings(w)
    n_rings, n_ring = ring_pts.shape[:2]
    verts = ring_pts.reshape(-1, 3)
    # The parent's last ring and the winglet's first are the same section, kept twice so each
    # body has its own colour; normals average both sides there, so shading stays seamless.
    i, j, k = triangles(n_rings, n_ring, skip={n_parent - 1})
    normals = vertex_normals(verts, *triangles(n_rings, n_ring))
    angle = np.degrees(np.arccos(np.clip(normals @ (LIGHT / np.linalg.norm(LIGHT)), -1, 1)))
    xyz = {"x": verts[:, 0].astype(np.float32), "y": verts[:, 1].astype(np.float32),
           "z": verts[:, 2].astype(np.float32), "i": i, "j": j, "k": k}
    colours = np.where(np.arange(n_rings).repeat(n_ring) < n_parent, WING, WINGLET)
    shaded = go.Mesh3d(**xyz, vertexcolor=colours, hoverinfo="skip", showlegend=False, visible=False,
                       lighting={"ambient": 0.62, "diffuse": 0.6, "specular": 0.25, "roughness": 0.5,
                                 "fresnel": 0.05},
                       lightposition={"x": 1500, "y": -2500, "z": 4000},
                       meta={"variant": name, "kind": "shaded"})
    stripes = go.Mesh3d(**xyz, intensity=np.sin(2 * math.pi * angle / STRIPE_PERIOD_DEG).astype(np.float32),
                        intensitymode="vertex", cmin=-1, cmax=1, showscale=False,
                        colorscale=[[0, "#111"], [0.5, "#111"], [0.5, "#f4f4f4"], [1, "#f4f4f4"]],
                        lighting={"ambient": 1.0, "diffuse": 0.0, "specular": 0.0},
                        hoverinfo="skip", showlegend=False, visible=False,
                        meta={"variant": name, "kind": "stripes"})

    s = np.linspace(0, w.s_tip, 160)
    fr = np.linspace(0, w.f0, 60)
    le = np.array([w.parent_frame(float(f)).le for f in fr] + [w.frame(float(t)).le for t in s])
    te = np.array([w.parent_frame(float(f)).te for f in fr] + [w.frame(float(t)).te for t in s])
    edges = np.vstack([display(le), [[np.nan] * 3], display(te)]).astype(np.float32)
    edge_trace = go.Scatter3d(x=edges[:, 0], y=edges[:, 1], z=edges[:, 2], mode="lines",
                              line={"color": "#b03a2e", "width": 5}, hoverinfo="skip", showlegend=False,
                              visible=False, meta={"variant": name, "kind": "edges"})
    outlines = []
    for t in np.linspace(0, w.s_tip, 12):
        outlines += [display(w.frame(float(t)).section()), [[np.nan] * 3]]
    outlines = np.vstack(outlines).astype(np.float32)
    section_trace = go.Scatter3d(x=outlines[:, 0], y=outlines[:, 1], z=outlines[:, 2], mode="lines",
                                 line={"color": "#1f4e79", "width": 3}, hoverinfo="skip", showlegend=False,
                                 visible=False, meta={"variant": name, "kind": "sections"})
    return [shaded, stripes, edge_trace, section_trace], ring_pts


PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Wing B winglet · Stage 1 bend variants</title>
<style>
:root { font-family: system-ui, sans-serif; color: #223646; background: #eef2f5; }
* { box-sizing: border-box; }
html, body { margin: 0; height: 100%; }
body { display: grid; grid-template-columns: 300px 1fr; height: 100vh; }
aside { padding: 16px 16px 12px; overflow-y: auto; border-right: 1px solid #d3dde5; background: #f7f9fb; }
h1 { font-size: 17px; margin: 0 0 4px; letter-spacing: -.01em; }
.note { font-size: 12px; color: #52687a; line-height: 1.45; margin: 0 0 12px; }
h2 { font-size: 11px; text-transform: uppercase; letter-spacing: .06em; color: #5d7385; margin: 14px 0 6px; }
.group { display: flex; flex-wrap: wrap; gap: 5px; }
.group.column { flex-direction: column; }
button { font: inherit; font-size: 12.5px; color: #233749; background: #fff; border: 1px solid #c3d1dc;
         border-radius: 6px; padding: 6px 9px; cursor: pointer; text-align: left; }
button:hover { border-color: #5a7489; }
button[aria-pressed="true"] { background: #dfeaf3; border-color: #4f6f88; color: #10314a; }
button small { display: block; color: #5c7183; font-size: 11px; margin-top: 2px; }
dl { display: grid; grid-template-columns: auto 1fr; gap: 3px 10px; font-size: 12px; margin: 0;
     font-variant-numeric: tabular-nums; }
dt { color: #5d7385; } dd { margin: 0; text-align: right; }
#plot { width: 100%; height: 100vh; }
.legend { font-size: 12px; line-height: 1.6; }
.swatch { display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 6px; vertical-align: -1px; }
</style></head>
<body>
<aside>
  <h1>Wing B blended winglet</h1>
  <p class="note">Stage 1 law: review mesh from the generated sections. The winglet tip is still open (Stage 3 caps it).</p>
  <h2>Bend distribution</h2>
  <div class="group column" id="variants"></div>
  <h2>Selected variant</h2>
  <dl id="facts"></dl>
  <h2>Surface</h2>
  <div class="group" id="modes">
    <button data-mode="shaded" aria-pressed="true">Shaded</button>
    <button data-mode="stripes" aria-pressed="false">Stripes</button>
  </div>
  <p class="note" style="margin-top:6px">Stripes are isophotes, fixed to the surface. A kink or a jump in a stripe marks a break in the skin.</p>
  <h2>Show</h2>
  <div class="group" id="layers">
    <button data-layer="edges" aria-pressed="true">LE and TE</button>
    <button data-layer="sections" aria-pressed="true">Sections</button>
    <button data-layer="span" aria-pressed="true">600 mm span</button>
  </div>
  <h2>Camera</h2>
  <div class="group" id="cameras"></div>
  <h2>Frame</h2>
  <div class="group" id="frames">
    <button data-frame="tip" aria-pressed="true">Tip</button>
    <button data-frame="whole" aria-pressed="false">Whole half-wing</button>
  </div>
  <h2>Colours</h2>
  <div class="legend"><span class="swatch" style="background:WING_COLOUR"></span>Wing B<br>
  <span class="swatch" style="background:WINGLET_COLOUR"></span>Winglet<br>
  <span class="swatch" style="background:#b03a2e"></span>Leading and trailing edges<br>
  <span class="swatch" style="background:#7d3c98"></span>600 mm projected half-span</div>
</aside>
PLOT
<script>
const VARIANTS = VARIANTS_JSON;
const FRAMES = FRAMES_JSON;
const CAMERAS = {
  "Oblique": {eye: {x: 1.4, y: -1.6, z: 0.9}, up: {x: 0, y: 0, z: 1}},
  "Front": {eye: {x: 0, y: -2.5, z: 0}, up: {x: 0, y: 0, z: 1}},
  "Plan": {eye: {x: 0, y: 0, z: 2.5}, up: {x: 0, y: -1, z: 0}},
  "From behind": {eye: {x: 0, y: 2.5, z: 0.25}, up: {x: 0, y: 0, z: 1}},
  "From outboard": {eye: {x: 2.5, y: 0, z: 0.2}, up: {x: 0, y: 0, z: 1}},
  "Below": {eye: {x: 0.3, y: -0.6, z: -2.4}, up: {x: 0, y: -1, z: 0}},
};
const gd = document.getElementById("plot");
const state = {variant: "DEFAULT_VARIANT", mode: "shaded", edges: true, sections: true, span: true};

function apply() {
  const visible = gd.data.map(t => {
    const m = t.meta || {};
    if (m.kind === "span") return state.span;
    if (m.variant !== state.variant) return false;
    if (m.kind === "shaded" || m.kind === "stripes") return m.kind === state.mode;
    return state[m.kind];
  });
  Plotly.restyle(gd, {visible});
  const v = VARIANTS[state.variant];
  document.getElementById("facts").innerHTML = v.facts.map(([k, val]) => "<dt>" + k + "</dt><dd>" + val + "</dd>").join("");
  document.querySelectorAll("#variants button").forEach(b => b.setAttribute("aria-pressed", b.dataset.variant === state.variant));
  document.querySelectorAll("#modes button").forEach(b => b.setAttribute("aria-pressed", b.dataset.mode === state.mode));
}
for (const [name, v] of Object.entries(VARIANTS)) {
  const b = document.createElement("button");
  b.dataset.variant = name;
  b.innerHTML = v.label + "<small>" + v.detail + "</small>";
  b.onclick = () => { state.variant = name; apply(); };
  document.getElementById("variants").append(b);
}
document.querySelectorAll("#modes button").forEach(b => b.onclick = () => { state.mode = b.dataset.mode; apply(); });
document.querySelectorAll("#layers button").forEach(b => b.onclick = () => {
  state[b.dataset.layer] = !state[b.dataset.layer];
  b.setAttribute("aria-pressed", state[b.dataset.layer]);
  apply();
});
for (const [name, cam] of Object.entries(CAMERAS)) {
  const b = document.createElement("button");
  b.textContent = name;
  b.onclick = () => Plotly.relayout(gd, {"scene.camera": {...cam, center: {x: 0, y: 0, z: 0}, projection: {type: "orthographic"}}});
  document.getElementById("cameras").append(b);
}
document.querySelectorAll("#frames button").forEach(b => b.onclick = () => {
  const r = FRAMES[b.dataset.frame];
  Plotly.relayout(gd, {"scene.xaxis.range": r.x, "scene.yaxis.range": r.y, "scene.zaxis.range": r.z,
                       "scene.aspectratio": r.aspect});
  document.querySelectorAll("#frames button").forEach(x => x.setAttribute("aria-pressed", x === b));
});
apply();
</script>
</body></html>
"""


def frame_ranges(x, y, z):
    span = [x[1] - x[0], y[1] - y[0], z[1] - z[0]]
    largest = max(span)
    return {"x": x, "y": y, "z": z, "aspect": {a: s / largest * 1.55 for a, s in zip("xyz", span, strict=True)}}


def main():
    args = argparse.ArgumentParser()
    args.add_argument("--converter-src", required=True, type=Path)
    args.add_argument("--variant", default=build.DEFAULT_VARIANT, choices=list(build.VARIANTS))
    args.add_argument("--mode", default="shaded", choices=("shaded", "stripes"))
    args.add_argument("--output", type=Path, default=HERE / "viewer.html")
    opts = args.parse_args()
    sys.path.insert(0, str(opts.converter_src))
    from airfoil_converter import geometry as converter

    summary = json.loads((HERE / "variants" / "summary.json").read_text())
    traces, variants, tops = [], {}, []
    for name, (title, overrides) in build.VARIANTS.items():
        w = build.make(converter, **overrides)
        items, ring_pts = variant_traces(name, w)
        traces += items
        tops.append(float(ring_pts[..., 2].max()))
        m = summary[name]
        label, _, detail = title.partition(": ")
        variants[name] = {"label": label, "detail": detail, "facts": [
            ["Starts at", f"{m['attach_mm']:.1f} mm"],
            ["Height", f"{m['height_mm']:.1f} mm"],
            ["Tightest bend radius", f"{m['tightest_bend_radius_mm']:.1f} mm"],
            ["Bend ramps in / out", f"{overrides['ramp_in']:.0%} / {overrides['ramp_out']:.0%}"],
            ["Tip cant", f"{m['tip_cant_deg']:.1f}°"],
            ["Tip LE sweep", f"{m['tip_le_sweep_deg']:.1f}°"],
            ["Tip TE sweep", (f"{-m['trailing_edge_tip_sweep_deg']:.1f}° forward"
                              if m["trailing_edge_tip_sweep_deg"] < 0
                              else f"{m['trailing_edge_tip_sweep_deg']:.1f}° aft")],
            ["Toe-out", f"{m['toe_out_deg']:.1f}°"],
            ["Tip chord", f"{m['tip_chord_mm']:.1f} mm"],
            ["Smallest LE radius", f"{m['le_radius_min_mm']:.2f} mm"],
            ["Trailing edge", "0.4 mm"],
        ]}
    top = max(tops) + 12
    plane = np.array([[600, 55, 35], [600, 165, 35], [600, 165, top], [600, 55, top], [600, 55, 35]])
    traces.append(go.Scatter3d(x=plane[:, 0], y=plane[:, 1], z=plane[:, 2], mode="lines",
                               line={"color": "#7d3c98", "width": 4, "dash": "dash"}, hoverinfo="skip",
                               showlegend=False, meta={"kind": "span"}))
    frames = {"tip": frame_ranges([495, 620], [50, 190], [30, 135]),
              "whole": frame_ranges([-10, 620], [-10, 260], [-40, top])}
    tip = frames["tip"]
    fig = go.Figure(traces)
    fig.update_layout(
        margin={"l": 0, "r": 0, "t": 0, "b": 0}, paper_bgcolor="#eef2f5",
        scene={"xaxis": {"title": "outboard (mm)", "range": tip["x"]},
               "yaxis": {"title": "aft (mm)", "range": tip["y"]},
               "zaxis": {"title": "up (mm)", "range": tip["z"]},
               "aspectmode": "manual", "aspectratio": tip["aspect"], "bgcolor": "#eef2f5",
               "camera": {"eye": {"x": 1.4, "y": -1.6, "z": 0.9}, "projection": {"type": "orthographic"}}},
    )
    plot = fig.to_html(include_plotlyjs=True, full_html=False, div_id="plot", default_height="100vh",
                       config={"displaylogo": False, "responsive": True})
    page = (PAGE.replace("PLOT", plot).replace("VARIANTS_JSON", json.dumps(variants))
            .replace("FRAMES_JSON", json.dumps(frames)).replace("DEFAULT_VARIANT", opts.variant)
            .replace("WING_COLOUR", WING).replace("WINGLET_COLOUR", WINGLET)
            .replace('mode: "shaded"', f'mode: "{opts.mode}"'))
    opts.output.write_text(page, encoding="utf-8")
    print(opts.output, f"{opts.output.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
