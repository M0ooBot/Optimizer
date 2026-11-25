"""Offline 3D review of the Stage 1 winglet configurations, grouped by the setting they vary.

Run from the repository root (same Converter argument as build.py), after build.py:

    uv run --locked --extra phase0 --extra cad python runs/winglet-b/stage1/viewer.py \
        --converter-src <dir>/src

Writes viewer.html next to this file: one self-contained page, no network needed. Each skin is
the law's own sections joined by flat facets (a review mesh, not the CAD surface). The page
ships only the section rings of each configuration; it builds the meshes, normals and stripes
in the browser.
"""

import argparse
import base64
import json
import sys
from pathlib import Path

import numpy as np
from plotly.offline import get_plotlyjs

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build

PROFILE_STEP = 3  # every third profile point; TE corners kept
PARENT_RINGS, WINGLET_RINGS = 31, 101


def display(points):
    """Study axes (-X outboard, +Y up, +Z forward) to display axes (outboard, aft, up)."""
    p = np.asarray(points, float)
    return np.stack([-p[..., 0], -p[..., 2], p[..., 1]], axis=-1)


def packed(array):
    return base64.b64encode(np.ascontiguousarray(array, dtype=np.float32).tobytes()).decode("ascii")


def geometry(w):
    parent = [w.parent_frame(float(f)) for f in np.linspace(0.0, w.f0, PARENT_RINGS)]
    winglet = [w.frame(float(s)) for s in w.stations(WINGLET_RINGS)]
    n = len(parent[0].section())
    keep = np.unique(np.r_[np.arange(0, n - 1, PROFILE_STEP), n - 2, n - 1])
    rings = np.array([display(f.section()[keep]) for f in parent + winglet])
    s = np.linspace(0, w.s_tip, 160)
    fr = np.linspace(0, w.f0, 60)
    le = display([w.parent_frame(float(f)).le for f in fr] + [w.frame(float(t)).le for t in s])
    te = display([w.parent_frame(float(f)).te for f in fr] + [w.frame(float(t)).te for t in s])
    sections = np.array([display(w.frame(float(t)).section()) for t in np.linspace(0, w.s_tip, 12)])
    return {"rings": packed(rings), "ring": int(rings.shape[1]), "le": packed(le), "te": packed(te),
            "sections": packed(sections), "sectionPoints": int(sections.shape[1])}


def facts(row):
    return [
        ["Starts at", f"{row['attach_mm']:.1f} mm"],
        ["Height", f"{row['height_mm']:.1f} mm"],
        ["Tip cant", f"{row['tip_cant_deg']:.1f}°"],
        ["Bend", f"{row['bend_fraction']:.0%} of the winglet"
                 + ("" if row["bend_fraction"] >= 1 else ", then straight")],
        ["Tightest bend radius", f"{row['tightest_bend_radius_mm']:.1f} mm"],
        ["LE in planform", f"{row['le_root_angle_deg']:.0f}° → {row['le_tip_angle_deg']:.0f}° aft"],
        ["TE in planform", f"{row['te_root_angle_deg']:.0f}° → {row['te_tip_angle_deg']:.0f}° aft"],
        ["Guides bend aft progressively", "yes" if row["guides_bend_aft_progressively"] else "NO"],
        ["Toe-out", f"{row['toe_out_deg']:.1f}°"],
        ["Tip chord", f"{row['tip_chord_mm']:.1f} mm"],
        ["LE radius, smallest", (f"{row['le_radius_min_mm']:.2f} mm across the edge "
                                f"({row['le_radius_streamwise_min_mm']:.2f} streamwise)")],
        ["Droop below Wing B", f"{row['front_view_sag']['outline_top_below_tangent_mm']:.3f} mm"],
    ]


def main():
    args = argparse.ArgumentParser()
    args.add_argument("--converter-src", required=True, type=Path)
    args.add_argument("--select", default="bend-centred", help="configuration shown on opening")
    args.add_argument("--mode", default="shaded", choices=("shaded", "stripes"))
    args.add_argument("--camera", default="Oblique", help="camera preset on opening (for screenshots)")
    args.add_argument("--all", action="store_true",
                      help="view the whole configuration study (default: the baseline being tuned)")
    args.add_argument("--output", type=Path, default=None)
    opts = args.parse_args()
    source = HERE / ("configs" if opts.all else "tuning")
    if opts.output is None:
        opts.output = HERE / ("viewer.html" if opts.all else "tuning.html")
    if not opts.all and opts.select == "bend-centred":
        opts.select = "baseline"
    sys.path.insert(0, str(opts.converter_src))
    from airfoil_converter import geometry as converter

    summary = json.loads((source / "summary.json").read_text())
    groups = {g: {"title": t, "members": []} for g, t in summary["groups"].items()}
    configs, shapes = {}, {}
    for name, row in summary["configs"].items():
        groups[row["group"]]["members"].append(name)
        entry = {"label": row["label"], "fits": row["fits"], "cant": row.get("cant_deg", 75.0)}
        if not row["fits"]:
            entry["reason"] = row["reason"]
        else:
            origin = row.get("same_as") or name
            if origin not in shapes:
                shapes[origin] = geometry(build.make(converter, **summary["configs"][origin]["overrides"]))
            entry.update(shape=origin, facts=facts(row))
        configs[name] = entry

    page = (TEMPLATE.replace("/*PLOTLY*/", get_plotlyjs())
            .replace("/*DATA*/", json.dumps({"groups": groups, "configs": configs, "shapes": shapes,
                                             "baseline": summary["baseline_label"],
                                             "select": opts.select, "mode": opts.mode,
                                             "camera": opts.camera})))
    opts.output.write_text(page, encoding="utf-8")
    print(opts.output, f"{opts.output.stat().st_size / 1e6:.1f} MB")


TEMPLATE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Wing B winglet · Stage 1 configurations</title>
<style>
:root { font-family: system-ui, sans-serif; color: #223646; background: #eef2f5; }
* { box-sizing: border-box; }
html, body { margin: 0; height: 100%; }
body { display: grid; grid-template-columns: 330px 1fr; height: 100vh; }
aside { padding: 16px 16px 12px; overflow-y: auto; border-right: 1px solid #d3dde5; background: #f7f9fb; }
h1 { font-size: 17px; margin: 0 0 4px; letter-spacing: -.01em; }
.note { font-size: 12px; color: #52687a; line-height: 1.45; margin: 0 0 10px; }
h2 { font-size: 11px; text-transform: uppercase; letter-spacing: .06em; color: #5d7385; margin: 14px 0 6px; }
h3 { font-size: 12px; margin: 10px 0 5px; color: #2c4458; font-weight: 600; }
.group { display: flex; flex-wrap: wrap; gap: 5px; }
button { font: inherit; font-size: 12.5px; color: #233749; background: #fff; border: 1px solid #c3d1dc;
         border-radius: 6px; padding: 5px 8px; cursor: pointer; }
button:hover { border-color: #5a7489; }
button[aria-pressed="true"] { background: #dfeaf3; border-color: #4f6f88; color: #10314a; }
button:disabled { color: #9aa9b5; border-style: dashed; cursor: not-allowed; background: #f3f5f7; }
button.baseline::after { content: " ●"; color: #6d8aa0; font-size: 9px; vertical-align: 1px; }
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
  <p class="note">Stage 1 configurations, one setting at a time around the baseline (●). Review mesh from the generated sections; the winglet tip is still open.</p>
  <p class="note" id="baseline"></p>
  <div id="configs"></div>
  <h2>Selected</h2>
  <dl id="facts"></dl>
  <h2>Surface</h2>
  <div class="group" id="modes">
    <button data-mode="shaded">Shaded</button>
    <button data-mode="stripes">Stripes</button>
  </div>
  <p class="note" style="margin-top:6px">Stripes are isophotes, fixed to the surface: a kink or jump in a stripe marks a break in the skin.</p>
  <h2>Show</h2>
  <div class="group" id="layers">
    <button data-layer="edges">LE and TE</button>
    <button data-layer="sections">Sections</button>
    <button data-layer="span">600 mm span</button>
    <button data-layer="ghost">Baseline ghost</button>
    <button data-layer="mirror">Mirror</button>
  </div>
  <h2>Camera</h2>
  <div class="group" id="cameras"></div>
  <h2>Frame</h2>
  <div class="group" id="frames">
    <button data-frame="tip">Tip</button>
    <button data-frame="whole">Whole half-wing</button>
  </div>
  <h2>Colours</h2>
  <div class="legend"><span class="swatch" style="background:#7fa7c9"></span>Wing B
  &nbsp; <span class="swatch" style="background:#e8a868"></span>Winglet
  &nbsp; <span class="swatch" style="background:#b03a2e"></span>Edges<br>
  <span class="swatch" style="background:#7d3c98"></span>600 mm projected half-span
  &nbsp; <span class="swatch" style="background:#9fb3c4"></span>Baseline ghost</div>
</aside>
<div id="plot"></div>
<script>/*PLOTLY*/</script>
<script>
const DATA = /*DATA*/;
const WING = [0x7f / 255, 0xa7 / 255, 0xc9 / 255], WINGLET = [0xe8 / 255, 0xa8 / 255, 0x68 / 255];
const LIGHT = (() => { const v = [0.35, -0.55, 0.76], n = Math.hypot(...v); return v.map(x => x / n); })();
const PERIOD = 7, PARENT_RINGS = 31;
function planformCamera() {
  const g = (DATA.configs[state.config].cant || 75) * Math.PI / 180;
  return [[2.5 * Math.sin(g), 0, -2.5 * Math.cos(g)], [Math.cos(g), 0, Math.sin(g)]];
}
const CAMERAS = {"Oblique": [[1.4, -1.6, 0.9], [0, 0, 1]], "Front": [[0, -2.5, 0], [0, 0, 1]], "Plan": [[0, 0, 2.5], [0, -1, 0]],
  "From behind": [[0, 2.5, 0.25], [0, 0, 1]], "From outboard": [[2.5, 0, 0.2], [0, 0, 1]], "Below": [[0.3, -0.6, -2.4], [0, -1, 0]]};
const state = {config: DATA.select, mode: DATA.mode, edges: true, sections: true, span: true, ghost: false,
               mirror: false, frame: "tip"};
const unpack = b64 => { const s = atob(b64), u = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) u[i] = s.charCodeAt(i); return new Float32Array(u.buffer); };
const cache = {};

function build(shapeName) {
  if (cache[shapeName]) return cache[shapeName];
  const sh = DATA.shapes[shapeName], pts = unpack(sh.rings), n = sh.ring, rings = pts.length / 3 / n;
  const x = new Float32Array(rings * n), y = new Float32Array(rings * n), z = new Float32Array(rings * n);
  for (let v = 0; v < rings * n; v++) { x[v] = pts[3 * v]; y[v] = pts[3 * v + 1]; z[v] = pts[3 * v + 2]; }
  const I = [], J = [], K = [], nx = new Float32Array(rings * n), ny = new Float32Array(rings * n), nz = new Float32Array(rings * n);
  for (let r = 0; r < rings - 1; r++) {
    for (let p = 0; p < n - 1; p++) {
      const a = r * n + p, b = (r + 1) * n + p;
      for (const [i, j, k] of [[a, b, a + 1], [a + 1, b, b + 1]]) {
        const ux = x[j] - x[i], uy = y[j] - y[i], uz = z[j] - z[i], vx = x[k] - x[i], vy = y[k] - y[i], vz = z[k] - z[i];
        const fx = uy * vz - uz * vy, fy = uz * vx - ux * vz, fz = ux * vy - uy * vx;
        for (const q of [i, j, k]) { nx[q] += fx; ny[q] += fy; nz[q] += fz; }
        // the parent's last ring and the winglet's first are the same section: draw no faces between them
        if (r !== PARENT_RINGS - 1) { I.push(i); J.push(j); K.push(k); }
      }
    }
  }
  const stripes = new Float32Array(rings * n), colours = [];
  for (let v = 0; v < rings * n; v++) {
    const m = Math.hypot(nx[v], ny[v], nz[v]) || 1;
    const c = Math.max(-1, Math.min(1, (nx[v] * LIGHT[0] + ny[v] * LIGHT[1] + nz[v] * LIGHT[2]) / m));
    stripes[v] = Math.sin(2 * Math.PI * (Math.acos(c) * 180 / Math.PI) / PERIOD);
    const col = Math.floor(v / n) < PARENT_RINGS ? WING : WINGLET;
    colours.push("rgb(" + col.map(t => Math.round(t * 255)).join(",") + ")");
  }
  const line = b64 => { const p = unpack(b64), out = {x: [], y: [], z: []};
    for (let i = 0; i < p.length; i += 3) { out.x.push(p[i]); out.y.push(p[i + 1]); out.z.push(p[i + 2]); } return out; };
  const le = line(sh.le), te = line(sh.te), sec = line(sh.sections);
  const edges = {x: [...le.x, null, ...te.x], y: [...le.y, null, ...te.y], z: [...le.z, null, ...te.z]};
  const sections = {x: [], y: [], z: []};
  for (let i = 0; i < sec.x.length; i++) {
    if (i && i % sh.sectionPoints === 0) { sections.x.push(null); sections.y.push(null); sections.z.push(null); }
    sections.x.push(sec.x[i]); sections.y.push(sec.y[i]); sections.z.push(sec.z[i]);
  }
  return cache[shapeName] = {x, y, z, I, J, K, stripes, colours, edges, sections};
}

const flip = (arr, on) => on ? Array.from(arr, v => v === null ? null : -v) : arr;
function meshTrace(g, mirror, style) {
  const base = {type: "mesh3d", x: flip(g.x, mirror), y: g.y, z: g.z, i: g.I, j: mirror ? g.K : g.J, k: mirror ? g.J : g.K,
                hoverinfo: "skip", showlegend: false, flatshading: false};
  if (style === "ghost") return {...base, color: "#9fb3c4", opacity: 0.28, lighting: {ambient: 0.9, diffuse: 0.2}};
  if (style === "stripes") return {...base, intensity: g.stripes, intensitymode: "vertex", cmin: -1, cmax: 1, showscale: false,
                colorscale: [[0, "#111"], [0.5, "#111"], [0.5, "#f4f4f4"], [1, "#f4f4f4"]],
                lighting: {ambient: 1, diffuse: 0, specular: 0}};
  return {...base, vertexcolor: g.colours, lighting: {ambient: 0.62, diffuse: 0.6, specular: 0.25, roughness: 0.5, fresnel: 0.05},
          lightposition: {x: 1500, y: -2500, z: 4000}};
}
const lineTrace = (l, mirror, color, width, dash) => ({type: "scatter3d", mode: "lines", x: flip(l.x, mirror), y: l.y, z: l.z,
  line: {color, width, dash}, hoverinfo: "skip", showlegend: false});

const FRAMES = {tip: {x: [495, 620], y: [50, 190], z: [30, 140]}, whole: {x: [-10, 620], y: [-10, 260], z: [-40, 150]}};
function sceneRanges() {
  const r = FRAMES[state.frame], x = state.mirror ? [-r.x[1], r.x[1]] : r.x;
  const span = [x[1] - x[0], r.y[1] - r.y[0], r.z[1] - r.z[0]], big = Math.max(...span);
  return {x, y: r.y, z: r.z, aspect: {x: span[0] / big * 1.55, y: span[1] / big * 1.55, z: span[2] / big * 1.55}};
}

function draw() {
  const cfg = DATA.configs[state.config], g = build(cfg.shape), traces = [];
  const sides = state.mirror ? [false, true] : [false];
  if (state.ghost) {
    const base = build((DATA.configs["bend-centred"] || DATA.configs["baseline"]).shape);
    for (const m of sides) traces.push(meshTrace(base, m, "ghost"));
  }
  for (const m of sides) traces.push(meshTrace(g, m, state.mode));
  if (state.edges) for (const m of sides) traces.push(lineTrace(g.edges, m, "#b03a2e", 5));
  if (state.sections) for (const m of sides) traces.push(lineTrace(g.sections, m, "#1f4e79", 3));
  if (state.span) {
    const top = 150, plane = {x: [600, 600, 600, 600, 600], y: [55, 175, 175, 55, 55], z: [35, 35, top, top, 35]};
    for (const m of sides) traces.push(lineTrace(plane, m, "#7d3c98", 4, "dash"));
  }
  const r = sceneRanges();
  const layout = {margin: {l: 0, r: 0, t: 0, b: 0}, paper_bgcolor: "#eef2f5", uirevision: "keep",
    scene: {xaxis: {title: {text: "outboard (mm)"}, range: r.x}, yaxis: {title: {text: "aft (mm)"}, range: r.y},
            zaxis: {title: {text: "up (mm)"}, range: r.z}, aspectmode: "manual", aspectratio: r.aspect,
            bgcolor: "#eef2f5", uirevision: "keep"}};
  if (!document.getElementById("plot").layout) {
    const [eye, up] = DATA.camera === "Winglet planform" ? planformCamera() : (CAMERAS[DATA.camera] || CAMERAS["Oblique"]);
    layout.scene.camera = {eye: {x: eye[0], y: eye[1], z: eye[2]}, up: {x: up[0], y: up[1], z: up[2]},
                           projection: {type: "orthographic"}};
  }
  Plotly.react("plot", traces, layout, {displaylogo: false, responsive: true});
  document.getElementById("facts").innerHTML = cfg.facts.map(([k, v]) => "<dt>" + k + "</dt><dd>" + v + "</dd>").join("");
  document.querySelectorAll("#configs button").forEach(b => b.setAttribute("aria-pressed", b.dataset.config === state.config));
  document.querySelectorAll("#modes button").forEach(b => b.setAttribute("aria-pressed", b.dataset.mode === state.mode));
  document.querySelectorAll("#layers button").forEach(b => b.setAttribute("aria-pressed", !!state[b.dataset.layer]));
  document.querySelectorAll("#frames button").forEach(b => b.setAttribute("aria-pressed", b.dataset.frame === state.frame));
}

document.getElementById("baseline").textContent = DATA.baseline + ".";
const host = document.getElementById("configs");
for (const [g, group] of Object.entries(DATA.groups)) {
  const h = document.createElement("h3"); h.textContent = group.title; host.append(h);
  const row = document.createElement("div"); row.className = "group"; host.append(row);
  for (const name of group.members) {
    const cfg = DATA.configs[name], b = document.createElement("button");
    b.textContent = cfg.label; b.dataset.config = name;
    if (cfg.fits && (cfg.shape === "bend-centred" || cfg.shape === "baseline")) b.classList.add("baseline");
    if (!cfg.fits) { b.disabled = true; b.title = cfg.reason; }
    b.onclick = () => { state.config = name; draw(); };
    row.append(b);
  }
}
document.querySelectorAll("#modes button").forEach(b => b.onclick = () => { state.mode = b.dataset.mode; draw(); });
document.querySelectorAll("#layers button").forEach(b => b.onclick = () => { state[b.dataset.layer] = !state[b.dataset.layer]; draw(); });
document.querySelectorAll("#frames button").forEach(b => b.onclick = () => { state.frame = b.dataset.frame; draw(); });
for (const [name, [eye, up]] of Object.entries(CAMERAS)) {
  const b = document.createElement("button"); b.textContent = name;
  b.onclick = () => Plotly.relayout("plot", {"scene.camera": {eye: {x: eye[0], y: eye[1], z: eye[2]}, up: {x: up[0], y: up[1], z: up[2]},
                                               center: {x: 0, y: 0, z: 0}, projection: {type: "orthographic"}}});
  document.getElementById("cameras").append(b);
}
{ // looking straight at the winglet's outer face, its span direction up the screen
  const b = document.createElement("button"); b.textContent = "Winglet planform";
  b.onclick = () => { const [eye, up] = planformCamera(); Plotly.relayout("plot", {"scene.camera": {
    eye: {x: eye[0], y: eye[1], z: eye[2]}, up: {x: up[0], y: up[1], z: up[2]}, center: {x: 0, y: 0, z: 0},
    projection: {type: "orthographic"}}}); };
  document.getElementById("cameras").prepend(b);
}
draw();
</script>
</body></html>
"""


if __name__ == "__main__":
    main()
