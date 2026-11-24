"""Create a standalone interactive comparison of two half-wing STEP solids."""

import argparse
import json
from pathlib import Path

import cadquery as cq
import numpy as np
import plotly.graph_objects as go
from OCP.TopAbs import TopAbs_REVERSED


CAMERAS = {
    "front": {"eye": {"x": 0, "y": -2, "z": 0}, "up": {"x": 0, "y": 0, "z": 1}},
    "top": {"eye": {"x": 0, "y": 0, "z": 2}, "up": {"x": 0, "y": 1, "z": 0}},
    "orbit": {"eye": {"x": 1, "y": -1.7, "z": 1.15}, "up": {"x": 0, "y": 0, "z": 1}},
}
for camera in CAMERAS.values():
    camera["projection"] = {"type": "orthographic"}

COLORS = (("#3382c4", "#194f7c"), ("#e17835", "#a13e14"))


def display_xyz(points):
    return [(-point.x, -point.z, point.y) for point in points]


def mesh_wing(solid):
    """Sample the loft skin and mesh the planar caps against the same end rings."""
    solid_faces = solid.Faces()
    caps = [face for face in solid_faces if face.geomType() == "PLANE"]
    skins = [face for face in solid_faces if face.geomType() == "BSPLINE"]
    if len(skins) != 1 or len(caps) != 2 or len(solid_faces) != 3:
        raise ValueError("Expected one B-spline wing skin and two planar end caps")
    skin = skins[0]
    surface = skin._geomAdaptor()
    u_min, u_max, v_min, v_max = skin._uvBounds()
    chord = np.linspace(u_min, u_max, 241)
    span = np.linspace(v_min, v_max, 241)
    # A STEP loft parameter spans each complete airfoil around the chord. The
    # other parameter runs root to tip. Default OCCT triangles can connect
    # far-apart spans, visibly flattening the front silhouette between them.
    grid = np.array([
        [(-point.X(), -point.Z(), point.Y()) for point in surface_points]
        for surface_points in (
            (surface.Value(u, v) for u in chord) for v in span
        )
    ])
    if not np.allclose(grid[:, :, 0], grid[:, :1, 0], atol=1e-5, rtol=0):
        raise ValueError("Expected constant-span airfoil sections on the wing skin")
    if not np.allclose(grid[:, 0], grid[:, -1], atol=1e-5, rtol=0):
        raise ValueError("Expected the wing skin to close at the trailing edge")
    # Keep the seam vertices separate for the sharp trailing-edge normals.
    ids = np.arange(len(span) * len(chord), dtype=np.int32).reshape(grid.shape[:2])
    a, b = ids[:-1, :-1].ravel(), ids[:-1, 1:].ravel()
    c, d = ids[1:, :-1].ravel(), ids[1:, 1:].ravel()
    faces = np.concatenate([np.column_stack([a, b, c]), np.column_stack([b, d, c])])
    if skin.wrapped.Orientation() == TopAbs_REVERSED:
        faces = faces[:, [0, 2, 1]]
    points = [grid.reshape(-1, 3)]
    triangles = [faces]
    offset = len(points[0])
    end_rings = sorted((grid[0], grid[-1]), key=lambda ring: ring[0, 0])
    caps.sort(key=lambda face: -face.Center().x)
    for cap, ring in zip(caps, end_rings):
        # Use these skin boundary points on the cap as well. Separate indices
        # preserve the real hard edge where the planar face meets the skin.
        wire = cq.Wire.makePolygon([cq.Vector(-x, z, -y) for x, y, z in ring[:-1]], close=True)
        polygon = cq.Face.makeFromWires(wire)
        vertices, faces = polygon.tessellate(0.02, 0.05)
        faces = np.asarray(faces, dtype=np.int32)
        if polygon.normalAt().dot(cap.normalAt()) < 0:
            faces = faces[:, [0, 2, 1]]
        points.append(np.asarray(display_xyz(vertices)))
        triangles.append(faces + offset)
        offset += len(vertices)
    nose = grid[np.arange(len(span)), np.argmin(grid[:, :, 1], axis=1)]
    return np.concatenate(points).astype(np.float32), np.concatenate(triangles), nose


def add_wing(figure, path, label, colors):
    solids = cq.importers.importStep(str(path)).solids().vals()
    if len(solids) != 1 or not solids[0].isValid():
        raise ValueError(f"Expected one valid wing solid in {path}")
    solid = solids[0]
    points, faces, nose = mesh_wing(solid)
    x, y, z = points.T
    # The display axis transform reverses handedness, so reverse winding.
    i, k, j = faces.T
    figure.add_trace(go.Mesh3d(
        x=x, y=y, z=z, i=i, j=j, k=k,
        name=label, color=colors[0], opacity=.42, flatshading=False,
        lighting=dict(ambient=.75, diffuse=.7, specular=.14, roughness=.62),
        lightposition=dict(x=240, y=-180, z=400), hoverinfo="skip",
        showscale=False, showlegend=False,
    ))
    edges = []
    for edge in solid.Edges():
        count = max(8, min(250, int(edge.Length() / 3)))
        edges.extend(display_xyz(edge.sample(count)[0]))
        edges.append((None, None, None))
    x, y, z = zip(*edges)
    figure.add_trace(go.Scatter3d(
        x=x, y=y, z=z, mode="lines", name=label,
        line=dict(color=colors[1], width=3), hoverinfo="skip", showlegend=False,
    ))
    x, y, z = nose.T
    figure.add_trace(go.Scatter3d(
        x=x, y=y, z=z, mode="lines", name=f"{label} leading edge",
        line=dict(color=colors[1], width=5), showlegend=False,
        hovertemplate=f"{label}, sampled nose trace<br>Outboard %{{x:.1f}} mm<br>"
                      "Aft %{y:.1f} mm<br>Up %{z:.1f} mm<extra></extra>",
    ))
    low, high = points.min(axis=0), points.max(axis=0)
    print(f"{label}: {len(points):,} vertices, {len(faces):,} triangles, bounds {low} to {high}")
    return low, high


def make_html(wing_a, wing_b):
    figure = go.Figure()
    bounds_a = add_wing(figure, wing_a, "Wing A", COLORS[0])
    bounds_b = add_wing(figure, wing_b, "Wing B", COLORS[1])
    low = np.minimum(bounds_a[0], bounds_b[0])
    high = np.maximum(bounds_a[1], bounds_b[1])
    low -= 10
    high += 10
    ranges = high - low
    # Orthographic zoom uses the scene aspect ratio. Scale every axis by the
    # same amount so both wings keep their physical proportions in every view.
    aspect = dict(zip("xyz", (float(size) for size in 2.05 * ranges / max(ranges))))
    axes = dict(showbackground=True, backgroundcolor="#f7fafc", gridcolor="#c8d2dc",
                zerolinecolor="#8c9cb0", showspikes=False, showline=True, linecolor="#8999ab")
    figure.update_layout(
        paper_bgcolor="#f8fafc", plot_bgcolor="#f8fafc", font_color="#233c50",
        margin=dict(l=0, r=0, t=0, b=0), uirevision="shape-comparison",
        scene=dict(
            xaxis={**axes, "title": "Outboard (mm)", "range": [float(low[0]), float(high[0])]},
            yaxis={**axes, "title": "Aft (mm)", "range": [float(low[1]), float(high[1])]},
            zaxis={**axes, "title": "Up (mm)", "range": [float(low[2]), float(high[2])]},
            aspectmode="manual", aspectratio=aspect, dragmode="orbit", camera=CAMERAS["front"],
        ),
    )
    plot = figure.to_html(
        full_html=False, include_plotlyjs=True, div_id="wing-comparison",
        default_height="100%", default_width="100%",
        config={"responsive": True, "scrollZoom": True, "displaylogo": False,
                "modeBarButtonsToRemove": ["resetCameraLastSave3d", "resetCameraDefault3d",
                                          "orbitRotation", "tableRotation"]},
    )
    return (
        TEMPLATE.replace("{{PLOT}}", plot)
        .replace("{{CAMERAS}}", json.dumps(CAMERAS))
    )


TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Wing comparison</title>
<style>
:root { font-family: system-ui,sans-serif; color: #132d43; }
body { margin: 0; background: #f8fafc; height: 100dvh; display: grid;
       grid-template-rows: auto minmax(0,1fr); }
header { padding: 15px 24px; background: #fff; border-bottom: 1px solid #c8d5e1;
         display: flex; justify-content: space-between; flex-wrap: wrap; gap: 12px; }
h1 { margin: 0 0 5px; font-size: 20px; font-weight: 650; }
p { margin: 0; font-size: 13px; color: #476076; }
nav { display: flex; flex-wrap: wrap; align-items: center; gap: 6px; }
nav + nav { margin-left: 10px; }
.controls { display: flex; flex-wrap: wrap; gap: 6px; }
button { background: #fff; color: #233c50; font: inherit; font-size: 13px;
         padding: 7px 10px; border: 1px solid #b7c6d4; border-radius: 5px; cursor: pointer; }
button:hover, button:focus-visible { border-color: #39719d; }
button[aria-pressed="true"] { color: #11354f; background: #e4f0fa; border-color: #5486ad; }
.key { border-left: 3px solid; padding-left: 5px; }
.wing_a { border-color: #194f7c; }
.wing_b { border-color: #a13e14; margin-left: 8px; }
main { min-height: 0; overflow: hidden; }
main > div { height: 100%; }
</style>
</head>
<body>
<header>
<div>
  <h1>Wing comparison</h1>
  <p><span class="key wing_a">Wing A</span><span class="key wing_b">Wing B</span>
  · Drag to orbit · Scroll to zoom</p>
</div>
<div class="controls">
  <nav aria-label="Visible wings">
    <button data-layer="both" aria-pressed="true">Overlay</button>
    <button data-layer="wing_a" aria-pressed="false">Wing A</button>
    <button data-layer="wing_b" aria-pressed="false">Wing B</button>
  </nav>
  <nav aria-label="Camera">
    <button data-camera="front" aria-pressed="true">Front</button>
    <button data-camera="top" aria-pressed="false">Top</button>
    <button data-camera="orbit" aria-pressed="false">3/4</button>
  </nav>
</div>
</header>
<main>{{PLOT}}</main>
<script>
const wing = document.getElementById("wing-comparison");
const cameras = {{CAMERAS}};
const cameraButtons = document.querySelectorAll("[data-camera]");
const layerButtons = document.querySelectorAll("[data-layer]");
const layers = {both: [true, true], wing_a: [true, false], wing_b: [false, true]};
function markActive(buttons, active) {
    buttons.forEach(button => button.setAttribute("aria-pressed", String(button === active)));
}
cameraButtons.forEach(button => button.addEventListener("click", async () => {
    await Plotly.relayout(wing, {"scene.camera": cameras[button.dataset.camera]});
    markActive(cameraButtons, button);
}));
layerButtons.forEach(button => button.addEventListener("click", async () => {
    const layer = button.dataset.layer;
    const [showA, showB] = layers[layer];
    const opacity = layer === "both" ? .42 : 1;
    await Plotly.restyle(wing, {
        visible: [showA, showA, showA, showB, showB, showB],
        opacity: [opacity, 1, 1, opacity, 1, 1]
    });
    markActive(layerButtons, button);
}));
</script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wing_a", type=Path)
    parser.add_argument("wing_b", type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("compare-wings.html"))
    args = parser.parse_args()
    args.output.write_text(make_html(args.wing_a, args.wing_b), encoding="utf-8")
    print(args.output.resolve())


if __name__ == "__main__":
    main()
