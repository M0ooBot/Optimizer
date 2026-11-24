"""Render the mirrored full-span CAD fit study with spanwise tip dimensions."""

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import cadquery as cq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go


HERE = Path(__file__).resolve().parent
DIM = "#345269"
FAIL = "#a23d37"
PASS = "#146b5d"
WITNESS = "#992d62"
SKIN = "#afcadb"
BACKGROUND = "#f2f6f9"
MIRROR = np.array([-1, 1, 1])


def helpers():
    spec = importlib.util.spec_from_file_location("spar_view", HERE / "render-viewer.py")
    module = importlib.util.module_from_spec(spec)
    previous = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous
    return module


def fit_text(result, wall_digits):
    wall = result["min_wall_mm"]
    if wall is None:
        return f'Best tested axis fails fit · {result["failure_reason"]}'
    if result["geometry_pass"]:
        return f"CAD clearance passes · Minimum wall {wall:.{wall_digits}f} mm"
    if wall < 0:
        return f"Best tested axis fails fit · Skin breach {abs(wall):.2f} mm"
    return f"Best tested axis fails fit · Minimum wall {wall:.3f} mm"


def dimensions(points, end):
    tip = points[np.argmax(points[:, 0])]
    y = max(float(tip[1]), float(end[1])) + 49
    setback = float(tip[0] - end[0])
    guides = np.array([
        end, [end[0], y, end[2]],
        [np.nan] * 3,
        [tip[0], y, end[2]], tip,
    ])
    arrows = np.array([[end[0], y, end[2]], [tip[0], y, end[2]]])
    return guides, arrows, setback


def make_html(out, view, points, faces, outline, result, ends, guides, arrows, labels):
    fit_color = PASS if result["geometry_pass"] else FAIL
    fig = go.Figure()
    for side, winding in ((1, faces[:, [0, 2, 1]]), (MIRROR, faces)):
        x, y, z = (points * side).T
        i, j, k = winding.T
        fig.add_trace(go.Mesh3d(
            x=x, y=y, z=z, i=i, j=j, k=k, color=SKIN, opacity=.48,
            lighting={"ambient": .94, "diffuse": .45, "specular": .10},
            hoverinfo="skip", showlegend=False,
        ))
        fig.add_trace(view.line3d(outline * side, "Wing outline", DIM, mode="lines", hoverinfo="skip"))
    overlay_start = len(fig.data)
    channel = go.Figure()
    view.add_channel(channel, result, "Verified channel" if result["geometry_pass"] else "Proposed channel", True)
    for sign, side in ((1, 1), (-1, MIRROR)):
        clone = channel.data[0].to_plotly_json()
        clone.update(x=np.asarray(channel.data[0].x) * sign,
                     opacity=.48, color=fit_color)
        if sign == -1:
            clone["j"], clone["k"] = clone["k"], clone["j"]
        fig.add_trace(clone)
        fig.add_trace(view.line3d(
            ends * side, "Spar axis", fit_color, width=7, mode="lines+markers",
            marker={"size": 5, "color": fit_color}, text=["Root", "Channel cap center"],
            hovertemplate="%{text}<br>Span coordinate %{x:.1f} mm<extra></extra>",
        ))
        fig.add_trace(view.line3d(guides * side, "Spanwise guides", DIM, width=2, mode="lines",
                                  hoverinfo="skip"))
        fig.add_trace(view.line3d(
            arrows * side, "Spanwise setback", DIM, width=4, mode="lines+markers",
            marker={"size": 3, "color": DIM}, hoverinfo="skip",
        ))
        fig.add_trace(view.line3d(
            arrows.mean(axis=0, keepdims=True) * side + (0, 9, 0),
            "Spanwise setback", DIM, mode="text", text=[labels["setback"]],
            textfont={"color": DIM, "size": 14}, hoverinfo="skip",
        ))
        if labels["witness"]:
            witness = view.display([result["worst_channel_point_mm"], result["worst_skin_point_mm"]])
            fig.add_trace(view.line3d(
                witness * side, labels["witness"], labels["witness_color"], width=7,
                mode="lines+markers",
                marker={"size": 6, "symbol": "x", "color": labels["witness_color"]},
                hovertemplate=labels["witness"] + "<br>Span coordinate %{x:.1f} mm<extra></extra>",
            ))
    low = points.min(axis=0).astype(float) - 9
    high = points.max(axis=0).astype(float) + 9
    high[1] = max(high[1], arrows[0, 1] + 26)
    low[0] = -high[0]
    scene = view.axis_layout(low, high)
    for axis in "xyz":
        scene[f"{axis}axis"]["visible"] = False
    cameras = {
        "Top": view.camera(0, 0, 1.95, (0, -1, 0)),
        "Oblique": view.camera(-.10, -1.42, 1.30),
        "Front": view.camera(0, -1.95, 0),
    }
    fig.update_layout(
        margin=dict(l=0, r=0, t=0, b=0), paper_bgcolor="#ffffff",
        scene=dict(**scene, dragmode="orbit", camera=cameras["Top"]),
    )
    plot = fig.to_html(full_html=False, include_plotlyjs=True, div_id="wing-view",
                       config=dict(responsive=True, scrollZoom=True, displayModeBar=False))
    template = (HERE / "tip-setback.template.html").read_text(encoding="utf-8")
    replacements = {
        "PLOT": plot, "FIT": labels["fit"],
        "COLOR": fit_color, "TITLE": labels["title"], "LENGTH": labels["length"],
        "WITNESS": labels["witness_note"],
        "SPAN": f"{2 * np.max(points[:, 0]):,.0f}",
        "EDGE": f'{points[:, 0].max() - result["channel_farthest_span_mm"]:.1f}',
        "DATA": json.dumps(dict(cameras=cameras, overlayStart=overlay_start)),
    }
    for key, value in replacements.items():
        template = template.replace("{{" + key + "}}", value)
    (out / "full-wing.html").write_text(template, encoding="utf-8")


def make_png(out, points, result, ends, guides, arrows, labels):
    fit_color = PASS if result["geometry_pass"] else FAIL
    axis_style = "o-" if result["geometry_pass"] else "o--"
    grid = points[:241**2].reshape(241, 241, 3)
    rings = np.arange(len(grid))
    front = grid[rings, np.argmin(grid[:, :, 1], axis=1)]
    back = grid[rings, np.argmax(grid[:, :, 1], axis=1)]
    profile = np.concatenate((front, back[::-1]))

    fig = plt.figure(figsize=(12, 5.25), dpi=145, facecolor=BACKGROUND)
    axes = fig.add_axes((.035, .16, .93, .66), aspect="equal", facecolor=BACKGROUND)
    for sign, side in ((1, 1), (-1, MIRROR)):
        x, y, _ = (profile * side).T
        axes.fill(x, y, color=SKIN, edgecolor=DIM, linewidth=1.05)
        x, y, _ = (ends * side).T
        axes.plot(x, y, axis_style, color=fit_color, linewidth=2.4, markersize=5, zorder=4)
        if labels["witness"]:
            cad_x, _, cad_z = result["worst_skin_point_mm"]
            axes.scatter(-cad_x * sign, -cad_z, marker="x", s=66,
                         c=labels["witness_color"], zorder=5)
        x, y, _ = (guides * side).T
        axes.plot(x, y, color=DIM, linewidth=.8, alpha=.72, zorder=3)
        start, finish = (arrows * side)[:, :2]
        axes.annotate("", xy=finish, xytext=start, arrowprops={
            "arrowstyle": "<->", "color": DIM, "linewidth": 1.05, "shrinkA": 0, "shrinkB": 0,
        })
        label = (start + finish) / 2 + (0, 11)
        axes.text(*label, labels["setback"], color=DIM, fontsize=10, ha="center", va="top")
    wing_tip = float(points[:, 0].max())
    axes.set_xlim(-wing_tip - 14, wing_tip + 14)
    axes.set_ylim(float(arrows[0, 1]) + 31, float(points[:, 1].min()) - 13)
    axes.axis("off")
    fig.text(.043, .919, labels["title"], fontsize=23, color="#1e384c", fontweight="bold")
    fig.text(.044, .862, f"{2 * wing_tip:,.0f} mm tip to tip · Revised half-wing reflected at the root plane",
             fontsize=10.5, color=DIM)
    headline = " · ".join(part for part in (labels["fit"], labels["length"]) if part)
    fig.text(.044, .091, headline, fontsize=11, color=fit_color, fontweight="bold")
    witness_note = f' {labels["witness_note"]}' if labels["witness_note"] else ""
    axis_description = "Green lines" if result["geometry_pass"] else "Dashed lines"
    fig.text(.044, .052,
             f"{axis_description} show the channel axes.{witness_note} End dimensions are spanwise.",
             fontsize=9.5, color=DIM)
    fig.savefig(out / "full-wing.png", dpi=fig.dpi, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", type=Path)
    parser.add_argument("result", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--longest", action="store_true", help="Label the longest verified passing axis")
    parser.add_argument("--wall-digits", type=int, choices=range(1, 7), default=3)
    parser.add_argument("--length-digits", type=int, choices=range(1, 7), default=1)
    args = parser.parse_args()
    view = helpers()
    result = json.loads(args.result.read_text(encoding="utf-8"))
    if hashlib.sha256(args.step.read_bytes()).hexdigest() != result["input_sha256"]:
        raise ValueError("Fit result was measured on a different STEP")
    solid, = cq.importers.importStep(str(args.step)).solids().vals()
    renderer = view.load_renderer()
    points, faces, nose = renderer.mesh_wing(solid)
    outline = list(nose)
    for edge in solid.Edges():
        outline.append((np.nan,) * 3)
        outline.extend(renderer.display_xyz(edge.sample(181)[0]))
    ends = view.display([result["axis_root_mm"], result["axis_channel_end_mm"]])
    guides, arrows, setback = dimensions(points, ends[1])
    if args.longest and not result["geometry_pass"]:
        raise ValueError("The longest verified axis requires a passing CAD result")
    if not args.longest and abs(setback - 50) > .02:
        raise ValueError(f"Expected a 50 mm setback, measured {setback:.3f} mm")
    has_breach = result["min_wall_mm"] is not None and result["min_wall_mm"] < 0
    witness_note = ""
    if args.longest:
        witness_span = abs(result["worst_skin_point_mm"][0])
        witness_note = f"Crosses mark minimum wall at {witness_span:.1f} mm outboard each side."
    elif has_breach:
        witness_note = "Crosses mark the skin breach."
    labels = {
        "title": "Wing B · Longest verified spar" if args.longest else "Wing B · 50 mm tip setback",
        "setback": f"{setback:.{args.length_digits}f} mm" if args.longest else f"{setback:.0f} mm",
        "length": (
            f'{result["axis_length_mm"]:.{args.length_digits}f} mm root to cap per half'
            if args.longest else ""
        ),
        "fit": fit_text(result, args.wall_digits),
        "witness": "Minimum wall" if args.longest else "Skin breach" if has_breach else "",
        "witness_color": WITNESS if args.longest else FAIL,
        "witness_note": witness_note,
    }
    out = args.out or HERE / ("longest-spar" if args.longest else "tip-setback-50")
    out.mkdir(parents=True, exist_ok=True)
    make_html(out, view, points, faces, np.asarray(outline), result, ends, guides, arrows, labels)
    make_png(out, points, result, ends, guides, arrows, labels)


if __name__ == "__main__":
    main()
