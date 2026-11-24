"""Render the rounded-tip STEP and its baseline into standalone review files."""

import argparse
import base64
import json
from dataclasses import dataclass
from pathlib import Path

import cadquery as cq
import matplotlib
import numpy as np
import plotly.graph_objects as go
from OCP.BRepGProp import BRepGProp_Face
from OCP.TopAbs import TopAbs_REVERSED
from OCP.gp import gp_Pnt, gp_Vec

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


HERE = Path(__file__).resolve().parent
BASELINE = HERE.parent / "paired-wing-trial/a/measured/01-wing-occt.step"
REVISED = HERE / "wing-b-rounded.step"
TRANSITION = 540
DETAIL_START = 525
CUTS = (90, 110, 135)
DEFAULT_CUT = 110
PAPER = "#f4f6f5"
INK = "#132e34"
MUTED = "#50676c"
GRID = "#e4ebea"
BASE = "#62747c"
NEW = "#007b78"
GUIDE = "#b07045"
LABELS = ("Original", "Rounded tip")
COLORS = (BASE, NEW)
SCOPES = ("Full wing", "Full half-wing", "Tip detail")
ROOT_TOL = 1e-5


@dataclass
class Wing:
    points: np.ndarray
    triangles: np.ndarray
    edges: list
    profiles: dict
    zebra_mesh: tuple


def display_xyz(vectors):
    coordinates = np.asarray([point.toTuple() for point in vectors])
    return coordinates[:, [0, 2, 1]] * [-1, -1, 1]


def wire_points(edges, minimum=60, interval=1.5):
    lines = []
    for edge in edges:
        count = max(minimum, int(edge.Length() / interval) + 1)
        lines.append(display_xyz(edge.sample(count)[0]))
    return lines


def evaluate_surface(face, parameters):
    """Evaluate CAD position and oriented normal, including collapsed UV poles."""
    points, normals = [], []
    position, normal = gp_Pnt(), gp_Vec()
    evaluator = BRepGProp_Face(face.wrapped)
    for u, v in parameters:
        evaluator.Normal(u, v, position, normal)
        points.append((-position.X(), -position.Z(), position.Y()))
        normals.append((-normal.X(), -normal.Z(), normal.Y()))
    normals = np.asarray(normals)
    lengths = np.linalg.norm(normals, axis=1)
    regular = np.isfinite(lengths) & (lengths > 1e-14)
    if not np.any(regular):
        raise ValueError("STEP face has no regular surface normals")
    normals[regular] /= lengths[regular, None]
    if not np.all(regular):
        # A collapsed UV pole takes its limiting normal from the same face.
        scales = np.maximum(np.ptp(parameters, axis=0), 1e-12)
        for index in np.flatnonzero(~regular):
            distance = np.linalg.norm((parameters[regular] - parameters[index]) / scales, axis=1)
            normals[index] = normals[regular][distance.argmin()]
    return np.asarray(points, dtype=np.float32), normals.astype(np.float32)


def rectangular_face_mesh(face):
    """The Wing B spline faces have rectangular UV domains, including collapsed sides."""
    u0, u1, v0, v1 = face._uvBounds()
    is_body = face.BoundingBox().xmax > -1
    u_count, v_count = (321, 121) if is_body else (201, 101)
    cosine = (1 - np.cos(np.linspace(0, np.pi, u_count))) / 2
    u = u0 + (u1 - u0) * cosine
    v = np.linspace(v0, v1, v_count)
    if not is_body:
        v = v0 + (v1 - v0) * (1 - np.cos(np.linspace(0, np.pi, v_count))) / 2
    uu, vv = np.meshgrid(u, v)
    points, normals = evaluate_surface(face, np.column_stack((uu.ravel(), vv.ravel())))
    ids = np.arange(len(points), dtype=np.int32).reshape(v_count, u_count)
    a, b = ids[:-1, :-1].ravel(), ids[:-1, 1:].ravel()
    c, d = ids[1:, :-1].ravel(), ids[1:, 1:].ravel()
    triangles = np.concatenate((np.column_stack((a, b, c)), np.column_stack((b, d, c))))
    if face.wrapped.Orientation() != TopAbs_REVERSED:
        triangles = triangles[:, [0, 2, 1]]
    return points, normals, triangles


def zebra_mesh(solid):
    points, normals, triangles = [], [], []
    offset = 0
    for face in solid.Faces():
        if face.geomType() == "PLANE":
            vertices, facets = face.tessellate(0.06, 0.10)
            positions = display_xyz(vertices).astype(np.float32)
            direction = display_xyz([face.normalAt()])
            directions = np.repeat(direction, len(positions), axis=0).astype(np.float32)
            facets = np.asarray(facets, dtype=np.int32)[:, [0, 2, 1]]
        else:
            positions, directions, facets = rectangular_face_mesh(face)
        points.append(positions)
        normals.append(directions)
        triangles.append(facets + offset)
        offset += len(positions)
    return np.concatenate(points), np.concatenate(normals), np.concatenate(triangles)


def load_wing(path):
    solids = cq.importers.importStep(str(path)).solids().vals()
    if len(solids) != 1 or not solids[0].isValid():
        raise ValueError(f"Expected one valid STEP solid in {path}")
    solid = solids[0]
    vertices, triangles = solid.tessellate(0.06, 0.10)
    points = display_xyz(vertices).astype(np.float32)
    # The display coordinate transform reverses the OCCT triangle winding.
    triangles = np.asarray(triangles, dtype=np.int32)[:, [0, 2, 1]]
    profiles = {}
    for aft in CUTS:
        plane = cq.Face.makePlane(basePnt=(0, 0, -aft), dir=(0, 0, 1))
        edges = solid.intersect(plane).Edges()
        if not edges:
            raise ValueError(f"Empty CAD section at {aft} mm aft in {path}")
        profiles[aft] = wire_points(edges, minimum=170, interval=0.45)
    print(f"{path.name}: {len(points):,} vertices, {len(triangles):,} triangles")
    return Wing(points, triangles, wire_points(solid.Edges()), profiles, zebra_mesh(solid))


def envelope(wing, stations):
    """Intersect the tessellation with span planes for top and front silhouettes."""
    triangles = wing.triangles
    edges = np.sort(np.concatenate((triangles[:, :2], triangles[:, 1:],
                                    triangles[:, [0, 2]])), axis=1)
    segments = wing.points[np.unique(edges, axis=0)]
    segments = segments[segments[:, 0, 0] != segments[:, 1, 0]]
    first, second = segments.transpose(1, 0, 2)
    lower = np.minimum(first[:, 0], second[:, 0])
    upper = np.maximum(first[:, 0], second[:, 0])
    results = []
    for span in stations:
        mask = (lower <= span) & (span <= upper)
        delta = second[mask] - first[mask]
        fraction = (span - first[mask, 0]) / delta[:, 0]
        points = first[mask, 1:] + fraction[:, None] * delta[:, 1:]
        if not len(points):
            raise ValueError(f"Tessellation has no intersections at {span} mm outboard")
        results.append((points.min(axis=0), points.max(axis=0)))
    return np.asarray(results)


def break_lines(lines):
    separator = np.full((1, 3), np.nan, dtype=np.float32)
    return np.concatenate([np.vstack((line, separator)) for line in lines]).T


def clip_tip(wing):
    """Keep the outboard mesh and split triangles that cross the display plane."""
    geometry = wing.points[wing.triangles]
    inside = geometry[:, :, 0] >= DETAIL_START
    counts = inside.sum(axis=1)
    pieces = [geometry[counts == 3]]

    def intersections(inside_point, outside_point):
        delta = outside_point.astype(np.float64) - inside_point
        weight = (DETAIL_START - inside_point[:, 0].astype(np.float64)) / delta[:, 0]
        cut = inside_point + weight[:, None] * delta
        cut[:, 0] = DETAIL_START
        return cut.astype(np.float32)

    for count, isolated in ((1, True), (2, False)):
        selected = counts == count
        separate = np.argmax(inside[selected] == isolated, axis=1)
        order = (separate[:, None] + np.arange(3)) % 3
        single, next_point, last_point = geometry[selected][
            np.arange(len(order))[:, None], order
        ].transpose(1, 0, 2)
        if isolated:
            cut_next = intersections(single, next_point)
            cut_last = intersections(single, last_point)
            pieces.append(np.stack((single, cut_next, cut_last), axis=1))
        else:
            cut_next = intersections(next_point, single)
            cut_last = intersections(last_point, single)
            pieces.extend((np.stack((cut_next, next_point, last_point), axis=1),
                           np.stack((cut_next, last_point, cut_last), axis=1)))

    vertices, indices = np.unique(np.concatenate(pieces).reshape(-1, 3), axis=0, return_inverse=True)
    return vertices, indices.reshape(-1, 3).astype(np.int32)


def outside_triangles(wing):
    root = np.abs(wing.points[:, 0]) <= ROOT_TOL
    return wing.triangles[~np.all(root[wing.triangles], axis=1)]


def mirror_wing(wing):
    """Reflect the half-wing at its root and omit the internal cap and wire."""
    root = np.abs(wing.points[:, 0]) <= ROOT_TOL
    exterior = outside_triangles(wing)
    mirrored = wing.points[~root] * np.array([-1, 1, 1], dtype=np.float32)
    mirror_ids = np.arange(len(wing.points), dtype=np.int32)
    mirror_ids[~root] = len(wing.points) + np.arange(len(mirrored), dtype=np.int32)
    points = np.concatenate((wing.points, mirrored))
    points[root.nonzero()[0], 0] = 0
    triangles = np.concatenate((exterior, mirror_ids[exterior[:, [0, 2, 1]]]))
    edges = [line for line in wing.edges if not np.all(np.abs(line[:, 0]) <= ROOT_TOL)]
    return points, triangles, edges + [line * [-1, 1, 1] for line in edges]


def add_geometry(figure, wing, label, color):
    for scope in SCOPES:
        visible = scope == SCOPES[0] and label == LABELS[1]
        if scope == SCOPES[0]:
            points, triangles, lines = mirror_wing(wing)
        elif scope == "Tip detail":
            points, triangles = clip_tip(wing)
            lines = [line.copy() for line in wing.edges]
            for line in lines:
                line[line[:, 0] < DETAIL_START] = np.nan
        else:
            points, triangles, lines = wing.points, wing.triangles, wing.edges
        x, y, z = points.T
        i, j, k = triangles.T
        surface = "#9bafb5" if label == LABELS[0] else "#51b5ac"
        figure.add_trace(go.Mesh3d(
            x=x, y=y, z=z, i=i, j=j, k=k, name=label, color=surface,
            opacity=1, visible=visible, hoverinfo="skip",
            lighting=dict(ambient=.78, diffuse=.9, specular=.38, roughness=.44),
            lightposition=(dict(x=200, y=-390, z=550) if scope == SCOPES[0]
                           else dict(x=590, y=38, z=170)), showlegend=False,
        ))
        x, y, z = break_lines(lines)
        figure.add_trace(go.Scatter3d(
            x=x, y=y, z=z, mode="lines", name=label, visible=visible,
            line=dict(color=color, width=2.5), hoverinfo="skip", showlegend=False,
        ))


def section_ranges(wings, cut):
    lines = [line for wing in wings for line in wing.profiles[cut]]
    heights = np.concatenate([line[line[:, 0] >= DETAIL_START, 2] for line in lines])
    return [[DETAIL_START, 608], [float(heights.min() - 3), float(heights.max() + 3)]]


def camera(eye, up=(0, 0, 1)):
    return dict(eye=dict(zip("xyz", eye)), up=dict(zip("xyz", up)),
                projection=dict(type="orthographic"))


def scene_settings(ranges, *, full=False):
    axis_titles = ("Span (mm)" if full else "Outboard (mm)", "Aft (mm)", "Height (mm)")
    axis = dict(showbackground=False, gridcolor="#dbe5e4", zeroline=False, showspikes=False,
                tickfont=dict(size=10))
    settings = {}
    for name, limits, title in zip("xyz", ranges, axis_titles, strict=True):
        settings[f"scene.{name}axis"] = dict(axis, title=dict(text=title, font=dict(size=11)),
                                              range=limits)
    lengths = [abs(end - start) for start, end in ranges]
    scale = 5.0 if full else 1.9
    settings["scene.aspectratio"] = dict(zip("xyz", (scale * n / max(lengths) for n in lengths)))
    return settings


def zebra_data(wings):
    def encode(array, dtype):
        return base64.b64encode(array.astype(dtype).tobytes()).decode("ascii")

    models = {}
    for label, wing in zip(LABELS, wings, strict=True):
        points, normals, triangles = wing.zebra_mesh
        root = np.abs(points[:, 0]) <= ROOT_TOL
        skin = triangles[~np.all(root[triangles], axis=1)]
        models[label] = dict(positions=encode(points, "<f4"), normals=encode(normals, "<f4"),
                             triangles=encode(triangles, "<u4"), skin=encode(skin, "<u4"))
    return json.dumps(models)


def make_html(wings, output):
    view = go.Figure()
    section = go.Figure()
    for wing, label, color in zip(wings, LABELS, COLORS, strict=True):
        add_geometry(view, wing, label, color)
        for cut in CUTS:
            x, _, z = break_lines(wing.profiles[cut])
            section.add_trace(go.Scatter(
                x=x, y=z, name=label, visible=cut == DEFAULT_CUT, mode="lines",
                line=dict(color=color, width=2.7, dash="dash" if label == LABELS[0] else "solid",
                          simplify=False),
                hovertemplate=f"{label}<br>%{{x:.2f}} mm outboard<br>%{{y:.2f}} mm high<extra></extra>",
                showlegend=False,
            ))
    cameras = {
        "Angled": camera((1.25, -1.45, .79)),
        "Full angled": camera((0, -1.58, 1.52)),
        "Lip": camera((.26, -1.62, .19)),
        "Top": camera((0, 0, 2), up=(0, 1, 0)),
        "Front": camera((0, -2, 0)),
    }
    scopes = {
        "Full wing": scene_settings([[-618, 618], [-8, 208], [-9, 76]], full=True),
        "Full half-wing": scene_settings([[-10, 613], [-8, 208], [-9, 76]]),
        "Tip detail": scene_settings([[DETAIL_START, 608], [54, 159], [41, 73]]),
    }
    base = dict(paper_bgcolor="#ffffff", plot_bgcolor="#ffffff",
                font=dict(family="system-ui, sans-serif", size=12, color=INK))
    view.update_layout(**base, margin=dict(l=0, r=0, t=0, b=0), uirevision="rounded-tip",
                       scene=dict(aspectmode="manual", camera=cameras["Full angled"], dragmode="orbit"))
    for field, value in scopes["Full wing"].items():
        view.layout.scene[field.split(".")[1]] = value
    for name in "xyz":
        view.layout.scene[f"{name}axis"].visible = False
    ranges = {cut: section_ranges(wings, cut) for cut in CUTS}
    chart_axis = dict(gridcolor=GRID, linecolor="#b8c8c6", showline=True,
                      ticks="outside", tickcolor="#b8c8c6", zeroline=False)
    section.update_layout(
        **base, margin=dict(l=60, r=22, t=7, b=51), dragmode="pan",
        xaxis=dict(**chart_axis, title="Outboard (mm)", range=ranges[DEFAULT_CUT][0]),
        yaxis=dict(**chart_axis, title="Height (mm)", range=ranges[DEFAULT_CUT][1],
                   scaleanchor="x", scaleratio=1, constrain="domain"),
    )
    config = dict(responsive=True, scrollZoom=True, displayModeBar=False)
    replacements = {
        "VIEW": view.to_html(full_html=False, include_plotlyjs=True, div_id="wing", config=config,
                             default_height="100%"),
        "SECTION": section.to_html(full_html=False, include_plotlyjs=False, div_id="cut", config=config,
                                   default_height="100%"),
        "CAMERAS": json.dumps(cameras), "SCOPES": json.dumps(scopes), "CUTS": json.dumps(ranges),
        "ZEBRA_DATA": zebra_data(wings),
        "ZEBRA_RENDERER": Path(__file__).with_name("zebra.js").read_text(encoding="utf-8"),
    }
    page = Path(__file__).with_name("review.template.html").read_text(encoding="utf-8")
    for marker, value in replacements.items():
        page = page.replace("{{" + marker + "}}", value)
    output.write_text(page, encoding="utf-8")


def make_png(wings, output):
    # Include the final millimetre at higher density without evaluating beyond
    # the tessellated tip, where the geometric and mesh extrema can differ.
    limit = min(wing.points[:, 0].max() for wing in wings)
    stations = np.unique(np.concatenate((
        np.arange(0, 528, 2), np.arange(528, 599, .5), np.linspace(599, limit - .0002, 70)
    )))
    outlines = [envelope(wing, stations) for wing in wings]
    figure = plt.figure(figsize=(15.6, 9), facecolor=PAPER, dpi=150)
    figure.text(.052, .947, "Wing B rounded tip", fontsize=31, color=INK, weight="bold", va="top")
    figure.text(.054, .888, "STEP geometry comparison", fontsize=14, color=MUTED, va="top")
    handles = [Line2D([0], [0], color=color, lw=2.7, ls=style)
               for color, style in zip(COLORS, ("--", "-"), strict=True)]
    figure.legend(handles, LABELS, loc="upper right", bbox_to_anchor=(.951, .928), ncol=2,
                  frameon=False, labelcolor=INK, fontsize=12, columnspacing=1.6)
    full_plan = figure.add_axes((.075, .475, .50, .337))
    full_front = figure.add_axes((.075, .16, .50, .185))
    tip_plan = figure.add_axes((.66, .513, .282, .285))
    lip_cut = figure.add_axes((.66, .14, .282, .22))
    axes = (
        (full_plan, "Full planform", "Aft (mm)"),
        (full_front, "Full front view", "Height (mm)"),
        (tip_plan, "Rounded planform detail", "Aft (mm)"),
        (lip_cut, f"Lip section · {DEFAULT_CUT} mm aft", "Height (mm)"),
    )
    for ax, title, ylabel in axes:
        ax.set_facecolor("white")
        ax.set_title(title, loc="left", fontsize=15, color=INK, weight="bold", pad=11)
        ax.set_xlabel("Outboard (mm)", fontsize=11, color=MUTED, labelpad=6)
        ax.set_ylabel(ylabel, fontsize=11, color=MUTED, labelpad=6)
        ax.tick_params(labelsize=10, color="#a1b2b0", labelcolor=MUTED, length=3)
        ax.set_axisbelow(True)
        ax.grid(color=GRID, linewidth=.8)
        for spine in ax.spines.values():
            spine.set_color("#c9d7d5")
        ax.set_aspect("equal", adjustable="box")
        ax.axvline(TRANSITION, color=GUIDE, linewidth=1.1, linestyle=(0, (3, 3)), zorder=4)

    def silhouette(ax, index):
        lows, highs = outlines[1][:, :, index].T
        ax.fill_between(stations, lows, highs, color=NEW, alpha=.105, zorder=2)
        for case in (1, 0):
            low, high = outlines[case][:, :, index].T
            x = np.concatenate((stations, stations[::-1], stations[:1]))
            y = np.concatenate((low, high[::-1], low[:1]))
            ax.plot(x, y, color=COLORS[case], linestyle="-" if case else "--",
                    linewidth=2 if case else 1.45, zorder=5 + (1 - case))

    silhouette(full_plan, 0)
    silhouette(full_front, 1)
    silhouette(tip_plan, 0)
    full_plan.set(xlim=(-8, 614), ylim=(211, -11))
    full_front.set(xlim=(-8, 614), ylim=(-9, 78), yticks=(0, 25, 50, 75))
    tip_plan.set(xlim=(521, 609), ylim=(159, 59), xticks=(525, 550, 575, 600))
    tip_plan.axhline(DEFAULT_CUT, color=GUIDE, linewidth=1.2, linestyle=":", zorder=4)
    tip_plan.annotate("SECTION", (525, DEFAULT_CUT), xytext=(2, -11), textcoords="offset points",
                      color=GUIDE, fontsize=7.5, weight="bold")
    for wing, color, style in zip(wings, COLORS, ("--", "-"), strict=True):
        for line in wing.profiles[DEFAULT_CUT]:
            lip_cut.plot(line[:, 0], line[:, 2], color=color, linestyle=style,
                         linewidth=2.2, solid_capstyle="round", dash_capstyle="round")
    x_range, y_range = section_ranges(wings, DEFAULT_CUT)
    lip_cut.set(xlim=x_range, ylim=y_range, xticks=(525, 550, 575, 600))
    figure.text(.055, .043, "Dashed brown line marks the start of the tip revision at 540 mm outboard.",
                fontsize=10.5, color=MUTED)
    figure.savefig(output, dpi=150)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", type=Path, default=BASELINE)
    parser.add_argument("--rounded", type=Path, default=REVISED)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    wings = (load_wing(args.original), load_wing(args.rounded))
    make_html(wings, args.output_dir / "review.html")
    make_png(wings, args.output_dir / "review.png")
    print(f"Review files written to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
