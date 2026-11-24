"""Render an offline, equal-axis-scale inspection of the Wing B spar evidence."""

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import cadquery as cq
import numpy as np
import plotly.graph_objects as go


HERE = Path(__file__).resolve().parent
BASELINE = HERE.parent / "wing-b-twist/reversal-study/measured/02-wing-occt.step"
COLORS = ("#286d9f", "#a55425")
CHANNEL = "#147c68"
WALL = "#a22c50"
REACH = "#a45c13"
DIAMETER = 6.4


def load_renderer():
    spec = importlib.util.spec_from_file_location(
        "wing_renderer", HERE.parent / "paired-wing-trial/render-wing-comparison.py"
    )
    renderer = importlib.util.module_from_spec(spec)
    previous = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec.loader.exec_module(renderer)
    finally:
        sys.dont_write_bytecode = previous
    return renderer


def display(points):
    return np.asarray(points)[..., [0, 2, 1]] * (-1, -1, 1)


def line3d(points, name, color, width=4, **kwargs):
    x, y, z = np.asarray(points).T
    return go.Scatter3d(x=x, y=y, z=z, name=name, showlegend=False,
                        line={"color": color, "width": width}, **kwargs)


def add_wing(wing, renderer, solid, label, color):
    points, faces, nose = renderer.mesh_wing(solid)
    x, y, z = points.T
    i, k, j = faces.T
    wing.add_trace(go.Mesh3d(
        x=x, y=y, z=z, i=i, j=j, k=k, name=label, color=color, opacity=.24,
        lighting=dict(ambient=.82, diffuse=.58, specular=.08, roughness=.75),
        hoverinfo="skip", showlegend=False,
    ))
    edges = list(nose)
    for edge in solid.Edges():
        edges.append((float("nan"),) * 3)
        edges.extend(renderer.display_xyz(edge.sample(241)[0]))
    wing.add_trace(line3d(edges, label, color, mode="lines", hoverinfo="skip"))
    return points.min(axis=0), points.max(axis=0)


def profile_at(solid, span):
    plane = cq.Face.makePlane(basePnt=(-span, 0, 0), dir=(1, 0, 0))
    faces = solid.intersect(plane).Faces()
    if len(faces) != 1:
        raise ValueError(f"Expected one CAD profile at {span} mm outboard")
    outline = []
    for edge in faces[0].outerWire().Edges():
        outline.extend(display([(p.x, p.y, p.z) for p in edge.sample(701)[0]]))
        outline.append((float("nan"),) * 3)
    return np.asarray(outline)


def circle_basis(axis):
    first = np.cross(axis, (0, 0, 1))
    first /= np.linalg.norm(first)
    second = np.cross(axis, first)
    theta = np.linspace(0, 2 * np.pi, 97)
    return DIAMETER / 2 * (np.cos(theta[:, None]) * first + np.sin(theta[:, None]) * second)


def add_channel(wing, result, label, visible):
    ends = display([result["axis_root_mm"], result["axis_channel_end_mm"]])
    axis = ends[1] - ends[0]
    axis /= np.linalg.norm(axis)
    circle = circle_basis(axis)
    # The channel opens at the root mounting plane. Its blind end is normal to the axis.
    root_ring = ends[0] + circle - circle[:, :1] / axis[0] * axis
    end_ring = ends[1] + circle
    width = len(circle)
    near = np.arange(width - 1)
    far = near + width
    points = np.concatenate((root_ring, end_ring, ends[1:]))
    triangles = np.concatenate((
        np.column_stack((near, near + 1, far)),
        np.column_stack((near + 1, far + 1, far)),
        np.column_stack((far, far + 1, np.full(width - 1, 2 * width))),
    ))
    x, y, z = points.T
    i, j, k = triangles.T
    wing.add_trace(go.Mesh3d(
        x=x, y=y, z=z, i=i, j=j, k=k, name=label, color=CHANNEL, opacity=.78,
        visible=visible, hoverinfo="skip", showlegend=False,
    ))
    rings = np.concatenate((root_ring, [[np.nan] * 3], end_ring))
    wing.add_trace(line3d(rings, label, CHANNEL, mode="lines",
                          visible=visible, hoverinfo="skip"))
    wing.add_trace(line3d(
        ends, "Straight spar axis", "#17382f", width=6, mode="lines+markers",
        marker={"size": 3, "color": "#17382f"}, visible=visible,
        hovertemplate="Axis<br>Outboard %{x:.2f}<br>Aft %{y:.2f}<br>Up %{z:.2f} mm<extra></extra>",
    ))
    wall = display([result["worst_channel_point_mm"], result["worst_skin_point_mm"]])
    wing.add_trace(line3d(
        wall, "Minimum wall", WALL, width=8, mode="lines+markers",
        marker={"size": 4, "color": WALL}, visible=visible,
        text=["Channel surface", "Wing skin"],
        hovertemplate="%{text}<br>Outboard %{x:.3f}<br>Aft %{y:.3f}<br>Up %{z:.3f} mm<extra></extra>",
    ))
    reach = np.array([display(result["axis_last_joint_mm"]), ends[1]])
    wing.add_trace(line3d(
        reach, "Tip engagement along axis", REACH, width=8, mode="lines+markers",
        marker={"size": 6, "color": REACH, "line": {"color": "#fff8e8", "width": 1.5}},
        visible=visible,
        text=["Tip joint", "Channel end"],
        hovertemplate="%{text}<br>%{x:.2f} mm outboard<extra></extra>",
    ))
    return ends, axis, circle, wall, reach


def section_channel(ends, axis, circle, span):
    points = ends[0] + circle
    along = (span - points[:, :1]) / axis[0]
    # Keep only the part of the oblique plane intersection inside the blind cylinder.
    total = np.linalg.norm(ends[1] - ends[0])
    intersection = points + along * axis
    intersection[along[:, 0] > total + 1e-8] = np.nan
    return intersection, ends[0] + ((span - ends[0, 0]) / axis[0]) * axis


def axis_layout(low, high):
    size = high - low
    return dict(
        aspectmode="manual", aspectratio=dict(zip("xyz", (2.2 * size / max(size)).tolist())),
        **{
            f"{axis}axis": dict(title=dict(text=f"{name} (mm)"), range=[float(a), float(b)],
                                showbackground=False, gridcolor="#dce5eb",
                                zeroline=False, showspikes=False)
            for axis, name, a, b in zip("xyz", ("Outboard", "Aft", "Up"), low, high)
        },
    )


def camera(x, y, z, up=(0, 0, 1)):
    return dict(eye=dict(zip("xyz", (x, y, z))), up=dict(zip("xyz", up)),
                projection=dict(type="orthographic"))


CAMERAS = {
    "Tip reach": camera(.08, -1.8, .64),
    "Orbit": camera(1.18, -1.62, 1.05),
    "Front": camera(0, -2.2, 0),
    "Top": camera(0, 0, 2.2, (0, 1, 0)),
    "Minimum wall": camera(.95, -1.46, 1.16),
}


def make_html(groups):
    renderer = load_renderer()
    wing, section = go.Figure(), go.Figure()
    bounds, scenarios, models, profiles = [], [], {}, {}
    for group_index, cases in enumerate(groups):
        records = []
        selected = len(cases) - 1
        scenario_settings = None
        for index, (path, result_path, variant) in enumerate(cases):
            result = json.loads(result_path.read_text(encoding="utf-8"))
            if path not in models:
                solids = cq.importers.importStep(str(path)).solids().vals()
                if len(solids) != 1 or not solids[0].isValid():
                    raise ValueError(f"Expected one valid wing solid in {path}")
                color = COLORS[index]
                wing_start = len(wing.data)
                bounds.append(add_wing(wing, renderer, solids[0], variant, color))
                models[path] = (solids[0], wing_start, color, hashlib.sha256(path.read_bytes()).hexdigest())
            solid, wing_start, color, digest = models[path]
            if digest != result["input_sha256"]:
                raise ValueError(f"Result {result_path} measured a different STEP file from {path}")
            show = group_index == 0 and index == selected
            channel_start = len(wing.data)
            channel_label = "Spar channel" if result["geometry_pass"] else "Best tested channel"
            ends, axis, circle, wall, reach = add_channel(wing, result, channel_label, show)
            cut_views = {}
            cuts = [-int(x) if x == int(x) else -x for x in result["cut_plane_x_mm"]]
            settings = dict(
                cuts=cuts,
                tip=result.get("tip_piece_length_mm", 600 - cuts[-1]),
                engagement=result["required_full_diameter_engagement_mm"],
                height=result.get("max_print_height_mm", 250),
            )
            if scenario_settings is not None and settings != scenario_settings:
                raise ValueError(f"Result {result_path} uses different settings from the other wing")
            scenario_settings = settings
            for span in cuts:
                key = (path, span)
                if key not in profiles:
                    profiles[key] = profile_at(solid, span)
                outline = profiles[key]
                profile_start = len(section.data)
                section.add_trace(go.Scatter(
                    x=outline[:, 1], y=outline[:, 2], mode="lines", name=variant,
                    line=dict(color=color, width=2.5, simplify=False), showlegend=False,
                    visible=group_index == 0 and span == cuts[0], hovertemplate="%{x:.2f} aft, %{y:.2f} up mm<extra></extra>",
                ))
                intersection, center = section_channel(ends, axis, circle, span)
                section.add_trace(go.Scatter(
                    x=intersection[:, 1], y=intersection[:, 2], mode="lines",
                    name=channel_label, line=dict(color=CHANNEL, width=3, simplify=False),
                    showlegend=False, visible=show and span == cuts[0], hoverinfo="skip",
                ))
                section.add_trace(go.Scatter(
                    x=[center[1]], y=[center[2]], mode="markers",
                    marker=dict(size=8, symbol="cross", color="#17382f"), name="Axis center",
                    showlegend=False, visible=show and span == cuts[0],
                    hovertemplate="Axis<br>%{x:.2f} aft, %{y:.2f} up mm<extra></extra>",
                ))
                cut_views[span] = dict(
                    trace=profile_start, center=center[1:].tolist(),
                    low=np.nanmin(outline[:, 1:], axis=0).tolist(),
                    high=np.nanmax(outline[:, 1:], axis=0).tolist(),
                )
            records.append(dict(
                variant=variant, color=color, wing=wing_start, channel=channel_start, cuts=cut_views,
                wall=float(result["min_wall_mm"]),
                wallSpan=float(wall[1, 0]), passGeometry=bool(result["geometry_pass"]),
                engagement=result["full_diameter_engagement_mm"],
                reach=reach.tolist(),
                wallView=axis_layout(wall.mean(axis=0) - 10, wall.mean(axis=0) + 10),
            ))
        scenarios.append(dict(records=records, **scenario_settings))
    # Make the close-up scale fixed for cases that share the same tip joint.
    for scenario in scenarios:
        relatives = [record["reach"] for other in scenarios
                     if other["cuts"][-1] == scenario["cuts"][-1] for record in other["records"]]
        points = np.asarray(relatives).reshape(-1, 3)
        pad = np.array([11., 16., 13.])
        scenario["reachView"] = axis_layout(points.min(axis=0) - pad, points.max(axis=0) + pad)
    low = np.min([item[0] for item in bounds], axis=0)
    high = np.max([item[1] for item in bounds], axis=0)
    planes = {}
    for span in sorted({span for scenario in scenarios for span in scenario["cuts"]}):
        planes[span] = len(wing.data)
        aft0, up0 = low[1:] - 4
        aft1, up1 = high[1:] + 4
        wing.add_trace(go.Mesh3d(
            x=[span] * 4, y=[aft0, aft1, aft1, aft0], z=[up0, up0, up1, up1],
            i=[0, 0], j=[1, 2], k=[2, 3], opacity=.13, color="#6758a5",
            name=f"{span} mm print cut", hoverinfo="name", showlegend=False,
        ))
        wing.add_trace(line3d(
            [(span, aft0, up0), (span, aft1, up0), (span, aft1, up1),
             (span, aft0, up1), (span, aft0, up0)],
            f"{span} mm print cut", "#7264a1", width=2, mode="lines", hoverinfo="name",
        ))
    whole_view = axis_layout(low - 10, high + 10)
    base = dict(paper_bgcolor="white", plot_bgcolor="white",
                font=dict(family="system-ui, sans-serif", color="#233749", size=12))
    wing.update_layout(
        **base, margin=dict(l=4, r=4, t=4, b=4), uirevision="spar-wing",
        scene=dict(**scenarios[0]["reachView"], dragmode="orbit", camera=CAMERAS["Tip reach"]),
    )
    section_axis = dict(showline=True, gridcolor="#dce5eb", zeroline=False, constrain="domain")
    section.update_layout(
        **base, margin=dict(l=60, r=20, t=12, b=50), dragmode="pan",
        xaxis=dict(**section_axis, title="Aft (mm)"),
        yaxis=dict(**section_axis, title="Up (mm)", scaleanchor="x", scaleratio=1),
    )
    config = dict(responsive=True, scrollZoom=True, displayModeBar=False, doubleClick=False)
    source = (HERE / "spar-viewer.template.html").read_text(encoding="utf-8")
    replacements = {
        "WING": wing.to_html(full_html=False, include_plotlyjs=True, div_id="wing-view", config=config),
        "SECTION": section.to_html(full_html=False, include_plotlyjs=False, div_id="section-view", config=config),
        "DATA": json.dumps(dict(scenarios=scenarios, cameras=CAMERAS,
                                wholeView=whole_view, planes=planes)),
    }
    for key, value in replacements.items():
        source = source.replace("{{" + key + "}}", value)
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=BASELINE)
    parser.add_argument("--result", type=Path, default=HERE / "approved/result.json")
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--candidate-result", type=Path)
    parser.add_argument(
        "--scenario", type=Path, action="append", nargs=2,
        metavar=("APPROVED_RESULT", "CANDIDATE_RESULT"),
        help="Additional measured pair, using the same two wing STEP files",
    )
    parser.add_argument("--output", type=Path, default=HERE / "wing-b-spar.html")
    args = parser.parse_args()
    if bool(args.candidate) != bool(args.candidate_result):
        parser.error("--candidate and --candidate-result must be supplied together")
    cases = [(args.baseline, args.result, "Approved Wing B")]
    if args.candidate:
        cases.append((args.candidate, args.candidate_result, "Constrained Wing B"))
    groups = [cases]
    for approved_result, candidate_result in args.scenario or ():
        if not args.candidate:
            parser.error("--scenario requires --candidate")
        groups.append([(args.baseline, approved_result, "Approved Wing B"),
                       (args.candidate, candidate_result, "Constrained Wing B")])
    args.output.write_text(make_html(groups), encoding="utf-8")
    print(args.output.resolve())


if __name__ == "__main__":
    main()
