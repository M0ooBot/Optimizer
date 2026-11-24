"""Phase 0 evidence from the saved geometry definition and optional STEP readback."""

from __future__ import annotations

from hashlib import file_digest
import json
from pathlib import Path

import numpy as np

from .phase0_geometry import build_curves, hinge_fit, hinge_points, projected_skin_span, section_frame, section_points


_OVERLAY_STATIONS = (0.75, 0.9)


def _stations(surface, subdivisions=4):
    ribs = sorted({0.0, 1.0, *surface.station_fractions})
    return np.concatenate(
        [np.linspace(a, b, subdivisions + 1)[:-1] for a, b in zip(ribs, ribs[1:])] + [[1.0]]
    )


def _samples(surface, baseline_dir, fractions):
    return np.asarray([section_points(surface, float(f), baseline_dir) for f in fractions])


def _smoothness(fractions, centers):
    """Finite-difference geometry of the sampled centerline in the X–Y view."""
    span = -centers[:, 0]
    slope = np.gradient(centers[:, 1], span, edge_order=2)
    angle = np.degrees(np.arctan(slope))
    curvature = np.gradient(slope, span, edge_order=2) / (1 + slope**2)**1.5
    k = int(np.argmax(np.abs(curvature)))
    decrease = int(np.argmin(np.diff(angle)))
    return {
        "minimum_dihedral_increment_deg": float(np.diff(angle)[decrease]),
        "dihedral_decrease_from_fraction": float(fractions[decrease]),
        "dihedral_decrease_to_fraction": float(fractions[decrease + 1]),
        "tip_dihedral_deg": float(angle[-1]),
        "minimum_sampled_bend_radius_mm": (
            1 / float(abs(curvature[k])) if abs(curvature[k]) > 1e-10 else None
        ),
        "peak_curvature_fraction": float(fractions[k]),
    }


def _generated(surface, baseline, baseline_dir, tolerance):
    stations = _stations(surface)
    sections = _samples(surface, baseline_dir, stations)
    frames = [section_frame(surface, float(f)) for f in stations]
    centers = np.asarray([f.quarter_chord for f in frames])
    le = np.asarray([f.le for f in frames])
    te = np.asarray([f.te for f in frames])
    root_displacement = np.linalg.norm(sections[0] - section_points(baseline, 0, baseline_dir), axis=1)
    physical_span = projected_skin_span(surface, baseline_dir)
    span_change = abs(physical_span - baseline.projected_span_mm)
    guide_angle = []
    plane_residual = []
    for frame, section in zip(frames, sections):
        _, _, vh = np.linalg.svd(section - section.mean(axis=0), full_matrices=False)
        normal = vh[-1]
        tangent = np.asarray(frame.span_tangent)
        tangent = tangent / np.linalg.norm(tangent)
        guide_angle.append(np.degrees(np.arccos(np.clip(abs(normal @ tangent), 0, 1))))
        plane_residual.append(np.max(np.abs((section - section.mean(axis=0)) @ normal)))
    worst = int(np.argmax(guide_angle))
    hinge_margin, hinge_station = hinge_fit(surface, baseline_dir)
    return {
        "root_profile_max_displacement_mm": float(root_displacement.max()),
        "projected_span_mm": surface.projected_span_mm,
        "projected_span_change_mm": span_change,
        "leading_edge_datum_projected_span_mm": float(le[0, 0] - le[-1, 0]),
        "baseline_root_attachment_within_tolerance": bool(root_displacement.max() <= tolerance),
        "baseline_projected_span_within_tolerance": bool(span_change <= tolerance),
        "hinge": {
            "minimum_sampled_section_height_margin_mm": hinge_margin,
            "limiting_fraction": hinge_station,
        },
        "leading_edge_tip_rise_mm": float(le[-1, 1] - le[0, 1]),
        "physical_projected_extent_mm": physical_span,
        "span_curve": surface.span_curve,
        "rise_onset_fraction": surface.rise_start,
        "rake_onset_fraction": surface.rake_start,
        "leading_edge_midspan_rise_mm": float(section_frame(surface, .5).le[1] - le[0, 1]),
        "leading_edge_tip_aft_mm": float(le[0, 2] - le[-1, 2]),
        "quarter_chord_tip_rise_mm": float(centers[-1, 1] - centers[0, 1]),
        "largest_guide_normal_mismatch_deg": float(guide_angle[worst]),
        "largest_guide_normal_mismatch_fraction": float(stations[worst]),
        "largest_section_plane_residual_mm": float(max(plane_residual)),
        "orientation_guide_tip_dihedral_deg": float(np.degrees(np.arctan2(
            frames[-1].span_tangent[1], abs(frames[-1].span_tangent[0])
        ))),
        "section_orientation": surface.section_orientation,
        "tip_section_front_tilt_deg": float(np.degrees(np.arctan2(
            abs(_normal(frames[-1])[1]), abs(_normal(frames[-1])[0])
        ))),
        "tip_twist_deg": surface.tip_twist_deg,
        "quarter_chord_smoothness": _smoothness(stations, centers),
        "leading_edge_smoothness": _smoothness(stations, le),
        "trailing_edge_smoothness": _smoothness(stations, te),
    }


def _profile_wire(cq, section):
    points = [cq.Vector(*point) for point in section]
    if (points[0] - points[-1]).Length > 1e-6:
        # Keep the declared sharp trailing edge when a profile has a finite opening.
        return cq.Wire.assembleEdges([cq.Edge.makeSpline(points), cq.Edge.makeLine(points[-1], points[0])])
    return cq.Wire.assembleEdges([cq.Edge.makeSpline(points)])


def _fallback_step(cq, surface, baseline_dir, path):
    # Add three derived profiles per rib gap to follow the continuous definition.
    # Eight measurement intervals per gap retain independent held-out stations.
    profiles = _samples(surface, baseline_dir, _stations(surface, 4))
    solid = cq.Solid.makeLoft([_profile_wire(cq, section) for section in profiles])
    if not solid.isValid():
        raise ValueError("OCCT loft produced an invalid solid")
    cq.exporters.export(solid, str(path))


def _load_solid(cq, path):
    solids = cq.importers.importStep(str(path)).solids().vals()
    if len(solids) != 1 or not solids[0].isValid():
        raise ValueError("STEP must contain one valid solid for this half-surface")
    return solids[0]


def _cut(cq, solid, frame):
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
    from OCP.gp import gp_Dir, gp_Pln, gp_Pnt

    plane = gp_Pln(gp_Pnt(*frame.le), gp_Dir(*_normal(frame)))
    section = BRepAlgoAPI_Section(solid.wrapped, plane, True)
    if not section.IsDone():
        return None
    return cq.Shape.cast(section.Shape())


def _normal(frame):
    return np.cross(frame.chord_axis, frame.up_axis)


def _cap(solid, frame):
    normal = _normal(frame)
    caps = [
        face for face in solid.Faces()
        if face.geomType() == "PLANE" and abs(np.dot(face.normalAt().toTuple(), normal)) > 0.98
    ]
    return min(caps, key=lambda face: abs(np.dot(np.asarray(face.Center().toTuple()) - frame.le, normal))) if caps else None


def _edges_extrema(edges, frame):
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from scipy.optimize import minimize_scalar

    chord_axis = np.asarray(frame.chord_axis)
    extremes = []
    for sign in [1, -1]:
        candidates = []
        for edge in edges:
            curve = BRepAdaptor_Curve(edge.wrapped)
            params = np.linspace(curve.FirstParameter(), curve.LastParameter(), 31)
            def position(parameter):
                p = curve.Value(parameter)
                return p.X(), p.Y(), p.Z()
            def objective(parameter):
                return sign * np.dot(position(parameter), chord_axis)
            best = int(np.argmin([objective(f) for f in params]))
            lo, hi = params[max(best - 1, 0)], params[min(best + 1, len(params) - 1)]
            optimum = minimize_scalar(objective, bounds=(lo, hi), method="bounded")
            at = min([lo, hi, params[best], optimum.x], key=objective)
            candidates.append((objective(at), position(at)))
        extremes.append(min(candidates)[1])
    return extremes


def _worst(errors, sections, stations, ribs):
    i, j = np.unravel_index(int(np.argmax(errors)), errors.shape)
    return {
        "distance_mm": float(errors[i, j]),
        "xyz_mm": sections[i, j].tolist(),
        "station_fraction": float(stations[i]),
        "between_defining_ribs": bool(not np.any(np.isclose(ribs, stations[i], atol=1e-9, rtol=0))),
    }


def _measure(cq, solid, surface, baseline, baseline_dir, tolerance, fallback=False):
    from OCP.BRepBndLib import BRepBndLib
    from OCP.Bnd import Bnd_Box

    stations = _stations(surface, 8)
    keep = ~np.isclose(stations[:, None], _OVERLAY_STATIONS, atol=1e-10, rtol=0).any(axis=1)
    stations = np.sort([*stations[keep], *_OVERLAY_STATIONS])
    sections = _samples(surface, baseline_dir, stations)
    ribs = sorted({0.0, 1.0, *surface.station_fractions})
    skin = cq.Compound.makeCompound(solid.Faces())
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(solid.wrapped, box, False, False)
    xmin, _, _, xmax, _, _ = box.Get()
    root_frame = section_frame(baseline, 0)
    tip_frame = section_frame(surface, 1)
    cap_sections = [_cap(solid, root_frame), _cap(solid, tip_frame)]
    cap_extrema = [_edges_extrema(cap.Edges(), frame) if cap else None
                   for cap, frame in zip(cap_sections, [root_frame, tip_frame])]
    span = cap_extrema[0][0][0] - cap_extrema[1][0][0] if all(cap_extrema) else None
    root_plane_error = abs(np.dot(
        np.asarray(cap_sections[0].Center().toTuple()) - root_frame.le, _normal(root_frame)
    )) if cap_sections[0] else None
    root_samples = np.asarray(section_points(baseline, 0, baseline_dir))
    cap_cuts = [cq.Compound.makeCompound(cap.Edges()) if cap else None for cap in cap_sections]
    root_cut = cap_cuts[0]
    root_error = None
    if root_cut is not None and root_cut.Edges():
        root_wire = _profile_wire(cq, root_samples)
        expected_to_cut = [root_cut.distance(cq.Vertex.makeVertex(*p)) for p in root_samples]
        cut_to_expected = [
            root_wire.distance(cq.Vertex.makeVertex(*p.toTuple()))
            for edge in root_cut.Edges() for p in edge.sample(61)[0]
        ]
        root_error = max(*expected_to_cut, *cut_to_expected)
    leading_edge_span_error = abs(span - (root_frame.le[0] - tip_frame.le[0])) if span is not None else None

    # Plane cuts at intermediate stations also support the reverse comparison.
    cut_samples = {}
    cut_deviation = []
    cap_deviation = []
    measured_le = []
    measured_te = []
    sample_indices = np.rint(np.linspace(0, sections.shape[1] - 1, min(81, sections.shape[1]))).astype(int)
    forward_sections = sections[:, sample_indices]
    errors = np.full(forward_sections.shape[:2], np.nan)
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape

    for i, (station, profile) in enumerate(zip(stations, sections)):
        frame = section_frame(surface, float(station))
        end_index = i in (0, len(stations) - 1)
        cut = cap_cuts[0 if i == 0 else 1] if end_index else _cut(cq, solid, frame)
        if cut is None or not cut.Edges():
            continue
        errors[i] = [cut.distance(cq.Vertex.makeVertex(*p)) for p in forward_sections[i]]
        wire = _profile_wire(cq, profile)
        samples = []
        for edge in cut.Edges():
            points = np.asarray([p.toTuple() for p in edge.sample(61)[0]])
            samples.append(points)
            for point in points:
                distance = wire.distance(cq.Vertex.makeVertex(*point))
                (cap_deviation if end_index else cut_deviation).append(
                    (float(distance), float(station), point.tolist())
                )
        if end_index:
            continue
        actual_le, actual_te = _edges_extrema(cut.Edges(), frame)
        measured_le.append(actual_le)
        measured_te.append(actual_te)
        cut_samples[float(station)] = samples

    cuts_present = np.isfinite(errors[:, 0])
    between_ribs = ~np.isclose(stations[:, None], ribs, atol=1e-10, rtol=0).any(axis=1) & cuts_present
    fallback_ribs = _stations(surface, 4)
    independent_cuts = (
        ~np.isclose(stations[:, None], fallback_ribs, atol=1e-10, rtol=0).any(axis=1)
        & cuts_present
    )
    complete_cuts = bool(cuts_present.all())
    worst = _worst(errors[cuts_present], forward_sections[cuts_present], stations[cuts_present], ribs) if cuts_present.any() else None
    reverse = max(cut_deviation, default=None)
    reverse_cap = max(cap_deviation, default=None)
    hinge = None
    if surface.hinge is not None:
        start, end = np.asarray(hinge_points(surface, baseline_dir))
        axis = cq.Edge.makeLine(cq.Vector(*start), cq.Vector(*end))
        extrema = BRepExtrema_DistShapeShape(axis.wrapped, skin.wrapped)
        if not extrema.IsDone() or not extrema.NbSolution():
            raise RuntimeError("Could not measure the hinge axis against the STEP skin")
        point = extrema.PointOnShape1(1)
        hinge = {
            "start_station_fraction": surface.hinge.start,
            "end_station_fraction": surface.hinge.end,
            "inboard_chord_fraction": surface.hinge.inboard_chord_fraction,
            "outboard_chord_fraction": surface.hinge.outboard_chord_fraction,
            "axis_endpoints_xyz_mm": [start.tolist(), end.tolist()],
            "axis_length_mm": float(np.linalg.norm(end - start)),
            "axis_inside_solid": bool(
                all(solid.isInside(cq.Vector(*p), 1e-7) for p in [start, (start + end) / 2, end])
                and extrema.Value() > 1e-7
            ),
            "minimum_axis_skin_distance_mm": float(extrema.Value()),
            "minimum_axis_skin_distance_xyz_mm": [point.X(), point.Y(), point.Z()],
            "control_travel_clearance": "unsupported on the continuous solid",
        }
    center_smoothness = None
    if len(cut_samples) == len(stations) - 2:
        centers = 0.75 * np.asarray(measured_le) + 0.25 * np.asarray(measured_te)
        center_smoothness = _smoothness(stations[1:-1], centers)

    mid_frame = section_frame(surface, .5)
    mid_cut = _cut(cq, solid, mid_frame)
    if mid_cut is None or not mid_cut.Edges():
        mid_cut_le = None
    else:
        mid_cut_le, _ = _edges_extrema(mid_cut.Edges(), mid_frame)

    return {
        "valid_solid": True,
        "tolerance_mm": tolerance,
        "sample_count": int(np.isfinite(errors).sum()),
        "station_fractions": stations.tolist(),
        "defining_rib_fractions": ribs,
        "maximum_definition_to_step_deviation": worst,
        "maximum_between_ribs_definition_to_step_deviation": _worst(
            errors[between_ribs], forward_sections[between_ribs], stations[between_ribs], ribs
        ) if between_ribs.any() else None,
        "fallback_loft_input_fractions": fallback_ribs.tolist() if fallback else None,
        "maximum_held_out_definition_to_step_deviation": _worst(
            errors[independent_cuts], forward_sections[independent_cuts], stations[independent_cuts], fallback_ribs
        ) if fallback and independent_cuts.any() else None,
        "maximum_intermediate_cut_to_definition_deviation": (
            {"distance_mm": reverse[0], "station_fraction": reverse[1], "xyz_mm": reverse[2]}
            if reverse else None
        ),
        "maximum_end_cap_perimeter_to_definition_deviation": (
            {"distance_mm": reverse_cap[0], "station_fraction": reverse_cap[1], "xyz_mm": reverse_cap[2]}
            if reverse_cap else None
        ),
        "sampled_surface_within_tolerance": bool(
            complete_cuts and worst and reverse and reverse_cap
            and max(worst["distance_mm"], reverse[0], reverse_cap[0]) <= tolerance
        ),
        "leading_edge_projected_span_mm": float(span) if span is not None else None,
        "bounding_projected_extent_mm": float(xmax - xmin),
        "baseline_bounding_projected_extent_within_tolerance": bool(
            abs(xmax - xmin - baseline.projected_span_mm) <= tolerance
        ),
        "leading_edge_tip_rise_mm": (
            float(cap_extrema[1][0][1] - cap_extrema[0][0][1]) if all(cap_extrema) else None
        ),
        "leading_edge_midspan_rise_mm": (
            float(mid_cut_le[1] - cap_extrema[0][0][1])
            if mid_cut_le is not None and cap_extrema[0] is not None else None
        ),
        "leading_edge_tip_aft_mm": (
            float(cap_extrema[0][0][2] - cap_extrema[1][0][2]) if all(cap_extrema) else None
        ),
        "leading_edge_projected_span_error_mm": float(leading_edge_span_error) if leading_edge_span_error is not None else None,
        "root_plane_displacement_mm": float(root_plane_error) if root_plane_error is not None else None,
        "root_profile_max_sampled_displacement_mm": root_error,
        "baseline_root_attachment_within_tolerance": bool(
            root_error is not None and root_plane_error is not None
            and max(root_error, root_plane_error) <= tolerance
        ),
        "leading_edge_projected_span_within_tolerance": bool(leading_edge_span_error is not None and leading_edge_span_error <= tolerance),
        "native_internal_section_orientation": "unsupported from STEP topology alone",
        "tip_cap_front_tilt_deg": (
            float(np.degrees(np.arctan2(
                abs(cap_sections[1].normalAt().y), abs(cap_sections[1].normalAt().x)
            ))) if cap_sections[1] else None
        ),
        "end_cap_plane_deviation_deg": [
            float(np.degrees(np.arccos(np.clip(abs(np.dot(cap.normalAt().toTuple(), _normal(frame))), 0, 1))))
            if cap else None for cap, frame in zip(cap_sections, [root_frame, tip_frame])
        ],
        "actual_smoothness": {
            "skin_continuity": "differential continuity unverified; see the sampled cut deviations",
            "method": "finite differences of interior section cuts using installed chord extrema",
            "cuts_expected": len(stations) - 2,
            "cuts_measured": len(cut_samples),
            "quarter_chord": center_smoothness,
            "leading_edge": _smoothness(stations[1:-1], np.asarray(measured_le)) if center_smoothness is not None else None,
            "trailing_edge": _smoothness(stations[1:-1], np.asarray(measured_te)) if center_smoothness is not None else None,
        },
        "hinge": hinge,
    }, {"sections": cut_samples, "le": measured_le, "te": measured_te}


def _overlay(outdir, revisions, baseline_dir, measured):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D

    fig = Figure(figsize=(16, 9), layout="constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(2, 4, width_ratios=[1.2, 1.2, 1, 1])
    colors = ("#667085", "#176b92", "#c15b28", "#267b52", "#7751a2")
    curves_by_label = {label: build_curves(definition, baseline_dir) for label, definition in revisions.items()}
    for j, surface_name in enumerate(("wing", "tail")):
        plan, front, *section_axes = axes[j]
        for i, (label, definition) in enumerate(revisions.items()):
            color = colors[i % len(colors)]
            surface = getattr(definition, surface_name)
            if surface is None:
                continue
            stations = _stations(surface)
            frames = [section_frame(surface, float(f)) for f in stations]
            le = np.asarray([f.le for f in frames])
            te = np.asarray([f.te for f in frames])
            quarter = np.asarray([f.quarter_chord for f in frames])
            plan.plot(-le[:, 0], le[:, 2], color=color, label=label)
            plan.plot(-te[:, 0], te[:, 2], color=color)
            plan.plot(-quarter[:, 0], quarter[:, 2], color=color, alpha=.35, linewidth=.8)
            front.plot(-quarter[:, 0], quarter[:, 1], color=color)
            readback = measured.get((label, surface_name))
            if readback and readback["le"]:
                actual_le = np.asarray(readback["le"])
                actual_te = np.asarray(readback["te"])
                actual_quarter = .75 * actual_le + .25 * actual_te
                for edge in (actual_le, actual_te):
                    plan.plot(-edge[:, 0], edge[:, 2], color=color, linestyle="--", linewidth=1)
                front.plot(
                    -actual_quarter[:, 0], actual_quarter[:, 1], color=color, linestyle="--", linewidth=1
                )
            if surface.hinge is not None:
                hinge = np.asarray(hinge_points(surface, baseline_dir))
                plan.plot(-hinge[:, 0], hinge[:, 2], color=color, linestyle=":")
                front.plot(-hinge[:, 0], hinge[:, 1], color=color, linestyle=":")
                outline = np.asarray(curves_by_label[label].get(f"{surface_name}_control_upper", []))
                if outline.size:
                    plan.plot(-outline[:, 0], outline[:, 2], color=color, alpha=.55, linewidth=.8)
            for target, section_ax in zip(_OVERLAY_STATIONS, section_axes):
                profile = np.asarray(section_points(surface, target, baseline_dir))
                section_ax.plot(profile[:, 2], profile[:, 1], color=color)
                if readback:
                    for points in readback["sections"].get(target, []):
                        section_ax.plot(points[:, 2], points[:, 1], color=color, linestyle="--", linewidth=1)
        plan.set(title=f"{surface_name.title()} plan", xlabel="Outboard (mm)", ylabel="Noseward +Z (mm)")
        front.set(title=f"{surface_name.title()} front, quarter chord", xlabel="Outboard (mm)", ylabel="Up +Y (mm)")
        for target, section_ax in zip(_OVERLAY_STATIONS, section_axes):
            section_ax.set(
                title=f"Section at {target:.0%} span",
                xlabel="Projected noseward +Z (mm)", ylabel="Projected up +Y (mm)",
            )
        for ax in (plan, front, *section_axes):
            ax.set_aspect("equal", adjustable="datalim")
            ax.grid(alpha=.25)
        handles, labels = plan.get_legend_handles_labels()
        if handles:
            plan.legend(handles, labels, loc="best")
    fig.legend(
        handles=[Line2D([], [], color="black", linestyle=style, label=label) for style, label in (
            ("-", "Generated definition"), ("--", "STEP readback"), (":", "Proposed hinge axis")
        )], loc="outside lower center", ncols=3
    )
    path = outdir / "overlays.png"
    fig.savefig(path, dpi=150)
    return path


def _step_pair(cq, a, b, tolerance, fractions, surface):
    """Symmetric sampled distance at matching STEP section planes and end caps."""
    if a.samefile(b):
        return {
            "status": "shared_step_input",
            "sampled_reproduction_within_tolerance": False,
        }
    solids = [_load_solid(cq, p) for p in [a, b]]
    paired_cuts = []
    for station in fractions:
        frame = section_frame(surface, float(station))
        if station in (0, 1):
            caps = [_cap(solid, frame) for solid in solids]
            paired_cuts.append([cq.Compound.makeCompound(cap.Edges()) if cap else None for cap in caps])
        else:
            paired_cuts.append([_cut(cq, solid, frame) for solid in solids])
    directions = []
    for i in range(2):
        worst = None
        n = 0
        missing_cuts = []
        for station, cuts in zip(fractions, paired_cuts):
            cut, other = cuts[i], cuts[1 - i]
            if any(part is None or not part.Edges() for part in cuts):
                missing_cuts.append(float(station))
                continue
            for edge in cut.Edges():
                for p in edge.sample(31)[0]:
                    point = p.toTuple()
                    distance = other.distance(cq.Vertex.makeVertex(*point))
                    sample = (float(distance), float(station), list(point))
                    worst = max(worst, sample) if worst else sample
                    n += 1
        directions.append({
            "sample_count": n,
            "empty_cuts_at_fractions": missing_cuts,
            "worst": {"distance_mm": worst[0], "station_fraction": worst[1], "xyz_mm": worst[2]} if worst else None,
        })
    return {
        "status": "measured",
        "original_to_rebuild": directions[0],
        "rebuild_to_original": directions[1],
        "tolerance_mm": tolerance,
        "sampled_reproduction_within_tolerance": all(
            d["worst"] and not d["empty_cuts_at_fractions"]
            and d["worst"]["distance_mm"] <= tolerance for d in directions
        ),
    }


def write_evidence(
    outdir, revisions, baseline_dir, step_paths=None, fallback=False, *, step_sources=None,
    rebuild_of=None,
):
    """Write JSON and PNG evidence for named definitions in comparison order.

    CAD inputs are keyed by (revision label, "wing" or "tail"). Sources are
    declared by the caller. OCCT fallback files use their own explicit source.
    The first revision supplies the locked baseline attachments.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if not revisions:
        raise ValueError("provide a baseline definition before candidate definitions")
    step_paths = {k: Path(v) for k, v in (step_paths or {}).items()}
    step_sources = step_sources or {}
    cq = None
    cad_errors = (ImportError, OSError, RuntimeError, ValueError)
    if step_paths or fallback:
        try:
            import cadquery as cq
            from OCP.Standard import Standard_Failure
        except ImportError:
            cq = None
        else:
            cad_errors += (Standard_Failure,)

    baseline_label, baseline = next(iter(revisions.items()))
    report = {
        "baseline_revision": baseline_label,
        "frame": {"origin": "wing root leading edge", "outboard": "-X", "up": "+Y", "nose": "+Z"},
        "units": "mm, degrees",
        "method": (
            "Eight readback station intervals per defining rib gap plus common overlay stations. "
            "Forward distances use up to 81 definition profile points against the STEP section curve. "
            "Reverse profile deviation samples 61 points per edge of each CAD cut and end cap. "
            "Sampling does not bound unsampled deviations or prove differential continuity."
        ),
        "revisions": {},
        "rebuild_comparisons": {},
        "revision_changes": {},
    }
    measured = {}
    available_steps = {}
    for i, (label, definition) in enumerate(revisions.items()):
        entry = {}
        report["revisions"][label] = entry
        for surface_name in ("wing", "tail"):
            surface = getattr(definition, surface_name)
            if surface is None:
                entry[surface_name] = {"status": "undefined"}
                continue
            reference = getattr(baseline, surface_name)
            if reference is None:
                raise ValueError(f"no baseline definition for {surface_name}")
            tolerance = definition.tolerance_mm
            generated = _generated(surface, reference, baseline_dir, tolerance)
            generated["agreed_locks_within_tolerance"] = (
                generated["baseline_root_attachment_within_tolerance"]
                and (surface_name == "tail" or (
                    generated["baseline_projected_span_within_tolerance"]
                    and abs(surface.projected_span_mm - reference.projected_span_mm) <= tolerance
                ))
            )
            evidence = {
                "generated": generated,
                "step_readback": {"status": "unavailable", "reason": "no STEP supplied"},
            }
            entry[surface_name] = evidence
            key = (label, surface_name)
            path = step_paths.get(key)
            source = step_sources.get(key, "supplied_step")
            if path is None and fallback:
                path = outdir / f"{i:02d}-{surface_name}-occt.step"
                source = "occt_fallback"
            if path is None:
                continue
            try:
                if cq is None:
                    raise ImportError("STEP readback requires the optional CAD dependencies")
                if source == "occt_fallback" and key not in step_paths:
                    _fallback_step(cq, surface, baseline_dir, path)
                solid = _load_solid(cq, path)
                readback, measured[key] = _measure(
                    cq, solid, surface, reference, baseline_dir, tolerance, fallback=source == "occt_fallback"
                )
                readback["agreed_locks_within_tolerance"] = (
                    readback["baseline_root_attachment_within_tolerance"]
                    and (surface_name == "tail" or (
                        readback["leading_edge_projected_span_within_tolerance"]
                        and readback["baseline_bounding_projected_extent_within_tolerance"]
                    ))
                )
                with path.open("rb") as file:
                    digest = file_digest(file, "sha256").hexdigest()
                evidence["step_readback"] = {
                    "status": "measured", "source": source, "path": str(path), "sha256": digest,
                    **readback,
                }
                available_steps[key] = path
            except cad_errors as exc:
                evidence["step_readback"] = {
                    "status": "unavailable", "source": source, "path": str(path), "reason": str(exc)
                }

    labels = [label for label in revisions if label not in (rebuild_of or {})]
    for old_label, new_label in zip(labels, labels[1:]):
        old, new = revisions[old_label], revisions[new_label]
        old_onset, new_onset = old.wing.rise_start, new.wing.rise_start
        change = {
            "previous": old_label,
            "wing_leading_edge_rise_change_mm": new.wing.rise_mm - old.wing.rise_mm,
            "wing_span_curve_previous": old.wing.span_curve,
            "wing_span_curve": new.wing.span_curve,
            "wing_rise_onset_change_fraction": (
                new_onset - old_onset if new_onset is not None and old_onset is not None else None
            ),
            "tail_saved_definition_matches": old.tail == new.tail,
        }
        old_wing = report["revisions"][old_label]["wing"]["step_readback"].get("leading_edge_tip_rise_mm")
        new_wing = report["revisions"][new_label]["wing"]["step_readback"].get("leading_edge_tip_rise_mm")
        if old_wing is not None and new_wing is not None:
            change["measured_wing_leading_edge_tip_rise_change_mm"] = new_wing - old_wing
            change["measured_wing_leading_edge_tip_rise_unchanged_within_tolerance"] = (
                abs(new_wing - old_wing) <= new.tolerance_mm
            )
        old_tail = available_steps.get((old_label, "tail"))
        new_tail = available_steps.get((new_label, "tail"))
        if change["tail_saved_definition_matches"] and new_tail and old_tail:
            change["measured_unchanged_tail"] = _step_pair(
                cq, old_tail, new_tail, new.tolerance_mm, _stations(new.tail), new.tail
            )
        report["revision_changes"][new_label] = change

    for rebuild_label, original_label in (rebuild_of or {}).items():
        if rebuild_label not in revisions or original_label not in revisions:
            raise ValueError("rebuild comparisons require both definitions in revisions")
        definition = revisions[rebuild_label]
        original_definition = revisions[original_label]
        entry = {}
        report["rebuild_comparisons"][rebuild_label] = {"original": original_label, "surfaces": entry}
        for name in ("wing", "tail"):
            surface = getattr(definition, name)
            if surface is None:
                continue
            original_surface = getattr(original_definition, name)
            path = available_steps.get((rebuild_label, name))
            original = available_steps.get((original_label, name))
            entry[name] = {"saved_definition_matches": surface == original_surface}
            if path and original:
                entry[name]["measured"] = _step_pair(
                    cq, original, path, definition.tolerance_mm, _stations(surface), surface
                )

    report["overlay"] = str(_overlay(outdir, revisions, baseline_dir, measured))
    path = outdir / "evidence.json"
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report
