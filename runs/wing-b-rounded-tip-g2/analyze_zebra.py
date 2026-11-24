"""Compare the inline Wing B viewer's analytical zebra field across the STEP seam.

Run with Python 3.12+ and CadQuery, NumPy, and SciPy in the environment.
"""

import argparse
import json
from dataclasses import dataclass
from math import cos, degrees, pi, radians, sin
from pathlib import Path

import cadquery as cq
import numpy as np
from scipy.optimize import brentq

from check_g2 import (AFT, CUT, HERE, SPAN, Boundary, Jet, Surface, angle_deg,
                      attached_faces, box_of, digest_of, edge_curve, norm,
                      to_design, unit)


PREVIOUS = HERE.parent / "wing-b-rounded-tip" / "wing-b-rounded.step"
REVISED = HERE / "wing-b-rounded.step"
PRESETS = {"Angled": (pi, 1.10), "Top": (0.0, pi / 2), "Front": (pi, 0.0)}
CHORD_TRIM = 0.001
SAMPLE_COUNT = 1001
AFT_SECTION = 150.0
CONDITION_MIN = 1e-7
DIRECTION_MIN_CYCLES_PER_MM = 1e-4
TRACE_SPAN_OFFSETS = (-10, -5, -2, -1, -.25, 0, .25, 1, 2, 5, 10, 15)


@dataclass
class View:
    label: str
    azimuth: float
    elevation: float
    stripe_deg: float
    side: int

    def vectors(self):
        az, el, stripe = self.azimuth, self.elevation, radians(self.stripe_deg)
        right = np.array([cos(az), 0.0, -sin(az)])
        eye = np.array([cos(el) * sin(az), sin(el), cos(el) * cos(az)])
        up = np.cross(eye, right)
        stripe_axis = cos(stripe) * right + sin(stripe) * up
        # The shader reflects each half's eye x only. Normals remain in the
        # positive-span mesh frame, and stripeAxis remains in the camera frame.
        return eye * [self.side, 1, 1], stripe_axis

    def describe(self):
        eye, axis = self.vectors()
        return dict(view=self.label, azimuth_deg=degrees(self.azimuth),
                    elevation_deg=degrees(self.elevation), stripe_deg=self.stripe_deg,
                    displayed_half="starboard" if self.side == 1 else "port",
                    diagnostic_eye=eye.tolist(), stripe_axis=axis.tolist())


def make_views():
    views = [View("Angled", *PRESETS["Angled"], 51, 1)]
    cameras = [(name, *angles) for name, angles in PRESETS.items()]
    cameras.extend(("Orbit", radians(az), radians(el)) for az in range(0, 360, 45)
                   for el in (-60, 0, 60) if (az, el) != (180, 0))
    views.extend(View(name, az, el, stripe, side)
                 for name, az, el in cameras for stripe in (0, 30, 51, 60, 90, 120, 150)
                 for side in (1, -1) if (name, stripe, side) != ("Angled", 51, 1))
    return views


def normal_terms(normal, eyes, axes):
    n_axis = axes @ normal
    n_eye = eyes @ normal
    # GLSL reflect(-eye, n) = -eye + 2 * dot(eye, n) * n.
    # Shader phase is 6 * dot(reflected, axis), with sin(pi * phase).
    cycles = 3 * (2 * n_eye * n_axis - np.einsum("ij,ij->i", eyes, axes))
    normal_gradient = 6 * (n_axis[:, None] * eyes + n_eye[:, None] * axes)
    return cycles, normal_gradient


def field(jet, frame, eyes, axes):
    """Return unwrapped cycles and d(cycles)/d(distance) in the supplied unit frame."""
    cycles, normal_gradient = normal_terms(jet.normal, eyes, axes)
    # S has the convention dn/ds = -frame @ S.
    gradient = -(normal_gradient @ frame) @ jet.shape_operator(frame)
    return cycles, gradient


def transport(frame, source, target):
    """Rotate an orthonormal frame by the shortest rotation source -> target."""
    axis = np.cross(source, target)
    cross_frame = np.cross(axis, frame.T).T
    return frame + cross_frame + np.cross(axis, cross_frame.T).T / (1 + source @ target)


def join_frames(body, cap):
    # Average the tangent planes with consistent normal hemispheres. Compare
    # the two gradients by rotationally transporting them into this plane.
    cap_normal = cap.normal if body.normal @ cap.normal >= 0 else -cap.normal
    normal = unit(body.normal + cap_normal)
    across = unit(SPAN - (SPAN @ normal) * normal)
    along = unit(np.cross(across, normal))
    if along @ AFT < 0:
        along = -along
    common = np.column_stack((across, along))
    return (transport(common, normal, body.normal),
            transport(common, normal, cap_normal))


@dataclass
class Join:
    body: Boundary
    cap: Boundary
    nose: float
    tail: float

    def read_pairs(self, aft):
        pairs = [(boundary.surface, boundary.at_aft(aft)) for boundary in (self.body, self.cap)]
        for surface, uv in pairs:
            if surface.condition_ratio(uv) < CONDITION_MIN:
                raise ValueError(f"Ill-conditioned STEP derivatives at {aft:.9f} mm aft")
        return pairs

    def jets(self, aft):
        return [surface.jet(uv) for surface, uv in self.read_pairs(aft)]


def read_joins(path):
    digest = digest_of(path)
    shape = cq.importers.importStep(str(path)).val()
    if digest_of(path) != digest:
        raise RuntimeError(f"STEP changed during import: {path}")
    faces = shape.Faces()
    cap = [face for face in faces if box_of(face)[0] < -CUT - .001]
    body = [face for face in faces if face not in cap and face.geomType() != "PLANE"]
    joins = {}
    for edge in shape.Edges():
        bounds = box_of(edge)
        if max(abs(bounds[i] + CUT) for i in (0, 3)) >= .001:
            continue
        inner, outer = attached_faces(edge, body), attached_faces(edge, cap)
        if len(inner) != 1 or len(outer) != 1:
            continue
        curve, start, end = edge_curve(edge)
        nose, tail = sorted(to_design(curve.Value(p))[2] for p in (start, end))
        inside, outside = (Boundary(edge, Surface(faces[0])) for faces in (inner, outer))
        side = "upper" if inside.surface.jet(inside.at_aft((nose + tail) / 2)).normal[1] > 0 else "lower"
        if side in joins:
            raise ValueError(f"Multiple {side} body/cap seams in {path}")
        joins[side] = Join(inside, outside, nose, tail)
    if set(joins) != {"upper", "lower"}:
        raise ValueError(f"Expected upper and lower body/cap seams in {path}")
    return joins, digest


def gradient_angle(body, cap):
    magnitude = np.stack((np.linalg.norm(body, axis=1),
                          np.linalg.norm(cap, axis=1)), axis=1)
    cross = body[:, 0] * cap[:, 1] - body[:, 1] * cap[:, 0]
    dot = np.einsum("ij,ij->i", body, cap)
    angle = np.degrees(np.arctan2(np.abs(cross), np.abs(dot)))
    angle[np.min(magnitude, axis=1) < DIRECTION_MIN_CYCLES_PER_MM] = np.nan
    return magnitude, angle


def read_field(joins, stations, eyes, axes):
    result = {}
    for side, join in joins.items():
        rows = []
        for aft in stations:
            body, cap = join.jets(aft)
            body_frame, cap_frame = join_frames(body, cap)
            body_phase, body_gradient = field(body, body_frame, eyes, axes)
            cap_phase, cap_gradient = field(cap, cap_frame, eyes, axes)
            magnitude, angle = gradient_angle(body_gradient, cap_gradient)
            rows.append({
                "phase_gap_cycles": np.abs(cap_phase - body_phase),
                "gradient_gap_cycles_per_mm": np.linalg.norm(cap_gradient - body_gradient, axis=1),
                "level_tangent_angle_deg": angle,
                "front_facing": (eyes @ body.normal > 0) & (eyes @ cap.normal > 0),
                "body_cycles": body_phase,
                "cap_cycles": cap_phase,
                "body_gradient_cycles_per_mm": body_gradient,
                "cap_gradient_cycles_per_mm": cap_gradient,
                "gradient_magnitudes_cycles_per_mm": magnitude,
                "seam_gap_mm": norm(body.position - cap.position),
                "normal_angle_deg": angle_deg(body.normal, cap.normal),
            })
        result[side] = {key: np.array([row[key] for row in rows]) for key in rows[0]}
    return result


METRICS = ("phase_gap_cycles", "gradient_gap_cycles_per_mm", "level_tangent_angle_deg")


def metrics_for(rows, stations, views, column_slice, front_facing=False):
    report = {}
    for key in METRICS:
        values = np.stack([row[key][:, column_slice] for row in rows.values()])
        valid = np.isfinite(values)
        if front_facing:
            valid &= np.stack([row["front_facing"][:, column_slice] for row in rows.values()])
        samples = values[valid]
        if not samples.size:
            report[key] = dict(sample_count=0, p95=None, max=None, worst_at=None)
            continue
        i_side, i_station, i_view = np.unravel_index(np.argmax(np.where(valid, values, -np.inf)), values.shape)
        worst_at = dict(skin=list(rows)[i_side], aft_mm=float(stations[i_station]),
                        **views[column_slice][i_view].describe())
        report[key] = dict(sample_count=samples.size, p95=float(np.percentile(samples, 95)),
                           max=float(values[i_side, i_station, i_view]), worst_at=worst_at)
    return report


def summarize(rows, stations, views, column_slice):
    groups = {**{side: {side: values} for side, values in rows.items()}, "combined": rows}
    report = {}
    for label, group in groups.items():
        seam = {side: {key: values[:-1] for key, values in items.items()} for side, items in group.items()}
        section = {side: {key: values[-1:] for key, values in items.items()} for side, items in group.items()}
        report[label] = dict(
            seam=metrics_for(seam, stations[:-1], views, column_slice),
            at_150_mm=metrics_for(section, stations[-1:], views, column_slice),
            front_facing_seam=metrics_for(seam, stations[:-1], views, column_slice, front_facing=True),
            front_facing_at_150_mm=metrics_for(section, stations[-1:], views, column_slice, front_facing=True),
        )
    return report


def default_profiles(rows, stations):
    profiles = {}
    for side, values in rows.items():
        data = {}
        for key, value in values.items():
            raw = value if value.ndim == 1 else value[:, 0]
            data[key] = np.where(np.isfinite(raw), raw.astype(object), None).tolist()
        profiles[side] = {"aft_mm": stations.tolist(), **data}
    return profiles


def stripe_crossings(join, stations, body_cycles, cap_cycles, eyes, axes):
    """Solve default dark/light transition crossings on each true STEP pcurve."""
    extrema = np.min((body_cycles, cap_cycles)), np.max((body_cycles, cap_cycles))
    half_cycles = range(int(np.ceil(2 * extrema[0])), int(np.floor(2 * extrema[1])) + 1)

    def roots(boundary, sampled, level):
        def phase_error(aft):
            jet = boundary.surface.jet(boundary.at_aft(aft))
            return float(normal_terms(jet.normal, eyes, axes)[0][0]) - level

        brackets = np.flatnonzero((sampled[:-1] > level) != (sampled[1:] > level))
        return [(brentq(phase_error, stations[i], stations[i + 1]), bool(sampled[i + 1] > level))
                for i in brackets]

    transitions = []
    unmatched = []
    for half_cycle in half_cycles:
        level = half_cycle / 2
        body_roots = roots(join.body, body_cycles, level)
        cap_roots = roots(join.cap, cap_cycles, level)
        if [up for _, up in body_roots] != [up for _, up in cap_roots]:
            unmatched.append(dict(level_cycles=level, body_roots=body_roots, cap_roots=cap_roots))
            continue
        for (body_aft, rising), (cap_aft, _) in zip(body_roots, cap_roots):
            body = join.body.surface.jet(join.body.at_aft(body_aft))
            cap = join.cap.surface.jet(join.cap.at_aft(cap_aft))
            body_frame, cap_frame = join_frames(body, cap)
            _, body_grad = field(body, body_frame, eyes, axes)
            _, cap_grad = field(cap, cap_frame, eyes, axes)
            magnitudes, angles = gradient_angle(body_grad, cap_grad)
            transitions.append(dict(level_cycles=level, phase_increases_with_aft=rising,
                                    body_aft_mm=body_aft, cap_aft_mm=cap_aft,
                                    aft_offset_mm=abs(cap_aft - body_aft),
                                    point_gap_mm=norm(cap.position - body.position),
                                    level_tangent_angle_deg=float(angles[0]) if np.isfinite(angles[0]) else None,
                                    gradient_magnitudes_cycles_per_mm=magnitudes[0].tolist(),
                                    front_facing=bool(eyes[0] @ body.normal > 0 and eyes[0] @ cap.normal > 0)))
    return dict(paired=sorted(transitions, key=lambda row: row["body_aft_mm"]), unmatched=unmatched)


def local_trace(join, eyes, axes):
    traces = []
    for boundary, offsets in ((join.body, [x for x in TRACE_SPAN_OFFSETS if x <= 0]),
                              (join.cap, [x for x in TRACE_SPAN_OFFSETS if x >= 0])):
        surface = boundary.surface
        guess = boundary.at_aft(AFT_SECTION)
        for offset in sorted(offsets, key=abs):
            uv = surface.invert_plan((CUT + offset, AFT_SECTION), guess) if offset else guess
            if uv is None or surface.condition_ratio(uv) < CONDITION_MIN:
                break
            guess = uv
            jet = surface.jet(uv)
            fixed_aft = unit(np.cross(AFT, jet.normal))
            if fixed_aft @ SPAN < 0:
                fixed_aft = -fixed_aft
            frame = np.column_stack((fixed_aft, np.cross(jet.normal, fixed_aft)))
            phase, gradient = field(jet, frame, eyes, axes)
            traces.append(dict(span_mm=float(jet.position[0]), offset_from_seam_mm=offset,
                               face="body" if boundary is join.body else "cap",
                               cycles=float(phase[0]), gradient_cycles_per_mm=norm(gradient[0]),
                               fixed_aft_slope_cycles_per_mm=float(gradient[0, 0])))
    return sorted(traces, key=lambda row: row["offset_from_seam_mm"])


def verify(joins, eyes, axes):
    phase_error = gradient_error = flip_error = 0.0
    step_mm = .001
    checks = 0
    for join in joins.values():
        for fraction in (.25, .5, .75, .9):
            aft = join.nose + fraction * (join.tail - join.nose)
            for boundary, span in ((join.body, CUT - 1), (join.cap, CUT + 1)):
                surface = boundary.surface
                uv = surface.invert_plan((span, aft), boundary.at_aft(aft))
                if uv is None:
                    raise ValueError(f"Failed to sample an interior derivative check at {span}, {aft}")
                jet = surface.jet(uv)
                local_frame = np.linalg.qr(jet.first)[0]
                cycles, gradient = field(jet, local_frame, eyes, axes)
                flipped = Jet(jet.position, jet.first, jet.second, -jet.normal)
                flip_cycles, flip_gradient = field(flipped, local_frame, eyes, axes)
                flip_error = max(flip_error, np.max(np.abs(cycles - flip_cycles)),
                                 np.max(np.abs(gradient - flip_gradient)))

                def reflected_cycles(at_uv):
                    _, first = surface.d1(at_uv)
                    normal = surface.sign * unit(np.cross(first[:, 0], first[:, 1]))
                    incident = -eyes
                    reflected = incident - 2 * ((incident @ normal)[:, None] * normal)
                    return 3 * np.einsum("ij,ij->i", reflected, axes)

                speeds = np.linalg.lstsq(jet.first, local_frame, rcond=None)[0]
                finite = np.stack([(reflected_cycles(uv + step_mm * speed)
                                    - reflected_cycles(uv - step_mm * speed)) / (2 * step_mm)
                                   for speed in speeds.T], axis=1)
                phase_error = max(phase_error, np.max(np.abs(cycles - reflected_cycles(uv))))
                gradient_error = max(gradient_error, np.max(np.abs(gradient - finite)))
                checks += 1
    if phase_error > 1e-11 or gradient_error > 1e-5 or flip_error > 1e-11:
        raise ValueError(f"Zebra verification failed: phase {phase_error}, gradient {gradient_error}, flip {flip_error}")
    return dict(interior_surface_points=checks, views_per_point=len(eyes),
                finite_difference_half_step_mm=step_mm,
                max_direct_reflection_disagreement_cycles=float(phase_error),
                max_finite_difference_component_error_cycles_per_mm=float(gradient_error),
                max_normal_flip_phase_or_gradient_component_error=float(flip_error))


def analyze(previous, revised):
    models = {}
    for label, path in (("previous", previous), ("revised", revised)):
        joins, digest = read_joins(path)
        models[label] = dict(joins=joins, digest=digest, path=path)
    all_joins = [join for model in models.values() for join in model["joins"].values()]
    nose = max(join.nose for join in all_joins)
    tail = min(join.tail for join in all_joins)
    trim = CHORD_TRIM * (tail - nose)
    stations = np.r_[np.linspace(nose + trim, tail - trim, SAMPLE_COUNT), AFT_SECTION]
    if not stations[0] <= AFT_SECTION <= stations[-2]:
        raise ValueError(f"{AFT_SECTION} mm is outside the common regular sampling range")
    views = make_views()
    eyes, axes = (np.asarray(vectors) for vectors in zip(*(view.vectors() for view in views)))
    report = {
        "method": {
            "source": "Inline full-wing-g2.html mirrored zebra shader on exact STEP surface derivatives",
            "field": "shader_phase = 6 * dot(reflect(-diagnostic_eye, unit_normal), stripe_axis); "
                     "cycles = shader_phase / 2; band = sin(2*pi*cycles). Phase gap is the absolute "
                     "difference in unwrapped cycles.",
            "coordinate_axes": ["positive span", "height up", "aft"],
            "gradient": "Analytic STEP D1/D2 derivatives, per mm of surface distance. Gradients from each "
                        "face are rotated by the shortest normal alignment into an orthonormal common "
                        "tangent frame with axes toward positive span and increasing aft.",
            "direction": "Acute angle between unoriented iso-phase tangents after the same rotation. "
                         "Undefined below the minimum gradient magnitude on either side.",
            "scope": "Fixed-aft matches on the 540 mm span body/cap seam only, with an additional "
                     "fixed-aft span trace at 150 mm. No raster sampling, occlusion, mesh interpolation or "
                     "antialiasing in the field calculation. Front-facing stats require both "
                     "oriented face normals to face the mirrored diagnostic eye.",
            "trace": "Fixed-aft slope is the signed derivative in the positive-span direction along "
                     "the actual 150 mm aft surface section. The trace includes both seam limits.",
            "stripe_edges": "Default dark/light transitions occur where cycles is a multiple of 0.5. "
                            "Bracketed crossings are solved separately on the body and cap pcurves. "
                            "Angle and point gap compare corresponding roots; tangential touches or "
                            "multiple roots within one sampling interval can escape bracketing.",
        },
        "sampling": {
            "seam_span_mm": CUT, "common_edge_aft_bounds_mm": [nose, tail],
            "included_aft_bounds_mm": [float(stations[0]), float(stations[-2])],
            "excluded_end_fraction_each": CHORD_TRIM,
            "seam_stations_per_skin": SAMPLE_COUNT, "extra_section_aft_mm": AFT_SECTION,
            "minimum_d1_singular_value_ratio": CONDITION_MIN,
            "endpoints": "The seam vertices are excluded because UV parameterizations can collapse "
                         "or become ill-conditioned there.",
            "direction_minimum_gradient_cycles_per_mm": DIRECTION_MIN_CYCLES_PER_MM,
            "default_view": views[0].describe(), "sweep_view_count": len(views),
            "sweep_views": [view.describe() for view in views],
            "aggregation": "Uniform in aft on each skin and equally weighted across sweep settings. "
                           "The exact 150 mm station is summarized separately. Maxima are sampled maxima.",
        },
        "models": {},
    }
    for label, model in models.items():
        joins = model["joins"]
        rows = read_field(joins, stations, eyes, axes)
        report["models"][label] = dict(
            step=model["path"].relative_to(HERE.parent).as_posix()
            if model["path"].is_relative_to(HERE.parent) else str(model["path"]),
            sha256=model["digest"],
            default=summarize(rows, stations, views, slice(0, 1)),
            sweep=summarize(rows, stations, views, slice(None)),
            default_profiles=default_profiles(rows, stations),
            default_stripe_edges={
                side: stripe_crossings(join, stations[:-1], rows[side]["body_cycles"][:-1, 0],
                                       rows[side]["cap_cycles"][:-1, 0], eyes[:1], axes[:1])
                for side, join in joins.items()
            },
            default_at_150_mm_span_traces={
                side: local_trace(join, eyes[:1], axes[:1]) for side, join in joins.items()
            },
            verification=verify(joins, eyes, axes),
        )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", type=Path, default=PREVIOUS)
    parser.add_argument("--revised", type=Path, default=REVISED)
    parser.add_argument("--output", type=Path, default=HERE / "zebra-analysis.json")
    args = parser.parse_args()
    report = analyze(args.previous.resolve(), args.revised.resolve())
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for version, model in report["models"].items():
        for case in ("default", "sweep"):
            seam = model[case]["combined"]["seam"]
            print(version, case, {key: (row["p95"], row["max"]) for key, row in seam.items()})
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
