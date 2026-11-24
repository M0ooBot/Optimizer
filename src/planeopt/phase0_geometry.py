"""Offline Phase 0 layouts in millimetres, based on the saved Converter rib.

The local study axes are -X outboard, +Y up and +Z forward. Parallel sections
use constant X. Front-normal sections are normal to a constant-Z bend guide
with the same rise and projected length as the leading edge. Tilted tips can
shift the skin's projected extent relative to that guide.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Iterator, Mapping

Point = tuple[float, float, float]
Curve = tuple[Point, ...]
Curves = dict[str, Curve]

CONVERTER_REVISION = "d77bbdf45a8c0aa92efbc2a3d3007cfe4f904105"
SKIN_FRACTIONS = (0.01, 0.02, 0.035, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25,
                  0.30, 0.40, 0.50, 0.60, 0.75, 0.90, 0.95)


@dataclass(frozen=True)
class HingeDefinition:
    start: float
    end: float
    inboard_chord_fraction: float
    outboard_chord_fraction: float

    def __post_init__(self) -> None:
        if not 0 <= self.start < self.end <= 1:
            raise ValueError("Hinge start and end must be ordered span fractions")
        if not all(0 < f < 1 for f in
                   (self.inboard_chord_fraction, self.outboard_chord_fraction)):
            raise ValueError("Hinge chord fractions must lie inside the section")


@dataclass(frozen=True)
class SurfaceDefinition:
    root_le_mm: Point
    root_chord_mm: float
    projected_span_mm: float
    tip_chord_mm: float
    sweep_mm: float
    rake_mm: float
    rake_start: float | None
    rise_mm: float
    rise_start: float | None
    tip_twist_deg: float
    station_fractions: tuple[float, ...]
    hinge: HingeDefinition
    section_orientation: str = "parallel"
    tip_span_trim_mm: float = 0.
    span_curve: str = "onset"
    midspan_twist_peak_deg: float | None = None
    full_span_linear_rise_fraction: float = 0.

    def __post_init__(self) -> None:
        if len(self.root_le_mm) != 3:
            raise ValueError("Root leading edge must have three coordinates")
        values = (*self.root_le_mm, self.root_chord_mm, self.projected_span_mm,
                  self.tip_chord_mm, self.sweep_mm, self.rake_mm, self.rise_mm,
                  self.tip_twist_deg, self.tip_span_trim_mm)
        if not all(math.isfinite(value) for value in values):
            raise ValueError("Surface dimensions must be finite")
        peak = self.midspan_twist_peak_deg
        if peak is not None and (not isinstance(peak, (int, float)) or isinstance(peak, bool)
                                 or not math.isfinite(peak)):
            raise ValueError("midspan_twist_peak_deg must be a finite angle or null")
        linear_rise = self.full_span_linear_rise_fraction
        if (not isinstance(linear_rise, (int, float)) or isinstance(linear_rise, bool)
                or not 0 <= linear_rise <= 1):
            raise ValueError("full_span_linear_rise_fraction must be in [0, 1]")
        if linear_rise and self.span_curve != "full_span":
            raise ValueError("A nonzero full_span_linear_rise_fraction requires full_span")
        if linear_rise and self.rise_mm and self.section_orientation == "front_normal":
            raise ValueError("full_span_linear_rise_fraction would rotate the locked front_normal root")
        if min(self.root_chord_mm, self.projected_span_mm, self.tip_chord_mm) <= 0:
            raise ValueError("Chords and projected span must be positive")
        if self.tip_span_trim_mm < 0 or self.tip_span_trim_mm >= self.projected_span_mm:
            raise ValueError("Tip span trim must be nonnegative and keep positive guide length")
        if self.span_curve == "onset":
            if not all(start is not None and 0 <= start < 1
                       for start in (self.rake_start, self.rise_start)):
                raise ValueError("Onset rake and rise start must be in [0, 1)")
        elif self.span_curve == "full_span":
            if self.rake_start is not None or self.rise_start is not None:
                raise ValueError("Full-span curves require null rake_start and rise_start")
            if self.sweep_mm < 0 or self.rake_mm < 0:
                raise ValueError("Full-span curves require nonnegative aft sweep and rake")
        else:
            raise ValueError("Span curve must be onset or full_span")
        stations = self.station_fractions
        if (len(stations) < 2 or stations[0] != 0 or stations[-1] != 1
                or not all(a < b for a, b in zip(stations, stations[1:]))):
            raise ValueError("Stations must increase strictly from 0 to 1")
        if self.section_orientation not in {"parallel", "front_normal"}:
            raise ValueError("Section orientation must be parallel or front_normal")


@dataclass(frozen=True)
class ShapeDefinition:
    name: str
    wing: SurfaceDefinition
    tail: SurfaceDefinition
    tolerance_mm: float = 0.25

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("Shape name must be nonempty")
        if not math.isfinite(self.tolerance_mm) or self.tolerance_mm <= 0:
            raise ValueError("CAD tolerance must be finite and positive")


@dataclass(frozen=True)
class SectionFrame:
    le: Point
    te: Point
    quarter_chord: Point
    span_tangent: Point
    leading_edge_tangent: Point
    chord_axis: Point
    up_axis: Point
    chord_mm: float


@dataclass(frozen=True)
class _ControlSample:
    fraction: float
    frame: SectionFrame
    chord_fraction: float
    margin_mm: float
    upper: float
    lower: float


def _surface_from_dict(values: Mapping) -> SurfaceDefinition:
    fields = dict(values)
    fields["root_le_mm"] = tuple(fields["root_le_mm"])
    fields["station_fractions"] = tuple(fields["station_fractions"])
    fields["hinge"] = HingeDefinition(**fields["hinge"])
    if fields.get("span_curve") == "full_span":
        fields.setdefault("rise_start", None)
        fields.setdefault("rake_start", None)
    return SurfaceDefinition(**fields)


def load_definition(path: str | Path) -> ShapeDefinition:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    data["wing"] = _surface_from_dict(data["wing"])
    data["tail"] = _surface_from_dict(data["tail"])
    return ShapeDefinition(**data)


def save_definition(definition: ShapeDefinition, path: str | Path) -> None:
    data = asdict(definition)
    for surface in (data["wing"], data["tail"]):
        if surface["midspan_twist_peak_deg"] is None:
            del surface["midspan_twist_peak_deg"]
        if surface["full_span_linear_rise_fraction"] == 0:
            del surface["full_span_linear_rise_fraction"]
    Path(path).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def validate_edit(
    previous: ShapeDefinition,
    revised: ShapeDefinition,
    *,
    lock_tail: bool = False,
    lock_tip_height: bool = False,
    lock_tip_aft: bool = False,
) -> None:
    conflicts = []
    for surface_name in ("wing", "tail"):
        old = getattr(previous, surface_name)
        new = getattr(revised, surface_name)
        for field in ("root_le_mm", "root_chord_mm"):
            if getattr(old, field) != getattr(new, field):
                conflicts.append(f"{surface_name}.{field}")
    if previous.wing.projected_span_mm != revised.wing.projected_span_mm:
        conflicts.append("wing.projected_span_mm")
    if lock_tail and previous.tail != revised.tail:
        conflicts.append("tail")
    if lock_tip_height and previous.wing.rise_mm != revised.wing.rise_mm:
        conflicts.append("wing.rise_mm")
    if lock_tip_aft and not math.isclose(
        previous.wing.sweep_mm + previous.wing.rake_mm,
        revised.wing.sweep_mm + revised.wing.rake_mm,
        rel_tol=0., abs_tol=1e-9,
    ):
        conflicts.append("wing.sweep_mm + wing.rake_mm")
    if previous.tolerance_mm != revised.tolerance_mm:
        conflicts.append("tolerance_mm")
    if conflicts:
        raise ValueError("Conflicts with locked values: " + ", ".join(conflicts))


def apply_edit(
    definition: ShapeDefinition,
    *,
    name: str,
    wing: Mapping | None = None,
    tail: Mapping | None = None,
    lock_tail: bool = False,
    lock_tip_height: bool = False,
    lock_tip_aft: bool = False,
) -> ShapeDefinition:
    """Merge a partial edit and enforce the study locks before generating geometry."""
    surfaces = {}
    for surface_name, changes in (("wing", wing), ("tail", tail)):
        if changes:
            if "tip_span_trim_mm" in changes:
                raise ValueError(f"{surface_name}.tip_span_trim_mm is derived from the physical span")
            fields = asdict(getattr(definition, surface_name))
            if changes.get("span_curve") == "full_span":
                fields.update(rise_start=None, rake_start=None)
            hinge_changes = changes.get("hinge")
            fields.update(changes)
            if hinge_changes is not None:
                fields["hinge"] = {**asdict(getattr(definition, surface_name).hinge),
                                   **hinge_changes}
            surfaces[surface_name] = _surface_from_dict(fields)
    revised = replace(definition, name=name, **surfaces)
    validate_edit(definition, revised, lock_tail=lock_tail, lock_tip_height=lock_tip_height,
                  lock_tip_aft=lock_tip_aft)
    return revised


def _unit(vector: Point) -> Point:
    length = math.hypot(*vector)
    return tuple(value / length for value in vector)


def _shift(origin: Point, axis: Point, distance: float) -> Point:
    return tuple(value + direction * distance for value, direction in zip(origin, axis))


def _twist_deg(surface: SurfaceDefinition, fraction: float) -> float:
    peak = surface.midspan_twist_peak_deg
    if peak is None:
        return fraction * surface.tip_twist_deg
    if fraction <= .5:
        start, end, t = 0., peak, 2 * fraction
    else:
        start, end, t = peak, surface.tip_twist_deg, 2 * fraction - 1
    # Quintic smoothstep has zero first and second derivatives at both ends.
    weight = t**3 * (10 - 15 * t + 6 * t**2)
    return (1 - weight) * start + weight * end


def section_frame(surface: SurfaceDefinition, fraction: float) -> SectionFrame:
    if not 0 <= fraction <= 1:
        raise ValueError("Section fraction must be in [0, 1]")
    if surface.span_curve == "full_span":
        linear_rise = surface.full_span_linear_rise_fraction
        rise = surface.rise_mm * ((1 - linear_rise) * fraction**2 + linear_rise * fraction)
        rise_rate = surface.rise_mm * (2 * (1 - linear_rise) * fraction + linear_rise)
        rake = surface.rake_mm * fraction**2
        rake_rate = 2 * surface.rake_mm * fraction
    else:
        rise_start, rake_start = surface.rise_start, surface.rake_start
        assert rise_start is not None and rake_start is not None
        rise_progress = max(0., (fraction - rise_start) / (1 - rise_start))
        rake_progress = max(0., (fraction - rake_start) / (1 - rake_start))
        # Cubic rise has zero slope and curvature at its onset. Quintic
        # rake transitions smoothly into and out of the added aft offset.
        rise = surface.rise_mm * rise_progress**3
        rise_rate = 3 * surface.rise_mm * rise_progress**2 / (1 - rise_start)
        rake = surface.rake_mm * rake_progress**3 * (10 - 15 * rake_progress + 6 * rake_progress**2)
        rake_rate = 30 * surface.rake_mm * rake_progress**2 * (1 - rake_progress)**2 / (1 - rake_start)
    root_x, root_y, root_z = surface.root_le_mm
    guide_span = surface.projected_span_mm - surface.tip_span_trim_mm
    le = (root_x - fraction * guide_span,
          root_y + rise,
          root_z - (fraction * surface.sweep_mm + rake))

    leading_edge_tangent = _unit((-guide_span, rise_rate, -surface.sweep_mm - rake_rate))
    span_tangent = (-1., 0., 0.)
    up = (0., 1., 0.)
    if surface.section_orientation == "front_normal":
        span_tangent = _unit((-guide_span, rise_rate, 0.))
        up = _unit((rise_rate, guide_span, 0.))
    angle = math.radians(_twist_deg(surface, fraction))
    sin_angle, cos_angle = math.sin(angle), math.cos(angle)
    chord_axis = (-up[0] * sin_angle, -up[1] * sin_angle, -cos_angle)
    up_axis = (up[0] * cos_angle, up[1] * cos_angle, -sin_angle)
    chord = surface.root_chord_mm + fraction * (surface.tip_chord_mm - surface.root_chord_mm)
    return SectionFrame(le, _shift(le, chord_axis, chord), _shift(le, chord_axis, .25 * chord),
                        span_tangent, leading_edge_tangent, chord_axis, up_axis, chord)


@lru_cache(maxsize=4)
def _profile(baseline_dir: Path) -> tuple[tuple[float, float], ...]:
    # Use the Converter's resampled wing section, including its point ordering,
    # continuous nose treatment and sharp trailing-edge closure.
    rows = (baseline_dir / "curves" / "baseline_s01.sldcrv").read_text(encoding="ascii")
    points = []
    for row in rows.splitlines():
        x, y, z = (float(value.removesuffix("mm")) for value in row.split())
        if x != 0:
            raise ValueError("Expected the baseline section at the root station")
        points.append((-z / 200., y / 200.))
    if len(points) != 321 or points[0] != (1., 0.) or points[-1] != (1., 0.):
        raise ValueError("Expected the saved 200 mm closed Converter SD7037 root curve")
    return tuple(points)


def _guide_anchor(baseline_dir: Path, name: str) -> tuple[float, float]:
    rows = (baseline_dir / "curves" / f"baseline_{name}.sldcrv").read_text(encoding="ascii")
    _, y, z = (float(value.removesuffix("mm")) for value in rows.splitlines()[0].split())
    return -z / 200., y / 200.


def _point(frame: SectionFrame, chord_fraction: float, height_fraction: float) -> Point:
    chord_point = _shift(frame.le, frame.chord_axis, frame.chord_mm * chord_fraction)
    return _shift(chord_point, frame.up_axis, frame.chord_mm * height_fraction)


def section_points(surface: SurfaceDefinition, fraction: float, baseline_dir: str | Path) -> Curve:
    frame = section_frame(surface, fraction)
    return tuple(_point(frame, x, y) for x, y in _profile(Path(baseline_dir).resolve()))


def _guide_stations(surface: SurfaceDefinition) -> list[float]:
    samples = {index / 120 for index in range(121)}
    samples.update(surface.station_fractions)
    if surface.rise_start is not None and surface.rake_start is not None:
        samples.update((surface.rise_start, surface.rake_start))
    return sorted(samples)


def projected_skin_span(surface: SurfaceDefinition, baseline_dir: str | Path) -> float:
    """Sampled maximum outboard distance in X from the root mounting plane."""
    profile = _profile(Path(baseline_dir).resolve())
    min_x = surface.root_le_mm[0]
    for fraction in _guide_stations(surface):
        frame = section_frame(surface, fraction)
        chord_x = frame.chord_mm * frame.chord_axis[0]
        up_x = frame.chord_mm * frame.up_axis[0]
        min_x = min(min_x, frame.le[0] + min(chord_x * x + up_x * y for x, y in profile))
    return surface.root_le_mm[0] - min_x


def calibrate_projected_span(surface: SurfaceDefinition, baseline_dir: str | Path) -> SurfaceDefinition:
    """Trim the LE reference inward to meet the sampled outer-skin span target."""
    target = surface.projected_span_mm
    lower, upper = 0., target / 2
    untrimmed = replace(surface, tip_span_trim_mm=lower)
    if abs(projected_skin_span(untrimmed, baseline_dir) - target) < 1e-6:
        return untrimmed
    if projected_skin_span(replace(surface, tip_span_trim_mm=upper), baseline_dir) > target:
        raise ValueError("Cannot match projected skin span within half the original guide length")
    for _ in range(32):
        trim = (lower + upper) / 2
        calibrated = replace(surface, tip_span_trim_mm=trim)
        deviation = projected_skin_span(calibrated, baseline_dir) - target
        if abs(deviation) < 1e-6:
            return calibrated
        if deviation > 0:
            lower = trim
        else:
            upper = trim
    raise ValueError("Could not converge projected skin span calibration")


def _skin_heights(profile: tuple[tuple[float, float], ...], x: float) -> tuple[float, float]:
    intersections = []
    for (x0, y0), (x1, y1) in zip(profile, profile[1:]):
        if x0 != x1 and min(x0, x1) <= x <= max(x0, x1):
            intersections.append(y0 + (y1 - y0) * (x - x0) / (x1 - x0))
    if len(intersections) < 2:
        raise ValueError(f"Profile does not give upper and lower skin at {x:g} chord")
    return max(intersections), min(intersections)


def hinge_points(surface: SurfaceDefinition, baseline_dir: str | Path | None = None) -> tuple[Point, Point]:
    """Straight hinge through mid-thickness at its two span cuts."""
    if baseline_dir is None:
        baseline_dir = Path(__file__).resolve().parents[2] / "aircraft" / "phase0_sd7037"
    profile = _profile(Path(baseline_dir).resolve())
    spec = surface.hinge
    start_camber = sum(_skin_heights(profile, spec.inboard_chord_fraction)) / 2
    end_camber = sum(_skin_heights(profile, spec.outboard_chord_fraction)) / 2
    return (
        _point(section_frame(surface, spec.start), spec.inboard_chord_fraction, start_camber),
        _point(section_frame(surface, spec.end), spec.outboard_chord_fraction, end_camber),
    )


def _control_samples(surface: SurfaceDefinition, baseline_dir: Path) -> Iterator[_ControlSample]:
    profile = _profile(baseline_dir)
    start, end = hinge_points(surface, baseline_dir)
    delta = tuple(b - a for a, b in zip(start, end))
    spec = surface.hinge
    previous_along = -1.
    for index in range(121):
        fraction = spec.start + (spec.end - spec.start) * index / 120
        frame = section_frame(surface, fraction)
        offset = tuple(a - b for a, b in zip(frame.le, start))
        projected_length = math.fsum(a * b for a, b in zip(delta, frame.span_tangent))
        if abs(projected_length) < 1e-9:
            raise ValueError(f"Hinge is parallel to the section plane at span fraction {fraction:.4f}")
        along = math.fsum(a * b for a, b in zip(offset, frame.span_tangent)) / projected_length
        if not (-1e-9 <= along <= 1 + 1e-9 and along > previous_along):
            raise ValueError(f"Hinge stations fall outside the ordered finite axis at {fraction:.4f}")
        previous_along = along
        intersection = _shift(start, delta, along)
        relative = tuple(a - b for a, b in zip(intersection, frame.le))
        chord_fraction = math.fsum(a * b for a, b in zip(relative, frame.chord_axis)) / frame.chord_mm
        height = math.fsum(a * b for a, b in zip(relative, frame.up_axis)) / frame.chord_mm
        if not 0 < chord_fraction < 1:
            raise ValueError(f"Hinge crosses outside the chord at span fraction {fraction:.4f}")
        upper, lower = _skin_heights(profile, chord_fraction)
        margin_mm = frame.chord_mm * min(upper - height, height - lower)
        yield _ControlSample(fraction, frame, chord_fraction, margin_mm, upper, lower)


def hinge_fit(surface: SurfaceDefinition, baseline_dir: str | Path) -> tuple[float, float]:
    """Minimum sampled section-normal height margin in mm and its span fraction."""
    samples = _control_samples(surface, Path(baseline_dir).resolve())
    worst = min(samples, key=lambda sample: sample.margin_mm)
    return worst.margin_mm, worst.fraction


def _control_outline(
    samples: list[_ControlSample], profile: tuple[tuple[float, float], ...], upper: bool
) -> Curve:
    front = [_point(sample.frame, sample.chord_fraction, sample.upper if upper else sample.lower)
             for sample in samples]
    nose_index = min(range(len(profile)), key=lambda index: profile[index][0])
    side = profile[:nose_index + 1] if upper else profile[nose_index:]

    def cut(sample: _ControlSample) -> list[Point]:
        skin = reversed(side) if upper else side
        return [_point(sample.frame, x, y) for x, y in skin if x > sample.chord_fraction]

    rear = [sample.frame.te for sample in reversed(samples[:-1])]
    return tuple(front + cut(samples[-1]) + rear + list(reversed(cut(samples[0])[:-1])) + front[:1])


def surface_section_names(surface: SurfaceDefinition, prefix: str = "wing") -> tuple[str, ...]:
    intermediates = tuple(f"{prefix}_s{index:02d}" for index in range(1, len(surface.station_fractions) - 1))
    return (f"{prefix}_root", *intermediates, f"{prefix}_tip")


def build_curves(definition: ShapeDefinition, baseline_dir: str | Path) -> Curves:
    baseline_dir = Path(baseline_dir).resolve()
    nose = _guide_anchor(baseline_dir, "le")
    profile = _profile(baseline_dir)
    curves = {}
    for prefix, surface in (("wing", definition.wing), ("tail", definition.tail)):
        for name, fraction in zip(surface_section_names(surface, prefix), surface.station_fractions):
            curves[name] = section_points(surface, fraction, baseline_dir)
        frames = [section_frame(surface, station) for station in _guide_stations(surface)]
        curves[f"{prefix}_le"] = tuple(_point(frame, *nose) for frame in frames)
        curves[f"{prefix}_te"] = tuple(frame.te for frame in frames)
        curves[f"{prefix}_bend_guide"] = tuple(
            (frame.le[0], frame.le[1], surface.root_le_mm[2]) for frame in frames
        )
        curves[f"{prefix}_hinge"] = hinge_points(surface, baseline_dir)
        control = list(_control_samples(surface, baseline_dir))
        curves[f"{prefix}_control_upper"] = _control_outline(control, profile, upper=True)
        curves[f"{prefix}_control_lower"] = _control_outline(control, profile, upper=False)
        for fraction in SKIN_FRACTIONS:
            percentage = 100 * fraction
            tag = (f"{int(percentage):02d}" if percentage.is_integer()
                   else f"{percentage:g}".replace(".", "p"))
            for side in ("upper", "lower"):
                name = f"{side}_{tag}"
                anchor = _guide_anchor(baseline_dir, name)
                curves[f"{prefix}_{name}"] = tuple(_point(frame, *anchor) for frame in frames)
    return curves


def write_curves(curves: Curves, output_dir: str | Path) -> None:
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    for name, curve in curves.items():
        values = ([f"{coordinate:.10f}" for coordinate in point] for point in curve)
        rows = ("\t".join(f"{s.removeprefix('-') if float(s) == 0 else s}mm" for s in row)
                for row in values)
        (directory / f"{name}.sldcrv").write_bytes(("\r\n".join(rows) + "\r\n").encode("ascii"))
