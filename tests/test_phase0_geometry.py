from pathlib import Path
from dataclasses import replace
import math

import pytest

from planeopt.phase0_geometry import (
    apply_edit,
    build_curves,
    calibrate_projected_span,
    hinge_fit,
    load_definition,
    projected_skin_span,
    save_definition,
    section_frame,
    section_points,
    surface_section_names,
    write_curves,
)


AIRCRAFT = Path(__file__).resolve().parents[1] / "aircraft"
BASELINE_DIR = AIRCRAFT / "phase0_sd7037"


@pytest.fixture
def baseline():
    return load_definition(AIRCRAFT / "phase0_shape_study" / "baseline.json")


def test_baseline_regenerates_pinned_converter_profiles_and_intersects_guides(baseline, tmp_path):
    curves = build_curves(baseline, BASELINE_DIR)
    write_curves(curves, tmp_path)
    for number, name in enumerate(surface_section_names(baseline.wing), start=1):
        expected = (BASELINE_DIR / "curves" / f"baseline_s{number:02d}.sldcrv").read_bytes()
        assert (tmp_path / f"{name}.sldcrv").read_bytes() == expected
        station_x = curves[name][0][0]
        for guide_name, points in curves.items():
            if guide_name in ("wing_le", "wing_te") or guide_name.startswith(("wing_upper", "wing_lower")):
                at_station = next(point for point in points if point[0] == station_x)
                assert at_station in curves[name]


def test_local_edit_preserves_locked_root_span_tip_and_tail(baseline, tmp_path):
    revised = apply_edit(
        baseline, name="swept",
        wing={
            "tip_chord_mm": 85, "sweep_mm": 90, "rake_mm": 30,
            "rise_mm": 45, "rise_start": .45, "tip_twist_deg": -2.,
            "section_orientation": "front_normal",
            "hinge": {"end": .73},
        },
        tail={"tip_chord_mm": 65, "sweep_mm": 22, "hinge": {"end": .80}},
    )
    follow_up = apply_edit(
        revised, name="late rise", wing={"rise_start": .65},
        lock_tip_height=True, lock_tail=True,
    )
    path = tmp_path / "definition.json"
    save_definition(follow_up, path)
    assert load_definition(path) == follow_up
    assert build_curves(load_definition(path), BASELINE_DIR) == build_curves(follow_up, BASELINE_DIR)
    assert section_points(revised.wing, 0, BASELINE_DIR) == section_points(baseline.wing, 0, BASELINE_DIR)
    assert section_frame(follow_up.wing, 1).le == section_frame(revised.wing, 1).le
    assert section_frame(follow_up.wing, .65).le[1] < section_frame(revised.wing, .65).le[1]
    assert follow_up.tail == revised.tail
    for fraction in (0., .45, .65, .8, 1.):
        frame = section_frame(follow_up.wing, fraction)
        assert math.hypot(*frame.span_tangent) == pytest.approx(1.)
        assert math.fsum(a * b for a, b in zip(frame.chord_axis, frame.span_tangent)) == pytest.approx(0.)
        assert math.fsum(a * b for a, b in zip(frame.up_axis, frame.span_tangent)) == pytest.approx(0.)
    assert section_frame(follow_up.wing, 1).span_tangent[1] > 0

    for definition in (revised, follow_up):
        surface = definition.wing
        assert projected_skin_span(surface, BASELINE_DIR) > 600
        calibrated = calibrate_projected_span(surface, BASELINE_DIR)
        assert projected_skin_span(calibrated, BASELINE_DIR) == pytest.approx(600, abs=1e-6)
        assert calibrated.tip_span_trim_mm > 0
        assert calibrated.projected_span_mm == 600
        assert section_frame(calibrated, 1).le[1] == section_frame(surface, 1).le[1]
        save_definition(replace(definition, wing=calibrated), path)
        assert load_definition(path).wing == calibrated


def test_incompatible_edits_name_each_conflict(baseline):
    with pytest.raises(ValueError, match=r"wing.root_chord_mm.*tail.root_le_mm.*wing.projected_span_mm"):
        apply_edit(baseline, name="broken",
                   wing={"root_chord_mm": 205, "projected_span_mm": 605},
                   tail={"root_le_mm": [0, 0, -700]})
    with pytest.raises(ValueError, match=r"tail, wing.rise_mm"):
        apply_edit(baseline, name="broken", wing={"rise_mm": 50}, tail={"tip_chord_mm": 50},
                   lock_tail=True, lock_tip_height=True)
    with pytest.raises(ValueError, match=r"wing.tip_span_trim_mm is derived"):
        apply_edit(baseline, name="manual trim", wing={"tip_span_trim_mm": 1})


def test_reversing_twist_rotates_about_leading_edge_with_smooth_midpoint(baseline, tmp_path):
    linear = apply_edit(baseline, name="linear", wing={"tip_twist_deg": -4})
    reversing = apply_edit(
        linear, name="reversing", wing={"midspan_twist_peak_deg": 4}, lock_tail=True,
    )

    def incidence(fraction):
        frame = section_frame(reversing.wing, fraction)
        return math.degrees(math.atan2(frame.le[1] - frame.te[1], frame.le[2] - frame.te[2]))

    stations = (0, .25, .5, .75, 1)
    assert [incidence(fraction) for fraction in stations] == pytest.approx([0, 2, 4, 0, -4])
    for fraction in stations:
        revised = section_frame(reversing.wing, fraction)
        prior = section_frame(linear.wing, fraction)
        assert revised.le == prior.le
        assert math.dist(revised.le, revised.te) == pytest.approx(prior.chord_mm)
    assert reversing.tail == baseline.tail

    # The residual from the peak scales cubically on either side, giving
    # zero one-sided slope and curvature at the join.
    for direction in (-1, 1):
        farther = 4 - incidence(.5 + direction * .002)
        nearer = 4 - incidence(.5 + direction * .001)
        assert farther / nearer == pytest.approx(8, rel=.004)

    path = tmp_path / "twist.json"
    save_definition(reversing, path)
    assert load_definition(path) == reversing
    reset = apply_edit(reversing, name=linear.name, wing={"midspan_twist_peak_deg": None})
    assert reset == linear
    save_definition(reset, path)
    assert "midspan_twist_peak_deg" not in path.read_text(encoding="utf-8")
    assert section_frame(load_definition(path).wing, .5).chord_axis == (
        0, -math.sin(math.radians(-2)), -math.cos(math.radians(-2)),
    )


@pytest.mark.parametrize("invalid", [float("nan"), float("inf"), "4", True])
def test_midspan_twist_peak_requires_finite_angle(baseline, invalid):
    with pytest.raises(ValueError, match="midspan_twist_peak_deg"):
        apply_edit(baseline, name="invalid", wing={"midspan_twist_peak_deg": invalid})


def test_full_span_law_distributes_rise_and_rake_with_tip_locks(baseline, tmp_path):
    onset = apply_edit(
        baseline, name="onset-based",
        wing={"tip_chord_mm": 22, "sweep_mm": 242, "rake_mm": 28, "rise_mm": 65,
              "section_orientation": "front_normal", "tip_twist_deg": -2,
              "hinge": {"start": .39, "end": .64}},
    )
    full = apply_edit(
        onset, name="full span",
        wing={"span_curve": "full_span", "sweep_mm": 205, "rake_mm": 65},
        lock_tip_height=True, lock_tip_aft=True, lock_tail=True,
    )
    assert full.wing.rise_start is None
    assert full.wing.rake_start is None
    calibrated = calibrate_projected_span(full.wing, BASELINE_DIR)
    assert projected_skin_span(calibrated, BASELINE_DIR) == pytest.approx(600, abs=1e-6)
    assert section_points(calibrated, 0, BASELINE_DIR) == section_points(onset.wing, 0, BASELINE_DIR)
    assert section_frame(calibrated, 1).le[1:] == section_frame(onset.wing, 1).le[1:]
    assert full.tail == onset.tail
    for fraction in (.1, .25, .5, .75, 1.):
        frame = section_frame(calibrated, fraction)
        assert frame.le[1] == pytest.approx(65 * fraction**2)
        assert frame.le[2] == pytest.approx(-205 * fraction - 65 * fraction**2)
        assert frame.span_tangent[1] / -frame.span_tangent[0] == pytest.approx(
            130 * fraction / (600 - calibrated.tip_span_trim_mm)
        )
        assert frame.leading_edge_tangent[2] / frame.leading_edge_tangent[0] == pytest.approx(
            (205 + 130 * fraction) / (600 - calibrated.tip_span_trim_mm)
        )
    path = tmp_path / "full-span.json"
    save_definition(replace(full, wing=calibrated), path)
    assert load_definition(path).wing == calibrated
    with pytest.raises(ValueError, match="wing.sweep_mm"):
        apply_edit(full, name="aft change", wing={"sweep_mm": 210}, lock_tip_aft=True)
    with pytest.raises(ValueError, match="require null"):
        apply_edit(full, name="onset conflict", wing={"rise_start": .2})


@pytest.mark.parametrize("linear_fraction", [.35, 1.])
def test_full_span_linear_rise_blends_height_and_tangent(baseline, tmp_path, linear_fraction):
    full = apply_edit(
        baseline, name="quadratic",
        wing={"span_curve": "full_span", "rise_mm": 65, "section_orientation": "parallel"},
    )
    revised = apply_edit(
        full, name="blended", wing={"full_span_linear_rise_fraction": linear_fraction},
        lock_tip_height=True, lock_tip_aft=True, lock_tail=True,
    )
    assert revised.tail == full.tail
    for endpoint in (0., 1.):
        assert section_frame(revised.wing, endpoint).le == section_frame(full.wing, endpoint).le
    assert section_points(revised.wing, 0, BASELINE_DIR) == section_points(full.wing, 0, BASELINE_DIR)
    for fraction in (0., .25, .5, 1.):
        frame = section_frame(revised.wing, fraction)
        quadratic_frame = section_frame(full.wing, fraction)
        height = 65 * ((1 - linear_fraction) * fraction**2 + linear_fraction * fraction)
        derivative = 65 * (2 * (1 - linear_fraction) * fraction + linear_fraction)
        assert frame.le[::2] == quadratic_frame.le[::2]
        assert frame.le[1] == pytest.approx(full.wing.root_le_mm[1] + height)
        assert frame.leading_edge_tangent[1] / -frame.leading_edge_tangent[0] == pytest.approx(
            derivative / revised.wing.projected_span_mm
        )
    path = tmp_path / "rise.json"
    save_definition(revised, path)
    assert load_definition(path) == revised
    save_definition(full, path)
    assert "full_span_linear_rise_fraction" not in path.read_text(encoding="utf-8")
    assert load_definition(path) == full


@pytest.mark.parametrize("value", [-.1, 1.1, float("nan"), float("inf"), True, None, ".4"])
def test_full_span_linear_rise_fraction_requires_number_in_unit_interval(baseline, value):
    with pytest.raises(ValueError, match="full_span_linear_rise_fraction"):
        apply_edit(baseline, name="invalid", wing={
            "span_curve": "full_span", "full_span_linear_rise_fraction": value,
        })


def test_nonzero_full_span_linear_rise_requires_full_span(baseline):
    with pytest.raises(ValueError, match="requires full_span"):
        apply_edit(baseline, name="invalid", wing={"full_span_linear_rise_fraction": .5})


def test_full_span_linear_rise_cannot_rotate_locked_root(baseline):
    with pytest.raises(ValueError, match="rotate the locked front_normal root"):
        apply_edit(baseline, name="invalid", wing={
            "span_curve": "full_span", "rise_mm": 65, "section_orientation": "front_normal",
            "full_span_linear_rise_fraction": .35,
        })


def test_both_control_outlines_close_around_the_trailing_edge(baseline):
    curves = build_curves(baseline, BASELINE_DIR)
    for name, surface in (("wing", baseline.wing), ("tail", baseline.tail)):
        hinge = curves[f"{name}_hinge"]
        assert len(hinge) == 2
        margin, fraction = hinge_fit(surface, BASELINE_DIR)
        assert margin > 0
        assert surface.hinge.start <= fraction <= surface.hinge.end
        for side in ("upper", "lower"):
            outline = curves[f"{name}_control_{side}"]
            assert outline[0] == outline[-1]
            assert section_frame(surface, surface.hinge.start).te in outline
            assert section_frame(surface, surface.hinge.end).te in outline
