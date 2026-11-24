from dataclasses import replace
from pathlib import Path

import pytest

from planeopt.phase0_evidence import _fallback_step, _load_solid, _measure, write_evidence
from planeopt.phase0_geometry import apply_edit, calibrate_projected_span, load_definition


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "aircraft" / "phase0_sd7037"
BASELINE = ROOT / "aircraft" / "phase0_shape_study" / "baseline.json"


def test_generated_evidence_does_not_claim_step_readback(tmp_path):
    baseline = load_definition(BASELINE)
    changed = apply_edit(baseline, name="outer rise", wing={"rise_mm": 25, "rise_start": 0.65})
    report = write_evidence(tmp_path, {"baseline": baseline, "candidate": changed}, FIXTURE)
    candidate = report["revisions"]["candidate"]["wing"]
    assert candidate["generated"]["leading_edge_tip_rise_mm"] == pytest.approx(25)
    assert candidate["generated"]["agreed_locks_within_tolerance"]
    assert candidate["step_readback"]["status"] == "unavailable"
    assert (tmp_path / "evidence.json").is_file()
    assert (tmp_path / "overlays.png").is_file()


def test_reference_step_is_measured_at_between_rib_stations(tmp_path):
    pytest.importorskip("cadquery")
    baseline = load_definition(BASELINE)
    report = write_evidence(
        tmp_path,
        {"baseline": baseline},
        FIXTURE,
        step_paths={("baseline", "wing"): FIXTURE / "SD7037-Baseline.step"},
        step_sources={("baseline", "wing"): "reference_step"},
    )
    readback = report["revisions"]["baseline"]["wing"]["step_readback"]
    assert readback["source"] == "reference_step"
    assert readback["status"] == "measured"
    assert readback["sampled_surface_within_tolerance"]
    assert readback["agreed_locks_within_tolerance"]
    assert readback["maximum_between_ribs_definition_to_step_deviation"]["between_defining_ribs"]
    assert readback["hinge"]["axis_inside_solid"]
    assert report["revisions"]["baseline"]["tail"]["step_readback"]["status"] == "unavailable"



def _bent_wing(baseline):
    return replace(
        baseline.wing,
        sweep_mm=100,
        rake_mm=25,
        rise_mm=60,
        tip_chord_mm=110,
        section_orientation="front_normal",
    )


def test_generated_uncalibrated_bend_exceeds_physical_span_lock(tmp_path):
    baseline = load_definition(BASELINE)
    changed = replace(baseline, name="uncalibrated bend", wing=_bent_wing(baseline))
    report = write_evidence(tmp_path, {"baseline": baseline, "candidate": changed}, FIXTURE)
    wing = report["revisions"]["candidate"]["wing"]["generated"]
    assert wing["projected_span_mm"] == 600
    assert wing["physical_projected_extent_mm"] > 600.25
    assert not wing["baseline_projected_span_within_tolerance"]
    assert not wing["agreed_locks_within_tolerance"]


@pytest.mark.slow
def test_fallback_bend_matches_independent_cad_cuts_and_physical_span(tmp_path):
    cq = pytest.importorskip("cadquery")
    baseline = load_definition(BASELINE)
    wing = calibrate_projected_span(_bent_wing(baseline), FIXTURE)
    path = tmp_path / "fallback.step"
    _fallback_step(cq, wing, FIXTURE, path)
    readback, _ = _measure(
        cq, _load_solid(cq, path), wing, baseline.wing, FIXTURE, baseline.tolerance_mm,
        fallback=True,
    )
    assert readback["sampled_surface_within_tolerance"]
    held_out = readback["maximum_held_out_definition_to_step_deviation"]
    assert held_out["station_fraction"] not in readback["fallback_loft_input_fractions"]
    assert held_out["distance_mm"] < baseline.tolerance_mm
    assert readback["baseline_root_attachment_within_tolerance"]
    assert readback["leading_edge_projected_span_within_tolerance"]
    assert readback["baseline_bounding_projected_extent_within_tolerance"]
    assert readback["tip_cap_front_tilt_deg"] > 15
