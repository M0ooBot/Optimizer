"""Blended winglet law (phase0_winglet): join continuity, locks, TE, toe sign, profiles."""

import math
from pathlib import Path

import numpy as np
import pytest

from planeopt.phase0_geometry import load_definition
from planeopt.phase0_winglet import (
    BlendedWinglet,
    ProfileFamily,
    WingletDefinition,
    decompose,
    naca_symmetric,
    nose_radius,
)

REPO = Path(__file__).resolve().parents[1]
WING_B = REPO / "runs/wing-b-spar/shape/constrained/redistributed/definition.json"
TE_MM = 0.4


def _blunt(loop, te_fraction):
    """Open a closed loop and add a linear thickness ramp: a stand-in for the Converter's cut
    that keeps these tests independent of it."""
    pts = np.asarray(loop, float)  # both ends at the TE, x = 1; the ramp opens them
    nose = int(np.argmin(pts[:, 0]))
    sign = np.where(np.arange(len(pts)) <= nose, 1.0, -1.0)
    return np.column_stack([pts[:, 0], pts[:, 1] + sign * pts[:, 0] * te_fraction / 2])


def _family(thickness, chords):
    return ProfileFamily(lambda c: _blunt(naca_symmetric(thickness), TE_MM / c), chords)


def _winglet(**overrides):
    return BlendedWinglet(load_definition(WING_B).wing, WingletDefinition(**overrides),
                          _family(0.12, np.linspace(55, 205, 31)), _family(0.10, np.linspace(18, 100, 42)))


@pytest.fixture(scope="module")
def winglet():
    return _winglet(attach_fraction=0.87)


def _one_sided(f, h, side):
    v = np.array([f(side * k * h) for k in range(6)])
    d1 = side * (-25 * v[0] + 48 * v[1] - 36 * v[2] + 16 * v[3] - 3 * v[4]) / (12 * h)
    d2 = (45 * v[0] - 154 * v[1] + 214 * v[2] - 156 * v[3] + 61 * v[4] - 10 * v[5]) / (12 * h * h)
    return v[0], d1, d2


def test_decompose_recovers_naca_nose_radius():
    profile = decompose(naca_symmetric(0.10))
    assert profile.le_radius == pytest.approx(1.1019 * 0.10**2, rel=2e-3)
    assert 2 * profile.half_thickness.max() == pytest.approx(0.10, rel=1e-3)
    assert nose_radius(naca_symmetric(0.12)) == pytest.approx(1.1019 * 0.12**2, rel=2e-3)


@pytest.mark.parametrize("index, upper", [(60, True), (160, False), (-1, True)])
def test_skin_is_curvature_continuous_with_the_parent(winglet, index, upper):
    def parent_side(s):
        frame = winglet.parent_frame(winglet.parent_fraction(s), s)
        x, y = (frame.profile.upper() if upper else frame.profile.lower())[index]
        return frame.point(x, y)

    def winglet_side(s):
        frame = winglet.frame(s)
        x, y = (frame.profile.upper() if upper else frame.profile.lower())[index]
        return frame.point(x, y)

    p0, a1, a2 = _one_sided(parent_side, 0.2, -1)
    q0, b1, b2 = _one_sided(winglet_side, 0.2, +1)
    assert np.linalg.norm(p0 - q0) < 1e-6
    assert np.linalg.norm(a1 - b1) < 1e-6
    assert np.linalg.norm(a2 - b2) < 1e-5 * max(np.linalg.norm(a2), 1e-3)


def test_bend_reaches_cant_with_rising_dihedral(winglet):
    s = np.linspace(0, winglet.s_tip, 400)
    gamma = np.array([winglet.gamma(float(t)) for t in s])
    assert np.all(np.diff(gamma) >= -1e-12)
    assert math.degrees(gamma[-1]) == pytest.approx(75.0, abs=1e-9)
    assert winglet.kappa(0.0) == pytest.approx(winglet.kappa0)
    assert winglet.kappa(winglet.s_tip) == pytest.approx(0.0, abs=1e-15)


def test_drawn_edges_give_a_tapering_chord_and_a_one_way_trailing_edge(winglet):
    s = np.linspace(0, winglet.s_tip, 200)
    chord = np.array([winglet.frame(float(t)).chord_mm for t in s])
    assert np.all(np.diff(chord) < 0)
    assert chord[-1] == pytest.approx(32.0, abs=1e-6)
    te_z = np.array([winglet.frame(float(t)).te[2] for t in s])
    assert np.all(np.diff(te_z) > 0) or np.all(np.diff(te_z) < 0)


@pytest.mark.parametrize("edge", ["le", "te"])
def test_edges_turn_one_way_seen_from_outboard(winglet, edge):
    """Each edge's direction seen from outboard (up against aft) changes one way only, and its
    overall sweep is the requested one; the tip direction is free."""
    s = winglet.stations(301)
    pts = np.array([winglet.frame(float(t)).le if edge == "le" else winglet.frame(float(t)).te for t in s])
    slope = np.gradient(pts[:, 2], s) / np.gradient(pts[:, 1], s)  # dZ/dy
    seen = np.gradient(pts[:, 1], s) > 0.3 * np.gradient(pts[:, 1], s).max()
    steps = np.diff(slope[seen])
    assert np.all(steps >= -1e-9) or np.all(steps <= 1e-9)
    sweep = winglet.d.sweep_deg if edge == "le" else None
    if sweep is not None:
        height = pts[-1, 1] - pts[0, 1]
        assert math.degrees(math.atan(-(pts[-1, 2] - pts[0, 2]) / height)) == pytest.approx(sweep, abs=1e-6)


def test_trailing_edge_holds_its_thickness_and_height_is_met(winglet):
    for s in winglet.stations(21):
        frame = winglet.frame(float(s))
        upper, lower = frame.section()[0], frame.section()[-2]
        assert np.linalg.norm(upper - lower) == pytest.approx(TE_MM, abs=1e-9)
    tip = winglet.frame(winglet.s_tip)
    assert tip.le[1] - winglet.le0[1] == pytest.approx(65.0, abs=1e-6)


def test_toe_out_points_the_nose_outboard(winglet):
    tip = winglet.frame(winglet.s_tip)
    away_from_upper_skin = -np.array([math.sin(tip.phi), math.cos(tip.phi), 0.0])
    nose_to_tail = tip.le - tip.te
    toe = math.degrees(math.asin(nose_to_tail @ away_from_upper_skin / np.linalg.norm(nose_to_tail)))
    assert toe == pytest.approx(2.0, abs=1e-9)


def test_span_lock_solves_the_attachment_outboard_of_the_hinge():
    solved = _winglet()
    assert solved.projected_span(True) == pytest.approx(600.0, abs=1e-5)
    assert solved.f0 >= solved.min_attach_fraction()


def _top_outline_sag(w):
    g0 = w.gamma0
    normal, along = np.array([math.sin(g0), math.cos(g0)]), np.array([math.cos(g0), math.sin(g0)])

    def top(frame):
        p = np.column_stack([-frame.section()[:, 0], frame.section()[:, 1]])
        k = np.argmax(p @ normal)
        return p[k] @ along, p[k] @ normal

    wing = np.array([top(w.parent_frame(w.parent_fraction(-h))) for h in (2.0, 1.0, 0.0)])
    winglet = np.array([top(w.frame(float(s))) for s in np.linspace(0, 0.6 * w.s_tip, 61)])
    slope = np.polyfit(wing[:, 0], wing[:, 1], 1)[0]
    return float((wing[-1, 1] + slope * (winglet[:, 0] - wing[-1, 0]) - winglet[:, 1]).max())


def test_shape_changes_that_follow_the_turn_sag_less_in_front_view():
    by_distance = _winglet(attach_fraction=0.87, shape_follows_turn=False, ramp_in=0.5, ramp_out=0.5)
    by_turn = _winglet(attach_fraction=0.87, shape_follows_turn=True, ramp_in=0.5, ramp_out=0.5)
    assert _top_outline_sag(by_turn) < 0.5 * _top_outline_sag(by_distance)


def test_a_bend_can_finish_early_and_leave_a_straight_leg():
    w = _winglet(attach_fraction=0.87, bend_fraction=0.7)
    assert w.s_bend == pytest.approx(0.7 * w.s_tip)
    assert w.kappa(0.8 * w.s_tip) == 0.0
    assert w.gamma(0.8 * w.s_tip) == pytest.approx(w.gamma_c, abs=1e-12)
    assert w.frame(w.s_tip).le[1] - w.le0[1] == pytest.approx(65.0, abs=1e-6)


def test_solving_the_bend_fits_a_winglet_that_a_full_length_bend_cannot():
    with pytest.raises(ValueError, match="hinge clearance"):
        _winglet(height_mm=70.0, cant_deg=70.0)
    w = _winglet(height_mm=70.0, cant_deg=70.0, bend_fraction=None)
    assert 0 < w.d.bend_fraction < 1
    assert w.f0 == pytest.approx(w.min_attach_fraction(), abs=1e-9)
    assert w.projected_span(True) == pytest.approx(600.0, abs=1e-4)


def test_a_bend_that_cannot_fit_names_the_way_out():
    with pytest.raises(ValueError, match="toward the winglet root"):
        _winglet(ramp_in=0.6, ramp_out=0.3)


def test_ramps_must_fit_in_the_winglet():
    with pytest.raises(ValueError, match="together at most 1"):
        WingletDefinition(ramp_in=0.7, ramp_out=0.4)


def test_out_of_family_chord_is_refused():
    family = _family(0.10, np.linspace(30, 60, 7))
    with pytest.raises(ValueError, match="outside the family"):
        family.at(20.0)
