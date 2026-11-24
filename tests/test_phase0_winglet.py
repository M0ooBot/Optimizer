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


@pytest.fixture(scope="module")
def winglet():
    parent = load_definition(WING_B).wing
    definition = WingletDefinition(bend_length_mm=70.0)
    return BlendedWinglet(parent, definition, _family(0.12, np.linspace(55, 205, 31)),
                          _family(0.10, np.linspace(18, 100, 42)))


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
    assert winglet.kappa(winglet.L) == 0.0
    chord = np.array([winglet.chord(float(t)) for t in s])
    assert np.all(np.diff(chord) < 0)
    assert chord[-1] == pytest.approx(32.0)


def test_trailing_edge_holds_its_thickness_and_height_is_met(winglet):
    for s in winglet.stations(21, 11):
        frame = winglet.frame(float(s))
        upper, lower = frame.section()[0], frame.section()[-2]
        assert np.linalg.norm(upper - lower) == pytest.approx(TE_MM, abs=1e-9)
    tip = winglet.frame(winglet.s_tip)
    assert tip.le[1] - winglet.le0[1] == pytest.approx(65.0, abs=1e-9)


def test_toe_out_points_the_nose_outboard(winglet):
    tip = winglet.frame(winglet.s_tip)
    away_from_upper_skin = -np.array([math.sin(tip.phi), math.cos(tip.phi), 0.0])
    nose_to_tail = tip.le - tip.te
    toe = math.degrees(math.asin(nose_to_tail @ away_from_upper_skin / np.linalg.norm(nose_to_tail)))
    assert toe == pytest.approx(2.0, abs=1e-9)


def test_span_lock_solves_the_bend_length():
    parent = load_definition(WING_B).wing
    solved = BlendedWinglet(parent, WingletDefinition(), _family(0.12, np.linspace(55, 205, 31)),
                            _family(0.10, np.linspace(18, 100, 42)))
    assert solved.projected_span(True) == pytest.approx(600.0, abs=1e-5)
    assert solved.projected_span(False) < 600.0


def test_out_of_family_chord_is_refused():
    family = _family(0.10, np.linspace(30, 60, 7))
    with pytest.raises(ValueError, match="outside the family"):
        family.at(20.0)
