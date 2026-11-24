"""Blended winglet that continues a Phase 0 wing past an attachment station.

The winglet is not a separate body blended onto the wing. Every quantity that
places a Phase 0 section (the leading-edge path, the section plane, chord,
twist and profile) is continued past the attachment by a law whose value,
first and second arc-length derivatives match the parent there. The surface is
therefore curvature-continuous with the parent by construction.

Coordinates follow phase0_geometry: -X outboard, +Y up, +Z forward, millimetres.
The common parameter s is the arc length of the leading edge in the front (X-Y)
view, measured outboard from the attachment.

Laws for s >= 0:

- Front-view bend. The leading edge turns from the parent's local dihedral to the
  cant angle with curvature kappa(s) = kappa0 (1 - W(x)) + K W'(x), x = s / L, where W
  is the quintic smoothstep. Curvature starts at the parent's value, rises to one
  peak and falls to zero with zero slope at s = L (a railway-style transition, not
  a circular arc), then the leg runs straight at the cant angle.
- Section plane. phi(s) = W(s / L) gamma(s): the parent's spanwise-parallel planes
  (phi = 0) rotate until they are normal to the front-view leading edge.
- Everything else (aft rate dZ/ds, chord, twist, profile morph) blends the parent's
  second-order Taylor expansion into a smooth target law with weight W(s / L).

The profile is split into camber and thickness on a fixed cosine grid. The
morph fades camber and blends the thickness without its trailing-edge ramp, then
adds back a ramp that holds the trailing edge at its specified thickness in mm.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass, replace

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from planeopt.phase0_geometry import SurfaceDefinition, _twist_deg, section_frame

N_BETA = 241
_BETA = np.linspace(0.0, math.pi, N_BETA)
XI = (1.0 - np.cos(_BETA)) / 2.0  # normalized chord position, 0 at the nose, 1 at the TE
_GAUSS_X, _GAUSS_W = np.polynomial.legendre.leggauss(48)


def smoothstep(x):
    """Quintic smoothstep W: W, W', W'' are 0, 0, 0 at x = 0 and 1, 0, 0 at x = 1."""
    x = np.clip(x, 0.0, 1.0)
    return x**3 * (10 - 15 * x + 6 * x**2)


def _smoothstep_integral(x):
    x = np.clip(x, 0.0, 1.0)
    return 2.5 * x**4 - 3 * x**5 + x**6


# ----------------------------------------------------------------- profiles


@dataclass(frozen=True)
class Profile:
    """One section in chord units, about the 2D origin the frame's LE point uses."""

    x0: float  # nose
    x1: float  # trailing edge
    camber: np.ndarray  # on XI
    half_thickness: np.ndarray  # on XI, includes the TE ramp
    le_radius: float  # nose radius, chord units

    def x(self):
        return self.x0 + (self.x1 - self.x0) * XI

    def upper(self):
        return np.column_stack([self.x(), self.camber + self.half_thickness])

    def lower(self):
        return np.column_stack([self.x(), self.camber - self.half_thickness])


def decompose(points: Sequence[Sequence[float]]) -> Profile:
    """Split an airfoil polyline (closed, or open at a blunt TE) into camber and thickness.

    The outline is a parametric cubic spline through its points in chord-length
    order, the way a CAD spline through the same points runs. The nose is the
    spline's foremost point (a sampled outline rarely has a vertex exactly there)
    and its radius is the spline's curvature there. Thickness is measured
    vertically on the cosine grid XI, so the surfaces are camber +/- half
    thickness exactly at the grid points.
    """
    t, xs, ys, t_nose, radius = _outline(points)
    x_nose, y_nose = float(xs(t_nose)), float(ys(t_nose))
    x_te = min(float(xs(t[0])), float(xs(t[-1])))
    grid = x_nose + (x_te - x_nose) * XI

    def side(t_from, t_to):
        # Knots walking away from the nose; each grid station is bracketed between two of
        # them and solved with brentq (spline root finders can miss a root on a knot).
        inner = t[(t > min(t_from, t_to)) & (t < max(t_from, t_to))]
        knots = np.concatenate([[t_from], inner if t_to > t_from else inner[::-1], [t_to]])
        x_knots = xs(knots) - x_nose
        values = np.empty_like(grid)
        values[0] = y_nose
        for k, x in enumerate(grid[1:], start=1):
            hits = np.nonzero((x_knots[:-1] - (x - x_nose)) * (x_knots[1:] - (x - x_nose)) <= 0)[0]
            if not len(hits):
                raise ValueError(f"The outline does not reach x = {x:.9f} on one side")
            a, b = knots[hits[0]], knots[hits[0] + 1]
            root = a if xs(a) == x else brentq(lambda r, x=x: float(xs(r)) - x, a, b, xtol=1e-15)
            values[k] = ys(root)
        return values

    first, second = side(t_nose, t[0]), side(t_nose, t[-1])
    upper, lower = (first, second) if first.mean() > second.mean() else (second, first)
    return Profile(x_nose, x_te, (upper + lower) / 2, (upper - lower) / 2, radius)


def _outline(points):
    """Parametric cubic spline through an outline, its foremost point and the radius there."""
    pts = np.asarray(points, float)
    if np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]
    step = np.hypot(*np.diff(pts, axis=0).T)
    pts = pts[np.concatenate([[True], step > 1e-12])]
    t = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(pts, axis=0).T))])
    xs, ys = CubicSpline(t, pts[:, 0]), CubicSpline(t, pts[:, 1])
    vertex = int(np.argmin(pts[:, 0]))
    lo, hi = t[max(vertex - 3, 0)], t[min(vertex + 3, len(t) - 1)]
    candidates = [r for r in xs.derivative().roots() if lo <= r <= hi]
    t_nose = min(candidates, key=lambda r: xs(r)) if candidates else t[vertex]
    dx, dy, ddx, ddy = xs(t_nose, 1), ys(t_nose, 1), xs(t_nose, 2), ys(t_nose, 2)
    radius = float((dx * dx + dy * dy) ** 1.5 / abs(dx * ddy - dy * ddx))
    return t, xs, ys, t_nose, radius


def nose_radius(points) -> float:
    """Radius at the foremost point of the spline through an outline, in its units."""
    return _outline(points)[4]


class ProfileFamily:
    """Profiles that depend smoothly on chord: a blunt TE in mm is a larger fraction of
    a smaller chord. Built from exact sections at a few chords and splined in chord,
    so the family is C2 in chord even when the source tool cuts on a polyline."""

    def __init__(self, build: Callable[[float], Sequence[Sequence[float]]],
                 chords: Sequence[float]):
        self.chords = np.asarray(sorted(chords), float)
        profiles = [decompose(build(float(c))) for c in self.chords]
        self.source = dict(zip(self.chords.tolist(), profiles))

        def stack(get):
            return np.array([get(p) for p in profiles])

        self._x0 = CubicSpline(self.chords, stack(lambda p: p.x0))
        self._x1 = CubicSpline(self.chords, stack(lambda p: p.x1))
        self._camber = CubicSpline(self.chords, stack(lambda p: p.camber), axis=0)
        self._thick = CubicSpline(self.chords, stack(lambda p: p.half_thickness), axis=0)
        self._radius = CubicSpline(self.chords, stack(lambda p: p.le_radius))

    def at(self, chord_mm: float) -> Profile:
        if not self.chords[0] - 1e-9 <= chord_mm <= self.chords[-1] + 1e-9:
            raise ValueError(f"Chord {chord_mm:.3f} mm is outside the family "
                             f"{self.chords[0]:g}-{self.chords[-1]:g} mm")
        return Profile(float(self._x0(chord_mm)), float(self._x1(chord_mm)),
                       self._camber(chord_mm), self._thick(chord_mm),
                       float(self._radius(chord_mm)))


def naca_symmetric(thickness_ratio: float, points_per_side: int = 401) -> np.ndarray:
    """Closed-trailing-edge NACA 00xx loop, upper TE -> nose -> lower TE, chord 1."""
    beta = np.linspace(0.0, math.pi, points_per_side)
    x = (1 - np.cos(beta)) / 2
    t = 5 * thickness_ratio * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2
                               + 0.2843 * x**3 - 0.1036 * x**4)
    upper = np.column_stack([x, t])[::-1]
    lower = np.column_stack([x, -t])[1:]
    return np.vstack([upper, lower])


# ----------------------------------------------------------------- definition


@dataclass(frozen=True)
class WingletDefinition:
    cant_deg: float = 75.0  # front-view angle of the leading edge at the winglet tip, above horizontal
    height_mm: float = 65.0  # LE rise from the attachment to the winglet tip station
    sweep_deg: float = 35.0  # LE sweep at the winglet tip, against its front-view direction
    tip_chord_mm: float = 32.0  # the trailing-edge law is solved to give this chord at the tip
    toe_out_deg: float = 2.0  # nose outboard is positive
    target_thickness_ratio: float = 0.10
    trailing_edge_mm: float = 0.4
    projected_span_mm: float = 600.0  # outboard skin extent, winglet and cap allowance included
    cap_allowance_mm: float = 12.0  # tip extended this far along its tangent for the span lock
    # Bend distribution, as fractions of the winglet's leading-edge length: curvature ramps up
    # from Wing B's over ramp_in, holds its peak, and ramps down to zero over the last ramp_out.
    # 0.5 / 0.5 is a pure bell; short ramps spread the turn evenly; unequal ramps move the peak.
    ramp_in: float = 0.35
    ramp_out: float = 0.35
    hinge_clearance_mm: float = 10.0  # the winglet starts at least this far outboard of the hinge end
    # True: sweep, taper, twist, morph and plane rotation advance with the turn already made,
    # so each section changes while the skin is already rising. False: with distance along it.
    shape_follows_turn: bool = True
    attach_fraction: float | None = None  # solved from the span lock when None

    def __post_init__(self):
        if not 0 < self.cant_deg < 90:
            raise ValueError("cant_deg must be between 0 and 90")
        if min(self.height_mm, self.tip_chord_mm, self.target_thickness_ratio,
               self.projected_span_mm) <= 0:
            raise ValueError("Lengths and thickness must be positive")
        if not (0 < self.ramp_in and 0 < self.ramp_out and self.ramp_in + self.ramp_out <= 1):
            raise ValueError("ramp_in and ramp_out must be positive and together at most 1")
        if self.attach_fraction is not None and not 0 < self.attach_fraction < 1:
            raise ValueError("attach_fraction must lie inside the parent span")


@dataclass(frozen=True)
class Frame:
    s: float
    le: np.ndarray
    chord_axis: np.ndarray
    up_axis: np.ndarray
    plane_normal: np.ndarray  # outboard-pointing section-plane normal
    chord_mm: float
    profile: Profile
    gamma: float  # front-view LE angle above horizontal, rad
    phi: float  # section-plane tilt, rad
    twist_deg: float

    def point(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        return (self.le + self.chord_mm * (np.multiply.outer(x, self.chord_axis)
                                           + np.multiply.outer(y, self.up_axis)))

    def section(self):
        """Closed outline: upper TE -> nose -> lower TE, then the straight TE line back."""
        up, lo = self.profile.upper(), self.profile.lower()
        loop = np.vstack([up[::-1], lo[1:], up[-1:]])
        return self.point(loop[:, 0], loop[:, 1])

    @property
    def te(self):
        return self.point(self.profile.x1, self.profile.camber[-1])


@dataclass
class _Taylor:
    value: float
    d1: float
    d2: float

    def __call__(self, s):
        return self.value + self.d1 * s + 0.5 * self.d2 * s * s


def _quad(f, a, b):
    nodes = a + 0.5 * (b - a) * (_GAUSS_X + 1)
    return 0.5 * (b - a) * float(np.dot(_GAUSS_W, [f(t) for t in nodes]))


class BlendedWinglet:
    """The parent wing continued into a winglet. frame(s) gives the section frame at s, the
    front-view leading-edge arc length from the attachment; negative s returns the parent.

    The whole winglet is the transition: the bend, section-plane rotation, sweep, twist and
    profile morph all run from the attachment to the tip. Solving order: the winglet length
    meets the height; the attachment station meets the projected span; the trailing-edge sweep
    meets the tip chord.
    """

    def __init__(self, parent: SurfaceDefinition, definition: WingletDefinition,
                 parent_profiles: ProfileFamily, target_profiles: ProfileFamily):
        if parent.span_curve != "full_span" or parent.section_orientation != "parallel":
            raise ValueError("The winglet law continues a full_span wing with parallel sections")
        for name, family in (("parent", parent_profiles), ("target", target_profiles)):
            for chord, profile in family.source.items():
                expected = definition.trailing_edge_mm / 2 / chord
                if abs(profile.half_thickness[-1] - expected) > 1e-9 * max(expected, 1e-3):
                    raise ValueError(f"The {name} profile at {chord:g} mm ends "
                                     f"{2 * profile.half_thickness[-1] * chord:.9f} mm thick, not "
                                     f"the {definition.trailing_edge_mm:g} mm trailing edge")
        self.parent = parent
        self.d = definition
        self.parent_profiles = parent_profiles
        self.target_profiles = target_profiles
        self.guide_span = parent.projected_span_mm - parent.tip_span_trim_mm
        if definition.attach_fraction is None:
            self.d = replace(definition, attach_fraction=self._solve_attachment())
        self._configure(self.d.attach_fraction)

    # -------------------------------------------------------------- parent side

    def _parent_q(self, f):
        frame = section_frame(self.parent, f)
        return np.array([*frame.le, frame.chord_mm, _twist_deg(self.parent, f), frame.te[2]])

    def _join(self, f0):
        h = 1e-3
        stencil = np.array([self._parent_q(f0 + k * h) for k in (-2, -1, 0, 1, 2)])
        q = stencil[2]
        q_f = (stencil[0] - 8 * stencil[1] + 8 * stencil[3] - stencil[4]) / (12 * h)
        q_ff = (-stencil[0] + 16 * stencil[1] - 30 * stencil[2] + 16 * stencil[3]
                - stencil[4]) / (12 * h * h)
        x_f, y_f, x_ff, y_ff = q_f[0], q_f[1], q_ff[0], q_ff[1]
        g = math.hypot(x_f, y_f)  # ds/df
        g_f = (x_f * x_ff + y_f * y_ff) / g
        q_s = q_f / g
        q_ss = (q_ff - q_s * g_f) / (g * g)
        self.f0, self.ds_df = f0, g
        self.le0 = q[:3]
        self.gamma0 = math.atan2(y_f, -x_f)
        self.kappa0 = (-x_f * y_ff + y_f * x_ff) / g**3
        self.le_aft = _Taylor(q_s[2], q_ss[2], 0.0)  # dZ/ds of the leading edge and its slope
        self.te_aft = _Taylor(q_s[5], q_ss[5], 0.0)  # the same for the trailing edge
        self.te_z0 = q[5]
        self.chord_p = _Taylor(q[3], q_s[3], q_ss[3])
        self.twist_p = _Taylor(q[4], q_s[4], q_ss[4])

    def _parent_speed(self, f):
        frame = section_frame(self.parent, f)
        slope = frame.leading_edge_tangent[1] / -frame.leading_edge_tangent[0]
        return self.guide_span * math.hypot(1.0, slope)

    def parent_fraction(self, s):
        """Parent span fraction at front-view arc length s <= 0 (Newton on the exact arc)."""
        f = self.f0 + s / self.ds_df
        for _ in range(30):
            step = (_quad(self._parent_speed, self.f0, f) - s) / self._parent_speed(f)
            f -= step
            if abs(step) < 1e-15:
                break
        return f

    # -------------------------------------------------------------- winglet side

    def _configure(self, f0):
        self._join(f0)
        self.gamma_c = math.radians(self.d.cant_deg)
        lo, hi = 1.0, 600.0
        for _ in range(100):  # the winglet length that meets the height
            mid = 0.5 * (lo + hi)
            self._set_length(mid)
            lo, hi = (lo, mid) if self._front(mid)[1] - self.le0[1] > self.d.height_mm else (mid, hi)
            if hi - lo < 1e-10:
                break
        self._set_length(0.5 * (lo + hi))
        if self.kappa_peak <= 0:
            raise ValueError("The winglet turns less than Wing B already bends; raise the cant")
        # The trailing-edge sweep that gives the tip chord: Z_te(S) is linear in tan(sweep).
        S = self.s_tip
        base = self.te_z0 + _quad(lambda t: (1 - self._weight(t)) * self.te_aft(t), 0, S)
        weight = _quad(self._weight, 0, S)
        cos_toe = math.cos(math.radians(self.d.toe_out_deg))
        self.te_tan = (base - self._z(S) + self.d.tip_chord_mm * cos_toe) / weight

    def _set_length(self, S):
        d = self.d
        self.s_tip = float(S)
        turn = self.gamma_c - self.gamma0
        self.kappa_peak = ((turn / S - self.kappa0 * d.ramp_in / 2)
                           / (1 - (d.ramp_in + d.ramp_out) / 2))

    def kappa(self, s):
        d, x = self.d, s / self.s_tip
        rise = float(smoothstep(x / d.ramp_in))
        fall = float(smoothstep((1 - x) / d.ramp_out))
        return self.kappa0 * (1 - rise) + self.kappa_peak * min(rise, fall)

    def gamma(self, s):
        """Front-view LE angle: the closed-form integral of kappa."""
        d, S = self.d, self.s_tip
        x = min(max(s / S, 0.0), 1.0)
        a, b = d.ramp_in, d.ramp_out
        xa = min(x, a)
        k0 = xa - a * _smoothstep_integral(xa / a)
        peak = a * _smoothstep_integral(xa / a) + max(0.0, min(x, 1 - b) - a)
        if x > 1 - b:
            peak += b * (0.5 - _smoothstep_integral((1 - x) / b))
        return float(self.gamma0 + S * (self.kappa0 * k0 + self.kappa_peak * peak))

    def _front(self, s):
        """Front-view LE (X, Y) at s, by Gauss quadrature of the tangent on each piece of the
        curvature law (its ramps join with a jump in the third derivative of the angle)."""
        S, d = self.s_tip, self.d
        cuts = [0.0] + [c for c in (d.ramp_in * S, (1 - d.ramp_out) * S) if 0 < c < s] + [s]
        x, y = self.le0[0], self.le0[1]
        for a, b in itertools.pairwise(cuts):
            nodes = a + 0.5 * (b - a) * (_GAUSS_X + 1)
            angles = np.array([self.gamma(t) for t in nodes])
            x -= 0.5 * (b - a) * np.dot(_GAUSS_W, np.cos(angles))
            y += 0.5 * (b - a) * np.dot(_GAUSS_W, np.sin(angles))
        return np.array([x, y])

    def _weight(self, s):
        """Blend weight for everything but the bend. Following the turn, it is W of the fraction
        of the turn made so far; W' and W'' vanish at both ends either way, so the join keeps
        its curvature continuity."""
        if not self.d.shape_follows_turn:
            return float(smoothstep(s / self.s_tip))
        progress = (self.gamma(s) - self.gamma0) / (self.gamma_c - self.gamma0)
        return float(smoothstep(progress))

    def le_aft_rate(self, s):
        w = self._weight(s)
        return (1 - w) * self.le_aft(s) - w * math.tan(math.radians(self.d.sweep_deg))

    def te_aft_rate(self, s):
        w = self._weight(s)
        return (1 - w) * self.te_aft(s) - w * self.te_tan

    def _z(self, s):
        return self.le0[2] + _quad(self.le_aft_rate, 0, s)

    def _te_z(self, s):
        return self.te_z0 + _quad(self.te_aft_rate, 0, s)

    def twist_deg(self, s):
        w = self._weight(s)
        return (1 - w) * self.twist_p(s) - w * self.d.toe_out_deg

    def chord(self, s):
        """From the two drawn edges: the chord spans their aft separation at the section's twist."""
        return (self._z(s) - self._te_z(s)) / math.cos(math.radians(self.twist_deg(s)))

    def morph(self, s):
        return self._weight(s)

    def profile(self, s, chord=None):
        """Camber fades and thickness blends; the TE ramp is re-added at the specified mm."""
        m = self.morph(s)
        chord = self.chord(s) if chord is None else chord
        par = self.parent_profiles.at(self.chord_p(s))
        tgt = self.target_profiles.at(chord)
        te_half = self.d.trailing_edge_mm / 2 / chord
        par_sharp = par.half_thickness - XI * par.half_thickness[-1]
        tgt_sharp = tgt.half_thickness - XI * tgt.half_thickness[-1]
        return Profile((1 - m) * par.x0 + m * tgt.x0, (1 - m) * par.x1 + m * tgt.x1,
                       (1 - m) * par.camber + m * tgt.camber,
                       (1 - m) * par_sharp + m * tgt_sharp + XI * te_half, float("nan"))

    def frame(self, s) -> Frame:
        if s < 0:
            return self.parent_frame(self.parent_fraction(s), s)
        front = self._front(s)
        le = np.array([front[0], front[1], self._z(s)])
        gamma = self.gamma(s)
        phi = self._weight(s) * gamma
        theta_deg = self.twist_deg(s)
        theta = math.radians(theta_deg)
        up = np.array([math.sin(phi), math.cos(phi), 0.0])
        chord_axis = np.array([-up[0] * math.sin(theta), -up[1] * math.sin(theta), -math.cos(theta)])
        up_axis = np.array([up[0] * math.cos(theta), up[1] * math.cos(theta), -math.sin(theta)])
        normal = np.array([-math.cos(phi), math.sin(phi), 0.0])
        chord = (le[2] - self._te_z(s)) / math.cos(theta)
        return Frame(float(s), le, chord_axis, up_axis, normal, chord, self.profile(s, chord),
                     gamma, phi, theta_deg)

    def parent_frame(self, f, s=float("nan")) -> Frame:
        pf = section_frame(self.parent, f)
        return Frame(float(s), np.array(pf.le), np.array(pf.chord_axis), np.array(pf.up_axis),
                     np.array(pf.span_tangent), pf.chord_mm, self.parent_profiles.at(pf.chord_mm),
                     math.atan2(pf.leading_edge_tangent[1], -pf.leading_edge_tangent[0]), 0.0,
                     _twist_deg(self.parent, f))

    def stations(self, n=161, extend_mm=0.0):
        s = self.s_tip * (1 - np.cos(np.linspace(0, math.pi, n))) / 2
        if extend_mm:
            s = np.concatenate([s, self.s_tip + np.linspace(0, extend_mm, 9)[1:]])
        return s

    def extended_frame(self, s):
        """Past the tip station the winglet continues along its tip tangent (cap allowance)."""
        if s <= self.s_tip:
            return self.frame(s)
        tip = self.frame(self.s_tip)
        direction = np.array([-math.cos(self.gamma_c), math.sin(self.gamma_c),
                              -math.tan(math.radians(self.d.sweep_deg))])
        return replace(tip, s=float(s), le=tip.le + direction * (s - self.s_tip))

    def projected_span(self, include_allowance=True):
        extra = self.d.cap_allowance_mm if include_allowance else 0.0
        outboard = max(float(-self.extended_frame(float(s)).section()[:, 0].min())
                       for s in self.stations(41, extra))
        return outboard + self.parent.root_le_mm[0]

    def min_attach_fraction(self):
        hinge_end = self.parent.hinge.end * self.guide_span
        return (hinge_end + self.d.hinge_clearance_mm) / self.guide_span

    def _solve_attachment(self):
        """Inboard starts leave more room, so the span falls as the attachment moves inboard."""
        target = self.d.projected_span_mm
        lo, hi = self.min_attach_fraction(), 0.98
        self._configure(lo)
        if self.projected_span() > target:
            raise ValueError(
                f"Starting at the hinge clearance ({lo * self.guide_span:.1f} mm) the winglet still "
                f"reaches {self.projected_span():.2f} mm; lower the height, steepen the cant or "
                "move the peak of the bend toward the winglet root (shorter ramp_in)")
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            self._configure(mid)
            lo, hi = (lo, mid) if self.projected_span() > target else (mid, hi)
            if hi - lo < 1e-11:
                break
        return 0.5 * (lo + hi)

    def to_dict(self):
        out = asdict(self.d)
        s = np.linspace(0, self.s_tip, 2001)
        kappa = np.array([self.kappa(float(t)) for t in s])
        out.update(attach_mm=self.f0 * self.guide_span, length_mm=self.s_tip,
                   peak_curvature_per_mm=float(kappa.max()), peak_at_fraction=float(s[kappa.argmax()] / self.s_tip),
                   tightest_bend_radius_mm=float(1 / kappa.max()),
                   trailing_edge_tip_sweep_deg=math.degrees(math.atan(self.te_tan)),
                   parent_dihedral_deg=math.degrees(self.gamma0), parent_curvature_per_mm=self.kappa0)
        return out
