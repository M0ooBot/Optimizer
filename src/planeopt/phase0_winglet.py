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
    attach_fraction: float = 0.9
    cant_deg: float = 75.0  # front-view angle of the straight leg above horizontal
    height_mm: float = 65.0  # LE rise from the attachment to the winglet tip station
    sweep_deg: float = 35.0  # LE sweep of the leg, against its front-view direction
    tip_chord_mm: float = 32.0  # the target chord tapers linearly from the join chord to this
    toe_out_deg: float = 2.0  # nose outboard is positive
    target_thickness_ratio: float = 0.10
    trailing_edge_mm: float = 0.4
    projected_span_mm: float = 600.0  # outboard skin extent, winglet and cap allowance included
    cap_allowance_mm: float = 12.0  # leg extended this far when checking the span lock
    bend_length_mm: float | None = None  # solved from the span lock when None

    def __post_init__(self):
        if not 0 < self.attach_fraction < 1:
            raise ValueError("attach_fraction must lie inside the parent span")
        if not 0 < self.cant_deg < 90:
            raise ValueError("cant_deg must be between 0 and 90")
        if min(self.height_mm, self.tip_chord_mm,
               self.target_thickness_ratio, self.projected_span_mm) <= 0:
            raise ValueError("Lengths and thickness must be positive")


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


class BlendedWinglet:
    """The parent wing continued into a winglet. frame(s) gives the section frame at s;
    negative s returns the parent itself."""

    def __init__(self, parent: SurfaceDefinition, definition: WingletDefinition,
                 parent_profiles: ProfileFamily, target_profiles: ProfileFamily):
        if parent.span_curve != "full_span" or parent.section_orientation != "parallel":
            raise ValueError("The winglet law continues a full_span wing with parallel sections")
        self.parent = parent
        self.d = definition
        self.parent_profiles = parent_profiles
        self.target_profiles = target_profiles
        for name, family in (("parent", parent_profiles), ("target", target_profiles)):
            for chord, profile in family.source.items():
                expected = definition.trailing_edge_mm / 2 / chord
                if abs(profile.half_thickness[-1] - expected) > 1e-9 * max(expected, 1e-3):
                    raise ValueError(f"The {name} profile at {chord:g} mm ends "
                                     f"{2 * profile.half_thickness[-1] * chord:.9f} mm thick, not "
                                     f"the {definition.trailing_edge_mm:g} mm trailing edge")
        self._join()
        if definition.bend_length_mm is None:
            self.d = replace(definition, bend_length_mm=self._solve_bend_length())
        self._configure(self.d.bend_length_mm)

    # -------------------------------------------------------------- parent side

    def _parent_q(self, f):
        frame = section_frame(self.parent, f)
        return np.array([*frame.le, frame.chord_mm, _twist_deg(self.parent, f)])

    def _join(self):
        f0, h = self.d.attach_fraction, 1e-3
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
        self.aft = _Taylor(q_s[2], q_ss[2], 0.0)  # dZ/ds and its slope
        self.chord_p = _Taylor(q[3], q_s[3], q_ss[3])
        self.twist_p = _Taylor(q[4], q_s[4], q_ss[4])

    def _parent_speed(self, f):
        frame = section_frame(self.parent, f)
        span = self.parent.projected_span_mm - self.parent.tip_span_trim_mm
        slope = frame.leading_edge_tangent[1] / -frame.leading_edge_tangent[0]
        return span * math.hypot(1.0, slope)

    def _parent_arc(self, f_from, f_to):
        nodes = f_from + 0.5 * (f_to - f_from) * (_GAUSS_X + 1)
        return 0.5 * (f_to - f_from) * sum(w * self._parent_speed(t) for w, t in zip(_GAUSS_W, nodes))

    def parent_fraction(self, s):
        """Parent span fraction at front-view arc length s <= 0 (Newton on the exact arc)."""
        f = self.f0 + s / self.ds_df
        for _ in range(30):
            step = (self._parent_arc(self.f0, f) - s) / self._parent_speed(f)
            f -= step
            if abs(step) < 1e-15:
                break
        return f

    # -------------------------------------------------------------- winglet side

    def _configure(self, bend_length):
        d = self.d
        self.L = float(bend_length)
        self.gamma_c = math.radians(d.cant_deg)
        self.K = (self.gamma_c - self.gamma0 - self.kappa0 * self.L / 2) / self.L
        self.bend_end = self._front(self.L)
        rise = self.bend_end[1] - self.le0[1]
        leg_rise = d.height_mm - rise
        if leg_rise <= 0:
            raise ValueError(f"A {self.L:.1f} mm bend already rises {rise:.1f} mm, "
                             f"more than the {d.height_mm:g} mm winglet height")
        self.leg_length = leg_rise / math.sin(self.gamma_c)
        self.s_tip = self.L + self.leg_length

    def kappa(self, s):
        if s >= self.L:
            return 0.0
        x = s / self.L
        return self.kappa0 * (1 - smoothstep(x)) + self.K * 30 * x * x * (1 - x) ** 2

    def gamma(self, s):
        x = min(s / self.L, 1.0)
        return float(self.gamma0 + self.kappa0 * self.L * (x - _smoothstep_integral(x))
                     + self.K * self.L * smoothstep(x))

    def _front(self, s):
        """Front-view LE (X, Y) at s, by Gauss quadrature of the tangent."""
        s_bend = min(s, self.L)
        nodes = 0.5 * s_bend * (_GAUSS_X + 1)
        angles = np.array([self.gamma(t) for t in nodes])
        dx = -0.5 * s_bend * np.dot(_GAUSS_W, np.cos(angles))
        dy = 0.5 * s_bend * np.dot(_GAUSS_W, np.sin(angles))
        x, y = self.le0[0] + dx, self.le0[1] + dy
        if s > self.L:
            x -= math.cos(self.gamma_c) * (s - self.L)
            y += math.sin(self.gamma_c) * (s - self.L)
        return np.array([x, y])

    def _blend(self, s, taylor, target):
        w = float(smoothstep(s / self.L))
        return (1 - w) * taylor(s) + w * target(s)

    def aft_rate(self, s):
        return self._blend(s, self.aft, lambda _: -math.tan(math.radians(self.d.sweep_deg)))

    def _z(self, s):
        if s > self.L:
            return self._z(self.L) - math.tan(math.radians(self.d.sweep_deg)) * (s - self.L)
        nodes = 0.5 * s * (_GAUSS_X + 1)
        return self.le0[2] + 0.5 * s * np.dot(_GAUSS_W, [self.aft_rate(t) for t in nodes])

    def chord(self, s):
        """Blend of the parent's Taylor chord into a linear taper from the join chord to the
        tip chord. Both decrease and the taper falls faster, so the blend only decreases."""
        start, tip = self.chord_p.value, self.d.tip_chord_mm

        def taper(t):
            return start + (tip - start) * t / self.s_tip

        return self._blend(s, self.chord_p, taper)

    def twist_deg(self, s):
        return self._blend(s, self.twist_p, lambda _: -self.d.toe_out_deg)

    def morph(self, s):
        return float(smoothstep(s / self.L))

    def profile(self, s):
        """Camber fades and thickness blends; the TE ramp is re-added at the specified mm."""
        m = self.morph(s)
        chord = self.chord(s)
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
        phi = float(smoothstep(s / self.L)) * gamma
        theta_deg = self.twist_deg(s)
        theta = math.radians(theta_deg)
        up = np.array([math.sin(phi), math.cos(phi), 0.0])
        chord_axis = np.array([-up[0] * math.sin(theta), -up[1] * math.sin(theta), -math.cos(theta)])
        up_axis = np.array([up[0] * math.cos(theta), up[1] * math.cos(theta), -math.sin(theta)])
        normal = np.array([-math.cos(phi), math.sin(phi), 0.0])
        return Frame(float(s), le, chord_axis, up_axis, normal, self.chord(s), self.profile(s),
                     gamma, phi, theta_deg)

    def parent_frame(self, f, s=float("nan")) -> Frame:
        pf = section_frame(self.parent, f)
        return Frame(float(s), np.array(pf.le), np.array(pf.chord_axis), np.array(pf.up_axis),
                     np.array(pf.span_tangent), pf.chord_mm, self.parent_profiles.at(pf.chord_mm),
                     math.atan2(pf.leading_edge_tangent[1], -pf.leading_edge_tangent[0]), 0.0,
                     _twist_deg(self.parent, f))

    def stations(self, bend=81, leg=41, extend_mm=0.0):
        s_bend = self.L * (1 - np.cos(np.linspace(0, math.pi, bend))) / 2
        s_leg = np.linspace(self.L, self.s_tip + extend_mm, leg + 1)[1:]
        return np.concatenate([s_bend, s_leg])

    def extended_frame(self, s):
        """Past the tip station the leg continues unchanged (used for the cap allowance)."""
        if s <= self.s_tip:
            return self.frame(s)
        tip = self.frame(self.s_tip)
        direction = np.array([-math.cos(self.gamma_c), math.sin(self.gamma_c),
                              -math.tan(math.radians(self.d.sweep_deg))])
        return replace(tip, s=float(s), le=tip.le + direction * (s - self.s_tip))

    def projected_span(self, include_allowance=True):
        extra = self.d.cap_allowance_mm if include_allowance else 0.0
        outboard = max(float(-self.extended_frame(float(s)).section()[:, 0].min())
                       for s in self.stations(41, 21, extra))
        return outboard + self.parent.root_le_mm[0]

    def _solve_bend_length(self):
        target = self.d.projected_span_mm

        def excess(L):
            self._configure(L)
            return self.projected_span() - target

        # A longer, gentler bend reaches farther outboard before the leg turns up.
        lo, hi = 5.0, 400.0
        if excess(lo) > 0:
            raise ValueError("Even the tightest bend reaches past the projected span")
        reached = False
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            try:
                over = excess(mid) > 0
            except ValueError:  # the bend alone rises past the winglet height
                hi = mid
                continue
            reached |= over
            lo, hi = (lo, mid) if over else (mid, hi)
            if hi - lo < 1e-7:
                break
        if not reached:
            raise ValueError("No bend within the winglet height reaches the projected span")
        return 0.5 * (lo + hi)

    def to_dict(self):
        out = asdict(self.d)
        out.update(bend_length_mm=self.L, leg_length_mm=self.leg_length, tip_s_mm=self.s_tip,
                   chord_at_bend_end_mm=self.chord(self.L),
                   parent_dihedral_deg=math.degrees(self.gamma0),
                   parent_curvature_per_mm=self.kappa0,
                   peak_curvature_per_mm=max(self.kappa(s) for s in np.linspace(0, self.L, 2001)))
        return out
