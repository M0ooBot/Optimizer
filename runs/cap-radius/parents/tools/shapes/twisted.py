"""Swept, tapered, dihedral wing with washout and a cambered NACA 2412 closed-TE section (accepted twisted trial).

Exact transcription of the locked analytic parent definition supplied with the accepted model.
Section in native parameter u in [0, 1]: x = sin^2(pi u / 2); NACA 2412 camber
with the two parabolas blended by a quintic smoothstep over 0.36-0.44 chord;
thickness 0.6 (0.2969 sqrt(x) - 0.1260 x - 0.3516 x^2 + 0.2843 x^3 - 0.1036 x^4)
applied along the camber-line normal. At station |Y|: chord 180 - 80 |Y|/440,
leading edge at (90 |Y|/440, Y, |Y| tan 5 deg), true twist (washout)
3 deg |Y|/440 about the leading edge. The left half mirrors the right in XZ.
"""

SHAPE = "twisted"
VARIANT = "sharp"
TE_THICKNESS = 0.0
UNITS = "mm"
AXES = {"X": "aft, along the root chord", "Y": "span (+Y is the right tip)", "Z": "up"}
SPAN_AXIS = "Y"
SIDES = {+1: "upper", -1: "lower"}
CUT = 440.0
TIP = 500.0
ATTACHMENT = CUT
SPAN_LIMITS = (-CUT, CUT)
FULL_SPAN = True
DIHEDRAL = float(np.tan(np.deg2rad(5)))
WASHOUT = float(np.deg2rad(3))
REACH = TIP - CUT


def _airfoil(u, side):
    """Normalised section point (X, Z) on the given side (+1 upper, -1 lower, 0 camber line)."""
    half = jsin(u * (np.pi / 2))
    x = half * half
    forward = (0.02 / 0.4**2) * (0.8 * x - x * x)
    aft = (0.02 / 0.6**2) * (0.2 + 0.8 * x - x * x)
    raw = (x - 0.36) / 0.08
    fraction = jwhere(raw.v <= 0, _const(0.0, raw), jwhere(raw.v >= 1, _const(1.0, raw), raw))
    blend = fraction**3 * (10 - 15 * fraction + 6 * fraction * fraction)
    blend_derivative = 30 * fraction**2 * (1 - fraction)**2 / 0.08
    camber = forward + (aft - forward) * blend
    slope = (0.8 - 2 * x) * ((0.02 / 0.4**2) * (1 - blend) + (0.02 / 0.6**2) * blend) \
        + (aft - forward) * blend_derivative
    thickness = 0.6 * (0.2969 * half - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x
                       - 0.1036 * x * x * x * x)
    norm = jsqrt(1 + slope * slope)
    return x - side * thickness * slope / norm, camber + side * thickness / norm


def _real(u, side):
    j = Jet.constant(u, 0)
    x, z = _airfoil(j, side)
    return x.v, z.v


def _nose():
    # The upper surface reaches slightly ahead of the u = 0 point; the Converter
    # measures chord from the most forward point of the loop.
    slope = lambda u: _airfoil(Jet.variable(u, 0, 1), +1)[0].d[..., 0]
    u = _bisect(slope, np.array(1e-9), np.array(0.2))
    return float(u), float(_real(u, +1)[0])


NOSE_U, NOSE = _nose()
CHORD_N = 1.0 - NOSE  # normalised chord the Converter keeps (max X is 1 at the TE)


def chord(span):
    return 180.0 - 80.0 * np.abs(np.asarray(span, float)) / CUT


def _side_u(x, side, lo=None):
    lo = np.full(np.shape(x), 0.5 if lo is None else lo)
    return _bisect(lambda u: _real(u, side)[0] - x, lo, np.ones(np.shape(x)), iterations=80)


def _cut(span):
    """Jets of the upper and lower native end parameters and the keep-chord scale."""
    if TE_THICKNESS == 0:
        one = _const(1.0, span)
        return one, one, one
    # The cut depends on |Y| only: solve once per distinct station.
    stations, inverse = np.unique(np.abs(span.v), return_inverse=True)
    target = TE_THICKNESS / (180.0 - 80.0 * stations / CUT)

    def residual(unknowns, outer):
        c = 180.0 - 80.0 * jabs(outer) / CUT
        x1, z1 = _airfoil(unknowns[0], +1)
        x2, z2 = _airfoil(unknowns[1], -1)
        return [x1 - x2, c * (z1 - z2) * CHORD_N / (x1 - NOSE) - TE_THICKNESS]

    def same_u_gap(u):
        x1, z1 = _real(u, +1)
        return (z1 - _real(u, -1)[1]) * CHORD_N / (x1 - NOSE) - target

    guess = _bisect(same_u_gap, np.full(stations.shape, 0.6), np.ones(stations.shape), iterations=60)
    roots = _newton_polish(residual, [guess, guess], stations, steps=8)
    check = residual([Jet.constant(r, 0) for r in roots], Jet.constant(stations, 0))
    if max(np.abs(r.v).max() for r in check) > 1e-12:
        def gap(xc):
            z_up = _real(_side_u(xc, +1), +1)[1]
            z_lo = _real(_side_u(xc, -1), -1)[1]
            return (z_up - z_lo) * CHORD_N / (xc - NOSE) - target

        xc = _bisect(gap, np.full(stations.shape, 0.6), np.ones(stations.shape), iterations=80)
        roots = _newton_polish(residual, [_side_u(xc, +1), _side_u(xc, -1)], stations)
    roots = [np.reshape(r[inverse], span.v.shape) for r in roots]
    u_up, u_lo = _implicit(residual, roots, span)
    return u_up, u_lo, CHORD_N / (_airfoil(u_up, +1)[0] - NOSE)


def cut(span):
    u_up, u_lo, scale = _cut(Jet.constant(span, 2))
    return {"u_end_upper": u_up.v, "u_end_lower": u_lo.v,
            "cut_chord_fraction": _airfoil(u_up, +1)[0].v, "scale": scale.v}


def _place(span, xn, zn):
    """Place a normalised section point at a span station (taper, twist about the LE, dihedral)."""
    ya = jabs(span)
    theta = ya * (WASHOUT / CUT)
    c = 180.0 - 80.0 * ya / CUT
    cos, sin = jcos(theta), jsin(theta)
    return [ya * (90.0 / CUT) + c * (xn * cos - zn * sin), span + 0 * xn,
            ya * DIHEDRAL + c * (xn * sin + zn * cos)]


def _native(s, u_end):
    """Native parameter u of chord parameter s: u = s + (u_end - 1) max(0, 2s - 1)^5.

    The forward half keeps u = s, so the camber blend lines (CHORD_BREAKS)
    stay at fixed s; only the aft half stretches to the cut.
    """
    w = jwhere(s.v > 0.5, (2 * s - 1)**5, _const(0.0, s))
    return s + (u_end - 1) * w


# The camber blend makes the skin only C1 across these chord-parameter lines.
CHORD_BREAKS = tuple(float(2 / np.pi * np.arcsin(np.sqrt(x))) for x in (0.36, 0.44))


def _skin(s, span, side):
    u_up, u_lo, scale = _cut(span)
    u = _native(s, u_up if side > 0 else u_lo)
    xn, zn = _airfoil(u, side)
    return _place(span, scale * xn, scale * zn)


def skin_jet(s, span, side):
    """Skin point and exact derivatives at chord parameter s (0 = LE, 1 = TE end) and Y = span.

    side +1 is the upper skin, -1 the lower. The native parameter is
    u = s + (u_end(side, Y) - 1) max(0, 2s - 1)^5 (u = s on the sharp parent). |Y| enters the
    definition, so at Y = 0 the span derivatives are those of the right half.
    Returns a dict of arrays (last axis XYZ): point, d_chord, d_span,
    d_chord_chord, d_chord_span, d_span_span.
    """
    sj, yj = _seed(s, span)
    return _pack(_skin(sj, yj, side))


def skin(s, span, side):
    return skin_jet(s, span, side)["point"]


def te_face_jet(t, span):
    """Blunt trailing-edge face: t = 0 on the lower TE corner, t = 1 on the upper corner."""
    if TE_THICKNESS == 0:
        raise ValueError("The sharp variant has no trailing-edge face")
    tj, yj = _seed(t, span)
    one = _const(1.0, tj)
    upper, lower = _skin(one, yj, +1), _skin(one, yj, -1)
    return _pack([lo + tj * (up - lo) for up, lo in zip(upper, lower)], "t", "span")


def te_face(t, span):
    return te_face_jet(t, span)["point"]


def section_frame(span):
    """Section plane at Y = span: 2D origin (LE point), chord and up directions (twisted), chord."""
    span = float(span)
    ya = abs(span)
    theta = WASHOUT * ya / CUT
    return {"origin": np.array([90.0 * ya / CUT, span, ya * DIHEDRAL]),
            "chord_direction": np.array([np.cos(theta), 0.0, np.sin(theta)]),
            "up_direction": np.array([-np.sin(theta), 0.0, np.cos(theta)]),
            "chord": float(chord(span))}


def section_thickness(xi, span):
    """Thickness in mm across the finished section at chord fraction xi.

    Measured as the Converter does: perpendicular to the chord line (vertical
    in the section frame), between the outermost upper and lower crossings.
    """
    xi, span = np.broadcast_arrays(np.asarray(xi, float), np.asarray(span, float))
    u_up, u_lo, scale = _cut(Jet.constant(span, 0))
    x = xi / scale.v
    x_end = _airfoil(u_up, +1)[0].v
    at_end = x >= x_end
    u1 = np.where(at_end, u_up.v, _side_u(np.minimum(x, x_end), +1, lo=NOSE_U))
    u2 = np.where(at_end, u_lo.v, _bisect(lambda u: _real(u, -1)[0] - np.minimum(x, x_end),
                                          np.zeros_like(x), np.ones_like(x), 80))
    thickness = chord(span) * scale.v * (_real(u1, +1)[1] - _real(u2, -1)[1])
    return np.where(x <= x_end * (1 + 1e-12), thickness, np.nan)


END_FAIRING = 0.18  # fraction of p over which an end correction fades out


def _end_shift():
    """Offsets that move the rim's ends onto this parent's attachment LE point and TE.

    Zero on the sharp parent. On the blunt parent the Converter keeps the chord
    measured from the most forward point of the section, which lies slightly
    ahead of the LE point, so the TE moves forward by about 0.0001 mm; the rim's
    TE end is moved by that amount onto the blunt face (nearest point).
    """
    pj = Jet.constant(np.array([0.0, 1.0]), 0)
    xn, zn = _airfoil(pj, 0)
    rim = np.stack([c.v for c in _place(Jet.constant(np.array([CUT, CUT]), 0), xn, zn)], -1)
    le = skin(0.0, CUT, +1)
    a, b = skin(1.0, CUT, -1), skin(1.0, CUT, +1)
    t = np.clip((rim[1] - a) @ (b - a) / max((b - a) @ (b - a), 1e-300), 0, 1)
    return le - rim[0], a + t * (b - a) - rim[1]


def _fade(d):
    x = jwhere(d.v < END_FAIRING, d * (1 / END_FAIRING), _const(1.0, d))
    return 1 - x**3 * (10 - 15 * x + 6 * x * x)


def tip_outline_jet(p, side=+1):
    """Locked 3D rim of the tip, p = 0 at the LE end, 1 at the TE end.

    As in the accepted brief: for each chord guide p the rim lies on the parent's
    camber surface continued outboard to |Y| = 440 + 60 sin(pi p), at native
    chord parameter u = p (taper, twist, dihedral and sweep continued). On the
    blunt parent, sub-micron end corrections faded over 18 % of p make the rim
    start exactly at the attachment LE point and end on the blunt TE face.
    """
    pj = Jet.variable(np.asarray(p, float), 0, 1)
    station = CUT + REACH * jsin(pj * np.pi)
    xn, zn = _airfoil(pj, 0)
    xyz = _place(station, xn, zn)
    if TE_THICKNESS > 0:
        le_shift, te_shift = _end_shift()
        w_le, w_te = _fade(pj), _fade(1 - pj)
        xyz = [c + w_le * float(a) + w_te * float(b) for c, a, b in zip(xyz, le_shift, te_shift)]
    x, y, z = xyz
    return _pack1([x, side * y, z])


def tip_outline(p, side=+1):
    return tip_outline_jet(p, side)["point"]


OUTLINE = {
    "kind": "3D rim on the continued camber surface (its plan view is the locked outline), mirrored for the left tip",
    "formula": "|Y| = 440 + 60 sin(pi p); point = camber-line point of native parameter u = p placed with the parent's taper, twist, dihedral and sweep at that |Y|; p in [0, 1] from LE to TE. Blunt variant only: plus end corrections (le_shift, te_shift) times 1 - x^3(10 - 15x + 6x^2), x = min(p/0.18, 1) and min((1-p)/0.18, 1), which move the ends onto the blunt attachment's LE point and TE face (each under 0.0002 mm)",
    "reach_max_mm": 60.0,
    "reach_at": "p = 0.5, tip at |Y| = 500 mm (1000 mm projected span)",
}
