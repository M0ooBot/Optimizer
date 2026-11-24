"""Straight 800 mm wing with a symmetric polynomial section (accepted flat-wing trial).

Every section is the same: x = 120 u^3 (10 - 15u + 6u^2), half thickness
h = u (1 - u) (34 - 20u) mm, for the native chord parameter u in [0, 1].
The leading-edge point is at X = 0, Z = 0 and the chord is 120 mm.
"""

SHAPE = "flat"
VARIANT = "sharp"
TE_THICKNESS = 0.0
UNITS = "mm"
AXES = {"X": "aft, along the chord", "Y": "span (+Y is the right tip)", "Z": "up"}
SPAN_AXIS = "Y"
SIDES = {+1: "upper (+Z)", -1: "lower (-Z)"}
CHORD = 120.0
ATTACHMENT = 350.0
SPAN_LIMITS = (-ATTACHMENT, ATTACHMENT)
FULL_SPAN = True
REACH = 50.0


def _x(u):
    """Normalised chordwise position of the native parameter u."""
    return u * u * u * (10 - 15 * u + 6 * u * u)


def _h(u):
    """Normalised half thickness (fraction of the 120 mm chord)."""
    return u * (1 - u) * (34 - 20 * u) / CHORD


def chord(span):
    """Chord in mm. Keep-chord blunting leaves it unchanged."""
    return np.full(np.shape(span), CHORD)


def _cut(span):
    """Jets of the native end parameter u_end and the keep-chord scale f."""
    if TE_THICKNESS == 0:
        return _const(1.0, span), _const(1.0, span)

    def residual(unknowns, outer):
        u = unknowns[0]
        return [2 * CHORD * _h(u) / _x(u) - TE_THICKNESS + 0 * outer]

    def real(u):
        j = Jet.constant(u, 0)
        return (2 * CHORD * _h(j) / _x(j)).v - TE_THICKNESS

    root = _bisect(real, np.full(span.v.shape, 0.6), np.full(span.v.shape, 1.0))
    root = _newton_polish(residual, [root], span.v)
    u_end, = _implicit(residual, root, span)
    return u_end, 1 / _x(u_end)


def cut(span):
    """Where the keep-chord cut lands: native parameter, original chord fraction, scale."""
    u_end, scale = _cut(Jet.constant(span, 2))
    return {"u_end": u_end.v, "cut_chord_fraction": _x(u_end).v, "scale": scale.v}


def _skin(s, span, side):
    u_end, scale = _cut(span)
    u = s * u_end
    return [CHORD * scale * _x(u), span + 0 * u, side * CHORD * scale * _h(u)]


def skin_jet(s, span, side):
    """Skin point and exact derivatives at chord parameter s (0 = LE, 1 = TE end) and Y = span.

    side +1 is the upper skin, -1 the lower skin. Returns a dict of arrays
    (last axis XYZ): point, d_chord, d_span, d_chord_chord, d_chord_span, d_span_span.
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
    """Section plane at Y = span: 2D origin (the LE point), chord and up directions, chord."""
    return {"origin": np.array([0.0, float(span), 0.0]), "chord_direction": np.array([1.0, 0, 0]),
            "up_direction": np.array([0, 0, 1.0]), "chord": CHORD}


def section_thickness(xi, span):
    """Thickness in mm across the finished section at chord fraction xi (vertical in the section)."""
    xi, span = np.broadcast_arrays(np.asarray(xi, float), np.asarray(span, float))
    u_end, scale = _cut(Jet.constant(span, 0))
    target = xi / scale.v
    real = lambda u: _x(Jet.constant(u, 0)).v - target
    u = _bisect(real, np.zeros_like(target), np.ones_like(target))
    u = np.where(target >= _x(u_end).v, u_end.v, u)
    thickness = 2 * CHORD * scale.v * _h(Jet.constant(u, 0)).v
    return np.where(target <= _x(u_end).v * (1 + 1e-12), thickness, np.nan)


def tip_outline_jet(p, side=+1):
    """Locked plan outline of the tip at the mid-plane Z = 0, p = 0 at the LE end, 1 at the TE end.

    X = 120 x(p), |Y| = 350 + 200 p (1 - p): the outline uses the section's own
    chordwise polynomial, so it meets the LE and TE of the attachment.
    """
    pj = Jet.variable(np.asarray(p, float), 0, 1)
    xyz = [CHORD * _x(pj), side * (ATTACHMENT + 200 * pj * (1 - pj)), 0 * pj]
    return _pack1(xyz)


def tip_outline(p, side=+1):
    return tip_outline_jet(p, side)["point"]


OUTLINE = {
    "kind": "plan curve in the mid-plane Z = 0, mirrored for the left tip",
    "formula": "X = 120 p^3 (10 - 15p + 6p^2), |Y| = 350 + 200 p (1 - p), Z = 0, p in [0, 1] from LE to TE",
    "reach_max_mm": 50.0,
    "reach_at": "p = 0.5 (X = 60 mm), tip at |Y| = 400 mm",
}
