"""Symmetric 55-degree swept, tapered wing with a NACA 0009 closed-TE section (accepted swept55 trial).

Native chord parameter q in [0, 1], chord fraction U = 3q^2 - 2q^3, half
thickness per chord t = 0.45 (0.2969 q sqrt(3 - 2q) - 0.1260 U - 0.3516 U^2
+ 0.2843 U^3 - 0.1036 U^4) (q sqrt(3 - 2q) = sqrt(U)). Leading edge at
X = tan(55 deg) |Y|, Z = 0; chord 210 - (115/390) |Y| mm. Level and untwisted.
"""

SHAPE = "swept55"
VARIANT = "sharp"
TE_THICKNESS = 0.0
UNITS = "mm"
AXES = {"X": "aft, along the chord", "Y": "span (+Y is the right tip)", "Z": "up"}
SPAN_AXIS = "Y"
SIDES = {+1: "upper (+Z)", -1: "lower (-Z)"}
ATTACHMENT = 390.0
SPAN_LIMITS = (-ATTACHMENT, ATTACHMENT)
FULL_SPAN = True
ROOT_CHORD = 210.0
TAPER = 115.0 / ATTACHMENT
SWEEP = float(np.tan(np.radians(55.0)))
REACH = 60.0


def _U(q):
    return q * q * (3 - 2 * q)


def _t(q):
    u = _U(q)
    return 0.45 * (0.2969 * q * jsqrt(3 - 2 * q) - 0.1260 * u - 0.3516 * u * u
                   + 0.2843 * u * u * u - 0.1036 * u * u * u * u)


def chord(span):
    return ROOT_CHORD - TAPER * np.abs(np.asarray(span, float))


def _cut(span):
    if TE_THICKNESS == 0:
        return _const(1.0, span), _const(1.0, span)

    def residual(unknowns, outer):
        q = unknowns[0]
        c = ROOT_CHORD - TAPER * jabs(outer)
        return [2 * c * _t(q) / _U(q) - TE_THICKNESS]

    c = ROOT_CHORD - TAPER * np.abs(span.v)

    def real(q):
        j = Jet.constant(q, 0)
        return 2 * c * (_t(j) / _U(j)).v - TE_THICKNESS

    root = _bisect(real, np.full(span.v.shape, 0.6), np.full(span.v.shape, 1.0))
    root = _newton_polish(residual, [root], span.v)
    q_end, = _implicit(residual, root, span)
    return q_end, 1 / _U(q_end)


def cut(span):
    q_end, scale = _cut(Jet.constant(span, 2))
    return {"q_end": q_end.v, "cut_chord_fraction": _U(q_end).v, "scale": scale.v}


def _skin(s, span, side):
    q_end, scale = _cut(span)
    q = s * q_end
    ya = jabs(span)
    c = ROOT_CHORD - TAPER * ya
    return [SWEEP * ya + c * scale * _U(q), span + 0 * q, side * c * scale * _t(q)]


def skin_jet(s, span, side):
    """Skin point and exact derivatives at chord parameter s (0 = LE, 1 = TE end) and Y = span.

    side +1 is the upper skin, -1 the lower. |Y| enters the planform, so at
    Y = 0 the span derivatives are those of the right (+Y) half. Returns a dict
    of arrays (last axis XYZ): point, d_chord, d_span, d_chord_chord,
    d_chord_span, d_span_span.
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
    span = float(span)
    return {"origin": np.array([SWEEP * abs(span), span, 0.0]), "chord_direction": np.array([1.0, 0, 0]),
            "up_direction": np.array([0, 0, 1.0]), "chord": float(chord(span))}


def section_thickness(xi, span):
    """Thickness in mm across the finished section at chord fraction xi (vertical in the section)."""
    xi, span = np.broadcast_arrays(np.asarray(xi, float), np.asarray(span, float))
    q_end, scale = _cut(Jet.constant(span, 0))
    target = xi / scale.v
    q = _bisect(lambda q: _U(Jet.constant(q, 0)).v - target, np.zeros_like(target), np.ones_like(target))
    q = np.where(target >= _U(q_end).v, q_end.v, q)
    thickness = 2 * chord(span) * scale.v * _t(Jet.constant(q, 0)).v
    return np.where(target <= _U(q_end).v * (1 + 1e-12), thickness, np.nan)


def tip_outline_jet(p, side=+1):
    """Locked plan outline (lip) at Z = 0, p = 0 at the LE end, 1 at the TE end.

    The lip lies on the parent planform continued outboard: at station
    r = 390 + 240 p (1 - p) it sits at chord fraction U(p) = 3p^2 - 2p^3,
    X = tan(55 deg) r + (210 - (115/390) r) U(p).
    """
    pj = Jet.variable(np.asarray(p, float), 0, 1)
    r = ATTACHMENT + 240 * pj * (1 - pj)
    xyz = [SWEEP * r + (ROOT_CHORD - TAPER * r) * _U(pj), side * r, 0 * pj]
    return _pack1(xyz)


def tip_outline(p, side=+1):
    return tip_outline_jet(p, side)["point"]


OUTLINE = {
    "kind": "plan curve in the chord plane Z = 0 (the lip), mirrored for the left tip",
    "formula": "r = 390 + 240 p (1 - p); X = tan(55 deg) r + (210 - (115/390) r)(3p^2 - 2p^3); |Y| = r; Z = 0; p in [0, 1] from LE to TE",
    "reach_max_mm": 60.0,
    "reach_at": "p = 0.5, tip at |Y| = 450 mm (900 mm projected span)",
}
