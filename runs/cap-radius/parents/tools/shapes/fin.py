"""Swept, tapered vertical fin with exact 12 % elliptical sections (accepted fin trial).

At height Z the section is an ellipse with its leading edge at X = 65 Z/240,
chord c = 160 - 75 Z/240 and half thickness 0.06 c (12 % thickness-to-chord).
Native chord parameter: angle theta = pi s in [0, pi], chord fraction
(1 - cos theta)/2, half thickness 0.06 c sin theta. Flat root at Z = 0.
"""

SHAPE = "fin"
VARIANT = "sharp"
TE_THICKNESS = 0.0
REACH = 40.0
UNITS = "mm"
AXES = {"X": "aft, along the chord", "Y": "transverse (thickness)", "Z": "up (span)"}
SPAN_AXIS = "Z"
SIDES = {+1: "+Y skin", -1: "-Y skin"}
ATTACHMENT = 240.0
SPAN_LIMITS = (0.0, ATTACHMENT)
FULL_SPAN = False
ROOT_CHORD = 160.0
THICKNESS_RATIO = 0.12


def chord(z):
    return ROOT_CHORD - 75.0 * np.asarray(z, float) / ATTACHMENT


def _xi(theta):
    return 0.5 * (1 - jcos(theta))


def _cut(z):
    """Jets of the end angle theta_end and the keep-chord scale."""
    if TE_THICKNESS == 0:
        return _const(np.pi, z), _const(1.0, z)

    def residual(unknowns, outer):
        theta = unknowns[0]
        c = ROOT_CHORD - (75.0 / ATTACHMENT) * outer
        return [THICKNESS_RATIO * c * jsin(theta) / _xi(theta) - TE_THICKNESS]

    c = ROOT_CHORD - 75.0 * z.v / ATTACHMENT

    def real(theta):
        j = Jet.constant(theta, 0)
        return THICKNESS_RATIO * c * (jsin(j) / _xi(j)).v - TE_THICKNESS

    root = _bisect(real, np.full(z.v.shape, 0.6 * np.pi), np.full(z.v.shape, np.pi))
    root = _newton_polish(residual, [root], z.v)
    theta_end, = _implicit(residual, root, z)
    return theta_end, 1 / _xi(theta_end)


def cut(z):
    theta_end, scale = _cut(Jet.constant(z, 2))
    return {"theta_end": theta_end.v, "cut_chord_fraction": _xi(theta_end).v, "scale": scale.v}


def _skin(s, z, side):
    theta_end, scale = _cut(z)
    theta = s * theta_end
    c = ROOT_CHORD - (75.0 / ATTACHMENT) * z
    return [(65.0 / ATTACHMENT) * z + c * scale * _xi(theta),
            side * 0.5 * THICKNESS_RATIO * c * scale * jsin(theta), z + 0 * theta]


def skin_jet(s, z, side):
    """Skin point and exact derivatives at chord parameter s (0 = LE, 1 = TE end) and height z.

    side +1 is the +Y skin, -1 the -Y skin. theta = s * theta_end(z)
    (theta_end = pi on the sharp parent). Returns a dict of arrays (last axis
    XYZ): point, d_chord, d_span, d_chord_chord, d_chord_span, d_span_span,
    where span means the height Z.
    """
    sj, zj = _seed(s, z)
    return _pack(_skin(sj, zj, side))


def skin(s, z, side):
    return skin_jet(s, z, side)["point"]


def te_face_jet(t, z):
    """Blunt trailing-edge face: t = 0 on the -Y corner, t = 1 on the +Y corner."""
    if TE_THICKNESS == 0:
        raise ValueError("The sharp variant has no trailing-edge face")
    tj, zj = _seed(t, z)
    one = _const(1.0, tj)
    plus, minus = _skin(one, zj, +1), _skin(one, zj, -1)
    return _pack([lo + tj * (up - lo) for up, lo in zip(plus, minus)], "t", "span")


def te_face(t, z):
    return te_face_jet(t, z)["point"]


def section_frame(z):
    z = float(z)
    return {"origin": np.array([65.0 * z / ATTACHMENT, 0.0, z]), "chord_direction": np.array([1.0, 0, 0]),
            "up_direction": np.array([0, 1.0, 0]), "chord": float(chord(z))}


def section_thickness(xi, z):
    """Thickness in mm (along Y) across the finished section at chord fraction xi."""
    xi, z = np.broadcast_arrays(np.asarray(xi, float), np.asarray(z, float))
    theta_end, scale = _cut(Jet.constant(z, 0))
    target = xi / scale.v
    theta = np.where(target >= _xi(theta_end).v, theta_end.v, np.arccos(np.clip(1 - 2 * target, -1, 1)))
    thickness = THICKNESS_RATIO * chord(z) * scale.v * np.sin(theta)
    return np.where(target <= _xi(theta_end).v * (1 + 1e-12), thickness, np.nan)


def tip_outline(h, branch):
    """Locked side-view outline at Y = 0: X = m + branch * a sqrt(1 - h^3), Z = 240 + REACH h.

    h in [0, 1] is the normalised height above the attachment; branch -1 is the
    leading-edge side, +1 the trailing-edge side; they meet at the crown h = 1.
    m = 80 + (27.5/240) Z and a = 80 - (37.5/240) Z continue the parent's
    linear sweep and taper.
    """
    h = np.asarray(h, float)
    z = ATTACHMENT + REACH * h
    m = 80 + 27.5 * z / ATTACHMENT
    a = 80 - 37.5 * z / ATTACHMENT
    return np.stack([m + branch * a * np.sqrt(np.clip(1 - h**3, 0, None)), np.zeros_like(h), z], -1)


def tip_outline_jet(h, branch):
    """Outline point and derivatives with respect to h (singular at the crown h = 1)."""
    hj = Jet.variable(np.asarray(h, float), 0, 1)
    z = ATTACHMENT + REACH * hj
    m = 80 + (27.5 / ATTACHMENT) * z
    a = 80 - (37.5 / ATTACHMENT) * z
    return _pack1([m + branch * a * jsqrt(1 - hj**3), 0 * hj, z])


OUTLINE = {
    "kind": "side-view outline in the plane Y = 0, two branches meeting at the crown",
    "formula": "h = (Z - 240)/REACH in [0, 1]; X = m ± a sqrt(1 - h^3); m = 80 + (27.5/240) Z; a = 80 - (37.5/240) Z; minus = LE branch, plus = TE branch",
    "reach_max_mm": None,
    "reach_at": "crown at h = 1",
}
