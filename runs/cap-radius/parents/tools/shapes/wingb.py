"""Wing B parent: the bowed, untwisted SD7037 half-wing of the paired silhouette trial.

The skin is the Open CASCADE loft of the Phase 0 sections (exactly as the
accepted parent was built) and is only available as a B-spline surface, so
this module exposes it from the faces of parent.step (needs OCP, as shipped
with CadQuery) rather than from a closed-form formula.

Axes are the STEP's own: X = -span (outboard is -X), Y = up, Z = forward
(the nose points to +Z; aft is -Z). The root plane is X = 0 and the
attachment plane is X = -540. Sections are parallel (constant X), untwisted.
Section frame of the definition at projected span S, f = S/600: leading-edge
reference (X, Y, Z) = (-S, 60 f^2, -(35 f + 45 f^2)), chord 200 - 125 f, chord
direction -Z, up +Y. The 2D origin of every Converter section is that
leading-edge reference.

Skin parameters: u runs once around the section (native B-spline parameter;
see nose_u and side_range: lower skin from u = 0 at the trailing edge to the
nose, upper skin from the nose to u = 1), v runs from the root (v = 0) to the
attachment (v = SURFACE_DATA["v_attachment"]). skin_jet(u, span) re-expresses
v as the projected span. On the blunt variant the trailing-edge face is
te_face(t, span), t = 0 at the lower corner and 1 at the upper corner.
"""

from functools import lru_cache
from pathlib import Path

SHAPE = "wingb"
VARIANT = "sharp"
TE_THICKNESS = 0.0
UNITS = "mm"
AXES = {"X": "-span (outboard is -X; root plane X = 0)", "Y": "up", "Z": "forward (nose +Z)"}
SPAN_AXIS = "-X"
ATTACHMENT = 540.0
SPAN_LIMITS = (0.0, ATTACHMENT)
FULL_SPAN = False
PROJECTED_SPAN = 600.0
REACH = PROJECTED_SPAN - ATTACHMENT
STEP_PATH = Path(__file__).with_name("parent.step")
# Filled in by tools/build_wingb.py from the built surface.
SURFACE_DATA = {}
OUTLINE_DATA = {}


def section_frame(span):
    """Definition frame at projected span: origin (LE reference), chord and up directions, chord."""
    f = float(span) / PROJECTED_SPAN
    return {"origin": np.array([-float(span), 60.0 * f * f, -(35.0 * f + 45.0 * f * f)]),
            "chord_direction": np.array([0.0, 0.0, -1.0]), "up_direction": np.array([0.0, 1.0, 0.0]),
            "chord": 200.0 - 125.0 * f}


def chord(span):
    return 200.0 - 125.0 * np.asarray(span, float) / PROJECTED_SPAN


@lru_cache(maxsize=1)
def surfaces():
    """(skin, trailing-edge face or None) as OCP Geom surfaces read from parent.step."""
    from OCP.BRep import BRep_Tool
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from OCP.GeomAbs import GeomAbs_BSplineSurface
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_Reader
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    reader = STEPControl_Reader()
    if reader.ReadFile(str(STEP_PATH)) != IFSelect_RetDone:
        raise IOError(f"Cannot read {STEP_PATH}")
    reader.TransferRoots()
    shape = reader.OneShape()
    found = []
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = TopoDS.Face_s(explorer.Current())
        if BRepAdaptor_Surface(face).GetType() == GeomAbs_BSplineSurface:
            props = GProp_GProps()
            BRepGProp.SurfaceProperties_s(face, props)
            found.append((props.Mass(), BRep_Tool.Surface_s(face)))
        explorer.Next()
    found.sort(key=lambda item: -item[0])
    skin = found[0][1]
    te = found[1][1] if len(found) > 1 else None
    return skin, te


def _d2(surface, u, v):
    from OCP.gp import gp_Pnt, gp_Vec
    u, v = np.broadcast_arrays(np.asarray(u, float), np.asarray(v, float))
    out = np.zeros(u.shape + (6, 3))
    for index in np.ndindex(u.shape):
        p, du, dv, duu, dvv, duv = gp_Pnt(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec()
        surface.D2(float(u[index]), float(v[index]), p, du, dv, duu, dvv, duv)
        out[index] = [p.Coord(), du.Coord(), dv.Coord(), duu.Coord(), duv.Coord(), dvv.Coord()]
    return out


def skin_jet_uv(u, v):
    """Skin point and derivatives in the native (u, v) parameters."""
    d = _d2(surfaces()[0], u, v)
    keys = ("point", "d_u", "d_v", "d_u_u", "d_u_v", "d_v_v")
    return {key: d[..., i, :] for i, key in enumerate(keys)}


def v_at_span(span):
    """Native v of the section plane X = -span (X depends on v only on this loft)."""
    span = np.asarray(span, float)
    skin = surfaces()[0]
    lo = np.zeros(span.shape)
    hi = np.full(span.shape, SURFACE_DATA["v_attachment"])
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        x = np.vectorize(lambda v: skin.Value(SURFACE_DATA["u_reference"], float(v)).X())(mid)
        outboard = -x < span
        lo, hi = np.where(outboard, mid, lo), np.where(outboard, hi, mid)
    return 0.5 * (lo + hi)


def skin_jet(u, span):
    """Skin point and exact derivatives at native u and projected span (chain rule through v).

    Keys: point, d_chord (= d/du), d_span, d_chord_chord, d_chord_span,
    d_span_span; span derivatives are per mm of projected span (towards -X).
    """
    u, span = np.broadcast_arrays(np.asarray(u, float), np.asarray(span, float))
    v = v_at_span(span)
    j = skin_jet_uv(u, v)
    xv, xvv = j["d_v"][..., 0], j["d_v_v"][..., 0]
    v1 = -1.0 / xv
    v2 = -xvv * v1**2 / xv
    return {"point": j["point"], "d_chord": j["d_u"], "d_span": j["d_v"] * v1[..., None],
            "d_chord_chord": j["d_u_u"], "d_chord_span": j["d_u_v"] * v1[..., None],
            "d_span_span": j["d_v_v"] * (v1**2)[..., None] + j["d_v"] * v2[..., None]}


def skin(u, span):
    return skin_jet(u, span)["point"]


def nose_u(span):
    """Native u of the most forward skin point (largest Z) in the section at this span."""
    span = np.asarray(span, float)
    v = v_at_span(span)
    lo = np.full(span.shape, SURFACE_DATA["nose_bracket"][0])
    hi = np.full(span.shape, SURFACE_DATA["nose_bracket"][1])
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        dz = skin_jet_uv(mid, v)["d_u"][..., 2]
        rising = dz * SURFACE_DATA["nose_sign"] > 0
        lo, hi = np.where(rising, mid, lo), np.where(rising, hi, mid)
    return 0.5 * (lo + hi)


def side_range(span, side):
    """Native u interval of one skin, ordered from the nose to the trailing-edge end.

    side +1 is the upper skin (larger Y), -1 the lower skin.
    """
    n = nose_u(span)
    end = SURFACE_DATA["u_end_upper"] if side > 0 else SURFACE_DATA["u_end_lower"]
    return n, np.full(np.shape(n), end)


def te_face_jet(t, span):
    """Blunt trailing-edge face, t = 0 on the lower corner and 1 on the upper, at a projected span."""
    te = surfaces()[1]
    if te is None:
        raise ValueError("The sharp variant has no trailing-edge face")
    t, span = np.broadcast_arrays(np.asarray(t, float), np.asarray(span, float))
    u0, u1 = SURFACE_DATA["te_u_range"]
    if not SURFACE_DATA["te_u_upward"]:
        u0, u1 = u1, u0
    u = u0 + t * (u1 - u0)
    sign = u1 - u0  # d u / d t
    v = v_at_span(span)
    d = _d2(te, u, v)
    xv, xvv = d[..., 2, 0], d[..., 5, 0]
    v1 = -1.0 / xv
    v2 = -xvv * v1**2 / xv
    return {"point": d[..., 0, :], "d_t": sign * d[..., 1, :], "d_span": d[..., 2, :] * v1[..., None],
            "d_t_t": sign**2 * d[..., 3, :], "d_t_span": sign * d[..., 4, :] * v1[..., None],
            "d_span_span": d[..., 5, :] * (v1**2)[..., None] + d[..., 2, :] * v2[..., None]}


def te_face(t, span):
    return te_face_jet(t, span)["point"]


def section_thickness(xi, span):
    """Thickness in mm (along Y) of the section at chord fraction xi of the definition frame.

    The station is Z = Z_LE - xi * chord, measured from the frame's leading-edge
    reference (the Converter's 2D origin). NaN outside the section.
    """
    xi = np.atleast_1d(np.asarray(xi, float))
    span = float(span)
    frame = section_frame(span)
    v = float(v_at_span(span))
    skin = surfaces()[0]
    nose = float(nose_u(span))
    z_target = frame["origin"][2] - xi * frame["chord"]
    z_nose = skin.Value(nose, v).Z()
    for side in (+1, -1):
        end = SURFACE_DATA["u_end_upper"] if side > 0 else SURFACE_DATA["u_end_lower"]
        lo, hi = np.full(xi.shape, nose), np.full(xi.shape, end)
        z_end = skin.Value(end, v).Z()
        for _ in range(70):
            mid = 0.5 * (lo + hi)
            z = np.array([skin.Value(float(m), v).Z() for m in mid])
            ahead = z > z_target
            lo, hi = np.where(ahead, mid, lo), np.where(ahead, hi, mid)
        y = np.array([skin.Value(float(m), v).Y() for m in 0.5 * (lo + hi)])
        if side > 0:
            upper = np.where(z_target >= z_end - 1e-9, y, np.nan)
        else:
            lower = np.where(z_target >= z_end - 1e-9, y, np.nan)
    return np.where(z_target <= z_nose + 1e-12, upper - lower, np.nan)


def _transition_weight(distance, extent):
    x = np.minimum(distance / extent, 1.0)
    return 1 - x**3 * (10 - 15 * x + 6 * x**2)


def tip_outline(p):
    """Locked (faired) plan outline of the tip, p = 0 at the leading-edge end, 1 at the trailing-edge end.

    Returns an (..., 2) array of STEP (X, Z): X = -span, Z = forward. A cubic
    Bezier in (span, aft) with 80 mm control arms along the parent's plan
    tangents at the attachment, whose aft coordinate is faired over the end
    18 % of p into the quadratic continuation of the parent's own leading and
    trailing edges.
    """
    d = OUTLINE_DATA
    p = np.asarray(p, float)
    start, stop = np.array(d["nose_span_aft"]), np.array(d["tail_span_aft"])
    arm = d["control_arm_mm"]
    controls = [start, start + [arm, arm * d["nose_plan_slope"]],
                stop + [arm, arm * d["tail_plan_slope"]], stop]
    weights = [(1 - p)**3, 3 * p * (1 - p)**2, 3 * p**2 * (1 - p), p**3]
    span = sum(w * c[0] for w, c in zip(weights, controls))
    aft = sum(w * c[1] for w, c in zip(weights, controls))
    for side, distance in enumerate((p, 1 - p)):
        root = (start, stop)[side][1]
        slope, second = d["edge_aft_taylor"][side]
        delta = span - ATTACHMENT
        continued = root + delta * (slope + delta * second)
        fair = distance < d["edge_fairing_fraction"]
        aft = np.where(fair, aft + (continued - aft) * _transition_weight(distance, d["edge_fairing_fraction"]), aft)
    plan = np.stack([-span, -aft], -1)
    if TE_THICKNESS > 0:
        le_shift, te_shift = _end_shift()
        fade = d["edge_fairing_fraction"]
        plan = plan + _transition_weight(p, fade)[..., None] * le_shift \
            + _transition_weight(1 - p, fade)[..., None] * te_shift
    return plan


@lru_cache(maxsize=1)
def _end_shift():
    """Plan offsets moving the outline ends onto the blunt attachment's nose and TE face.

    The Converter scales the section about the frame origin, which sits a hair
    ahead of the SD7037 nose, so the blunt nose and TE move by about 0.0004 mm.
    """
    d = OUTLINE_DATA
    ends = np.array([[-d["nose_span_aft"][0], -d["nose_span_aft"][1]],
                     [-d["tail_span_aft"][0], -d["tail_span_aft"][1]]])
    nose = skin(nose_u(ATTACHMENT), ATTACHMENT)[[0, 2]]
    te = te_face(0.5, ATTACHMENT)[[0, 2]]
    return nose - ends[0], te - ends[1]


OUTLINE = {
    "kind": "plan curve (STEP X, Z): the faired tip perimeter locked from the accepted Wing B rounded-tip model",
    "formula": "cubic Bezier in (span, aft) from the attachment LE point to the attachment TE point, control arms 80 mm outboard along the parent's plan tangents; aft coordinate blended over the end 18% of p into the parent's LE/TE quadratic continuation with weight 1 - x^3(10 - 15x + 6x^2). Blunt variant only: plus plan end corrections (under 0.0004 mm) faded with the same weight over 18% of p, moving the ends onto the blunt attachment's nose and TE face",
    "reach_max_mm": 60.0,
    "reach_at": "p = 0.5, tip at 600 mm projected half-span",
}
