"""Write spec.json for every parent from its parent.py (and parent.step for Wing B)."""

import json

import numpy as np

from common import PARENTS, SHAPES, VARIANTS, load_parent

MINIMUM_RADIUS = 0.4  # printable minimum radius of the cap rule
THICK_RULE = 2 * MINIMUM_RADIUS  # 0.8 mm

SPAN_TEXT = {
    "flat": {"parent": "Y from -350 to +350 mm (full span between the attachments)",
             "finished_projected_span_mm": 800.0},
    "twisted": {"parent": "Y from -440 to +440 mm (full span between the attachments; left half mirrors the right in XZ)",
                "finished_projected_span_mm": 1000.0},
    "swept55": {"parent": "Y from -390 to +390 mm (full span between the attachments; left half mirrors the right in XZ)",
                "finished_projected_span_mm": 900.0},
    "fin40": {"parent": "Z from 0 (flat root) to 240 mm", "finished_height_mm": 280.0},
    "fin5": {"parent": "Z from 0 (flat root) to 240 mm", "finished_height_mm": 245.0},
    "wingb": {"parent": "projected span from the root plane X = 0 to the attachment plane X = -540 mm (half-wing)",
              "finished_projected_half_span_mm": 600.0},
}

CHORD_PARAMETER = {
    "flat": "s in [0, 1] from the leading edge to the trailing-edge end of each skin; native u = s * u_end (u_end = 1 sharp)",
    "twisted": "s in [0, 1] from the leading edge to the trailing-edge end of each skin; native u = s + (u_end - 1) max(0, 2s - 1)^5 (u = s sharp)",
    "swept55": "s in [0, 1] from the leading edge to the trailing-edge end of each skin; native q = s * q_end (q_end = 1 sharp)",
    "fin40": "s in [0, 1] from the leading edge to the trailing-edge end of each skin; ellipse angle theta = pi s * theta_end/pi",
    "fin5": "s in [0, 1] from the leading edge to the trailing-edge end of each skin; ellipse angle theta = pi s * theta_end/pi",
    "wingb": "native B-spline u of the skin (lower skin from u = 0 at the trailing edge to the nose near u = 0.5, upper skin from the nose to u = 1), v from the root (0) to the attachment; see parent.py nose_u and side_range",
}


def table_fractions(start, end):
    base = start + (end - start) * (1 - np.cos(np.linspace(0, np.pi, 161))) / 2
    offsets = np.array([1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3])
    return np.unique(np.r_[base, start + offsets, end - offsets])


def crossings(thickness, start, end, level=THICK_RULE):
    xs = np.unique(np.r_[start + (end - start) * (1 - np.cos(np.linspace(0, np.pi, 4001))) / 2,
                         start + np.geomspace(1e-9, 1e-2, 200), end - np.geomspace(1e-9, 1e-2, 200)])
    xs = xs[(xs >= start) & (xs <= end)]
    t = thickness(xs) - level
    found = []
    for i in np.flatnonzero(np.isfinite(t[:-1]) & np.isfinite(t[1:]) & (np.sign(t[:-1]) != np.sign(t[1:]))):
        lo, hi = xs[i], xs[i + 1]
        f_lo = t[i]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            f_mid = thickness(np.array([mid]))[0] - level
            if np.sign(f_mid) == np.sign(f_lo):
                lo, f_lo = mid, f_mid
            else:
                hi = mid
        found.append(float(0.5 * (lo + hi)))
    return found


def section_ends(module, span):
    """Chord fractions of the most forward point and of the trailing edge in the finished section."""
    frame = module.section_frame(span)
    o, c, ch = frame["origin"], frame["chord_direction"], frame["chord"]
    if module.SHAPE == "wingb":
        v = float(module.v_at_span(span))
        nose = float(module.nose_u(span))
        skin = module.surfaces()[0]
        pts = np.array([skin.Value(u, v).Coord() for u in (nose, 0.0, 1.0)])
        x = (pts - o) @ c / ch
        return max(0.0, float(x[0])), float(min(x[1], x[2]))
    corners = np.array([module.skin(1.0, span, side) for side in (+1, -1)])
    te = float(((corners - o) @ c / ch).min())
    return 0.0, te


def wingb_cut():
    """Keep-chord cut of the Converter on the Wing B definition sections at the root and the attachment."""
    import sys
    from common import REPO, converter
    sys.path.insert(0, str(REPO / "src"))
    from planeopt.phase0_geometry import _profile, load_definition, section_frame
    geometry = converter()
    wing = load_definition(REPO / "runs/paired-wing-trial/a/bowed/definition.json").wing
    profile = np.array(_profile(REPO / "aircraft/phase0_sd7037"))
    out = {}
    for name, fraction in (("root", 0.0), ("attachment", 0.9)):
        chord = section_frame(wing, fraction).chord_mm
        flat = np.array(geometry.blunt_trailing_edge([(x * chord, y * chord) for x, y in profile], 0.4, keep_chord=True))
        scale = flat[:, 1].max() / (profile[:, 1].max() * chord)
        out[name] = {"scale": float(scale), "cut_chord_fraction": float(flat[0, 0] / scale / chord),
                     "te_chord_position_mm": float(flat[0, 0]), "te_thickness_mm": float(abs(flat[0, 1] - flat[-1, 1])),
                     "note": "from the Converter's own output on the Phase 0 definition section"}
    return out


def spec(shape, variant):
    m = load_parent(shape, variant)
    attachment = m.ATTACHMENT
    frame = m.section_frame(attachment)
    start, end = section_ends(m, attachment)
    thickness = lambda xi: m.section_thickness(xi, attachment)
    xi = table_fractions(start, end)
    table = [{"chord_fraction": float(a), "thickness_mm": float(b)} for a, b in zip(xi, thickness(xi))]
    if shape == "wingb":
        stations = {"attachment_plane": "X = -540 mm", "projected_span_mm": attachment}
        if variant == "blunt":
            ends = m.te_face(np.array([0.0, 1.0]), attachment)
            te_thickness = float(np.linalg.norm(ends[1] - ends[0]))
        else:
            pts = m.skin(np.array([0.0, 1.0]), attachment)
            te_thickness = float(np.linalg.norm(pts[0] - pts[1]))
        cut = wingb_cut() if variant == "blunt" else None
        outline_parameters = m.OUTLINE_DATA
        p = np.linspace(0, 1, 21)
        samples = {"p": p.tolist(), "X_Z": m.tip_outline(p).tolist()}
    else:
        stations = ({"Y_mm": [-attachment, attachment]} if m.FULL_SPAN else {"Z_mm": [attachment]})
        corners = np.array([m.skin(1.0, attachment, side) for side in (+1, -1)])
        te_thickness = float(np.linalg.norm(corners[0] - corners[1]))
        cut = {"root": {k: float(np.ravel(v)[0]) for k, v in m.cut(np.array([0.0])).items()},
               "attachment": {k: float(np.ravel(v)[0]) for k, v in m.cut(np.array([attachment])).items()}}
        if shape.startswith("fin"):
            h = np.linspace(0, 1, 11)
            samples = {"h": h.tolist(), "leading_edge_branch_XYZ": m.tip_outline(h, -1).tolist(),
                       "trailing_edge_branch_XYZ": m.tip_outline(h, +1).tolist()}
        else:
            p = np.linspace(0, 1, 21)
            samples = {"p": p.tolist(), "right_tip_XYZ": m.tip_outline(p, +1).tolist()}
        outline_parameters = {"reach_mm": m.REACH}
    outline = dict(m.OUTLINE)
    outline.update({"parameters": outline_parameters, "samples": samples,
                    "function": "parent.py tip_outline / tip_outline_jet"})
    return {
        "shape": shape,
        "variant": variant,
        "units": "mm",
        "axes": m.AXES,
        "span_axis": m.SPAN_AXIS,
        "extent": SPAN_TEXT[shape],
        "attachment_stations": stations,
        "cap_reach_mm": m.REACH,
        "locked_outline": outline,
        "chord_parameter": CHORD_PARAMETER[shape],
        "attachment_section": {
            "frame_origin_mm": frame["origin"].tolist(),
            "chord_direction": frame["chord_direction"].tolist(),
            "thickness_direction": frame["up_direction"].tolist(),
            "chord_mm": float(frame["chord"]),
            "leading_point_chord_fraction": start,
            "trailing_edge_chord_fraction": end,
            "trailing_edge_thickness_mm": te_thickness,
        },
        "trailing_edge": ("sharp (closed to a point)" if variant == "sharp" else
                          "blunt: Airfoil Converter 'TE thickness' 0.4 mm with 'Keep chord after the cut', "
                          "applied at every span station; flat trailing-edge face 0.4 mm thick"),
        "keep_chord_cut": cut,
        "attachment_thickness": {
            "definition": "thickness across the finished attachment section, perpendicular to the chord line, "
                          "at chord fraction measured from the section's 2D origin (leading-edge point) along the chord",
            "rule_threshold_mm": THICK_RULE,
            "chord_fractions_where_thickness_is_0p8_mm": crossings(thickness, start, end),
            "table": table,
        },
        "files": {"parent.py": "exact definition (see its docstring)", "parent.step": "parent alone, one closed solid"},
    }


def main():
    for shape in SHAPES:
        for variant in VARIANTS:
            data = spec(shape, variant)
            (PARENTS / shape / variant / "spec.json").write_text(json.dumps(data, indent=1) + "\n")
            a = data["attachment_section"]
            print(shape, variant, "chord", round(a["chord_mm"], 4), "TE", round(a["trailing_edge_thickness_mm"], 6),
                  "0.8 crossings", [round(x, 6) for x in data["attachment_thickness"]["chord_fractions_where_thickness_is_0p8_mm"]],
                  "rows", len(data["attachment_thickness"]["table"]), flush=True)


if __name__ == "__main__":
    main()
