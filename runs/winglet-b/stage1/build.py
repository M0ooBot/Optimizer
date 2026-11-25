"""Stage 1 of the Wing B blended winglet: the law, its outline and the checks on it.

Builds each bend distribution in VARIANTS into variants/<name>/ (definition, report,
curves, review PNGs) and a comparison of their bends. Run from the repository root:

    uv run --locked --extra phase0 --extra cad python runs/winglet-b/stage1/build.py \
        --converter-src <Airfoil-Converter checkout at d77bbdf>/src
"""

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "src"))
from planeopt.phase0_geometry import _profile, load_definition
from planeopt.phase0_winglet import (
    XI,
    BlendedWinglet,
    ProfileFamily,
    WingletDefinition,
    naca_symmetric,
    nose_radius,
)

PARENT = REPO / "runs/wing-b-spar/shape/constrained/redistributed/definition.json"
FIXTURE = REPO / "aircraft/phase0_sd7037"
CONVERTER_REVISION = "d77bbdf45a8c0aa92efbc2a3d3007cfe4f904105"
RULE_MM = 0.4

# The configuration study: one setting at a time around the baseline. bend_fraction=None lets
# the law finish the bend early (a straight leg after it) only when the full-length bend does
# not fit the span; it then picks the gentlest bend that does.
BASELINE = {"height_mm": 60.0, "te_sweep_deg": 10.0, "tip_chord_mm": None, "bend_fraction": None,
            "ramp_in": 0.35, "ramp_out": 0.35}
BASELINE_LABEL = "Baseline: 60 mm, 75° cant, overall sweep 35° LE / 10° TE aft, 2° toe-out, centred bend"
GROUPS = {
    "bend": ("Bend distribution (ramps in / out)", [
        ("even", "Even 20 / 20", {"ramp_in": 0.2, "ramp_out": 0.2}),
        ("centred", "Centred 35 / 35", {}),
        ("bell", "Bell 50 / 50", {"ramp_in": 0.5, "ramp_out": 0.5}),
        ("early", "Early peak 30 / 60", {"ramp_in": 0.3, "ramp_out": 0.6}),
        ("late", "Late peak 50 / 25", {"ramp_in": 0.5, "ramp_out": 0.25}),
    ]),
    "te": ("Trailing-edge sweep, overall", [
        (f"te{v:+d}", f"{abs(v)}° {'forward' if v < 0 else 'aft'}" if v else "0°", {"te_sweep_deg": float(v)})
        for v in (-10, 0, 10, 15, 20, 25)]),
    "le": ("Leading-edge sweep, overall", [
        (f"le{v}", f"{v}°", {"sweep_deg": float(v)}) for v in (20, 25, 35, 45, 50)]),
    "height": ("Height", [(f"h{v}", f"{v} mm", {"height_mm": float(v)}) for v in (45, 50, 55, 60, 65, 70)]),
    "cant": ("Cant at the tip", [(f"cant{v}", f"{v}°", {"cant_deg": float(v)}) for v in (60, 65, 70, 75, 80, 85)]),
    "toe": ("Toe-out", [(f"toe{v}", f"{v}°", {"toe_out_deg": float(v)}) for v in (0, 2, 4, 6)]),
    "combined": ("Combinations", [
        ("tall-canted", "70 mm and 70° cant", {"height_mm": 70.0, "cant_deg": 70.0}),
    ]),
}


def families(converter, definition):
    sd7037 = [tuple(p) for p in _profile(FIXTURE.resolve())]
    naca = [tuple(p) for p in naca_symmetric(definition.target_thickness_ratio)]
    te = definition.trailing_edge_mm

    def blunt(points):
        return lambda chord: converter.blunt_trailing_edge(points, te / chord, keep_chord=True)

    parent = ProfileFamily(blunt(sd7037), np.linspace(55.0, 205.0, 51))
    target = ProfileFamily(blunt(naca), np.linspace(18.0, 100.0, 42))
    return parent, target


def make(converter, **overrides):
    definition = WingletDefinition(**{**BASELINE, **overrides})
    parent_family, target_family = families(converter, definition)
    return BlendedWinglet(load_definition(PARENT).wing, definition, parent_family, target_family)


def one_sided(f, s0, h, side):
    """First and second derivatives from one side, 4th-order stencils."""
    k = np.arange(6) * side * h
    v = np.array([f(s0 + kk) for kk in k])
    d1 = side * (-25 * v[0] + 48 * v[1] - 36 * v[2] + 16 * v[3] - 3 * v[4]) / (12 * h)
    d2 = (45 * v[0] - 154 * v[1] + 214 * v[2] - 156 * v[3] + 61 * v[4] - 10 * v[5]) / (12 * h * h)
    return v[0], d1, d2


def skin_point(w, xi_index, upper, parent=False):
    """A skin point at fixed profile index. With parent=True it is evaluated on the parent
    wing (s <= 0), so a join check compares the two independent definitions."""
    def f(s):
        frame = w.parent_frame(w.parent_fraction(s), s) if parent else w.frame(s)
        side = frame.profile.upper() if upper else frame.profile.lower()
        x, y = side[xi_index]
        return frame.point(x, y)
    return f


def angle_along(points):
    """Front-view direction of a polyline above horizontal, per segment, degrees."""
    seg = np.diff(np.asarray(points)[:, :2], axis=0)
    return np.degrees(np.arctan2(seg[:, 1], -seg[:, 0]))


def front_view_sag(w, upto=0.6):
    """How far the winglet falls below the wing's own front-view tangent lines at the join.

    Outline: each section's highest and lowest point along the wing's up-inboard normal at
    the join, compared with the outline's tangent from the wing side. Skin lines: every
    sixth profile point on each surface, compared with that line's own tangent. Also the
    wing's own outline curvature, which the join carries on for curvature continuity."""
    g0 = w.gamma0
    normal, along = np.array([math.sin(g0), math.cos(g0)]), np.array([math.cos(g0), math.sin(g0)])

    def front(points):
        return np.column_stack([-points[..., 0], points[..., 1]])

    wing_side = [w.parent_frame(w.parent_fraction(-h)) for h in np.linspace(6, 0, 13)]
    winglet_side = [w.frame(float(s)) for s in np.linspace(0, upto * w.s_tip, 121)]
    out = {}
    for name, pick in (("top", np.argmax), ("bottom", np.argmin)):
        def outline(frames, pick=pick):
            rows = []
            for f in frames:
                p = front(f.section())
                k = pick(p @ normal)
                rows.append((p[k] @ along, p[k] @ normal))
            return np.array(rows)

        wing, winglet = outline(wing_side), outline(winglet_side)
        slope = np.polyfit(wing[-5:, 0], wing[-5:, 1], 1)[0]
        tangent = wing[-1, 1] + slope * (winglet[:, 0] - wing[-1, 0])
        out[f"outline_{name}_below_tangent_mm"] = float(max(0.0, (tangent - winglet[:, 1]).max()))
        out[f"wing_outline_{name}_curvature_per_mm"] = float(2 * np.polyfit(wing[:, 0], wing[:, 1], 2)[0])
    for side in ("upper", "lower"):
        deepest = 0.0
        for index in range(2, len(winglet_side[0].profile.camber), 6):
            def line(frame, index=index, side=side):
                x, y = (frame.profile.upper() if side == "upper" else frame.profile.lower())[index]
                point = frame.point(x, y)
                return np.array([-point[0], point[1]])
            pts = np.array([line(f) for f in winglet_side])
            before = line(w.parent_frame(w.parent_fraction(-0.5)))
            t = (pts[0] - before) / np.linalg.norm(pts[0] - before)
            deepest = max(deepest, float(-(((pts - pts[0]) @ np.array([-t[1], t[0]])).min())))
        out[f"{side}_skin_lines_below_tangent_mm"] = deepest
    return out


def side_view_turns(w):
    """How each guide turns seen from outboard (up against aft), the view the edge laws are
    drawn in. Counter-turn is the total turning, in degrees, against the guide's own overall
    direction of turn. The first stretch, where the view sees the wing nearly end-on (side-view
    speed under 30% of its largest), is skipped: directions there are not visible."""
    s = w.stations(401)
    frames = [w.frame(float(t)) for t in s]

    def turn(points):
        h, v = -points[:, 2], points[:, 1]
        dh, dv = np.gradient(h, s), np.gradient(v, s)
        speed = np.hypot(dh, dv)
        seen = speed >= 0.3 * speed.max()
        angle = np.degrees(np.unwrap(np.arctan2(dh, dv)))[seen]
        steps = np.diff(angle)
        total = float(steps.sum())
        counter = float(np.clip(-np.sign(total) * steps, 0, None).sum())
        return {"first_seen_deg": float(angle[0]), "tip_deg": float(angle[-1]), "turn_deg": total,
                "counter_turn_deg": counter}

    out = {"leading_edge": turn(np.array([f.le for f in frames])),
           "trailing_edge": turn(np.array([f.te for f in frames]))}
    for x in (0.25, 0.5, 0.75):
        pts = np.array([f.point(f.profile.x0 + x * (f.profile.x1 - f.profile.x0),
                                float(np.interp(x, XI, f.profile.camber))) for f in frames])
        out[f"chord_line_{int(100 * x)}pct"] = turn(pts)
    skin = 0.0
    for index in range(4, len(XI) - 1, 12):
        for upper in (True, False):
            pts = np.array([f.point(*(f.profile.upper() if upper else f.profile.lower())[index]) for f in frames])
            skin = max(skin, turn(pts)["counter_turn_deg"])
    out["skin_guides_worst_counter_turn_deg"] = skin
    return out


def measure(w):
    stations = w.stations(161)
    frames = [w.frame(float(s)) for s in stations]
    tip = frames[-1]
    le = np.array([f.le for f in frames])
    te = np.array([f.te for f in frames])
    last = le[-1] - le[-2]
    chord_vec = tip.le - tip.te
    outboard = -np.array([math.sin(tip.phi), math.cos(tip.phi), 0.0])  # away from the upper skin
    te_gap = [float(np.linalg.norm(
        f.point(f.profile.x1, f.profile.camber[-1] + f.profile.half_thickness[-1])
        - f.point(f.profile.x1, f.profile.camber[-1] - f.profile.half_thickness[-1]))) for f in frames]
    # The trailing edge's forward rate along the winglet; its slope is the edge's bend in side view.
    te_rate = np.diff(te[:, 2]) / np.diff(stations)
    te_bend = np.diff(te_rate) / np.diff(stations)[1:]
    main_sign = np.sign(te_bend[np.argmax(np.abs(te_bend))])
    te_reversal = float(np.max(np.clip(-main_sign * te_bend, 0, None)) / np.abs(te_bend).max())
    kappa = np.array([w.kappa(float(t)) for t in np.linspace(0, w.s_tip, 4001)])

    le_rows = []
    for i, (s, f) in enumerate(zip(stations, frames, strict=True)):
        up, lo = f.profile.upper(), f.profile.lower()
        r_s = nose_radius(np.vstack([up[::-1], lo[1:]])) * f.chord_mm
        j = min(max(i, 1), len(frames) - 2)
        tangent = (le[j + 1] - le[j - 1]) / np.linalg.norm(le[j + 1] - le[j - 1])
        cos_sweep = abs(float(np.dot(tangent, f.plane_normal)))
        le_rows.append({"s_mm": float(s), "chord_mm": float(f.chord_mm), "streamwise_mm": float(r_s),
                        "normal_to_le_mm": float(r_s / cos_sweep), "morph": w.morph(float(s)),
                        "thickness_ratio": float(2 * f.profile.half_thickness.max())})
    worst_le = min(le_rows, key=lambda r: r["normal_to_le_mm"])

    join = {}
    checks = {
        "leading_edge": (lambda s: w.parent_frame(w.parent_fraction(s), s).le, lambda s: w.frame(s).le),
        "upper_skin_30pct": (skin_point(w, 60, True, True), skin_point(w, 60, True)),
        "lower_skin_70pct": (skin_point(w, 160, False, True), skin_point(w, 160, False)),
        "trailing_edge_upper": (skin_point(w, -1, True, True), skin_point(w, -1, True)),
    }
    for name, (left, right) in checks.items():
        p0, a1, a2 = one_sided(left, 0.0, 0.2, -1)
        q0, b1, b2 = one_sided(right, 0.0, 0.2, +1)
        join[name] = {"position_gap_mm": float(np.linalg.norm(p0 - q0)),
                      "first_derivative_jump": float(np.linalg.norm(a1 - b1)),
                      "second_derivative_jump_per_mm": float(np.linalg.norm(a2 - b2)),
                      "second_derivative_scale_per_mm": float(max(np.linalg.norm(a2), np.linalg.norm(b2)))}

    chords = np.array([f.chord_mm for f in frames])
    side = side_view_turns(w)
    dihedral = angle_along(le)
    return frames, {
        "law": w.to_dict(),
        "span": {"with_cap_allowance_mm": w.projected_span(True),
                 "without_allowance_mm": w.projected_span(False)},
        "measured": {
            "attach_mm": float(-w.le0[0]),
            "height_mm": float(tip.le[1] - w.le0[1]),
            "tip_cant_deg": math.degrees(math.atan2(last[1], -last[0])),
            "tip_le_sweep_deg": math.degrees(math.atan2(-last[2], math.hypot(last[0], last[1]))),
            "toe_out_deg": math.degrees(math.asin(np.dot(chord_vec, outboard) / np.linalg.norm(chord_vec))),
            "tip_chord_mm": float(tip.chord_mm),
            "chord_decreases_throughout": bool(np.all(np.diff(chords) < 0)),
            "bend_curvature_never_negative": bool(kappa.min() >= 0),
            "dihedral_rises_throughout": bool(np.all(np.diff(dihedral) >= -1e-7)),
            # Informational: an aft-swept trailing edge must turn Wing B's slightly forward-running
            # edge around once, so it cannot also keep one direction.
            "trailing_edge_keeps_one_direction": bool(np.all(te_rate > 0) or np.all(te_rate < 0)),
            # Seen from outboard, the leading and trailing edges each turn one way only.
            "edges_turn_one_way": bool(side["leading_edge"]["counter_turn_deg"] < 0.05
                                       and side["trailing_edge"]["counter_turn_deg"] < 0.05),
            # Wing B's own trailing-edge bend is carried through the join for smoothness, so a
            # small opposite bend can remain near it; this is its size against the main bend.
            "trailing_edge_counter_bend_ratio": te_reversal,
            "trailing_edge_mm": {"min": min(te_gap), "max": max(te_gap)},
            "every_section_has_positive_thickness": bool(all((f.profile.half_thickness[1:] > 0).all()
                                                             for f in frames)),
            "front_view_sag": front_view_sag(w),
            "side_view": side,
        },
        "le_radius": {"rule_mm": RULE_MM, "worst": worst_le,
                      "passes": bool(worst_le["normal_to_le_mm"] >= RULE_MM), "stations": le_rows},
        "join_continuity": join,
    }


def write_variant(w, name, title, out):
    out.mkdir(parents=True, exist_ok=True)
    frames, report = measure(w)
    report.update(variant=name, title=title, parent=str(PARENT.relative_to(REPO)),
                  converter_revision=CONVERTER_REVISION)
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    (out / "definition.json").write_text(json.dumps(w.to_dict(), indent=2) + "\n")
    pframes = [w.parent_frame(float(f)) for f in np.linspace(0.0, w.f0, 91)]
    curves = {
        "winglet_le": [f.le.tolist() for f in frames], "winglet_te": [f.te.tolist() for f in frames],
        "parent_le": [f.le.tolist() for f in pframes], "parent_te": [f.te.tolist() for f in pframes],
        "sections": {f"s{float(s):07.2f}": w.frame(float(s)).section().tolist()
                     for s in np.linspace(0, w.s_tip, 13)},
    }
    if name == "bend-centred":
        (out / "curves.json").write_text(json.dumps(curves) + "\n")
    plot_views(w, frames, pframes, title, out / "views.png")
    plot_laws(w, frames, report, title, out / "laws.png")
    return report


def plot_views(w, frames, pframes, title, path):
    sections = [w.frame(float(s)).section() for s in np.linspace(0, w.s_tip, 11)]
    psections = [pf.section() for pf in pframes[::10]] + [pframes[-1].section()]
    le, te = np.array([f.le for f in frames]), np.array([f.te for f in frames])
    ple, pte = np.array([f.le for f in pframes]), np.array([f.te for f in pframes])

    def view(ax, b, zoom, label):
        for sec in psections:
            ax.plot(-sec[:, 0], sec[:, b], color="0.75", lw=0.6)
        for sec in sections:
            ax.plot(-sec[:, 0], sec[:, b], color="tab:blue", lw=0.8)
        for curve, colour in ((ple, "0.3"), (pte, "0.3"), (le, "tab:red"), (te, "tab:green")):
            ax.plot(-curve[:, 0], curve[:, b], color=colour, lw=1.2)
        ax.axvline(600, color="tab:purple", ls="--", lw=0.8)
        ax.axvline(-w.le0[0], color="0.5", ls=":", lw=0.8)
        ax.set_aspect("equal")
        ax.set_title(label, fontsize=10)
        if zoom:
            ax.set_xlim(470, 615)
        ax.grid(alpha=0.3)
        ax.set_xlabel("outboard (mm)")
        ax.set_ylabel("up (mm)" if b == 1 else "forward (mm)")

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), gridspec_kw={"height_ratios": [1, 1.6]})
    view(axes[0, 0], 1, False, "Front view, whole half-wing")
    view(axes[0, 1], 2, False, "Plan view, whole half-wing")
    view(axes[1, 0], 1, True, "Front view, tip (dashed: 600 mm span; dotted: winglet start)")
    view(axes[1, 1], 2, True, "Plan view, tip")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def plot_laws(w, frames, report, title, path):
    s = np.array([f.s for f in frames])
    rows = report["le_radius"]["stations"]
    kappa = np.array([w.kappa(float(t)) for t in s])
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes[0, 0].plot(s, [math.degrees(f.gamma) for f in frames], label="LE dihedral")
    axes[0, 0].plot(s, [math.degrees(f.phi) for f in frames], label="section-plane tilt")
    axes[0, 0].set_title("Angles in the front view (deg)")
    axes[0, 1].plot(s, kappa * 1000)
    axes[0, 1].set_title(f"Bend curvature (1/m); tightest radius {1 / kappa.max():.1f} mm")
    axes[0, 2].plot(s, [f.chord_mm for f in frames])
    axes[0, 2].set_title("Chord (mm), from the drawn LE and TE")
    axes[1, 0].plot(s, [f.twist_deg for f in frames])
    axes[1, 0].set_title("Twist about the LE (deg, negative = toe-out at the tip)")
    axes[1, 1].plot(s, [r["morph"] for r in rows], label="morph (0 SD7037, 1 target)")
    axes[1, 1].plot(s, [10 * r["thickness_ratio"] for r in rows], label="max t/c x10")
    axes[1, 1].set_title("Profile morph")
    axes[1, 2].plot(s, [r["streamwise_mm"] for r in rows], label="streamwise")
    axes[1, 2].plot(s, [r["normal_to_le_mm"] for r in rows], label="normal to LE")
    axes[1, 2].axhline(RULE_MM, color="tab:red", ls="--", label="0.4 mm rule")
    axes[1, 2].set_title("Leading-edge radius (mm)")
    for ax in axes.flat:
        ax.set_xlabel("s: winglet LE length in the front view (mm)")
        ax.grid(alpha=0.3)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_group(title, members, path):
    """Front and plan views of one group's leading and trailing edges, overlaid."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    for label, w in members:
        s = np.linspace(0, w.s_tip, 161)
        frames = [w.frame(float(t)) for t in s]
        le, te = np.array([f.le for f in frames]), np.array([f.te for f in frames])
        line = axes[0].plot(-le[:, 0], le[:, 1], label=label)[0]
        axes[0].plot(-te[:, 0], te[:, 1], color=line.get_color(), lw=0.8, ls="--")
        axes[1].plot(-le[:, 0], le[:, 2], color=line.get_color())
        axes[1].plot(-te[:, 0], te[:, 2], color=line.get_color(), lw=0.8, ls="--")
        axes[2].plot(s / w.s_tip, [w.kappa(float(t)) * 1000 for t in s], color=line.get_color())
    axes[0].set_title("Front view: LE solid, TE dashed (mm)")
    axes[1].set_title("Plan view: LE solid, TE dashed (mm)")
    axes[2].set_title("Bend curvature (1/m) along the winglet")
    for ax in axes[:2]:
        ax.axvline(600, color="tab:purple", ls="--", lw=0.8)
        ax.set_aspect("equal")
        ax.set_xlabel("outboard (mm)")
    axes[2].set_xlabel("fraction of the winglet's LE length")
    for ax in axes:
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=105)
    plt.close(fig)


def summarise(name, group, label, overrides, w, report):
    m, law = report["measured"], report["law"]
    return {
        "group": group, "label": label, "overrides": overrides, "fits": True,
        "attach_mm": m["attach_mm"], "length_mm": law["length_mm"], "bend_fraction": law["bend_fraction"],
        "tightest_bend_radius_mm": law["tightest_bend_radius_mm"], "height_mm": m["height_mm"],
        "tip_cant_deg": m["tip_cant_deg"], "tip_le_sweep_deg": m["tip_le_sweep_deg"],
        "le_root_angle_deg": law["leading_edge_root_angle_deg"], "le_tip_angle_deg": law["leading_edge_tip_angle_deg"],
        "te_root_angle_deg": law["trailing_edge_root_angle_deg"], "te_tip_angle_deg": law["trailing_edge_tip_angle_deg"],
        "te_overall_sweep_deg": law["trailing_edge_overall_sweep_deg"],
        "side_view": m["side_view"],
        "edges_turn_one_way": m["edges_turn_one_way"],
        "chord_line_guides_worst_counter_turn_deg": max(m["side_view"][f"chord_line_{p}pct"]["counter_turn_deg"]
                                                        for p in (25, 50, 75)),
        "toe_out_deg": m["toe_out_deg"], "tip_chord_mm": m["tip_chord_mm"],
        "span_mm": report["span"]["with_cap_allowance_mm"],
        "le_radius_min_mm": report["le_radius"]["worst"]["normal_to_le_mm"],
        "le_radius_streamwise_min_mm": min(r["streamwise_mm"] for r in report["le_radius"]["stations"]),
        "le_radius_passes": report["le_radius"]["passes"],
        "trailing_edge_mm": m["trailing_edge_mm"],
        "checks": {k: m[k] for k in ("chord_decreases_throughout", "bend_curvature_never_negative",
                                     "dihedral_rises_throughout", "edges_turn_one_way",
                                     "trailing_edge_counter_bend_ratio",
                                     "every_section_has_positive_thickness")},
        "front_view_sag": m["front_view_sag"],
        "join_max_second_derivative_jump": max(j["second_derivative_jump_per_mm"]
                                               for j in report["join_continuity"].values()),
        "join_max_position_gap_mm": max(j["position_gap_mm"] for j in report["join_continuity"].values()),
    }


def main():
    args = argparse.ArgumentParser()
    args.add_argument("--converter-src", required=True, type=Path)
    opts = args.parse_args()
    sys.path.insert(0, str(opts.converter_src))
    from airfoil_converter import geometry as converter

    out = HERE / "configs"
    summary, built = {}, {}
    for group, (group_title, members) in GROUPS.items():
        drawn = []
        for name, label, overrides in members:
            key = json.dumps({**BASELINE, **overrides}, sort_keys=True)
            config = f"{group}-{name}"
            try:
                if key not in built:
                    w = make(converter, **overrides)
                    title = f"{group_title}: {label}"
                    built[key] = (config, w, write_variant(w, config, title, out / config))
                first, w, report = built[key]
                summary[config] = summarise(config, group, label, overrides, w, report)
                summary[config]["same_as"] = first if first != config else None
                drawn.append((label, w))
            except ValueError as error:
                summary[config] = {"group": group, "label": label, "overrides": overrides,
                                   "fits": False, "reason": str(error)}
        plot_group(group_title, drawn, out / f"group-{group}.png")
    (out / "summary.json").write_text(json.dumps({"baseline": BASELINE, "baseline_label": BASELINE_LABEL,
                                                  "groups": {g: t for g, (t, _) in GROUPS.items()},
                                                  "configs": summary}, indent=2) + "\n")
    for config, row in summary.items():
        if row["fits"]:
            print(f"{config:22s} start {row['attach_mm']:6.1f}  bend {row['bend_fraction']:.2f}  "
                  f"R {row['tightest_bend_radius_mm']:5.1f}  tip chord {row['tip_chord_mm']:5.1f}  "
                  f"LE r {row['le_radius_min_mm']:.2f}/{row['le_radius_streamwise_min_mm']:.2f}  "
                  f"sag {row['front_view_sag']['outline_top_below_tangent_mm']:.3f}  "
                  f"LE {row['le_root_angle_deg']:.0f}>{row['le_tip_angle_deg']:.0f} TE {row['te_root_angle_deg']:.0f}>{row['te_tip_angle_deg']:.0f} "
                  f"one-way {row['edges_turn_one_way']} chord-line {row['chord_line_guides_worst_counter_turn_deg']:.1f} "
                  f"skin {row['side_view']['skin_guides_worst_counter_turn_deg']:.1f}  "
                  f"checks {all(v for k, v in row['checks'].items() if k != 'trailing_edge_counter_bend_ratio')}")
        else:
            print(f"{config:22s} does not fit: {row['reason'][:90]}")


if __name__ == "__main__":
    main()
