"""Stage 1 of the Wing B blended winglet: the law, its outline and the checks on it.

Run from the repository root:

    uv run --locked --extra phase0 --extra cad python runs/winglet-b/stage1/build.py \
        --converter-src <Airfoil-Converter checkout at d77bbdf>/src

Writes definition.json, report.json, curves.json and the review PNGs next to this file.
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
    BlendedWinglet,
    ProfileFamily,
    WingletDefinition,
    decompose,
    naca_symmetric,
    nose_radius,
)

PARENT = REPO / "runs/wing-b-spar/shape/constrained/redistributed/definition.json"
FIXTURE = REPO / "aircraft/phase0_sd7037"
CONVERTER_REVISION = "d77bbdf45a8c0aa92efbc2a3d3007cfe4f904105"
RULE_MM = 0.4


def families(converter, definition):
    sd7037 = [tuple(p) for p in _profile(FIXTURE.resolve())]
    naca = [tuple(p) for p in naca_symmetric(definition.target_thickness_ratio)]
    te = definition.trailing_edge_mm

    def blunt(points):
        return lambda chord: converter.blunt_trailing_edge(points, te / chord, keep_chord=True)

    parent = ProfileFamily(blunt(sd7037), np.linspace(55.0, 205.0, 51))
    target = ProfileFamily(blunt(naca), np.linspace(18.0, 100.0, 42))
    return parent, target


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


def main():
    args = argparse.ArgumentParser()
    args.add_argument("--converter-src", required=True, type=Path)
    opts = args.parse_args()
    sys.path.insert(0, str(opts.converter_src))
    from airfoil_converter import geometry as converter

    parent_def = load_definition(PARENT)
    wing = parent_def.wing
    definition = WingletDefinition()
    parent_family, target_family = families(converter, definition)
    w = BlendedWinglet(wing, definition, parent_family, target_family)
    stations = w.stations(121, 61)
    frames = [w.frame(float(s)) for s in stations]
    tip = frames[-1]

    # ------------------------------------------------------------ measurements
    le = np.array([f.le for f in frames])
    te = np.array([f.te for f in frames])
    front = le[:, :2]
    seg = np.diff(front, axis=0)
    dihedral = np.degrees(np.arctan2(seg[:, 1], -seg[:, 0]))
    leg = stations >= w.L - 1e-9
    leg_dir = le[leg][-1] - le[leg][0]
    cant = math.degrees(math.atan2(leg_dir[1], -leg_dir[0]))
    sweep = math.degrees(math.atan2(-leg_dir[2], math.hypot(leg_dir[0], leg_dir[1])))
    chord_vec = tip.le - tip.te
    outboard = -np.array([math.sin(tip.phi), math.cos(tip.phi), 0.0])  # away from the upper skin
    toe_out = math.degrees(math.asin(np.dot(chord_vec, outboard) / np.linalg.norm(chord_vec)))
    te_gap = [float(np.linalg.norm(f.point(f.profile.x1, f.profile.camber[-1] + f.profile.half_thickness[-1])
                                   - f.point(f.profile.x1, f.profile.camber[-1] - f.profile.half_thickness[-1])))
              for f in frames]

    le_rows = []
    for s, f, i in zip(stations, frames, range(len(frames))):
        up, lo = f.profile.upper(), f.profile.lower()
        r_s = nose_radius(np.vstack([up[::-1], lo[1:]])) * f.chord_mm
        j = min(max(i, 1), len(frames) - 2)
        tangent = le[j + 1] - le[j - 1]
        tangent /= np.linalg.norm(tangent)
        cos_sweep = abs(float(np.dot(tangent, f.plane_normal)))
        le_rows.append({"s_mm": float(s), "chord_mm": f.chord_mm, "streamwise_mm": r_s,
                            "normal_to_le_mm": r_s / cos_sweep, "morph": w.morph(float(s)),
                            "thickness_ratio": float(2 * f.profile.half_thickness.max())})
    worst_le = min(le_rows, key=lambda r: r["normal_to_le_mm"])
    parent_le = []
    for fpar in np.linspace(0.0, w.f0, 19):
        pf = w.parent_frame(float(fpar))
        up, lo = pf.profile.upper(), pf.profile.lower()
        parent_le.append({"fraction": float(fpar), "chord_mm": pf.chord_mm,
                              "streamwise_mm": nose_radius(np.vstack([up[::-1], lo[1:]])) * pf.chord_mm})
    fitted = decompose(np.array(_profile(FIXTURE.resolve())))

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

    kappa = np.array([w.kappa(float(s)) for s in stations])
    report = {
        "parent": str(PARENT.relative_to(REPO)), "converter_revision": CONVERTER_REVISION,
        "law": w.to_dict(),
        "span": {"with_cap_allowance_mm": w.projected_span(True), "without_allowance_mm": w.projected_span(False)},
        "measured": {
            "height_mm": float(tip.le[1] - w.le0[1]), "cant_deg": cant, "leg_sweep_deg": sweep,
            "toe_out_deg": toe_out, "tip_chord_mm": tip.chord_mm,
            "dihedral_monotonic": bool(np.all(np.diff(dihedral) >= -1e-9)),
            "dihedral_at_join_deg": float(dihedral[0]), "dihedral_at_tip_deg": float(dihedral[-1]),
            "min_bend_radius_mm": float(1 / kappa.max()),
            "min_bend_radius_at_s_mm": float(stations[int(np.argmax(kappa))]),
            "trailing_edge_mm": {"min": min(te_gap), "max": max(te_gap)},
            "every_section_has_positive_thickness": bool(all((f.profile.half_thickness[1:] > 0).all() for f in frames)),
        },
        "le_radius": {"rule_mm": RULE_MM, "worst": worst_le,
                       "passes": bool(worst_le["normal_to_le_mm"] >= RULE_MM), "stations": le_rows,
                       "parent_stations": parent_le,
                       "sd7037_spline_nose_radius_per_chord": fitted.le_radius},
        "join_continuity": join,
    }
    (HERE / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    (HERE / "definition.json").write_text(json.dumps(w.to_dict(), indent=2) + "\n")

    parent_fracs = np.linspace(0.0, w.f0, 91)
    pframes = [w.parent_frame(float(f)) for f in parent_fracs]
    curves = {
        "winglet_le": le.tolist(), "winglet_te": te.tolist(),
        "winglet_quarter_chord": [f.point(f.profile.x0 + 0.25 * (f.profile.x1 - f.profile.x0),
                                       float(np.interp(0.25, (np.linspace(0, 1, len(f.profile.camber))), f.profile.camber))).tolist()
                               for f in frames],
        "parent_le": [f.le.tolist() for f in pframes], "parent_te": [f.te.tolist() for f in pframes],
        "sections": {f"s{float(s):07.2f}": w.frame(float(s)).section().tolist()
                  for s in np.concatenate([np.linspace(0, w.L, 9), np.linspace(w.L, w.s_tip, 6)[1:]])},
    }
    (HERE / "curves.json").write_text(json.dumps(curves) + "\n")
    plots(w, frames, pframes, stations, dihedral, kappa, le_rows, report)
    print(json.dumps({k: report[k] for k in ("law", "span", "measured")}, indent=2))
    print("LE radius worst", worst_le, "passes", report["le_radius"]["passes"])
    print("join", json.dumps(join, indent=2))


def plots(w, frames, pframes, stations, dihedral, kappa, le_rows, report):
    sect_s = np.concatenate([np.linspace(0, w.L, 7), np.linspace(w.L, w.s_tip, 5)[1:]])
    sections = [w.frame(float(s)).section() for s in sect_s]
    psections = [pf.section() for pf in pframes[::10]] + [pframes[-1].section()]
    le = np.array([f.le for f in frames]); te = np.array([f.te for f in frames])
    ple = np.array([f.le for f in pframes]); pte = np.array([f.te for f in pframes])

    def view(ax, a, b, zoom, title, flip_b=False):
        sign = -1 if flip_b else 1
        for sec in psections:
            ax.plot(-sec[:, 0], sign * sec[:, b], color="0.75", lw=0.6)
        for sec in sections:
            ax.plot(-sec[:, 0], sign * sec[:, b], color="tab:blue", lw=0.8)
        for curve, colour in ((ple, "0.3"), (pte, "0.3"), (le, "tab:red"), (te, "tab:green")):
            ax.plot(-curve[:, 0], sign * curve[:, b], color=colour, lw=1.2)
        ax.axvline(600, color="tab:purple", ls="--", lw=0.8)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=10)
        if zoom:
            ax.set_xlim(470, 615)
        ax.grid(alpha=0.3)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), gridspec_kw={"height_ratios": [1, 1.6]})
    view(axes[0, 0], 0, 1, False, "Front view, whole half-wing (outboard to the right, up)")
    view(axes[0, 1], 0, 2, False, "Plan view, whole half-wing (outboard right, aft down)")
    view(axes[1, 0], 0, 1, True, "Front view, tip (dashed: 600 mm projected half-span)")
    view(axes[1, 1], 0, 2, True, "Plan view, tip")
    for ax in axes.flat:
        ax.set_xlabel("outboard (mm)")
    axes[0, 0].set_ylabel("up (mm)"); axes[1, 0].set_ylabel("up (mm)")
    axes[0, 1].set_ylabel("forward (mm)"); axes[1, 1].set_ylabel("forward (mm)")
    fig.tight_layout()
    fig.savefig(HERE / "views.png", dpi=130)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    s = stations
    axes[0, 0].plot(s, [math.degrees(f.gamma) for f in frames], label="LE dihedral (law)")
    axes[0, 0].plot(s, [math.degrees(f.phi) for f in frames], label="section-plane tilt")
    axes[0, 0].set_title("Angles in the front view (deg)")
    axes[0, 1].plot(s, kappa * 1000)
    axes[0, 1].set_title(f"Bend curvature (1/m); tightest radius {1 / kappa.max():.1f} mm")
    axes[0, 2].plot(s, [f.chord_mm for f in frames])
    axes[0, 2].set_title("Chord (mm)")
    axes[1, 0].plot(s, [f.twist_deg for f in frames])
    axes[1, 0].set_title("Twist about the LE (deg, negative = toe-out on the leg)")
    axes[1, 1].plot(s, [r["morph"] for r in le_rows], label="morph (0 SD7037, 1 target)")
    axes[1, 1].plot(s, [10 * r["thickness_ratio"] for r in le_rows], label="max t/c x10")
    axes[1, 1].set_title("Profile morph")
    axes[1, 2].plot(s, [r["streamwise_mm"] for r in le_rows], label="streamwise")
    axes[1, 2].plot(s, [r["normal_to_le_mm"] for r in le_rows], label="normal to LE")
    axes[1, 2].axhline(RULE_MM, color="tab:red", ls="--", label="0.4 mm rule")
    axes[1, 2].set_title("Leading-edge radius (mm)")
    for ax in axes.flat:
        ax.axvline(w.L, color="0.5", ls=":", lw=0.8)
        ax.set_xlabel("s: front-view LE arc length from 540 mm (mm)")
        ax.grid(alpha=0.3)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "laws.png", dpi=110)
    plt.close(fig)

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(projection="3d")
    for sec in sections:
        ax.plot(-sec[:, 0], -sec[:, 2], sec[:, 1], color="tab:blue", lw=0.8)
    for sec in psections[-4:]:
        ax.plot(-sec[:, 0], -sec[:, 2], sec[:, 1], color="0.6", lw=0.6)
    for curve, colour in ((le, "tab:red"), (te, "tab:green")):
        ax.plot(-curve[:, 0], -curve[:, 2], curve[:, 1], color=colour)
    ax.set_xlabel("outboard"); ax.set_ylabel("aft"); ax.set_zlabel("up")
    ax.set_box_aspect((1, 1.4, 1))
    ax.set_xlim(470, 620); ax.set_ylim(40, 250); ax.set_zlim(20, 170)
    ax.view_init(elev=22, azim=-125)
    fig.tight_layout()
    fig.savefig(HERE / "oblique.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
