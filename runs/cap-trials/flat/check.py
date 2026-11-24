"""Measure STEP readback, plot centerline-symmetric zebra review and lip sections."""
import hashlib
import json
from pathlib import Path

import cadquery as cq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.gp import gp_Pnt, gp_Vec
from scipy.optimize import brentq

from build import CAP, CHORD, CUT, HEIGHT, OUTPUT, X

HERE = Path(__file__).resolve().parent


def unit(v):
    return v / np.linalg.norm(v)


def measure(surface, u, v, outward):
    point = gp_Pnt()
    vectors = [gp_Vec() for _ in range(5)]
    surface.D2(float(u), float(v), point, *vectors)
    du, dv, duu, dvv, duv = [np.array(vector.Coord()) for vector in vectors]
    normal = unit(np.cross(du, dv))
    if normal @ outward < 0:
        normal *= -1
    # Convert surface derivatives into the ambient tangent-frame normal derivative.
    jac = np.column_stack((du, dv))
    second_form = np.array([[normal @ duu, normal @ duv],
                            [normal @ duv, normal @ dvv]])
    inverse_jac = np.linalg.pinv(jac)
    normal_gradient = -inverse_jac.T @ second_form @ inverse_jac
    return np.array(point.Coord()), normal, normal_gradient


def zebra(measurement, eyes, axes):
    _, n, dn = measurement
    reflected = -eyes + 2 * (eyes @ n)[:, None] * n
    cycles = 3 * np.sum(reflected * axes, axis=1)
    gradient = 6 * ((axes @ n)[:, None] * eyes + (eyes @ n)[:, None] * axes) @ dn
    return cycles, gradient


def angle(a, b, unoriented=False):
    dot = np.sum(a * b, axis=-1)
    if unoriented:
        dot = np.abs(dot)
    return np.rad2deg(np.arctan2(np.linalg.norm(np.cross(a, b), axis=-1), dot))


def cameras(side, skin):
    for elevation in (15, 40, 65, 90):
        for azimuth in (-35, 0, 35):
            az, el = np.deg2rad([azimuth, elevation])
            eye = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
            right = np.array([-np.sin(az), np.cos(az), 0])
            up = np.cross(eye, right)
            for stripe in (0, 30, 51, 60, 90, 120, 150):
                theta = np.deg2rad(stripe)
                axis = np.cos(theta) * right + np.sin(theta) * up
                mirror = np.array([1, side, skin])
                yield [elevation, azimuth, stripe], eye * mirror, axis * mirror


def sample_u(aft):
    return brentq(lambda u: np.polynomial.polynomial.polyval(u, X) - aft, 0, 1)


def largest(values, current, location):
    best = np.argmax(values)
    value = float(values[best])
    if value > current["value"]:
        return {"value": value, **location(best)}
    return current


def verify(solid, parents, caps):
    summaries = {key: {"value": -1} for key in (
        "join_separation_mm", "join_angle_deg", "join_curvature_gap_per_mm",
        "lip_separation_mm", "lip_angle_deg", "lip_curvature_gap_per_mm",
        "parent_position_error_mm", "zebra_phase_gap_cycles",
        "zebra_gradient_gap_cycles_per_mm", "zebra_crossing_angle_deg")}
    chord_positions = np.linspace(0.12, 119.88, 1001)
    parameters = [sample_u(float(aft)) for aft in chord_positions]
    parent_error = 0.0
    mirror_error = 0.0
    outline_angle = 0.0
    for skin, parent in parents.items():
        for u in np.linspace(0, 1, 101):
            x = np.polynomial.polynomial.polyval(u, X)
            z = skin * np.polynomial.polynomial.polyval(u, HEIGHT)
            for v in np.linspace(0, 1, 15):
                expected = [x, CUT * (2 * v - 1), z]
                parent_error = max(parent_error, np.linalg.norm(
                    np.array(parent.Value(u, v).Coord()) - expected))
                right = np.array(caps[1, 1].Value(u, v).Coord())
                for side in (-1, 1):
                    actual = np.array(caps[side, skin].Value(u, v).Coord())
                    mirror_error = max(mirror_error, np.linalg.norm(actual - right * [1, side, skin]))
        for side in (-1, 1):
            for u in (0, 1):
                # The rim's chordwise derivative is defined at the collapsed corner.
                p, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
                caps[side, skin].D1(u, 1, p, du, dv)
                outline_angle = max(outline_angle, angle(np.array(du.Coord()), [0, 1, 0], unoriented=True))
    assert parent_error < 1e-8 and mirror_error < 1e-8 and outline_angle < 1e-8
    crossings = 0
    for (side, skin), cap in caps.items():
        parent = parents[skin]
        names, eyes, axes = zip(*cameras(side, skin))
        eyes, axes = np.array(eyes), np.array(axes)
        outward = np.array([0, 0, skin])
        parent_v = 1 if side == 1 else 0

        def join(u):
            return measure(parent, u, parent_v, outward), measure(cap, u, 0, outward)

        phases = []
        for u, aft in zip(parameters, chord_positions):
            body, tip = join(u)
            link = lambda _: {"side": side, "skin": skin, "aft_mm": float(aft)}
            with_view = lambda i: {**link(i), "view_elevation_azimuth_stripe_deg": names[i]}
            body_cycles, body_grad = zebra(body, eyes, axes)
            tip_cycles, tip_grad = zebra(tip, eyes, axes)
            phases.append(body_cycles)
            metrics = {
                "join_separation_mm": np.linalg.norm(body[0] - tip[0]),
                "join_angle_deg": angle(body[1], tip[1]),
                "join_curvature_gap_per_mm": np.linalg.norm(body[2] - tip[2]),
                "parent_position_error_mm": np.linalg.norm(body[0] - [
                    aft, side * CUT, skin * np.polynomial.polynomial.polyval(u, HEIGHT)]),
            }
            for key, value in metrics.items():
                summaries[key] = largest([value], summaries[key], link)
            for key, values in (
                ("zebra_phase_gap_cycles", np.abs(body_cycles - tip_cycles)),
                ("zebra_gradient_gap_cycles_per_mm", np.linalg.norm(body_grad - tip_grad, axis=1)),
            ):
                summaries[key] = largest(values, summaries[key], with_view)
            if skin == 1:
                rim_out = np.array([0, side, 0])
                top = measure(cap, u, 1, rim_out)
                bottom = measure(caps[side, -1], u, 1, rim_out)
                for key, value in (
                    ("lip_separation_mm", np.linalg.norm(top[0] - bottom[0])),
                    ("lip_angle_deg", angle(top[1], bottom[1])),
                    ("lip_curvature_gap_per_mm", np.linalg.norm(top[2] - bottom[2])),
                ):
                    summaries[key] = largest([value], summaries[key], link)

        for view, values in enumerate(np.array(phases).T):
            view_eye, view_axis = eyes[view:view+1], axes[view:view+1]
            for level in np.arange(-3, 3.1, 0.5):
                intervals = np.flatnonzero((values[:-1] - level) * (values[1:] - level) < 0)
                for i in intervals:
                    def residual(u):
                        return zebra(measure(parent, u, parent_v, outward), view_eye, view_axis)[0][0] - level

                    u = brentq(residual, parameters[i], parameters[i+1])
                    body, tip = join(u)
                    lines = [np.cross(jet[1], zebra(jet, view_eye, view_axis)[1][0]) for jet in (body, tip)]
                    if min(map(np.linalg.norm, lines)) < 1e-8:
                        continue
                    difference = angle(*lines, unoriented=True)
                    location = lambda _: {"side": side, "skin": skin, "aft_mm": float(body[0][0]),
                                          "view_elevation_azimuth_stripe_deg": names[view],
                                          "stripe_crossing_level_cycles": float(level)}
                    key = "zebra_crossing_angle_deg"
                    summaries[key] = largest([difference], summaries[key], location)
                    crossings += 1

    limits = {
        "join_separation_mm": 1e-6, "join_angle_deg": 1e-5, "join_curvature_gap_per_mm": 1e-5,
        "lip_separation_mm": 1e-6, "lip_angle_deg": 1e-5, "lip_curvature_gap_per_mm": 1e-5,
        "parent_position_error_mm": 1e-6, "zebra_phase_gap_cycles": 1e-6,
        "zebra_gradient_gap_cycles_per_mm": 1e-5, "zebra_crossing_angle_deg": 1e-4,
    }
    assert crossings > 0 and all(summaries[k]["value"] <= limit for k, limit in limits.items())
    box = solid.BoundingBox()
    assert abs(box.xlen - CHORD) < 1e-6 and abs(box.ylen - 2 * (CUT + CAP)) < 1e-6
    return {
        "step_sha256": hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
        "solid_count": 1, "shell_closed": True, "strict_occt_valid": True, "face_count": 6,
        "chord_mm": box.xlen, "span_mm": box.ylen, "thickness_mm": box.zlen,
        "volume_mm3": solid.Volume(), "cap_start_abs_span_mm": CUT,
        "parent_interior_readback_gap_mm": float(parent_error),
        "cap_reflection_symmetry_gap_mm": float(mirror_error),
        "outline_endpoint_tangent_gap_deg": float(outline_angle),
        "sampling": {"chord_range_mm": [0.12, 119.88], "stations_per_join": len(parameters),
                     "view_settings_per_join": len(names), "stripe_crossings": crossings,
                     "endpoints": "Cap surface parameters collapse at 0 and 120 mm aft."},
        "maxima": summaries, "study_acceptance_limits": limits,
    }


def sections(parents, caps):
    cuts = (0.12, 6, 18, 36, 60, 84, 102, 114, 118.8, 119.88)
    fig, (shape, curvature) = plt.subplots(2, 1, figsize=(10, 8), layout="constrained")
    stations = np.linspace(0, 1, 401)
    results = []
    for aft in cuts:
        u = sample_u(aft)
        cap = caps[1, 1]
        jets = [measure(cap, u, v, [0, 1, 1]) for v in stations]
        xyz = np.array([jet[0] for jet in jets])
        y, z = xyz[:, 1], xyz[:, 2]
        # Local curvature of the physical x=constant cut, not normal curvature
        # projected onto the varying cap normal.
        trace_curvature = []
        for v in stations:
            p, a, b, c, d, e = gp_Pnt(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec(), gp_Vec()
            cap.D2(u, float(v), p, a, b, c, d, e)
            trace_curvature.append(abs(b.Y() * d.Z() - b.Z() * d.Y()) / b.Magnitude()**3)
        parent_points = np.array([parents[1].Value(u, float(v)).Coord()
                                  for v in np.linspace(0.98, 1, 30)])
        plot_y = np.r_[parent_points[:, 1], y]
        plot_z = np.r_[parent_points[:, 2], z]
        line, = shape.plot(plot_y, plot_z, label=f"{aft:g}")
        shape.plot(plot_y, -plot_z, color=line.get_color())
        curvature.plot(plot_y, np.r_[np.zeros(len(parent_points)), trace_curvature], color=line.get_color())
        peak = int(np.argmax(trace_curvature))
        results.append({"aft_mm": aft, "rim_abs_span_mm": float(y[-1]),
                        "radius_at_rim_mm": 1 / trace_curvature[-1],
                        "peak_curvature_per_mm": trace_curvature[peak],
                        "peak_curvature_abs_span_mm": float(y[peak]),
                        "curvature_non_decreasing": bool(np.all(np.diff(trace_curvature) >= -1e-10))})
    shape.set(title="Fixed-aft STEP sections through the right cap", ylabel="Height (mm)",
              aspect="equal", xlim=(336, 402))
    shape.tick_params(labelbottom=False)
    shape.legend(title="Aft station (mm)", frameon=False, ncol=5, loc="lower center",
                 bbox_to_anchor=(0.5, -0.31))
    curvature.set(xlabel="Span station (mm)", ylabel="Absolute curvature (1/mm)", xlim=(336, 402),
                  yscale="symlog", ylim=(0, None))
    for axis in (shape, curvature):
        axis.axvline(CUT, color="#666666", lw=0.65, linestyle="--")
        axis.grid(alpha=0.2)
    fig.savefig(HERE / "sections.png", dpi=170)
    plt.close(fig)
    return results


def render(parents, caps):
    # Orthographic top and bottom. Since x is independent of v and each visible
    # skin is single valued over x,y, pixel centers can be evaluated directly.
    fig, plots = plt.subplots(2, 1, figsize=(13, 4.7), layout="constrained")
    pixel_span = np.linspace(-402, 402, 3861)
    pixel_aft = np.linspace(0.0125, 119.9875, 681)
    for skin, plot in zip((1, -1), plots):
        pixels = np.full((len(pixel_aft), len(pixel_span)), np.nan)
        for i, aft in enumerate(pixel_aft):
            u = sample_u(aft)
            for side in (-1, 1):
                parent = parents[skin]
                cap = caps[side, skin]
                _, _, axis = next((labels, eye, axis) for labels, eye, axis in cameras(side, skin)
                                  if labels == [90, 0, 51])
                center = (pixel_span * side >= 0) & (pixel_span * side <= CUT)
                _, n, _ = measure(parent, u, 0.5, [0, 0, skin])
                pixels[i, center] = np.sin(12 * np.pi * (n @ axis) * (n[2] * skin))
                rim = abs(cap.Value(u, 1).Y())
                mask = (pixel_span * side > CUT) & (pixel_span * side <= rim)
                roll = 1 - np.sqrt(1 - (pixel_span[mask] * side - CUT) / (rim - CUT))
                for j, v in zip(np.flatnonzero(mask), roll):
                    p, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
                    cap.D1(u, float(v), p, du, dv)
                    n = unit(np.array(du.Crossed(dv).Coord()))
                    pixels[i, j] = np.sin(12 * np.pi * (n @ axis) * (n[2] * skin))
        plot.set_facecolor("#bcc7d3")
        plot.imshow(pixels, cmap="gray", vmin=-0.2, vmax=0.2, origin="lower",
                    extent=(-402, 402, 0, 120), interpolation="nearest", aspect="equal")
        plot.set(ylabel="Aft (mm)", title="Upper skin" if skin == 1 else "Lower skin")
    plots[-1].set_xlabel("Span from centerline (mm)")
    fig.suptitle("Centered orthographic CAD-normal zebra • mirrored illumination • 51° stripes")
    fig.savefig(HERE / "zebra.png", dpi=180)
    plt.close(fig)


def main():
    solids = cq.importers.importStep(str(OUTPUT)).val().Solids()
    assert len(solids) == 1
    solid = solids[0]
    assert BRepCheck_Analyzer(solid.wrapped, True, False, True).IsValid()
    assert len(solid.Shells()) == 1 and BRep_Tool.IsClosed_s(solid.Shells()[0].wrapped)
    assert len(solid.Faces()) == 6
    parents, caps = {}, {}
    for face in solid.Faces():
        _, span, height = face.Center().toTuple()
        skin = 1 if height > 0 else -1
        surface = BRepAdaptor_Surface(face.wrapped)
        if abs(span) < 1:
            parents[skin] = surface
        else:
            caps[1 if span > 0 else -1, skin] = surface
    report = verify(solid, parents, caps)
    report["fixed_aft_sections"] = sections(parents, caps)
    render(parents, caps)
    (HERE / "readback.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
