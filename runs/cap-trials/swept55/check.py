"""Verify the STEP with OCCT surface derivatives and generate review images."""
import hashlib
import json
from math import atan2, degrees

import cadquery as cq
import matplotlib.pyplot as plt
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.gp import gp_Pnt, gp_Vec
from scipy.optimize import brentq

from build import HERE, STATION, SWEEP, THICKNESS

STEP = HERE / "swept55-wing.step"
FRACTIONS = [0.02, 0.05, 0.25, 0.40, 0.50, 0.65, 0.80, 0.90, 0.97]


def unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def xyz(patch, u, v):
    return np.array(patch.Value(float(u), float(v)).Coord())


def geometry(patch, u, v, orientation):
    p = gp_Pnt()
    uv, vv, uu, vvv, mixed = (gp_Vec() for _ in range(5))
    patch.D2(float(u), float(v), p, uv, vv, uu, vvv, mixed)
    a, b, aa, bb, ab = [np.array(vector.Coord()) for vector in (uv, vv, uu, vvv, mixed)]
    n = orientation * unit(np.cross(a, b))
    jacobian = np.column_stack((a, b))
    metric_inverse = np.linalg.solve(jacobian.T @ jacobian, jacobian.T)
    second_form = np.array([[n @ aa, n @ ab], [n @ ab, n @ bb]])
    curvature = metric_inverse.T @ second_form @ metric_inverse
    return np.array(p.Coord()), n, curvature


def load():
    solids = cq.importers.importStep(str(STEP)).solids().vals()
    assert len(solids) == 1
    solid = solids[0]
    assert BRepCheck_Analyzer(solid.wrapped, True, False, True).IsValid()
    assert len(solid.Shells()) == 1 and BRep_Tool.IsClosed_s(solid.Shells()[0].wrapped)
    patches = {}
    for face in solid.Faces():
        patch = BRep_Tool.Surface_s(face.wrapped)
        assert patch.Bounds() == (0., 1., 0., 1.)
        _, y, z = xyz(patch, 0.5, 0.5)
        patches[int(np.sign(y)), int(np.sign(z)), abs(y) > STATION] = patch
    assert len(patches) == 8
    return solid, patches


def readback(solid, patches):
    bounds = solid.BoundingBox()
    skin_error = 0.
    mirror_error = 0.
    for u in np.linspace(0, 1, 41):
        for t in np.linspace(0, 1, 21):
            for side, skin, cap in patches:
                point = xyz(patches[side, skin, cap], u, t)
                mirror = xyz(patches[-side, skin, cap], u, t) * [1, -1, 1]
                mirror_error = max(mirror_error, np.linalg.norm(point - mirror))
                if not cap:
                    chord = 210 - 115 * t
                    reference = [SWEEP * STATION * t + chord * u**2,
                                 side * STATION * t, skin * chord * THICKNESS(u)]
                    skin_error = max(skin_error, np.linalg.norm(point - reference))
    parent = patches[1, 1, False]
    leading = xyz(parent, 0., 1.) - xyz(parent, 0., 0.)
    sweep = degrees(atan2(leading[0], leading[1]))
    root_chord = np.linalg.norm(xyz(parent, 1., 0.) - xyz(parent, 0., 0.))
    cap_chord = np.linalg.norm(xyz(parent, 1., 1.) - xyz(parent, 0., 1.))
    lip_u = np.sqrt(0.5)
    cap = patches[1, 1, True]
    exact_tip = xyz(cap, lip_u, 1.)
    a, b, aa, bb, ab = (gp_Vec() for _ in range(5))
    cap.D2(0., 1., gp_Pnt(), a, b, aa, bb, ab)
    nose_tangent = unit(np.array(aa.Coord()))  # q = u² at the nose.
    cap.D1(1., 1., gp_Pnt(), a, b)
    tail_tangent = unit(np.array(a.Coord()))
    trailing = xyz(parent, 1., 0.) - xyz(parent, 1., 1.)
    endpoint_angles = [float(angle(nose_tangent, unit(leading))),
                       float(angle(tail_tangent, unit(trailing)))]
    assert abs(bounds.ylen - 900.) < 1e-5
    assert np.allclose([sweep, root_chord, cap_chord, exact_tip[1]], [55., 210., 95., 450.], atol=1e-8, rtol=0)
    assert skin_error < 1e-8 and mirror_error < 1e-8
    assert max(endpoint_angles) < 1e-6
    return dict(step_sha256=hashlib.sha256(STEP.read_bytes()).hexdigest(), solid_count=1,
                closed=True, strict_occt_check=True, faces=len(patches), volume_mm3=solid.Volume(),
                projected_span_mm=bounds.ylen, positive_tip_xyz_mm=exact_tip.tolist(),
                cap_cut_y_mm=leading[1], sweep_deg=sweep, root_chord_mm=root_chord,
                cap_chord_mm=cap_chord, lip_endpoint_tangent_angles_deg=endpoint_angles,
                sampled_parent_error_mm=skin_error,
                sampled_mirror_error_mm=mirror_error)


def angle(a, b):
    return np.degrees(np.arctan2(np.linalg.norm(np.cross(a, b), axis=-1),
                                 np.sum(a*b, axis=-1)))


def reflection(n, k, view, stripe):
    # Four cycles per unit component of the reflected unit-view vector.
    phase = 4 * (2 * (n @ view) * (n @ stripe) - view @ stripe)
    gradient = -8 * np.einsum("...ij,...j->...i", k,
                             (n @ stripe)[..., None]*view + (n @ view)[..., None]*stripe)
    return phase, gradient


def lighting():
    for azimuth in (0, 45, 90, 135, 180, 225, 270, 315):
        for elevation in (15, 40, 70):
            az, el = np.radians([azimuth, elevation])
            view = np.array([np.cos(az)*np.cos(el), np.sin(az)*np.cos(el), np.sin(el)])
            e1 = unit(np.cross(view, [0, 0, 1]))
            e2 = np.cross(view, e1)
            for stripe_angle in range(0, 180, 30):
                theta = np.radians(stripe_angle)
                yield (azimuth, elevation, stripe_angle), view, np.cos(theta)*e1 + np.sin(theta)*e2


def summarize(a, q):
    index = np.argmax(a)
    return {"max": float(a[index]), "chord_fraction": float(q[index])}


def joins(patches):
    q = np.linspace(0.001, 0.999, 701)
    report = {}
    for side in (-1, 1):
        for label, skin_a, cap_a, va, skin_b, cap_b, vb in (
                ("upper_cut", 1, False, 1., 1, True, 0.),
                ("lower_cut", -1, False, 1., -1, True, 0.),
                ("lip", 1, True, 1., -1, True, 1.)):
            a, b = patches[side, skin_a, cap_a], patches[side, skin_b, cap_b]
            orient_a, orient_b = side * skin_a, side * skin_b
            pos_a, n_a, k_a = map(np.array, zip(*(geometry(a, u, va, orient_a) for u in np.sqrt(q))))
            pos_b, n_b, k_b = map(np.array, zip(*(geometry(b, u, vb, orient_b) for u in np.sqrt(q))))
            maxima = dict(phase_cycles={"max": 0.}, gradient_cycles_per_mm={"max": 0.},
                          crossing_angle_deg={"max": 0.})
            crossings = 0
            for setup, view, stripe in lighting():
                phase_a, grad_a = reflection(n_a, k_a, view, stripe)
                phase_b, grad_b = reflection(n_b, k_b, view, stripe)
                differences = dict(phase_cycles=np.abs(phase_a-phase_b),
                                   gradient_cycles_per_mm=np.linalg.norm(grad_a-grad_b, axis=-1))
                for key, values in differences.items():
                    result = summarize(values, q)
                    if result["max"] > maxima[key]["max"]:
                        maxima[key] = dict(result, view_azimuth_elevation_stripe_deg=setup)
                if np.all(np.abs(grad_a) < 1e-13):
                    continue
                # The plotted light/dark boundary is phase = k/2 + 1/4.
                for target in np.arange(-3.75, 4., 0.5):
                    signed = phase_a - target
                    for j in np.flatnonzero(signed[:-1] * signed[1:] < 0):
                        def residual(fraction):
                            _, n, k = geometry(a, np.sqrt(fraction), va, orient_a)
                            return reflection(n, k, view, stripe)[0] - target
                        f = brentq(residual, q[j], q[j+1])
                        _, n, k = geometry(a, np.sqrt(f), va, orient_a)
                        _, m, h = geometry(b, np.sqrt(f), vb, orient_b)
                        ga = reflection(n, k, view, stripe)[1]
                        gb = reflection(m, h, view, stripe)[1]
                        if min(np.linalg.norm(ga), np.linalg.norm(gb)) < 1e-9:
                            continue
                        # Evaluate both stripe directions in the same tangent plane.
                        direction_a = unit(np.cross(n, ga))
                        direction_b = unit(np.cross(n, gb))
                        delta = angle(direction_a, direction_b)
                        delta = min(delta, 180-delta)
                        crossings += 1
                        if delta > maxima["crossing_angle_deg"]["max"]:
                            maxima["crossing_angle_deg"] = dict(max=float(delta), chord_fraction=f,
                                                                view_azimuth_elevation_stripe_deg=setup)
            report[f"{side:+d}_{label}"] = dict(
                gap_mm=summarize(np.linalg.norm(pos_a-pos_b, axis=-1), q),
                tangent_plane_angle_deg=summarize(angle(n_a, n_b), q),
                curvature_tensor_per_mm=summarize(np.linalg.norm(k_a-k_b, axis=(1, 2)), q),
                resolved_stripe_crossings=crossings, **maxima)
    return report


def sections(patches):
    cap, parent = patches[1, 1, True], patches[1, 1, False]
    figures, axes = plt.subplots(3, 3, figsize=(13, 9), layout="constrained")
    sections = []
    plane_normal = unit(np.array([1., -SWEEP, 0.]))
    for fraction, ax in zip(FRACTIONS, axes.flat):
        # Intersect each actual exported patch with x - sweep*y = 95*q.
        offset = 95 * fraction
        def intersection(patch, t):
            def residual(u):
                x, y, _ = xyz(patch, u, t)
                return x - SWEEP*y - offset
            u = brentq(residual, 0., 1.)
            return u, xyz(patch, u, t)

        cut = [intersection(parent, t)[1] for t in np.linspace(0.94, 1, 81)]
        rolled = [intersection(cap, t)[1] for t in np.linspace(0, 1, 241)]
        u = intersection(cap, 1.)[0]
        lip, n, k = geometry(cap, u, 1., 1)
        tangent = unit(np.cross(plane_normal, n))
        normal_in_plane = n - (n @ plane_normal)*plane_normal
        curvature = abs(tangent @ k @ tangent) / np.linalg.norm(normal_in_plane)
        radius = float(1 / curvature)
        def section_plot(points, color):
            _, y, z = np.asarray(points).T
            ax.plot(y, z, color=color)
            ax.plot(y, -z, color=color)
        section_plot(cut, "#8493a2")
        section_plot(rolled, "#166286")
        ax.axvline(390., color="#b8c0c5", lw=0.8)
        ax.set_title(f"{100*fraction:g}% chord at cut · lip R {radius:.3f} mm", fontsize=10)
        ax.set_xlabel("Spanwise y (mm)")
        ax.set_ylabel("z (mm)")
        ax.set_ylim(-5.1, 5.1)
        ax.grid(alpha=0.15)
        sections.append(dict(cut_chord_fraction=fraction, lip_xyz_mm=lip.tolist(),
                             lip_radius_mm=radius, lip_u=float(u)))
    figures.suptitle("Oblique spanwise sections through the exported right wing", fontsize=16)
    figures.savefig(HERE / "lip-sections.png", dpi=160)
    plt.close(figures)
    return sections


def main():
    solid, patches = load()
    results = readback(solid, patches)
    results["join_samples"] = {"count_per_join": 701, "u_squared_range": [0.001, 0.999],
                               "light_setups_per_join": 144}
    results["joins"] = joins(patches)
    limits = dict(gap_mm=1e-6, tangent_plane_angle_deg=1e-4,
                  curvature_tensor_per_mm=1e-6, phase_cycles=1e-6,
                  gradient_cycles_per_mm=1e-6, crossing_angle_deg=1e-3)
    for label, measurements in results["joins"].items():
        assert measurements["resolved_stripe_crossings"] > 0, label
        for metric, limit in limits.items():
            assert measurements[metric]["max"] <= limit, (label, metric, measurements[metric])
    results["study_limits"] = limits
    results["sections"] = sections(patches)
    (HERE / "readback.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
