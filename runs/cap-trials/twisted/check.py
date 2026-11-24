"""Check CAD derivatives and render reflection evidence from the reimported STEP."""

import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

plot_cache = TemporaryDirectory(dir=Path(__file__).resolve().parent)
os.environ.setdefault("MPLCONFIGDIR", plot_cache.name)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from OCP.BRep import BRep_Tool
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepTools import BRepTools
from OCP.Bnd import Bnd_Box
from OCP.gp import gp_Pnt, gp_Vec
from scipy.optimize import brentq, minimize_scalar

from build import CUT, HERE, read_solid, section


STEP = HERE / "twisted-wing.step"
CHORD_PARAMS = np.linspace(.001, .999, 1001)
VIEWS = [(0.36, 0, 0.93), (-0.58, 0, 0.82), (0.75, 0, 0.66),
         (0.97, 0, 0.24), (0.4, 0.37, 0.84), (0.4, -0.37, 0.84),
         (0.36, 0, -0.93)]
STRIPE_ANGLES = range(0, 180, 15)


def unit(array, axis=-1):
    return array / np.linalg.norm(array, axis=axis, keepdims=True)


def coordinate_tuple(item):
    return np.array(item.Coord())


@dataclass
class Patch:
    """Expose each cap patch from its attachment (0) to its perimeter (1)."""

    cad: object
    extent: float
    offset: float = 0
    direction: float = 1

    def Value(self, u, v):
        return self.cad.Value(u, self.offset + self.direction*v)

    def D1(self, u, v, point, du, dv):
        self.cad.D1(u, self.offset + self.direction*v, point, du, dv)
        dv.Multiply(self.direction)

    def D2(self, u, v, point, du, dv, uu, vv, uv):
        self.cad.D2(u, self.offset + self.direction*v, point, du, dv, uu, vv, uv)
        dv.Multiply(self.direction)
        uv.Multiply(self.direction)


def jet(surf, u, v):
    point = gp_Pnt()
    derivatives = [gp_Vec() for _ in range(5)]
    surf.D2(float(u), float(v), point, *derivatives)
    du, dv, uu, vv, uv = map(coordinate_tuple, derivatives)
    tangent = np.column_stack((du, dv))
    normal = unit(np.cross(du, dv))
    inverse = np.linalg.pinv(tangent)
    form = np.array([[normal@uu, normal@uv], [normal@uv, normal@vv]])
    normal_derivative = -inverse.T @ form @ inverse
    return coordinate_tuple(point), normal, normal_derivative


def classify(solid):
    surfaces = {}
    for face in solid.Faces():
        surface = BRep_Tool.Surface_s(face.wrapped)
        _, _, low, high = BRepTools.UVBounds_s(face.wrapped)
        midpoint = coordinate_tuple(surface.Value(.5, (low+high)/2))
        mirror = int(np.sign(midpoint[1]))
        span = abs(midpoint[1])
        stage = "cap" if span > CUT else "skin"
        middle = section(span, np.array([.5]), 0)[0, 2]
        side = "upper" if midpoint[2] > middle else "lower"
        if stage == "skin":
            patch = Patch(surface, CUT)
        elif side == "lower":
            patch = Patch(surface, 1., low, 1.)
        else:
            patch = Patch(surface, 1., high, -1.)
        surfaces[mirror, side, stage] = patch
    if len(surfaces) != 8:
        raise RuntimeError("Expected a skin and cap on each side of each wing half")
    return surfaces


def camera(view):
    viewing = unit(np.array(view))
    horizontal = unit(np.array([0., 1., 0.]) - viewing[1]*viewing)
    vertical = np.cross(viewing, horizontal)
    return viewing, horizontal, vertical


def stripe_axis(view, degrees, mirror):
    viewing, horizontal, vertical = camera(view)
    a = np.deg2rad(degrees)
    return viewing, np.cos(a)*mirror*horizontal + np.sin(a)*vertical


def phase(normal, view, axis):
    return 8*(2*(normal@view)*(normal@axis) - view@axis)


def phase_gradient(normal, derivative, view, axis):
    return 16*((normal@view)[..., None]*(axis@derivative)
               + (normal@axis)[..., None]*(view@derivative))


def mismatch_record(values, parameters, **details):
    index = int(np.argmax(values))
    return dict(max=float(values[index]), u=float(parameters[index]), **details)


def seam_metrics(a, va, b, vb, *, zebra, mirror):
    samples = [np.array(list(item)) for item in zip(*(jet(a, u, va) for u in CHORD_PARAMS))]
    second = [np.array(list(item)) for item in zip(*(jet(b, u, vb) for u in CHORD_PARAMS))]
    pa, na, da = samples
    pb, nb, db = second
    align = 1 if np.einsum("ij,ij->i", na, nb).mean() > 0 else -1
    nb *= align
    db *= align
    position = np.linalg.norm(pa-pb, axis=1)
    angles = np.rad2deg(np.arctan2(np.linalg.norm(np.cross(na, nb), axis=1),
                                 np.einsum("ij,ij->i", na, nb)))
    # Compare curvature in the shared tangent frame, represented here as
    # symmetric operators on XYZ projected onto the reference tangent plane.
    project = np.eye(3) - na[:, :, None]*na[:, None, :]
    curvature = np.linalg.norm(project@(da-db)@project, axis=(1, 2))
    result = {
        "position_mismatch_mm": mismatch_record(position, CHORD_PARAMS),
        "normal_angle_deg": mismatch_record(angles, CHORD_PARAMS),
        "curvature_tensor_mismatch_per_mm": mismatch_record(curvature, CHORD_PARAMS),
    }
    if zebra:
        max_phase, max_gradient, max_angle = ({"max": -1.} for _ in range(3))
        crossings = 0
        for view_tuple in VIEWS:
            for angle in STRIPE_ANGLES:
                view, axis = stripe_axis(view_tuple, angle, mirror)
                phase_a, phase_b = phase(na, view, axis), phase(nb, view, axis)
                grad_a = phase_gradient(na, da, view, axis)
                grad_b = phase_gradient(nb, db, view, axis)
                context = {"view_xyz": view_tuple, "angle_deg": angle}
                delta = abs(phase_a-phase_b)
                if delta.max() > max_phase["max"]:
                    max_phase = mismatch_record(delta, CHORD_PARAMS, **context)
                delta = np.linalg.norm(grad_a-grad_b, axis=1)
                if delta.max() > max_gradient["max"]:
                    max_gradient = mismatch_record(delta, CHORD_PARAMS, **context)
                steps = np.floor(2*phase_a).astype(int)
                for i in np.flatnonzero(np.diff(steps)):
                    for transition in range(min(steps[i:i+2])+1, max(steps[i:i+2])+1):
                        t = brentq(
                            lambda u: phase(jet(a, u, va)[1], view, axis)-transition/2,
                            CHORD_PARAMS[i], CHORD_PARAMS[i+1],
                        )
                        _, first_n, first_d = jet(a, t, va)
                        _, next_n, next_d = jet(b, t, vb)
                        ga = phase_gradient(first_n, first_d, view, axis)
                        gb = phase_gradient(next_n, next_d, view, axis)
                        if min(np.linalg.norm(ga), np.linalg.norm(gb)) < 1e-10:
                            continue
                        crossings += 1
                        dga = unit(np.cross(first_n, ga))
                        dgb = unit(np.cross(next_n, gb))
                        kink = np.rad2deg(np.arctan2(np.linalg.norm(np.cross(dga, dgb)),
                                                    abs(dga@dgb)))
                        if kink > max_angle["max"]:
                            max_angle = {"max": float(kink), "u": t, "phase_cycles": transition/2,
                                         "xyz_mm": coordinate_tuple(a.Value(t, va)).tolist(), **context}
        result.update({
            "stripe_phase_mismatch_cycles": max_phase,
            "stripe_gradient_mismatch_cycles_per_mm": max_gradient,
            "stripe_direction_at_light_dark_crossing_deg": max_angle,
            "stripe_crossing_count": crossings,
        })
    for metric in result.values():
        if isinstance(metric, dict) and "xyz_mm" not in metric:
            metric["xyz_mm"] = coordinate_tuple(a.Value(metric["u"], va)).tolist()
    return result


def guide_curvature(surface, u, v):
    p, d1, d2, uu, vv, uv = gp_Pnt(), *(gp_Vec() for _ in range(5))
    surface.D2(float(u), float(v), p, d1, d2, uu, vv, uv)
    t, a = coordinate_tuple(d2), coordinate_tuple(vv)
    bend = np.cross(t, a)
    return np.copysign(np.linalg.norm(bend), bend[0]) / np.linalg.norm(t)**3


def surface_picture(surfaces):
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), constrained_layout=True)
    for ax, view, degrees in zip(axes, [(0.36, 0, 0.93), (0.75, 0, 0.66)], [35, 70]):
        eye, horizontal, vertical = camera(view)
        quad_positions, colors = [], []
        for (mirror, _, stage), s in surfaces.items():
            extent = s.extent
            axis = stripe_axis(view, degrees, mirror)[1]
            # Finer spacing in chord and near the cap rim. Evaluate CAD
            # normals at cell centers rather than use triangulation normals.
            u = np.linspace(0, 1, 193)
            v = np.linspace(0, extent, 130 if stage == "cap" else 91)
            grid = np.array([[coordinate_tuple(s.Value(float(q), float(t))) for t in v] for q in u])
            for i, q in enumerate((u[:-1]+u[1:])/2):
                for j, t in enumerate((v[:-1]+v[1:])/2):
                    p, d1, d2 = gp_Pnt(), gp_Vec(), gp_Vec()
                    s.D1(float(q), float(t), p, d1, d2)
                    normal = unit(coordinate_tuple(d1.Crossed(d2)))
                    intensity = 0.18 + 0.68*(np.sin(2*np.pi*phase(normal, eye, axis)) >= 0)
                    quad_positions.append(grid[[i, i+1, i+1, i], [j, j, j+1, j+1]])
                    colors.append([intensity*.86, intensity*.96, intensity, 1])
        quads = np.array(quad_positions)
        order = np.argsort(quads.mean(axis=1)@eye)
        polygons = np.stack((quads@horizontal, quads@vertical), axis=-1)
        ax.add_collection(PolyCollection(polygons[order], facecolors=np.array(colors)[order],
                                         edgecolors="face", linewidths=.12, antialiased=False))
        ax.autoscale()
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"Whole wing · {degrees}° mirrored stripes")
    fig.savefig(HERE / "zebra.png", dpi=170, facecolor="#f4f6fa")
    plt.close(fig)


def trace_picture(surfaces):
    fig, axes = plt.subplots(2, 2, figsize=(12, 4.2), constrained_layout=True)
    for ax, x in zip(axes.flat, [110, 135, 165, 187]):
        for side in ("upper", "lower"):
            points = []
            for stage, params in [("skin", np.linspace(CUT-20, CUT, 150)),
                                  ("cap", np.linspace(0, 1, 450)[1:])]:
                surface = surfaces[1, side, stage]
                for v in params:
                    u = brentq(lambda q: surface.Value(q, v).X()-x, 0., 1.)
                    points.append(coordinate_tuple(surface.Value(u, v)))
            y, z = np.asarray(points)[:, 1:].T
            ax.plot(y, z, color="#11516e")
        ax.axvline(CUT, ls="--", color=".55", lw=.7)
        ax.set(xlim=(420, 502), ylim=(32, 53), xlabel="Projected span (mm)", ylabel="Height (mm)",
               title=f"Vertical section at {x} mm aft")
        ax.grid(alpha=.15)
        ax.set_aspect("equal")
    fig.savefig(HERE / "sections.png", dpi=160, facecolor="#f4f6fa")
    plt.close(fig)


def curvature_picture(surfaces):
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    for ax, fraction in zip(axes.flat, [.1, .25, .4, .6, .85, .97]):
        u = 2*np.arcsin(np.sqrt(fraction))/np.pi
        for side, color in [("upper", "#11516e"), ("lower", "#a84522")]:
            for stage, samples in [("skin", np.linspace(CUT-10, CUT, 90)),
                                   ("cap", np.linspace(0, 1, 300))]:
                surface = surfaces[1, side, stage]
                y = [surface.Value(u, v).Y() for v in samples]
                curvature = [guide_curvature(surface, u, v) for v in samples]
                ax.plot(y, curvature, label=side if stage == "cap" else None, color=color)
        ax.axvline(CUT, ls="--", color=".55", lw=.7)
        ax.set(xlabel="Projected span (mm)", ylabel="Signed guide curvature (1/mm)",
               title=f"{fraction:.0%} chord")
        ax.set_yscale("symlog", linthresh=1e-5)
        ax.grid(alpha=.15)
    axes[0, 0].legend()
    fig.savefig(HERE / "guide-curvature.png", dpi=160, facecolor="#f4f6fa")
    plt.close(fig)


def main():
    solid = read_solid(STEP)
    surfaces = classify(solid)
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(solid.wrapped, box, False, False)
    _, ymin, _, _, ymax, _ = box.Get()
    tip_extrema = {}
    for mirror in (1, -1):
        rim = surfaces[mirror, "upper", "cap"]
        tip = minimize_scalar(lambda u: -mirror*rim.Value(u, 1).Y(),
                              bounds=(.45, .55), method="bounded",
                              options={"xatol": 1e-13})
        tip_extrema[str(mirror)] = float(-mirror*tip.fun)
    skin = surfaces[1, "upper", "skin"]
    root_ends = np.array([coordinate_tuple(skin.Value(u, 0)) for u in (0, 1)])
    cut_ends = np.array([coordinate_tuple(skin.Value(u, CUT)) for u in (0, 1)])
    chord_vector = cut_ends[1] - cut_ends[0]
    metrics = {
        "step_sha256": hashlib.sha256(STEP.read_bytes()).hexdigest(),
        "closed_valid_solids": 1,
        "face_count": len(solid.Faces()),
        "volume_mm3": solid.Volume(),
        "non_tessellated_kernel_bounding_span_mm": ymax-ymin,
        "rim_extrema_y_mm": tip_extrema,
        "projected_span_mm": tip_extrema["1"]-tip_extrema["-1"],
        "root_chord_mm": float(np.linalg.norm(root_ends[1]-root_ends[0])),
        "cap_station_mm": float(cut_ends[0, 1]),
        "cap_station_chord_mm": float(np.linalg.norm(chord_vector)),
        "cap_station_projected_chord_mm": float(chord_vector[0]),
        "cap_station_le_sweep_mm": float(cut_ends[0, 0]-root_ends[0, 0]),
        "leading_reference_dihedral_deg": float(np.rad2deg(np.arctan2(cut_ends[0, 2], CUT))),
        "cap_station_washout_deg": float(np.rad2deg(np.arctan2(chord_vector[2], chord_vector[0]))),
        "sampled_u": [float(CHORD_PARAMS[0]), float(CHORD_PARAMS[-1]), len(CHORD_PARAMS)],
        "camera_count": len(VIEWS),
        "stripe_angles_per_camera": len(STRIPE_ANGLES),
        "parent_analytic_skin_deviation_mm": {},
        "attachment_seams": {},
        "perimeter_seams": {},
    }
    geometric_noses = {}
    for station in (0, CUT):
        geometric_noses[str(station)] = min(
            minimize_scalar(lambda u: surfaces[1, side, "skin"].Value(u, station).X(),
                            bounds=(0, .05), method="bounded", options={"xatol": 1e-13}).fun
            for side in ("upper", "lower")
        )
    metrics["geometric_le_sweep_mm"] = geometric_noses[str(CUT)]-geometric_noses["0"]
    for sign, side in [(1, "upper"), (-1, "lower")]:
        ref_u = np.linspace(0, 1, 401)
        surf = surfaces[1, side, "skin"]
        metrics["parent_analytic_skin_deviation_mm"][side] = max(
            float(np.linalg.norm(
                np.array([coordinate_tuple(surf.Value(float(u), float(y))) for u in ref_u])
                - section(y, ref_u, sign), axis=1).max())
            for y in np.linspace(0, CUT, 21)
        )
    for mirror in (1, -1):
        label = "right" if mirror > 0 else "left"
        for side in ("upper", "lower"):
            a, b = [surfaces[mirror, side, stage] for stage in ("skin", "cap")]
            metrics["attachment_seams"][f"{label}_{side}"] = seam_metrics(
                a, CUT, b, 0, zebra=True, mirror=mirror,
            )
        metrics["perimeter_seams"][label] = seam_metrics(
            surfaces[mirror, "upper", "cap"], 1., surfaces[mirror, "lower", "cap"], 1.,
            zebra=False, mirror=mirror,
        )
    metrics["right_lip_wrap_guides"] = []
    for fraction in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.75, 0.85, 0.9, 0.95, 0.97, 0.99]:
        u = 2*np.arcsin(np.sqrt(fraction))/np.pi
        upper, lower = [surfaces[1, side, "cap"] for side in ("upper", "lower")]
        metrics["right_lip_wrap_guides"].append({
            "chord_fraction": fraction,
            "xyz_mm": coordinate_tuple(upper.Value(u, 1)).tolist(),
            "upper_radius_mm": 1/abs(guide_curvature(upper, u, 1)),
            "lower_radius_mm": 1/abs(guide_curvature(lower, u, 1)),
        })
    (HERE / "checks.json").write_text(json.dumps(metrics, indent=2)+"\n", encoding="utf-8")
    trace_picture(surfaces)
    curvature_picture(surfaces)
    surface_picture(surfaces)


if __name__ == "__main__":
    main()
