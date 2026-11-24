"""Measure joins, reflection stripes and oblique spanwise sections on the STEP."""

from collections import namedtuple
from hashlib import sha256
import json
from pathlib import Path

import cadquery as cq
from OCP.Bnd import Bnd_Box
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.BRepCheck import BRepCheck_Analyzer, BRepCheck_NoError, BRepCheck_Shell
from OCP.GeomAbs import GeomAbs_BSplineSurface
from OCP.gp import gp_Pnt, gp_Vec
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from build_fin import CUT, HEIGHT, HERE, OUTPUT, THICKNESS_RATIO, center, radius

Jet = namedtuple("Jet", ["point", "normal", "frame", "normal_derivative"])


def unit(vector):
    return vector / np.linalg.norm(vector)


def jet(surface, u, v, side):
    point = gp_Pnt()
    vectors = [gp_Vec() for _ in range(5)]
    surface.D2(u, v, point, *vectors)
    du, dv, duu, dvv, duv = [np.array(vector.Coord()) for vector in vectors]
    # U runs from aft to forward on both skins and V points upward.
    normal = side * unit(np.cross(du, dv))
    frame = np.column_stack((unit(du), np.cross(normal, unit(du))))
    j = frame.T @ np.column_stack((du, dv))
    inv = np.linalg.inv(j)
    normal_curvature = np.array([[normal @ duu, normal @ duv], [normal @ duv, normal @ dvv]])
    return Jet(np.array(point.Coord()), normal, frame, -frame @ inv.T @ normal_curvature @ inv)


def at(face, u, v):
    return np.array(face.Value(u, v).Coord())


def patches(path):
    imported = cq.importers.importStep(str(path)).val()
    solids = imported.Solids()
    if len(solids) != 1:
        raise ValueError("Expected one solid after reimport")
    solid = solids[0]
    strict = BRepCheck_Analyzer(solid.wrapped, True, False, True).IsValid()
    closed = all(BRepCheck_Shell(sh.wrapped).Closed() == BRepCheck_NoError for sh in solid.Shells())
    if not strict or not closed or len(solid.Faces()) != 5:
        raise ValueError("Unexpected STEP topology")
    result = {}
    for face in solid.Faces():
        adaptor = BRepAdaptor_Surface(face.wrapped)
        if adaptor.GetType() != GeomAbs_BSplineSurface:
            continue
        region = "cap" if at(adaptor, 0.5, 0)[2] > CUT/2 else "parent"
        sign = 1 if at(adaptor, 0.5, 0.5)[1] > 0 else -1
        result[region, sign] = adaptor
    if len(result) != 4:
        raise ValueError("Expected four B-spline skin faces")
    box = Bnd_Box()
    BRepBndLib.AddOptimal_s(solid.wrapped, box, False, False)
    bbox = box.Get()
    body = result["parent", 1]
    root_le, root_te = (at(body, u, 0)[0] for u in (1, 0))
    cut_le, cut_te = (at(body, u, 1)[0] for u in (1, 0))
    props = GProp_GProps()
    volume_error = BRepGProp.VolumeProperties_s(solid.wrapped, props, 1e-9, True)
    design_volume = quad(lambda z: np.pi*THICKNESS_RATIO*radius(z)**2, 0, CUT)[0]
    design_volume += quad(lambda z: np.pi*THICKNESS_RATIO*radius(z)**2
                          * (1 - ((z-CUT)/(HEIGHT-CUT))**3), CUT, HEIGHT)[0]
    summary = {
        "step_sha256": sha256(Path(path).read_bytes()).hexdigest(),
        "solids": len(solids), "faces": len(solid.Faces()), "shells": len(solid.Shells()),
        "closed": bool(closed), "strict_occt_valid": bool(strict), "volume_mm3": props.Mass(),
        "adaptive_integration_error_estimate": volume_error,
        "analytic_design_volume_mm3": design_volume,
        "height_mm": float(bbox[5] - bbox[2]), "root_chord_mm": root_te-root_le,
        "cut_height_mm": at(body, 0, 1)[2], "cut_chord_mm": cut_te-cut_le,
        "leading_edge_sweep_at_cut_mm": cut_le-root_le,
    }
    actual = [summary[key] for key in (
        "height_mm", "root_chord_mm", "cut_height_mm", "cut_chord_mm",
        "leading_edge_sweep_at_cut_mm",
    )]
    if not np.allclose(actual, [HEIGHT, 160, CUT, 85, 65], rtol=0, atol=1e-6):
        raise ValueError(f"Locked dimensions changed: {actual}")
    return result, summary


def stripe_views():
    cameras = [(az, el) for az in range(0, 360, 45) for el in (-45, 0, 45)]
    angles = [0, 30, 51, 75, 90, 120, 150]
    views, eyes, axes = [], [], []
    for az, el in cameras:
        a, e = np.radians([az, el])
        eye = np.array([np.cos(e)*np.sin(a), np.cos(e)*np.cos(a), np.sin(e)])
        right = np.array([np.cos(a), -np.sin(a), 0])
        up = np.cross(right, eye)
        for angle in angles:
            theta = np.radians(angle)
            views.append([az, el, angle])
            eyes.append(eye)
            axes.append(np.cos(theta)*right + np.sin(theta)*up)
    return views, np.array(eyes), np.array(axes)


def zebra(j, eyes, axes, frame):
    ne, na = eyes @ j.normal, axes @ j.normal
    cycle = 3 * (2 * ne * na - np.einsum("ij,ij->i", eyes, axes))
    by_normal = 6 * (na[:, None] * eyes + ne[:, None] * axes)
    gradient = by_normal @ j.normal_derivative @ (j.frame.T @ frame)
    return cycle, gradient


def angle(a, b):
    return np.degrees(np.arctan2(np.linalg.norm(np.cross(a, b), axis=-1), np.sum(a*b, axis=-1)))


def peak(record, name, value, location):
    if name not in record or value > record[name]["value"]:
        record[name] = {"value": float(value), "location": location}


def compare(face_1, face_2, coordinates, label):
    views, eyes, axes = stripe_views()
    record = {"sample_count": len(coordinates), "view_count": len(views), "crossings": 0}
    phases = []
    pairs = []
    for parameters, first, second in coordinates:
        a = jet(face_1, *first)
        b = jet(face_2, *second)
        location = {"xyz_mm": a.point.tolist(), "parameter": float(parameters)}
        peak(record, "gap_mm", np.linalg.norm(a.point - b.point), location)
        peak(record, "normal_angle_deg", angle(a.normal, b.normal), location)
        peak(record, "normal_derivative_difference_per_mm",
             np.linalg.norm(a.normal_derivative - b.normal_derivative @ (b.frame.T @ a.frame), 2),
             location)
        phase_a, gradient_a = zebra(a, eyes, axes, a.frame)
        phase_b, gradient_b = zebra(b, eyes, axes, a.frame)
        for name, values in [("stripe_phase_difference_cycles", np.abs(phase_a-phase_b)),
                             ("stripe_gradient_difference_cycles_per_mm",
                              np.linalg.norm(gradient_a-gradient_b, axis=1))]:
            i = int(np.argmax(values))
            peak(record, name, values[i], dict(location, view_az_el_stripe_deg=views[i]))
        phases.append(phase_a)
        pairs.append((first, second))
    for index, (phase_1, phase_2) in enumerate(zip(phases, phases[1:])):
        candidates = np.flatnonzero(np.floor(2*phase_1) != np.floor(2*phase_2))
        for camera in candidates:
            lo, hi = sorted((2*phase_1[camera], 2*phase_2[camera]))
            eye, axis = eyes[[camera]], axes[[camera]]
            f1, s1 = pairs[index]
            f2, s2 = pairs[index+1]

            def crossing(t):
                first = np.array(f1) + t*(np.array(f2) - f1)
                return jet(face_1, *first)

            for edge in range(int(np.floor(lo))+1, int(np.floor(hi))+1):
                t = brentq(lambda t: zebra(crossing(t), eye, axis, np.eye(3))[0][0] - edge/2,
                           0, 1)
                a = crossing(t)
                second = np.array(s1) + t*(np.array(s2) - s1)
                b = jet(face_2, *second)
                _, ga = zebra(a, eye, axis, a.frame)
                _, gb = zebra(b, eye, axis, a.frame)
                if min(np.linalg.norm(ga), np.linalg.norm(gb)) < 1e-8:
                    continue
                determinant = ga[0, 0]*gb[0, 1] - ga[0, 1]*gb[0, 0]
                turn = np.degrees(np.arctan2(abs(determinant), abs(ga[0] @ gb[0])))
                record["crossings"] += 1
                peak(record, "stripe_crossing_turn_deg", turn,
                     {"xyz_mm": a.point.tolist(), "view_az_el_stripe_deg": views[camera]})
    print(f"{label}: {record['crossings']} stripe crossings")
    return record


def section(plus, parent, q):
    side_edge = 1 if q < 0.5 else 0

    def section_x(z):
        return center(z) + (2*q-1)*radius(z)

    def rim_error(v):
        p = at(plus, side_edge, v)
        return p[0] - section_x(p[2])

    lip_v = brentq(rim_error, 0, 1, xtol=1e-14)
    lip = jet(plus, side_edge, lip_v, 1)
    x, _, z = lip.point
    slope = (section_x(HEIGHT) - section_x(CUT)) / (HEIGHT-CUT)
    section_tangent = unit([slope, 0, 1])
    # At the rim the cross-sectional tangent is parallel to the thickness axis.
    thickness_dir = np.array([0, 1, 0])
    normal_curvature = abs(thickness_dir @ lip.normal_derivative @ (lip.frame.T @ thickness_dir))
    lip_radius = abs(lip.normal @ section_tangent) / normal_curvature

    def sample(v, surface):
        z = at(surface, 0, v)[2]
        x = section_x(z)
        u = brentq(lambda u: at(surface, u, v)[0] - x, 0, 1, xtol=1e-14)
        j = jet(surface, u, v, 1)
        plane_normal = unit([1, 0, -slope])
        tangent = unit(np.cross(j.normal, plane_normal))
        curve_normal = np.cross(tangent, plane_normal)
        k = abs(tangent @ j.normal_derivative @ (j.frame.T @ tangent))
        return z, j.point[1], k/abs(j.normal @ curve_normal)

    root = [sample(v, parent) for v in np.linspace(220/CUT, 1, 31)]
    crown = [sample(v, plus) for v in np.linspace(0, lip_v, 181)[:-1]]
    trace = np.array([*root, *crown[1:], (z, 0, 1/lip_radius)])
    return {"cut_chord_fraction_from_le": float(q), "lip_xyz_mm": [x, 0, z],
            "lip_radius_mm": float(lip_radius), "trace_z_y_curvature": trace.tolist()}


def design_checks(faces):
    profile_error = {"parent": 0.0, "cap": 0.0}
    for (region, side), face in faces.items():
        for u in np.linspace(0, 1, 31):
            for v in np.linspace(0, 1, 137):
                x, y, z = at(face, u, v)
                implicit = ((x-center(z))/radius(z))**2 + (y/(THICKNESS_RATIO*radius(z)))**2
                if region == "cap":
                    implicit += ((z-CUT)/(HEIGHT-CUT))**3
                profile_error[region] = max(profile_error[region], abs(1-implicit))

    perimeter_tangents = {}
    for edge, u in [("leading", 1), ("trailing", 0)]:
        tangents = []
        for region, v in [("parent", 1), ("cap", 0)]:
            p, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
            faces[region, 1].D1(u, v, p, du, dv)
            tangents.append(unit(dv.Coord()))
        perimeter_tangents[edge] = angle(*tangents)
    near_apex = []
    for (region, side), face in faces.items():
        if region == "cap":
            for u in np.linspace(0, 1, 9):
                j = jet(face, u, 0.99999, side)
                eigenvalues = np.linalg.eigvalsh(j.frame.T @ j.normal_derivative)
                near_apex.append([*map(float, sorted(1/abs(eigenvalues))),
                                  float(angle(j.normal, np.array([0, 0, 1])))])
    apex = np.array(near_apex)
    return {
        "max_dimensionless_implicit_profile_error": profile_error,
        "perimeter_to_parent_edge_tangent_angle_deg": perimeter_tangents,
        "apex_xyz_mm": at(faces["cap", 1], 0.5, 1).tolist(),
        "near_apex_parameter": 0.99999,
        "near_apex_principal_radii_mm_range": np.array([apex.min(axis=0)[:2],
                                                       apex.max(axis=0)[:2]]).T.tolist(),
        "near_apex_max_normal_angle_to_vertical_deg": float(max(apex[:, 2])),
    }


def check(path):
    faces, summary = patches(path)
    coordinates = np.linspace(0, 1, 501)
    joins = {
        str(side): compare(faces["parent", side], faces["cap", side],
                           [(u, (u, 1, side), (u, 0, side)) for u in coordinates],
                           f"skin {side} attachment")
        for side in (-1, 1)
    }
    # Exclude only the collapsed apex parameter and retain both root endpoints.
    rims = {
        edge: compare(faces["cap", 1], faces["cap", -1],
                      [(v, (u, v, 1), (u, v, -1)) for v in np.linspace(0, 0.999, 501)], edge)
        for edge, u in [("leading", 1), ("trailing", 0)]
    }
    sections = [section(faces["cap", 1], faces["parent", 1], q)
                for q in (0.05, 0.15, 0.25, 0.36, 0.45, 0.55, 0.64, 0.75, 0.85, 0.95)]
    return dict(readback=summary, design=design_checks(faces),
                attachment=joins, perimeter=rims, sections=sections)


if __name__ == "__main__":
    evidence = check(OUTPUT)
    (HERE / "readback.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps(evidence["readback"], indent=2))
