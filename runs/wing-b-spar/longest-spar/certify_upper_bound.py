"""Verify upper bounds for a straight Ø6.4 mm channel with a 1 mm wall.

Uses necessary disk/ellipse support constraints from exact positive-weight
STEP section Bezier control hulls, independent of the candidate-axis search.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.GeomConvert import GeomConvert_BSplineCurveToBezierCurve
from scipy.optimize import brentq, linprog
from scipy.spatial import ConvexHull

BORE_RADIUS = 3.2
WALL = 1.0
SUPPORT_PAD = 1e-6
SLOPE_TO_LENGTH_SCALE = 430.0


def section_constraints(step):
    solid, = cq.importers.importStep(str(step)).solids().vals()
    if not solid.isValid():
        raise ValueError("Expected a single valid half-wing STEP solid")
    skin, = (f for f in solid.Faces() if f.geomType() == "BSPLINE")
    surface = skin._geomAdaptor()
    u0, _, v0, v1 = surface.Bounds()
    root = -surface.Value(u0, v0).X()
    span = -surface.Value(u0, v1).X()
    if abs(root) > SUPPORT_PAD:
        raise ValueError("Expected root at X=0")

    stations = np.unique(np.r_[np.arange(.03, 2, .06), np.arange(2, 400, 3.),
                              np.arange(400, 452.01, .10)])
    all_stations, all_normals, all_support = [], [], []
    directions = np.r_[np.arange(0, 360, 15), np.arange(-15, 15.1, .5),
                       np.arange(165, 195.1, .5)]
    directions = np.deg2rad(np.unique(directions % 360))
    directions = np.column_stack((np.cos(directions), np.sin(directions)))
    min_weight = math.inf
    max_span_error = 0.0
    for station in stations:
        v = brentq(lambda t: -surface.Value(u0, t).X() - station, v0, v1, xtol=1e-13)
        converter = GeomConvert_BSplineCurveToBezierCurve(surface.VIso(v))
        controls, weights = [], []
        for i in range(1, converter.NbArcs()+1):
            arc = converter.Arc(i)
            for j in range(1, arc.NbPoles()+1):
                controls.append(arc.Pole(j).Coord())
                weights.append(arc.Weight(j))
        controls = np.asarray(controls)
        min_weight = min(min_weight, *weights)
        max_span_error = max(max_span_error, np.max(abs(controls[:, 0]+station)))
        hull = ConvexHull(controls[:, 1:])
        normals = np.vstack((directions, hull.equations[:, :2][::5]))
        all_stations.extend([station]*len(normals))
        all_normals.extend(normals)
        all_support.extend((hull.points[hull.vertices] @ normals.T).max(axis=0) + SUPPORT_PAD)
    if min_weight <= 0 or max_span_error >= SUPPORT_PAD:
        raise ValueError("Expected positive Bezier weights on fixed-span sections")
    box = solid.BoundingBox()
    return (np.asarray(all_stations), np.asarray(all_normals), np.asarray(all_support),
            np.array([box.ymin, box.zmin])-SUPPORT_PAD,
            np.array([box.ymax, box.zmax])+SUPPORT_PAD,
            float(min_weight), float(max_span_error), len(stations))


def certificate(section, target, metric, limit=2500):
    stations, normals, support, root_lo, root_hi, *_ = section
    axis_travel = root_hi-root_lo
    least_projection = target if metric == "projected" else math.sqrt(target**2-axis_travel @ axis_travel)
    initial_q = axis_travel / least_projection
    pending = [(-initial_q, initial_q)]
    checked = 0
    exclusions = []
    # Include smaller outer radii close to the rounded clearance cap. An
    # R-slice point whose projection is past the cap by d is at most
    # hypot(R-BORE_RADIUS, d) from the finite cylinder.
    radii = np.r_[4.2, 4.2 - np.geomspace(.0001, .12, 13), 3.9, 3.6, 3.2]
    while pending and checked < limit:
        qlo, qhi = pending.pop()
        checked += 1
        qmid = (qlo+qhi)/2
        qrad = (qhi-qlo)/2
        qabs = np.maximum(abs(qlo), abs(qhi))
        hmax = math.sqrt(1 + qabs @ qabs)
        reach = target if metric == "projected" else target/hmax
        nmid = normals @ qmid
        nrad = abs(normals) @ qrad
        nlower = nmid-nrad
        nupper = nmid+nrad
        outer = np.zeros(len(stations))
        for radius in radii:
            wall_radial = radius-BORE_RADIUS
            left = stations + wall_radial * nlower / np.hypot(1, nlower)
            right_limit = reach + math.sqrt(max(0, WALL**2-wall_radial**2))/hmax
            right = stations + radius*nupper/np.hypot(1, nupper)
            take = (outer == 0) & (left >= SUPPORT_PAD) & (right <= right_limit-SUPPORT_PAD)
            outer[take] = radius
        use = outer > 0
        n = normals[use]
        s = stations[use]
        # The support of the oblique cylinder's plane section along n is
        # R*sqrt(1+(n.q)^2). Its tangent at n.qmid is a global lower bound.
        center = nmid[use]
        tangent_denom = np.hypot(1, center)
        radius = outer[use]
        A = np.column_stack((n, ((s+radius*center/tangent_denom)/SLOPE_TO_LENGTH_SCALE)[:, None]*n))
        B = support[use]-radius/tangent_denom
        lo = np.r_[root_lo, qlo*SLOPE_TO_LENGTH_SCALE]
        hi = np.r_[root_hi, qhi*SLOPE_TO_LENGTH_SCALE]
        # Maximize uniform support slack. LP dual weights provide explicit
        # positive contradictions for each independently excluded slope box.
        augmented = np.column_stack((A,np.ones(len(A))))
        active = np.arange(0, len(A), 83)
        for attempt in range(30):
            res = linprog([0,0,0,0,-1], A_ub=augmented[active], b_ub=B[active],
                          bounds=[*zip(lo,hi), (None,None)], method="highs")
            if not res.success:
                raise RuntimeError(res.message)
            if res.fun > 1e-8:
                break
            violation = augmented @ res.x - B
            if violation.max() < 1e-8 or attempt == 29:
                break
            extra = np.argpartition(violation, -24)[-24:]
            active = np.union1d(active, extra)
        weights = np.zeros(len(A))
        weights[active] = np.maximum(-res.ineqlin.marginals, 0)
        weights /= weights.sum()
        global_coeffs = weights @ A
        bound = float(weights @ B - np.minimum(global_coeffs*lo, global_coeffs*hi).sum())
        if bound < -1e-8:
            exclusions.append({
                "slope_low_yz": qlo.tolist(), "slope_high_yz": qhi.tolist(),
                "weighted_slack_ceiling_mm": bound,
                "support_rows": [
                    {"station_mm": float(s[i]), "normal_yz": n[i].tolist(),
                     "section_support_mm": float(support[use][i]),
                     "outer_radius_mm": float(radius[i]), "weight": float(weights[i])}
                    for i in np.flatnonzero(weights)
                ],
            })
            continue
        if np.max(qrad) < 1e-6:
            return False, checked, exclusions, qmid.tolist(), bound
        j = int(np.argmax(qrad))
        a = qhi.copy()
        b = qlo.copy()
        a[j] = b[j] = qmid[j]
        pending.extend([(qlo,a),(b,qhi)])
    return not pending, checked, exclusions, None, None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", type=Path)
    parser.add_argument("--length", type=float, required=True)
    parser.add_argument("--reach", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    section = section_constraints(args.step)
    result = {
        "input_step": args.step.as_posix(),
        "input_sha256": hashlib.sha256(args.step.read_bytes()).hexdigest(),
        "channel_radius_mm": BORE_RADIUS,
        "required_wall_mm": WALL,
        "support_padding_mm": SUPPORT_PAD,
        "slope_to_length_scale_mm": SLOPE_TO_LENGTH_SCALE,
        "minimum_bezier_weight": section[5],
        "maximum_bezier_span_error_mm": section[6],
        "tested_section_count": section[7],
        "wing_box_low_yz_mm": section[3].tolist(),
        "wing_box_high_yz_mm": section[4].tolist(),
        "notes": [
            "Wing cross-sections and their interiors are contained in the convex hull of positive-weight Bezier control points.",
            "The oblique infinite cylinder ellipse at every section has support R sqrt(1+(n.q)^2), lower-bounded by its tangent. The analytic root and cap guards only admit ellipse extremizers within the required distance of the finite, root-clipped bore.",
            "For 3D length at least L, cap span is at least L/sqrt(1+|q|^2). Each slope box uses the maximum denominator in that box. For projected reach at least S, cap span is at least S.",
            "Each excluded box has nonnegative LP weights summing to one. Minimizing the weighted linear coefficients over the root and slope box bounds gives a negative upper bound on required support slack. Dyadic boxes cover all possible slopes.",
            "Bounds concern the axis cap center and apply independently of print joints and search candidate."
        ],
    }
    for metric, target in (("physical_3d_length", args.length), ("projected", args.reach)):
        legal_excluded, count, boxes, unresolved, unresolved_slack = certificate(section, target, metric)
        result[metric] = {"excluded_at_or_above_mm": target if legal_excluded else None,
                          "tested_target_mm": target, "tested_slope_boxes": count,
                          "unresolved_slope": unresolved, "unresolved_slack_mm": unresolved_slack,
                          "excluded_slope_boxes": boxes}
        print(metric, target, legal_excluded, count, unresolved, unresolved_slack, flush=True)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
