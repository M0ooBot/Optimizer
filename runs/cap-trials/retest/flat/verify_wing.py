"""Measure physical sections and reflection continuity from the exported STEP."""
from hashlib import sha256
import json

import cadquery as cq
import numpy as np
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.GeomAdaptor import GeomAdaptor_Curve
from OCP.gp import gp_Dir, gp_Pln, gp_Pnt, gp_Vec

from build_wing import OUT, d_at, h_at, x_at


def vec(value):
    return np.array(value.Coord())


def unit(value):
    return value / np.linalg.norm(value)


def angle(a, b):
    return float(np.degrees(np.arctan2(np.linalg.norm(np.cross(a, b)), a @ b)))


def differential(surface, u, v, skin, side, lip=False):
    p, a, b, aa, bb, ab = gp_Pnt(), *(gp_Vec() for _ in range(5))
    surface.D2(float(u), float(v), p, a, b, aa, bb, ab)
    a, b, aa, bb, ab = map(vec, (a, b, aa, bb, ab))
    area = np.cross(a, b)
    coordinate = 1 if lip else 2
    target = side if lip else skin
    sign = np.sign(area[coordinate]) * target
    n = sign * unit(area)
    projection = np.eye(3) - np.outer(n, n)
    normal_uv = sign * projection @ np.column_stack((
        np.cross(aa, b) + np.cross(a, ab),
        np.cross(ab, b) + np.cross(a, bb),
    )) / np.linalg.norm(area)
    frame = np.column_stack((a, b))
    return vec(p), n, normal_uv @ np.linalg.pinv(frame)


def compare_geometry(rows):
    return {label: max(({'value': value(a,b), 'u': u, 'side': side, 'skin': skin}
                        for side, skin, u, a, b in rows), key=lambda row: row['value'])
            for label, value in (
                ('gap_mm', lambda a,b: float(np.linalg.norm(a[0]-b[0]))),
                ('tangent_angle_deg', lambda a,b: angle(a[1],b[1])),
                ('normal_derivative_jump_per_mm', lambda a,b: float(np.linalg.norm(a[2]-b[2], ord=2))),
            )}


def reflected(n, dn, view, light):
    # A sinusoidal environment with light/dark crossings at integer pi.
    phase_scale = 8 * np.pi
    phase = phase_scale * (2 * (view @ n) * (light @ n) - light @ view)
    gradient = phase_scale * 2 * ((light @ n)*view + (view @ n)*light) @ dn
    return float(phase), gradient


def maximum(records, key):
    return max(records, key=lambda row: row[key])


def verify_zebra(surfaces, grid):
    settings = [(a, b) for a in ((0,0,1), (.45,.30,1), (-.35,.2,1))
                for b in ((1,0,0), (.7,.7,0), (0,1,0), (.5,.4,.75))]
    joins, rims, stripes, crossings = [], [], [], []
    for side in (-1, 1):
        for u in grid:
            upper = differential(surfaces[side, 1], u, 1, 1, side, lip=True)
            lower = differential(surfaces[side, -1], u, 1, -1, side, lip=True)
            rims.append((side, 0, u, upper, lower))
        for skin in (-1, 1):
            def pair(u):
                parent = differential(surfaces[0, skin], u, (side+1)//2, skin, side)
                cap = differential(surfaces[side, skin], u, 0, skin, side)
                return parent, cap
            pairs = [pair(u) for u in grid]
            joins.extend((side, skin, u, *p) for u,p in zip(grid,pairs))
            for view, light in settings:
                mirror = np.array([1,side,skin])
                v, l = unit(view)*mirror, unit(light)*mirror

                def phase_and_gradient(sample):
                    _, n, dn = sample
                    return reflected(n, dn, v, l)

                def record(u, samples):
                    (phase_a, grad_a), (phase_b, grad_b) = map(phase_and_gradient, samples)
                    row = {'u': u, 'side': side, 'skin': skin, 'view': view, 'light': light,
                           'phase_jump_rad': abs(phase_a-phase_b),
                           'gradient_jump_rad_per_mm': float(np.linalg.norm(grad_a-grad_b))}
                    return phase_a, grad_a, grad_b, row

                phases = []
                for u, samples in zip(grid,pairs):
                    phase, _, _, row = record(u, samples)
                    phases.append(phase)
                    stripes.append(row)
                for k, (a,b) in enumerate(zip(phases, phases[1:])):
                    for band in range(int(np.ceil(min(a,b)/np.pi)), int(np.floor(max(a,b)/np.pi))+1):
                        target = band*np.pi
                        lo, hi = grid[k:k+2]
                        for _ in range(31):
                            mid = (lo+hi)/2
                            value = phase_and_gradient(pair(mid)[0])[0]
                            if (value-target)*(a-target) > 0:
                                lo = mid
                            else:
                                hi = mid
                        u = (lo+hi)/2
                        _, ga, gb, row = record(u, pair(u))
                        if min(np.linalg.norm(ga), np.linalg.norm(gb)) > 1e-9:
                            crossings.append({**row, 'band': band, 'stripe_direction_jump_deg': angle(ga,gb)})
    return {
        'join_geometry': compare_geometry(joins), 'rim_geometry': compare_geometry(rims),
        'u_min_max': [grid[0],grid[-1]], 'grid_samples_per_join': len(grid),
        'undefined_derivative_endpoints': [0,1],
        'view_stripe_angle_combinations_per_join': len(settings),
        'crossing_count': len(crossings),
        'phase_max': maximum(stripes, 'phase_jump_rad'),
        'gradient_max': maximum(stripes, 'gradient_jump_rad_per_mm'),
        'crossing_angle_max': maximum(crossings, 'stripe_direction_jump_deg'),
    }


def verify_sections(solid, surfaces, values):
    results = []
    for u in values:
        plane = gp_Pln(gp_Pnt(x_at(u),0,0), gp_Dir(1,0,0))
        cut = BRepAlgoAPI_Section(solid.wrapped, plane, True)
        assert cut.IsDone()
        edges = cq.Shape.cast(cut.Shape()).Edges()
        assert len(edges) == 6, (u, len(edges))
        for edge in edges:
            if abs(edge.Center().y) < 349:
                continue
            side = int(np.sign(edge.Center().y))
            skin = int(np.sign(edge.Center().z))
            # The OCCT section edges confirm the physical intersections. Measure
            # curvature on the exact STEP surface isocurve, avoiding the section
            # algorithm's polyline approximation and its zero second derivatives.
            curve = GeomAdaptor_Curve(surfaces[side,skin].BSpline().UIso(u))
            a, b = curve.FirstParameter(), curve.LastParameter()
            assert abs(curve.Value((a+b)/2).X() - x_at(u)) < 1e-8
            assert min(np.linalg.norm(vec(curve.Value(a))-np.array(tuple(v)))
                       for v in (edge.startPoint(), edge.endPoint())) < 1e-6
            assert min(np.linalg.norm(vec(curve.Value(b))-np.array(tuple(v)))
                       for v in (edge.startPoint(), edge.endPoint())) < 1e-6
            if abs(curve.Value(a).Y()) > abs(curve.Value(b).Y()):
                a,b = b,a
            direction = np.sign(b-a)
            samples = []
            for t in np.linspace(a,b,1001):
                point, first, second = gp_Pnt(), gp_Vec(), gp_Vec()
                curve.D2(float(t),point,first,second)
                # Trace each end from upper attachment through the lip to the lower attachment.
                kappa = side*skin*direction*(first.Y()*second.Z()-first.Z()*second.Y()) / first.Magnitude()**3
                samples.append((float(kappa), point.Coord()))
            peak = max(range(len(samples)), key=lambda i: abs(samples[i][0]))
            nonzero = [x for x,_ in samples if abs(x)>1e-8]
            results.append({
                'u': u, 'x_mm': x_at(u), 'side': side, 'skin': skin,
                'min_signed_curvature_per_mm': min(x for x,_ in samples),
                'max_signed_curvature_per_mm': max(x for x,_ in samples),
                'reversals_above_1e-8_per_mm': sum(a*b<0 for a,b in zip(nonzero,nonzero[1:])),
                'minimum_radius_at_xyz_mm': samples[peak][1],
                'minimum_radius_mm': 1/abs(samples[peak][0]),
                'lip_radius_mm': 1/abs(samples[-1][0]),
                'lip_to_minimum_radius_ratio': abs(samples[peak][0]/samples[-1][0]),
                'samples_attachment_to_lip': len(samples),
            })
    return results


def main():
    step = OUT / 'wing.step'
    imported = cq.importers.importStep(str(step)).val()
    solids = imported.Solids()
    assert len(solids) == 1
    solid = solids[0]
    assert solid.isValid() and len(solid.Shells()) == 1 and solid.Shells()[0].Closed()
    faces = solid.Faces()
    assert len(faces) == 6
    surfaces = {}
    for face in faces:
        center = face.Center()
        side = 0 if abs(center.y) < 350 else int(np.sign(center.y))
        skin = int(np.sign(center.z))
        surfaces[side,skin] = BRepAdaptor_Surface(face.wrapped)
    assert len(surfaces) == 6

    grid = np.linspace(.005,.995,199)
    rim_error = attachment_error = parent_error = 0
    for u in np.linspace(0,1,201):
        for side in (-1,1):
            for skin in (-1,1):
                cap = surfaces[side,skin]
                rim_error = max(rim_error, np.linalg.norm(vec(cap.Value(u,1)) -
                                                          (x_at(u),side*(350+d_at(u)),0)))
                attachment_error = max(attachment_error, np.linalg.norm(vec(cap.Value(u,0)) -
                                                                        (x_at(u),side*350,skin*h_at(u))))
                for v in np.linspace(0,1,11):
                    parent = surfaces[0,skin]
                    parent_error = max(parent_error, np.linalg.norm(vec(parent.Value(u,v)) -
                                                                    (x_at(u),700*v-350,skin*h_at(u))))
    box = solid.BoundingBox()
    spans = [surfaces[side,1].Value(.5,1).Y() for side in (-1,1)]
    report = {
        'step_sha256': sha256(step.read_bytes()).hexdigest(), 'valid': True,
        'closed_shell': True, 'solid_count': len(solids), 'face_count': len(faces),
        'volume_mm3': solid.Volume(), 'chord_mm': box.xlen,
        'span_mm': max(spans) - min(spans), 'tip_y_mm': spans,
        'parent_max_error_mm': parent_error, 'rim_max_error_mm': rim_error,
        'attachment_max_error_mm': attachment_error,
        'zebra': verify_zebra(surfaces, grid),
        'physical_sections': verify_sections(solid, surfaces, [.025,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,.975]),
    }
    assert abs(report['span_mm'] - 800) < 1e-6 and abs(report['chord_mm'] - 120) < 1e-6
    assert max(parent_error, rim_error, attachment_error) < 1e-6
    (OUT / 'readback.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf8')


if __name__ == '__main__':
    main()
