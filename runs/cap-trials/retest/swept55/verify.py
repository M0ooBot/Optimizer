"""Verify the exported wing with CAD derivatives and swept physical sections."""
from hashlib import sha256
import json
from math import pi
from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.gp import gp_Pnt, gp_Vec
from scipy.optimize import brentq, minimize_scalar

from build import SWEEP, TAPER, thickness

HERE = Path(__file__).resolve().parent
STEP = HERE/"wing.step"


def unit(a):
    return a/np.linalg.norm(a)


def angle(a, b):
    return np.degrees(np.arctan2(np.linalg.norm(np.cross(a,b)),np.dot(a,b)))


def derivatives(surface, q, v):
    xyz = gp_Pnt()
    vectors = [gp_Vec() for _ in range(5)]
    surface.D2(q,v,xyz,*vectors)
    return [np.array(p.Coord()) for p in [xyz,*vectors]]


def frame(surface, q, v, side, skin):
    p, u, v, uu, vv, uv = derivatives(surface,q,v)
    w = np.cross(u,v)
    multiplier = side*skin
    normal = multiplier*unit(w)
    normal_jac = multiplier*(np.eye(3)-np.outer(normal,normal)) @ np.column_stack((
        np.cross(uu,v)+np.cross(u,uv),
        np.cross(uv,v)+np.cross(u,vv),
    ))/np.linalg.norm(w)
    inverse = np.linalg.pinv(np.column_stack((u,v)))
    return p,normal,normal_jac@inverse


def phase(normal, normal_derivative, view, axis):
    reflect = 2*np.dot(normal,view)*normal-view
    gradient = 6*(np.dot(normal,view)*axis+np.dot(normal,axis)*view)@normal_derivative
    return 3*np.dot(reflect,axis), gradient


def maximum(record, key, value, **where):
    if key not in record or value > record[key]["value"]:
        record[key] = {"value":float(value),**where}


def verify_joins(faces):
    metrics = {}
    zebra = {}
    views = [unit([.28, .34, .90]),unit([-.48,-.15,.86]),unit([.15,-.60,.79])]
    angles = np.radians([0,35,70,105,140])
    qs = np.linspace(.002,.998,501)
    crossings = 0
    for side in (-1,1):
        for skin in (-1,1):
            parent, cap = faces[side,skin,False], faces[side,skin,True]

            def at(q):
                return frame(parent,q,1,side,skin),frame(cap,q,0,side,skin)

            pairs = [at(q) for q in qs]
            for q, (a,b) in zip(qs,pairs):
                where = dict(q=float(q),side=side,skin=skin)
                maximum(metrics,"gap_mm",np.linalg.norm(a[0]-b[0]),**where)
                maximum(metrics,"tangent_plane_deg",angle(a[1],b[1]),**where)
                maximum(metrics,"normal_derivative_per_mm",np.linalg.norm(a[2]-b[2],ord=2),**where)
                if skin == 1:
                    upper = frame(cap,q,1,side,1)
                    lower = frame(faces[side,-1,True],q,1,side,-1)
                    # The upper and lower parameter directions reverse at the lip.
                    maximum(metrics,"lip_gap_mm",np.linalg.norm(upper[0]-lower[0]),**where)
                    maximum(metrics,"lip_tangent_plane_deg",angle(upper[1],lower[1]),**where)
                    maximum(metrics,"lip_normal_derivative_per_mm",
                            np.linalg.norm(upper[2]-lower[2],ord=2),**where)

            for index, view0 in enumerate(views):
                view = view0*np.array([1,side,skin])
                right = unit(np.cross(view,[0,0,1]))
                up = np.cross(right,view)
                for theta in angles:
                    axis = np.cos(theta)*right+np.sin(theta)*up
                    stripes = []
                    for q,(a,b) in zip(qs,pairs):
                        phase_a,grad_a = phase(a[1],a[2],view,axis)
                        phase_b,grad_b = phase(b[1],b[2],view,axis)
                        stripes.append(phase_a)
                        where = dict(q=float(q),side=side,skin=skin,view=index,
                                     stripe_deg=float(np.degrees(theta)))
                        maximum(zebra,"phase_cycles",abs(phase_a-phase_b),**where)
                        maximum(zebra,"gradient_cycles_per_mm",np.linalg.norm(grad_a-grad_b),**where)
                    for i,(a,b) in enumerate(zip(stripes,stripes[1:])):
                        for k in range(-6,6):
                            level = .25+.5*k
                            if (a-level)*(b-level) >= 0:
                                continue
                            def crossing(q):
                                p = frame(parent,q,1,side,skin)
                                return phase(p[1],p[2],view,axis)[0]-level
                            q = brentq(crossing,qs[i],qs[i+1])
                            fa,fb = at(q)
                            ga = phase(fa[1],fa[2],view,axis)[1]
                            gb = phase(fb[1],fb[2],view,axis)[1]
                            # In the common tangent plane the perpendicular to
                            # a stripe gradient is its stripe tangent direction.
                            ta = np.cross(fa[1],ga)
                            tb = np.cross(fb[1],gb)
                            crossings += 1
                            maximum(zebra,"crossing_direction_deg",angle(ta,tb),q=float(q),
                                    side=side,skin=skin,view=index,
                                    stripe_deg=float(np.degrees(theta)),level_cycles=level)
    return metrics, {"worst":zebra,"crossings":crossings,"views":len(views),
                     "base_view_directions_upper_positive_y":[v.tolist() for v in views],
                     "stripe_angles_deg":np.degrees(angles).tolist(),
                     "samples_per_join":len(qs),"q_domain":[float(qs[0]),float(qs[-1])]}


def curvature(surface,q,v,side,slope,skin):
    p, _, dv, _, dvv, _ = derivatives(surface,q,v)
    span = unit(np.array([slope,side,0.]))
    speed = np.array([np.dot(span,dv),dv[2]])
    accel = np.array([np.dot(span,dvv),dvv[2]])
    value = (speed[0]*accel[1]-speed[1]*accel[0])/np.linalg.norm(speed)**3
    return skin*value,p


def verify_sections(faces):
    sections = []
    qs = [.005,.015,.03,.05,.075,.1,.15,.2,.25,.3,.35,.4,.45,.5,
          .55,.6,.65,.7,.75,.8,.85,.9,.925,.95,.97,.985,.995]
    vs = np.linspace(0,1,2001)
    for q in qs:
        u = 3*q*q-2*q*q*q
        slope = SWEEP-TAPER*u
        worst_plane = 0.
        worst_pair = 0.
        reversals = 0
        smallest = np.inf
        at_v = 0.
        at_half = None
        at_point = None
        lip_radii = []
        for side in (-1,1):
            for skin in (-1,1):
                surface = faces[side,skin,True]
                mate = faces[side,-skin,True]

                def point(v):
                    return curvature(surface,q,v,side,slope,skin)

                values = np.array([point(v)[0] for v in vs])
                reversals += int(np.count_nonzero(values > 1e-9))
                index = np.argmin(values)
                v = vs[index]
                if 0 < index < len(vs)-1:
                    v = minimize_scalar(lambda t: point(t)[0],
                                        bounds=(vs[index-1],vs[index+1]),method="bounded").x
                k, p = point(v)
                radius = -1/k
                if radius < smallest:
                    smallest,at_v,at_half,at_point = radius,v,[side,skin],p.tolist()
                lip_radii.append(-1/values[-1])
                for v in vs[::20]:
                    _,p = point(v)
                    worst_plane = max(worst_plane,abs(p[0]-slope*side*p[1]-210*u))
                    opposite = np.array(mate.Value(q,v).Coord())
                    worst_pair = max(worst_pair,np.linalg.norm(p-opposite*[1,1,-1]))
        sections.append({"q":q,"chord_fraction":u,"plane_x_minus_slope_abs_y_mm":210*u,
                         "plane_slope":slope,"max_plane_error_mm":worst_plane,
                         "max_top_bottom_mirror_error_mm":worst_pair,
                         "positive_curvature_samples":reversals,
                         "min_radius_mm":float(smallest),"min_radius_v":float(at_v),
                         "min_radius_half_and_skin":at_half,"min_radius_xyz_mm":at_point,
                         "lip_radius_mm_range":[min(lip_radii),max(lip_radii)]})
    return sections


def run():
    solid = cq.importers.importStep(str(STEP)).val()
    assert solid.ShapeType()=="Solid" and solid.isValid() and solid.Volume() > 0
    assert len(solid.Shells()) == 1 and BRep_Tool.IsClosed_s(solid.Shells()[0].wrapped)
    faces = {}
    for face in solid.Faces():
        center = face.Center()
        side = 1 if center.y > 0 else -1
        skin = 1 if center.z > 0 else -1
        faces[side,skin,abs(center.y)>390] = BRep_Tool.Surface_s(face.wrapped)
    assert len(faces) == 8
    dimensions = {}
    for (side,skin,cap),surface in faces.items():
        for q in np.linspace(0,1,301):
            u = 3*q*q-2*q*q*q
            if cap:
                s = 390+240*q*(1-q)
                expected = [SWEEP*s+(210-TAPER*s)*u,side*s,0.]
                actual = np.array(surface.Value(q,1).Coord())
                maximum(dimensions,"locked_lip_error_mm",np.linalg.norm(actual-expected),
                        side=side,skin=skin,q=float(q))
            else:
                for v in [0.,.15,.4,.65,.85,1.]:
                    s = 390*v
                    chord = 210-TAPER*s
                    expected = [SWEEP*s+chord*u,side*s,skin*chord*thickness(q)]
                    actual = np.array(surface.Value(q,v).Coord())
                    maximum(dimensions,"parent_naca_error_mm",np.linalg.norm(actual-expected),
                            side=side,skin=skin,q=float(q),v=v)
    lower = faces[-1,1,True].Value(.5,1)
    upper = faces[1,1,True].Value(.5,1)
    span = upper.Y()-lower.Y()
    assert abs(span-900) < 1e-8
    assert all(item["value"] < 1e-6 for item in dimensions.values())
    metrics,zebra = verify_joins(faces)
    sections = verify_sections(faces)
    assert metrics["gap_mm"]["value"] < 1e-7
    assert metrics["tangent_plane_deg"]["value"] < 1e-4
    assert metrics["lip_tangent_plane_deg"]["value"] < 1e-4
    assert zebra["worst"]["crossing_direction_deg"]["value"] < 1e-3
    assert all(s["positive_curvature_samples"] == 0 and s["min_radius_v"] > .9999
               and s["max_plane_error_mm"] < 1e-7 for s in sections)
    out = {"step_sha256":sha256(STEP.read_bytes()).hexdigest(),"solid_valid":True,
           "shell_closed":True,"face_count":len(faces),"volume_mm3":solid.Volume(),
           "span_mm":span,"geometry":dimensions,"joins_and_lips":metrics,
           "zebra":zebra,"swept_sections":sections,"wrap_samples_per_skin":2001}
    (HERE/"readback.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")


if __name__ == "__main__":
    run()
