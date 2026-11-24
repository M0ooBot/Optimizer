"""Optimize a finite straight round spar on redistributed Wing B."""
import json
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np
from scipy.optimize import LinearConstraint, minimize

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from planeopt import spar

BASE = Path(__file__).resolve().parent
WING = Path("runs/wing-b-spar/shape/constrained/redistributed/wing.step")
SCALE = 450.
CUTS = (175.,350.)
MIN_ENGAGEMENT = 20.
RADIUS = 3.2
WALL = 1.


def axis(p):
    origin = np.r_[0., p[:2]]
    vector = np.r_[-1., p[2:4]/SCALE]
    direction = vector/np.linalg.norm(vector)
    return origin, direction


def engagement(p):
    _, direction = axis(p)
    return p[4]-(CUTS[-1]+RADIUS*np.linalg.norm(direction[1:]))/-direction[0]


def clearance(points, p):
    root,direction=axis(p)
    delta=points-root
    axial=np.sum(delta*direction,axis=1)
    radial=np.linalg.norm(delta-axial[:,None]*direction,axis=1)-RADIUS
    cap=axial-p[4]
    return np.hypot(np.maximum(radial,0.),np.maximum(cap,0.))+np.minimum(np.maximum(radial,cap),0.)


def station_constraints(profiles, stations):
    angles=np.linspace(0,2*np.pi,25)[:-1]
    normals=np.column_stack((np.cos(angles),np.sin(angles)))
    rows, supports=[],[]
    for s,profile in zip(stations,profiles):
        if s%10==0 and s<=400:
            rows.extend(np.column_stack((normals,normals*s/SCALE,np.zeros(len(normals)))))
            supports.extend((profile[:,1:]@normals.T).max(axis=0))
    return LinearConstraint(rows,-np.inf,supports)


def refine(p, skin):
    return spar.refine_clearance(skin,p,SCALE,CUTS[-1],engagement(p))


def optimization(seed, objective, points, profiles, stations, skin, bounds, margin):
    rings = points.reshape(len(stations), len(profiles[0]), 3)
    chosen = []
    for up in (-3., 0., 3.):
        shifted = seed.copy()
        shifted[0] += up
        closest = np.argmin(clearance(points, shifted).reshape(rings.shape[:2]),axis=1)
        chosen.extend(rings[np.arange(len(rings)),closest])
    p = seed.copy()
    best = seed.copy()
    interior = station_constraints(profiles, stations)
    tested = []
    for index in range(45):
        chosen.extend(point for value,point,witness in refine(p,skin) if value<1.2)
        support=np.unique(chosen,axis=0)
        step_bounds = [(max(lo,v-0.85),min(hi,v+0.85)) for v,(lo,hi) in zip(best,bounds[:4])]
        step_bounds.append((max(420.,best[4]-.05),min(438.,best[4]+1.8)))
        result=minimize(objective,best,method="SLSQP",bounds=step_bounds,
                        constraints=[interior,{"type":"ineq","fun":lambda x: engagement(x)-MIN_ENGAGEMENT},
                                     {"type":"ineq","fun":lambda x:clearance(support,x)-margin}],
                        options={"ftol":1e-10,"maxiter":120})
        p=result.x
        values=clearance(points,p)
        minima=refine(p,skin)
        continuous,nearest,_=min(minima,key=lambda item:item[0])
        print("search",index,"length",p[4],"span",-axis(p)[1][0]*p[4],
              "wall",values.min(),continuous,"at",-nearest[0],"status",result.message,flush=True)
        if min(values.min(),continuous)>=margin-5e-8:
            if abs(objective(p)-objective(best))<1e-6:
                return p
            best=p.copy()
        chosen.extend(points[np.argpartition(values,95)[:96]])
        if any(np.linalg.norm(p-last)<1e-8 for last in tested):
            break
        tested.append(p.copy())
    return best


def rounded_envelope(tool, radius, length, point, direction):
    if abs(radius-(RADIUS+WALL)) > 1e-12:
        return tool(radius,length,point,direction)
    local = spar.cq.Workplane("XZ").moveTo(0.,0.).lineTo(radius,0.).lineTo(radius,length-WALL)
    mid = WALL / np.sqrt(2.)
    local = local.threePointArc((RADIUS+mid,length-WALL+mid),(RADIUS,length))
    shape = local.lineTo(0.,length).close().revolve(360.,(0.,0.),(0.,1.)).val()
    d = np.asarray(direction.toTuple())
    rotation = np.cross([0.,0.,1.],d)
    angle = np.degrees(np.arccos(d[2]))
    return shape.rotate((0.,0.,0.),tuple(rotation),angle).translate(point)


def run():
    solid=spar.cq.importers.importStep(str(WING)).solids().val()
    skin=next(face for face in solid.Faces() if face.geomType()=="BSPLINE")
    stations=np.arange(0.,601.)
    profiles=spar.sections(solid,stations)
    points=np.concatenate(profiles)
    box=solid.BoundingBox()
    limits=[(box.ymin,box.ymax),(box.zmin,box.zmax),(-99.,99.),(-299.,299.),(350.,630.)]
    seeds=[]
    for name in ["constrained","constrained-250mm-20mm","constrained-250mm-30mm"]:
        payload=json.loads((BASE.parent/name/"result.json").read_text(encoding="utf-8"))
        direction=np.asarray(payload["axis_direction"])
        seeds.append(np.r_[payload["axis_root_mm"][1:],SCALE*direction[1:]/-direction[0],payload["axis_length_mm"]])
    objectives={
        "physical_length":lambda p:-p[4],
        "projected_reach":lambda p:axis(p)[1][0]*p[4],
    }
    optima_path = BASE / "optima.json"
    if "--verify-only" in sys.argv:
        optimized = json.loads(optima_path.read_text(encoding="utf-8"))
    else:
        seed = json.loads((BASE/"seed-75p5.json").read_text(encoding="utf-8"))
        direction = np.asarray(seed["axis_direction"])
        near = np.r_[seed["axis_root_mm"][1:],SCALE*direction[1:]/-direction[0],seed["axis_length_mm"]]
        optimized = {}
        for name, objective in objectives.items():
            starts = [near, seeds[0]]
            if name == "projected_reach":
                starts.insert(0,np.asarray(optimized["physical_length"]["params"]))
            alternatives = [optimization(starts[0],objective,points,profiles,stations,skin,limits,1.001)]
            alternatives.extend(starts)
            p = min([a for a in alternatives if min(v for v,_,_ in refine(a,skin)) >= 1.-1e-7],key=objective)
            optimized[name] = {"params":p.tolist(),"alternatives":[p.tolist() for p in alternatives]}
        optima_path.write_text(json.dumps(optimized,indent=2)+"\n",encoding="utf-8")

    p=np.asarray(optimized["physical_length"]["params"])
    if "--near-boundary" in sys.argv:
        alternatives = [p.copy()]
        for name,objective in objectives.items():
            trial=optimization(p,objective,points,profiles,stations,skin,limits,1.0001)
            alternatives.append(trial)
        (BASE/"boundary-trials.json").write_text(json.dumps(
            [{"params":q.tolist(), "projected_reach_mm":float(-axis(q)[1][0]*q[4]),
              "continuous_wall_mm":min(v for v,_,_ in refine(q,skin))} for q in alternatives],indent=2)
            +"\n",encoding="utf-8")
        p=max(alternatives,key=lambda q:q[4])

    # Re-importable STEP solids are written only for a passing rounded-envelope check.
    if "--reverify-boundary" in sys.argv:
        trials = json.loads((BASE/"boundary-trials.json").read_text(encoding="utf-8"))
        p = np.asarray(max(trials,key=lambda trial:trial["params"][4])["params"])
    near_boundary = "--near-boundary" in sys.argv or "--reverify-boundary" in sys.argv
    test_path = BASE / ("boundary-verified" if near_boundary else "verified")
    test_path.mkdir(exist_ok=True)
    result=None
    tool = spar.cq.Solid.makeCylinder
    planned = [round(p[4],1),p[4],p[4]-.01,p[4]-.05,p[4]-.1,p[4]-.2]
    if near_boundary:
        planned = [p[4], *planned]
    for length in planned:
        trial=p.copy()
        trial[4]=length
        with patch.object(spar,"search",return_value=trial), patch.object(
                spar.cq.Solid,"makeCylinder",side_effect=lambda *args:rounded_envelope(tool,*args)):
            result=spar.verify_and_export(
                solid,trial,stations,profiles,limits[:4],WING,test_path,600.,CUTS,engagement(trial),
                250.,250.,SCALE)
        print("B-rep trial",length,result["axis_length_mm"],result["min_wall_mm"],
              result["clearance_envelope_contained"],flush=True)
        if result["geometry_pass"]:
            break
    if not result["geometry_pass"]:
        raise RuntimeError("No tested candidate passed CAD containment.")
    result["physical_length_optimization"]=optimized
    result["projected_reach_mm"]=-result["axis_channel_end_mm"][0]
    result["tip_setback_cap_center_mm"]=600.-result["projected_reach_mm"]
    result["tip_setback_farthest_channel_point_mm"]=600.-result["channel_farthest_span_mm"]
    result["required_full_diameter_engagement_mm"]=MIN_ENGAGEMENT
    result["clearance_proof"]="B-rep Boolean subtract using a one-millimetre offset envelope with rounded blind-cap edge."
    result["length_adjustment_from_searched_candidate_mm"]=float(length-p[4])
    result["search"]["method"]="SLSQP on four axis placement/direction variables and physical length with exchanged sampled skin constraints and continuous local B-spline minima. B-rep envelope containment determines feasibility."
    (test_path/"result.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":
    run()
