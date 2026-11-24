"""Independent, paired physical-section and seam audit of both generations of STEP."""
import argparse
import hashlib
import json
from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve2d, BRepAdaptor_Surface
from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepTools import BRepTools
from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector
from OCP.TopAbs import TopAbs_REVERSED
from OCP.gp import gp_Dir, gp_Lin, gp_Pln, gp_Pnt, gp_Vec
from scipy.optimize import brentq, minimize, minimize_scalar

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CASES = {
    "flat": ("flat/flat-wing.step", 350, 50),
    "swept55": ("swept55/swept55-wing.step", 390, 60),
    "twisted": ("twisted/twisted-wing.step", 440, 60),
    "fin/40": ("fin/vertical-fin.step", 240, 40),
    "fin/5": ("fin/cap-5.0mm/vertical-fin.step", 240, 5),
}
SWEEP = np.tan(np.deg2rad(55))


def unit(v):
    return v / np.linalg.norm(v)


def chord_station(case, s):
    if case == "flat":
        return 0, 120
    if case == "swept55":
        return s*SWEEP, 210-115*s/390
    if case == "twisted":
        return 90*s/440, (180-80*s/440)*np.cos(np.deg2rad(3*s/440))
    return 65*s/240, 160-75*s/240


def rim(case, u):
    _, cut, reach = CASES[case]
    if case.startswith("fin"):
        # First branch ascends the leading edge; second descends the trailing edge.
        h = 1-abs(2*u-1)
        z = cut+reach*h
        le, chord = chord_station(case,z)
        x = le+chord/2*(1+np.sign(u-.5)*np.sqrt(1-h**3))
        return np.array([x,0,z])
    if case == "flat":
        return np.array([120*u**3*(10-15*u+6*u*u), cut+4*reach*u*(1-u), 0])
    if case == "swept55":
        y=cut+4*reach*u*(1-u)
        le,chord=chord_station(case,y)
        return np.array([le+chord*(3*u*u-2*u**3),y,0])
    y=cut+reach*np.sin(np.pi*u)
    x=np.sin(np.pi*u/2)**2
    blend=np.clip((x-.36)/.08,0,1)
    blend=blend**3*(10-15*blend+6*blend**2)
    camber=(1-blend)*.02/.16*(.8*x-x*x)+blend*.02/.36*(.2+.8*x-x*x)
    theta=np.deg2rad(3*y/cut)
    c=180-80*y/cut
    return np.array([90*y/cut+c*(x*np.cos(theta)-camber*np.sin(theta)),y,
                     y*np.tan(np.deg2rad(5))+c*(x*np.sin(theta)+camber*np.cos(theta))])


def planes(case):
    fin=case.startswith("fin")
    _,cut,_=CASES[case]
    fractions=([.05,.15,.25,.36,.45,.55,.64,.75,.85,.95] if fin else
               [.02,.05,.25,.40,.50,.65,.80,.90,.97]+([.99] if case!="swept55" else []))
    for f in fractions:
        le,chord=chord_station(case,cut)
        root=chord_station(case,0)[1]*f
        slope=(le+chord*f-root)/cut
        if case=="swept55":
            slope,root=SWEEP,95*f
        n=unit(np.array([1.,0,-slope] if fin else [1.,-slope,0]))
        origin=np.array([root,0.,0.])
        p=rim(case,brentq(lambda t: (rim(case,t)-origin)@n,0,1))
        yield f,n,origin,p


def jet(surface,u,v,orientation):
    p=gp_Pnt()
    a,b,c,d,e=[gp_Vec() for _ in range(5)]
    surface.D2(u,v,p,a,b,c,d,e)
    du,dv,duu,dvv,duv=(np.array(w.Coord()) for w in (a,b,c,d,e))
    n=orientation*unit(np.cross(du,dv))
    inv=np.linalg.pinv(np.column_stack((du,dv)))
    h=inv.T@np.array([[n@duu,n@duv],[n@duv,n@dvv]])@inv
    return np.array(p.Coord()),n,h


def orient(face):
    return -1 if face.wrapped.Orientation()==TopAbs_REVERSED else 1


def section_trace(face,n,origin,lip):
    surface=BRepAdaptor_Surface(face.wrapped)
    sec=BRepAlgoAPI_Section(face.wrapped,gp_Pln(gp_Pnt(*origin),gp_Dir(*n)),False)
    sec.ComputePCurveOn1(True)
    sec.Approximation(True)
    sec.Build()
    if not sec.IsDone():
        raise RuntimeError("CAD section failed")
    result=[]
    for edge in cq.Shape.cast(sec.Shape()).Edges():
        pc=BRepAdaptor_Curve2d(edge.wrapped,face.wrapped)
        start,end=pc.FirstParameter(),pc.LastParameter()
        def uv(t):
            return pc.Value(start+(end-start)*t).Coord()
        def value(t):
            p,normal,h=jet(surface,*uv(t),orient(face))
            cross=np.cross(normal,n)
            speed=np.linalg.norm(cross)
            tangent=cross/speed
            return p,float(tangent@h@tangent/speed)
        def distance(t):
            return np.linalg.norm(np.array(surface.Value(*uv(t)).Coord())-lip)
        t_lip=min([0.,1.,minimize_scalar(distance,bounds=(0,1),method="bounded",
                                        options={"xatol":1e-14}).x],key=distance)
        ts=np.linspace(0.,1.,401)
        xyz,ks=zip(*(value(t) for t in ts))
        ks=np.asarray(ks)
        idx=np.argmax(np.abs(ks))
        witness=minimize_scalar(lambda t:-abs(value(t)[1]),bounds=(ts[max(0,idx-1)],ts[min(400,idx+1)]),
                                method="bounded",options={"xatol":1e-13}).x
        t=max([ts[idx],witness],key=lambda t:abs(value(t)[1]))
        pos,k=value(t)
        # One common numerical deadband so roundoff around a flat attachment
        # cannot invent reversals. Report crossings on individual CAD branches.
        nonzero=np.flatnonzero(abs(ks)>1e-8)
        reversals=[value(brentq(lambda t:value(t)[1],ts[i],ts[j]))[0].tolist()
                   for i,j in zip(nonzero,nonzero[1:]) if ks[i]*ks[j]<0]
        result.append(dict(radius=1/abs(k), location=pos.tolist(), lip_distance=distance(t_lip),
                           lip_radius=1/abs(value(t_lip)[1]) if distance(t_lip)<.01 else None,
                           reversals=reversals, cut_plane_error=max(abs((p-origin)@n) for p in xyz)))
    return result


def shared_edges(faces,cut,axis):
    groups={"attachment":[],"rim":[]}
    for i,a in enumerate(faces):
        for b in faces[i+1:]:
            for edge in a.Edges():
                if edge.Length()<1e-5 or not any(edge.isSame(other) for other in b.Edges()):
                    continue
                samples=np.array([edge.positionAt(f).toTuple() for f in (.1,.5,.9)])
                s=np.abs(samples[:,axis])
                if np.max(abs(s-cut))<1e-4:
                    group="attachment"
                elif s.min()>cut:
                    group="rim"
                else:
                    continue
                groups[group].append((a,b,edge))
    return groups


def seam_stats(edges):
    # Identical diagnostic fields/units for all models. Both opposing faces
    # use the same physical axes and eyes; the sign of the normal is immaterial.
    eyes=[]; axes=[]
    for vec in ([1,.2,1],[-1,.3,-.7],[.4,1,.8],[.5,-1,-.6]):
        eye=unit(np.array(vec));right=unit(np.cross([0,0,1],eye));up=np.cross(eye,right)
        for angle in np.deg2rad([0,30,60,90,120,150]):
            eyes.append(eye);axes.append(np.cos(angle)*right+np.sin(angle)*up)
    eyes,axes=np.array(eyes),np.array(axes)
    peaks=np.zeros(5)
    for a,b,edge in edges:
        pcs=[BRepAdaptor_Curve2d(edge.wrapped,f.wrapped) for f in (a,b)]
        surfaces=[BRepAdaptor_Surface(f.wrapped) for f in (a,b)]
        for t in np.linspace(.004,.996,129):
            samples=[]
            for f,pc,surf in zip((a,b),pcs,surfaces):
                uv=pc.Value(pc.FirstParameter()+t*(pc.LastParameter()-pc.FirstParameter()))
                p,n,h=jet(surf,*uv.Coord(),orient(f))
                phase=4*(2*(eyes@n)*(axes@n)-np.sum(eyes*axes,axis=1))
                grad=-8*((eyes@n)[:,None]*axes+(axes@n)[:,None]*eyes)@h
                samples.append((p,n,h,phase,grad))
            (p,n,h,pa,ga),(q,m,j,pb,gb)=samples
            peaks=np.maximum(peaks,[np.linalg.norm(p-q),np.rad2deg(np.arctan2(np.linalg.norm(np.cross(n,m)),n@m)),
                np.linalg.norm(h-j),np.max(abs(pa-pb)),np.linalg.norm(ga-gb,axis=1).max()])
    return dict(edges=len(edges),**dict(zip(("gap_mm","normal_deg","curvature_per_mm","phase_cycles",
                     "gradient_cycles_per_mm"),peaks.tolist()))) if edges else {"edges":0}


def parent_samples(case,solid):
    _,cut,_=CASES[case];fin=case.startswith("fin")
    hit=IntCurvesFace_ShapeIntersector();hit.Load(solid.wrapped,1e-8)
    rows=[]
    for fspan in (.03,.27,.58,.83,.997):
        s=cut*fspan
        le,chord=chord_station(case,s)
        for f in [.02,.05,.15,.3,.5,.75,.9,.97]:
            for side in ([1] if fin else [-1,1]):
                x=le+chord*f
                start,direction=([x,-1500,s],[0,1,0]) if fin else ([x,side*s,-1500],[0,0,1])
                hit.Perform(gp_Lin(gp_Pnt(*start),gp_Dir(*direction)),0,3000)
                heights=sorted(hit.Pnt(i).Coord()[1 if fin else 2] for i in range(1,hit.NbPnt()+1))
                if len(heights)!=2: raise ValueError((case,fspan,f,heights))
                rows.append(heights)
    return rows


def outer_extent(face,axis,side):
    # These cap faces are rectangular in UV, including the half-domain trims.
    # Optimizing the actual surface also captures extrema in the interior of
    # a single wrap patch that the kernel bounding box can underestimate.
    u0,u1,v0,v1=BRepTools.UVBounds_s(face.wrapped)
    surf=BRepAdaptor_Surface(face.wrapped)
    def objective(t):
        return -side*surf.Value(u0+(u1-u0)*t[0],v0+(v1-v0)*t[1]).Coord()[axis]
    seeds=[(u,v) for u in np.linspace(0,1,11) for v in np.linspace(0,1,11)]
    start=min(seeds,key=objective)
    fit=minimize(objective,start,bounds=[(0.,1.)]*2,method="L-BFGS-B",options={"ftol":1e-15})
    return -min(objective(start),fit.fun)


def mirror_check(case,solid):
    _,cut,_=CASES[case]
    fin=case.startswith("fin")
    axis=2 if fin else 1
    thickness_axis=1 if fin else 2
    ray=IntCurvesFace_ShapeIntersector()
    ray.Load(solid.wrapped,1e-8)
    gaps=[]
    for fraction,n,origin,lip in planes(case):
        station=origin[0]-n[axis]*cut/n[0]
        for t in (.20,.55,.85):
            s=cut+t*(lip[axis]-cut)
            x=station+t*(lip[0]-station)
            heights=[]
            for side in ([1] if fin else [1,-1]):
                start,direction=([x,-1500,s],[0,1,0]) if fin else ([x,side*s,-1500],[0,0,1])
                ray.Perform(gp_Lin(gp_Pnt(*start),gp_Dir(*direction)),0,3000)
                hits=sorted(ray.Pnt(i).Coord()[thickness_axis] for i in range(1,ray.NbPnt()+1))
                if len(hits)!=2: raise ValueError(("mirror",case,fraction,t,hits))
                heights.append(np.array(hits))
            if fin:
                gaps.append(abs(sum(heights[0])))
            else:
                gaps.append(float(max(abs(heights[0]-heights[1]))))
    return max(gaps)


def audit(case,path):
    _,cut,_=CASES[case];axis=2 if case.startswith("fin") else 1
    solids=cq.importers.importStep(str(path)).val().Solids()
    if len(solids)!=1: raise ValueError("Expected one solid")
    solid=solids[0]
    faces=solid.Faces()
    cap=[f for f in faces if f.Center().toTuple()[axis]>cut+1e-6]
    sections=[]
    for fraction,n,origin,p in planes(case):
        segments=[s for f in cap for s in section_trace(f,n,origin,p)]
        lips=[s for s in segments if s["lip_radius"] is not None]
        if not lips: raise ValueError((case,fraction,"prescribed lip missing in section",segments))
        tight=min(segments,key=lambda s:s['radius'])
        sections.append(dict(fraction=fraction,smallest_radius=tight["radius"],at=tight["location"],
                    lip_radius=min(s["lip_radius"] for s in lips),max_plane_error=max(s["cut_plane_error"] for s in segments),
                    rim_error=min(s["lip_distance"] for s in lips),branches=segments))
    extent=max(outer_extent(f,axis,1) for f in cap)
    if axis==1:
        extent+=max(outer_extent(f,axis,-1) for f in faces if f.Center().y < -cut)
    return dict(step_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),solid_count=len(solids),
                closed=len(solid.Shells())==1 and BRep_Tool.IsClosed_s(solid.Shells()[0].wrapped),
                valid=BRepCheck_Analyzer(solid.wrapped,True,False,True).IsValid(),faces=len(faces),
                extent_mm=extent,parent=parent_samples(case,solid),mirror_error_mm=mirror_check(case,solid),
                seams={k:seam_stats(v) for k,v in shared_edges(faces,cut,axis).items()},sections=sections)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",choices=CASES)
    parser.add_argument("version",choices=["control","rebuilt"])
    args=parser.parse_args()
    source=ROOT/CASES[args.case][0] if args.version=="control" else HERE/args.case/"wing.step"
    results=audit(args.case,source)
    out=HERE/"measurements"/args.case/(args.version+".json")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(results,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(args.case,args.version,results['seams'],flush=True)
