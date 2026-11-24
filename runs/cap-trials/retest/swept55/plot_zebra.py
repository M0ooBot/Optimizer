"""Draw a centered, edge-free zebra preview from the STEP surface normals."""
import cadquery as cq
import matplotlib.pyplot as plt
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.gp import gp_Pnt, gp_Vec

from verify import HERE, STEP


def plot():
    solid = cq.importers.importStep(str(STEP)).val()
    fig, (full, detail) = plt.subplots(
        1, 2, figsize=(12,5), width_ratios=(2.4,1), layout="constrained"
    )
    for face in solid.Faces():
        center = face.Center()
        if center.z < 0:
            continue
        surface = BRep_Tool.Surface_s(face.wrapped)
        mesh = np.empty((321,181,3))
        for i,q in enumerate((1-np.cos(np.linspace(.001,np.pi-.001,321)))/2):
            for j,v in enumerate((1-np.cos(np.linspace(0,np.pi,181)))/2):
                p,du,dv = gp_Pnt(),gp_Vec(),gp_Vec()
                surface.D1(q,v,p,du,dv)
                normal = du.Crossed(dv).Normalized()
                light = normal.X()*normal.Z()
                mesh[i,j] = p.Y(),p.X(),.5+.5*np.cos(2*np.pi*6*light)

        for ax in (full,detail):
            ax.pcolormesh(*np.moveaxis(mesh,-1,0), shading="gouraud",
                          cmap="gray", vmin=0, vmax=1, rasterized=True)
            ax.set_facecolor("#e7ecf3")
            ax.set_aspect("equal")
            ax.set_xlabel("Span coordinate Y (mm)")
    full.set(xlim=(-465,465), ylim=(740,-25), title="Complete wing")
    detail.set(xlim=(340,460), ylim=(720,480), title="Right tip")
    full.set_ylabel("Aft coordinate X (mm)")
    fig.savefig(HERE/"zebra.png",dpi=180)


if __name__ == "__main__":
    plot()
