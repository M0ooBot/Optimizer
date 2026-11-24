"""Render the imported wing with a centered, edge-free zebra view."""
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.gp import gp_Pnt, gp_Vec

from verify_wing import OUT, cq, unit, vec


def samples(surface, span, side):
    chord = np.linspace(.0001,.9999,241)
    horizontal = np.empty((len(chord),len(span)))
    vertical = np.empty_like(horizontal)
    stripe = np.empty_like(horizontal)
    light = unit([.75,.65*side,.2])
    p, a, b = gp_Pnt(), gp_Vec(), gp_Vec()
    for i,u in enumerate(chord):
        for j,v in enumerate(span):
            surface.D1(float(u),float(v),p,a,b)
            normal = unit(np.cross(vec(a),vec(b)))
            reflection = 2*normal[2]*normal - np.array([0,0,1])
            horizontal[i,j], vertical[i,j] = p.Y(), p.X()
            stripe[i,j] = .5+.5*np.sin(8*np.pi*(reflection @ light))
    return horizontal, vertical, stripe


def main():
    wing = cq.importers.importStep(str(OUT / 'wing.step')).val()
    faces = [f for f in wing.Faces() if f.Center().z > 0]
    faces.sort(key=lambda f: f.Center().y)
    surfaces = [BRepAdaptor_Surface(f.wrapped) for f in faces]
    panels = [
        samples(surfaces[0], np.linspace(0,1,251), -1),
        samples(surfaces[1], np.linspace(0,.5,101), -1),
        samples(surfaces[1], np.linspace(.5,1,101), 1),
        samples(surfaces[2], np.linspace(0,1,251), 1),
    ]
    fig, ax = plt.subplots(figsize=(12,2.25), dpi=180)
    fig.patch.set_facecolor('#cbd2dc')
    ax.set_facecolor('#cbd2dc')
    for y,x,c in panels:
        ax.pcolormesh(y,x,c,shading='gouraud',cmap='gray',vmin=0,vmax=1)
    ax.set_aspect('equal')
    ax.set(xlim=(-418,418), ylim=(129,-9))
    ax.axis('off')
    fig.subplots_adjust(left=.012,right=.988,bottom=.04,top=.96)
    fig.savefig(OUT / 'zebra.png')
    plt.close(fig)


if __name__ == '__main__':
    main()
