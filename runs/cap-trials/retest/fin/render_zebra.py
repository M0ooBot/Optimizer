from pathlib import Path

import cadquery as cq
import matplotlib
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.gp import gp_Pnt, gp_Vec

from verify import ANGLES, FREQUENCY, arr, unit

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def sample(face, rows):
    surface = BRep_Tool.Surface_s(face.wrapped)
    view = np.array([0,1,0])
    stripe = np.array([np.cos(ANGLES[2]), 0, np.sin(ANGLES[2])])
    points = np.empty((rows, 801, 3))
    levels = np.empty((rows, 801))
    us = (1 - np.cos(np.linspace(0, np.pi, 801))) / 2
    vs = np.linspace(0, 1, rows)
    for i, v in enumerate(vs):
        for j, u in enumerate(us):
            point, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
            surface.D1(u, v, point, du, dv)
            cross = du.Crossed(dv)
            n = unit(arr(cross)) if cross.Magnitude() > 1e-8 else np.array([0,0,1])
            phase = FREQUENCY * stripe @ (2 * (view @ n)*n - view)
            points[i,j] = arr(point)
            levels[i,j] = 0.5 + 0.5 * np.tanh(3.5 * np.cos(phase))
    return points[:,:,0], points[:,:,2], levels


def plot(reach):
    folder = HERE / str(reach)
    solid = cq.importers.importStep(str(folder/"wing.step")).val().Solids()[0]
    faces = sorted((face for face in solid.Faces() if face.Center().y > .1), key=lambda face: face.Center().z)
    data = [sample(face, rows) for face, rows in zip(faces, [201,241])]
    figure, axes = plt.subplots(1, 2, figsize=(10,6), gridspec_kw={"width_ratios":[1,1.4]})
    for axis in axes:
        for x,z,intensity in data:
            axis.pcolormesh(x,z,intensity,shading="gouraud",cmap="gray",vmin=0,vmax=1)
        axis.set_aspect("equal")
        axis.set_facecolor("#b0b0b0")
        axis.set_xlabel("X aft (mm)")
    axes[0].set_xlim(-8, 168)
    axes[0].set_ylim(-8, 288)
    axes[0].set_ylabel("Z up (mm)")
    axes[0].set_title("Whole fin")
    axes[1].set_xlim(61, 154)
    axes[1].set_ylim(231, 283)
    axes[1].set_title("Cap detail")
    figure.suptitle(f"{reach} mm fin cap | orthographic zebra")
    figure.tight_layout()
    figure.savefig(folder/"zebra.png",dpi=170)
    plt.close(figure)


if __name__ == "__main__":
    for case in (40,5):
        plot(case)
