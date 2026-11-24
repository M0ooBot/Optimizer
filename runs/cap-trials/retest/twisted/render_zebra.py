
"""Render mirrored zebra lighting on the complete STEP solid."""
import cadquery as cq
import matplotlib.pyplot as plt
from matplotlib.artist import Artist
import numpy as np
from OCP.gp import gp_Pnt, gp_Vec

from verify_wing import HERE, STEP, classify, unit, xyz


class ShadedMesh(Artist):
    def __init__(self, triangles, colors):
        super().__init__()
        self.triangles = triangles
        self.colors = colors

    def draw(self, renderer):
        gc = renderer.new_gc()
        renderer.draw_gouraud_triangles(gc, self.triangles, self.colors, self.axes.transData)
        gc.restore()


def mesh(faces):
    all_positions, all_normals, all_triangles, all_sides = [], [], [], []
    offset = 0
    for (side, name), (surface, orientation) in faces.items():
        us = np.linspace(0, 1, 221)
        vs = (np.unique(np.r_[np.linspace(-1, 1, 161), np.linspace(-.06, .06, 121)])
              if name == "cap" else np.linspace(0, 440, 101))
        for u in us:
            for v in vs:
                p, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
                surface.D1(u, v, p, du, dv)
                if name == "cap" and u in (0, 1):
                    surface.D1(np.clip(u, 1e-6, 1-1e-6), v, gp_Pnt(), du, dv)
                all_positions.append(xyz(p))
                all_normals.append(orientation * unit(np.cross(xyz(du), xyz(dv))))
                all_sides.append(side)
        stride = len(vs)
        for i in range(len(us) - 1):
            for j in range(stride - 1):
                a = offset + i*stride+j
                all_triangles.extend(([a, a+stride, a+stride+1], [a, a+stride+1, a+1]))
        offset += len(us)*stride
    return tuple(np.array(item) for item in (all_positions, all_normals, all_triangles, all_sides))


def draw(ax, camera, geometry):
    vertices, normals, triangles, sides = geometry
    forward = unit(camera - [125, 0, 25])
    vertical = unit(np.cross(forward, [0, 1, 0]))
    horizontal = np.array([0., 1., 0.])
    view = camera - vertices
    view /= np.linalg.norm(view, axis=1)[:, None]
    dot = (view * normals).sum(axis=1)
    reflection = 2 * dot[:, None] * normals - view
    phase = 4 * (.5 * reflection[:, 0] + .8660254 * sides * reflection[:, 1])
    shade = .20 + .73 * (1 + np.tanh(2.8*np.cos(2*np.pi*phase))) / 2
    colors = np.column_stack([.84*shade, .92*shade, shade, np.ones(len(shade))])
    front = triangles[(dot[triangles].mean(axis=1) > 0)]
    ordered = front[np.argsort((vertices @ forward)[front].mean(axis=1))]
    projected = (vertices - camera) @ np.column_stack([horizontal, vertical])
    depth = -((vertices - camera) @ forward)
    projected *= np.linalg.norm(camera - [125, 0, 25]) / depth[:,None]
    ax.add_artist(ShadedMesh(projected[ordered], colors[ordered]))
    ax.update_datalim(projected)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.set_axis_off()


def main():
    solid = cq.importers.importStep(str(STEP)).val()
    geometry = mesh(classify(solid))
    fig, axes = plt.subplots(2, 1, figsize=(13, 6), layout="constrained")
    for ax, camera, title in zip(axes, [[125, 0, 1100], [200, 0, -1000]], ["Upper skin", "Lower skin"]):
        draw(ax, np.array(camera, dtype=float), geometry)
        ax.set_title(title, loc="left")
    fig.savefig(HERE / "zebra.png", dpi=170, facecolor="#f5f7fa")


if __name__ == "__main__":
    main()
