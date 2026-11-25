"""Zebra renders of the finished wing (verification step 5: review in context).

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/twisted/blunt/render_zebra.py

Orthographic views, exact CAD normals at the mesh nodes (from the STEP surfaces),
stripes from the reflected view vector against a family of parallel light planes.
The left tip uses mirrored stripe lighting (stripe direction mirrored in XZ) so the
two tips must look identical; the camera is centred on the wing's plane of symmetry
for the full-span view.  Edge lines are not drawn; camera, zoom and lighting are
fixed per view so revisions can be compared image to image.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import verify_cap as vc  # noqa: E402
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRep import BRep_Tool
from OCP.TopLoc import TopLoc_Location
from OCP.GeomLProp import GeomLProp_SLProps
from OCP.TopAbs import TopAbs_REVERSED

OUTDIR = os.path.join(HERE, "renders")


def mesh_faces(shape, defl):
    BRepMesh_IncrementalMesh(shape, defl, False, 0.05, True)
    tris = []
    for f in vc.faces_of(shape):
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(f, loc)
        if tri is None:
            continue
        surf = BRep_Tool.Surface_s(f)
        rev = f.Orientation() == TopAbs_REVERSED
        n = tri.NbNodes()
        P = np.array([[tri.Node(i).X(), tri.Node(i).Y(), tri.Node(i).Z()] for i in range(1, n + 1)])
        N = np.zeros_like(P)
        for i in range(1, n + 1):
            uv = tri.UVNode(i)
            pr = GeomLProp_SLProps(surf, uv.X(), uv.Y(), 1, 1e-9)
            if pr.IsNormalDefined():
                d = pr.Normal()
                N[i - 1] = (d.X(), d.Y(), d.Z())
        if rev:
            N = -N
        T = np.array([[tri.Triangle(i).Value(k) - 1 for k in (1, 2, 3)] for i in range(1, tri.NbTriangles() + 1)])
        tris.append((P, N, T))
    return tris


def render(tris, view, up, center, half, stripe, fname, period=0.12, clip=None, res=1600, mirror_left=True):
    view = np.asarray(view, float); view /= np.linalg.norm(view)
    up = np.asarray(up, float); up -= (up @ view) * view; up /= np.linalg.norm(up)
    right = np.cross(view, up)
    polys, cols, depth = [], [], []
    for P, N, T in tris:
        if clip is not None:
            keep = clip(P[T].mean(axis=1))
            TT = T[keep]
        else:
            TT = T
        if len(TT) == 0:
            continue
        Nt = N[TT]                                   # (m, 3, 3)
        facing = (Nt.mean(axis=1) @ view) < 0
        TT, Nt = TT[facing], Nt[facing]
        r = view[None, None, :] - 2 * (Nt @ view)[..., None] * Nt
        d = np.where((P[TT][..., 1:2] < 0) & mirror_left, stripe * np.array([1, -1, 1]), stripe)
        s = np.sum(r * d, -1).mean(axis=1)
        c = 0.5 + 0.5 * np.sign(np.sin(2 * np.pi * s / period))
        rel = P[TT] - center
        xy = np.stack([rel @ right, rel @ up], -1)
        polys.append(xy); cols.append(c); depth.append((P[TT].mean(axis=1) - center) @ view)
    polys = np.concatenate(polys); cols = np.concatenate(cols); depth = np.concatenate(depth)
    o = np.argsort(-depth)
    fig = plt.figure(figsize=(8, 8), dpi=res // 8)
    ax = fig.add_axes([0, 0, 1, 1])
    pc = PolyCollection(polys[o], array=cols[o], cmap="gray", edgecolors="none", antialiased=False)
    pc.set_clim(0, 1)
    ax.add_collection(pc)
    ax.set_xlim(-half, half); ax.set_ylim(-half, half); ax.set_aspect("equal"); ax.axis("off")
    fig.savefig(fname)
    plt.close(fig)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    shape = vc.read_step(vc.STEP)
    tris = mesh_faces(shape, 0.01)
    stripe = np.array([1.0, 0.35, 0.0]); stripe /= np.linalg.norm(stripe)
    tipc = np.array([140.0, 470.0, 46.0])
    views = [
        ("full_span_top.png", [0, 0, -1], [1, 0, 0], np.array([95.0, 0.0, 30.0]), 520, None),
        ("right_tip_top.png", [0, 0, -1], [1, 0, 0], tipc, 70, lambda c: c[:, 1] > 380),
        ("left_tip_top_mirrored_light.png", [0, 0, -1], [1, 0, 0], tipc * [1, -1, 1], 70, lambda c: c[:, 1] < -380),
        ("right_tip_bottom.png", [0, 0, 1], [1, 0, 0], tipc, 70, lambda c: c[:, 1] > 380),
        ("right_tip_oblique.png", [0.35, -0.55, -0.76], [0.6, 0.2, 0.4], tipc, 70, lambda c: c[:, 1] > 380),
        ("te_corner_top_2mm.png", [0, 0, -1], [1, 0, 0], np.array([189.0, 441.0, 44.0]), 2.0,
         lambda c: (c[:, 1] > 438) & (c[:, 0] > 185)),
        ("te_corner_aft_1mm.png", [-1, 0, 0], [0, 0, 1], np.array([189.86, 440.6, 43.8]), 1.0,
         lambda c: (c[:, 1] > 438) & (c[:, 0] > 185)),
        ("le_corner_top_10mm.png", [0, 0, -1], [1, 0, 0], np.array([95.0, 447.0, 39.0]), 10.0,
         lambda c: (c[:, 1] > 430) & (c[:, 0] < 115)),
    ]
    for name, v, u, c, h, clip in views:
        render(tris, v, u, c, h, stripe, os.path.join(OUTDIR, name), period=0.12 if h > 20 else 0.05, clip=clip)
        print("wrote", name)


if __name__ == "__main__":
    main()
