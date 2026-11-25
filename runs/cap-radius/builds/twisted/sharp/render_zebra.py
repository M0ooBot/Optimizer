"""Render zebra (reflection-stripe) images of wing.step from CAD normals of the reimported STEP.

Run from the project root:
    .venv/bin/python runs/cap-radius/builds/twisted/sharp/render_zebra.py

Each pixel is rasterised from a fine triangulation only to find its face and (u, v); the shading
normal is then evaluated on the CAD face at that (u, v), so the stripes carry no facet artefacts.
Edge lines are not drawn.  Cameras are mirrored for the left tip and stripe lighting is mirrored
about the centre plane (the stripe normal's Y component flips with the half of the wing).
"""
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import step_geom as SG  # noqa: E402

from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRep import BRep_Tool
from OCP.TopLoc import TopLoc_Location

OUT = HERE / "zebra"
N_BANDS = 14


def mesh_faces(shape, faces, deflection=0.02):
    BRepMesh_IncrementalMesh(shape, deflection, False, 0.2, True)
    data = []
    for key, F in faces.items():
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(F.face, loc)
        n = tri.NbNodes()
        xyz = np.array([[tri.Node(i).X(), tri.Node(i).Y(), tri.Node(i).Z()] for i in range(1, n + 1)])
        uv = np.array([[tri.UVNode(i).X(), tri.UVNode(i).Y()] for i in range(1, n + 1)])
        tris = np.array([tri.Triangle(i).Get() for i in range(1, tri.NbTriangles() + 1)]) - 1
        data.append((key, F, xyz, uv, tris))
    return data


def render(data, eye, target, up, fov_deg, size, stripe_d, fname, title, mirror_stripes=True):
    W, H = size
    eye, target, up = map(lambda a: np.asarray(a, float), (eye, target, up))
    fwd = target - eye; fwd /= np.linalg.norm(fwd)
    right = np.cross(fwd, up); right /= np.linalg.norm(right)
    upv = np.cross(right, fwd)
    f = 0.5 * H / np.tan(np.radians(fov_deg) / 2)
    depth = np.full((H, W), np.inf)
    fid = np.full((H, W), -1, int)
    U = np.zeros((H, W)); V = np.zeros((H, W))
    for k, (key, F, xyz, uv, tris) in enumerate(data):
        rel = xyz - eye
        z = rel @ fwd
        px = W / 2 + f * (rel @ right) / z
        py = H / 2 - f * (rel @ upv) / z
        for t in tris:
            xs, ys, zs = px[t], py[t], z[t]
            x0, x1 = int(max(np.floor(xs.min()), 0)), int(min(np.ceil(xs.max()), W - 1))
            y0, y1 = int(max(np.floor(ys.min()), 0)), int(min(np.ceil(ys.max()), H - 1))
            if x0 > x1 or y0 > y1:
                continue
            gx, gy = np.meshgrid(np.arange(x0, x1 + 1) + 0.5, np.arange(y0, y1 + 1) + 0.5)
            d = (ys[1] - ys[2]) * (xs[0] - xs[2]) + (xs[2] - xs[1]) * (ys[0] - ys[2])
            if abs(d) < 1e-12:
                continue
            l0 = ((ys[1] - ys[2]) * (gx - xs[2]) + (xs[2] - xs[1]) * (gy - ys[2])) / d
            l1 = ((ys[2] - ys[0]) * (gx - xs[2]) + (xs[0] - xs[2]) * (gy - ys[2])) / d
            l2 = 1 - l0 - l1
            inside = (l0 >= -1e-9) & (l1 >= -1e-9) & (l2 >= -1e-9)
            if not inside.any():
                continue
            zz = l0 * zs[0] + l1 * zs[1] + l2 * zs[2]
            sub = depth[y0:y1 + 1, x0:x1 + 1]
            upd = inside & (zz < sub)
            if not upd.any():
                continue
            sub[upd] = zz[upd]
            fid[y0:y1 + 1, x0:x1 + 1][upd] = k
            U[y0:y1 + 1, x0:x1 + 1][upd] = (l0 * uv[t[0], 0] + l1 * uv[t[1], 0] + l2 * uv[t[2], 0])[upd]
            V[y0:y1 + 1, x0:x1 + 1][upd] = (l0 * uv[t[0], 1] + l1 * uv[t[1], 1] + l2 * uv[t[2], 1])[upd]
    img = np.ones((H, W, 3)) * np.array([0.93, 0.93, 0.95])
    ys_, xs_ = np.nonzero(fid >= 0)
    d0 = np.asarray(stripe_d, float); d0 /= np.linalg.norm(d0)
    for y, x in zip(ys_, xs_):
        key, F = data[fid[y, x]][0], data[fid[y, x]][1]
        u = min(max(U[y, x], 0.0), 1.0)
        v = min(max(V[y, x], F.v0), F.v1)
        if key[1] == "cap":
            u = min(max(u, 1e-6), 1 - 1e-6)
        X, Su, Sv, *_ = F.d2(u, v)
        N = np.cross(Su, Sv) * F.sign
        nn = np.linalg.norm(N)
        if nn == 0:
            continue
        N /= nn
        ray = X - eye; ray /= np.linalg.norm(ray)
        R = ray - 2 * (ray @ N) * N
        d = d0 * (np.array([1, -1, 1]) if (mirror_stripes and X[1] < 0) else 1)
        s = np.sin(np.pi * N_BANDS * (R @ d))
        g = 0.12 + 0.83 * (s > 0)
        shade = 0.85 + 0.15 * abs(ray @ N)
        img[y, x] = g * shade
    OUT.mkdir(exist_ok=True)
    fig = plt.figure(figsize=(W / 100, H / 100 + 0.35), dpi=100)
    ax = fig.add_axes([0, 0, 1, H / (H + 35)])
    ax.imshow(img, interpolation="nearest"); ax.axis("off")
    fig.text(0.01, 1 - 0.3 / (H / 100 + 0.35), title, fontsize=9, va="center")
    fig.savefig(OUT / fname, dpi=100)
    plt.close(fig)
    return OUT / fname


def main():
    shape = SG.read_step(HERE / "wing.step")
    faces = SG.wing_faces(shape)
    data = mesh_faces(shape, faces)
    views = [
        ("full_top.png", (95, 0, 1150), (95, 0, 20), (1, 0, 0), 18.6, (1500, 520), (0.25, 0.6, 0.75), "Whole wing from above, centred camera, mirrored stripes"),
        ("full_front_low.png", (-700, 0, -160), (110, 0, 30), (0, 0, 1), 25, (1500, 520), (0.0, 0.5, 0.86), "Whole wing from ahead and below, centred camera, mirrored stripes"),
    ]
    for side, sg in (("right", 1), ("left", -1)):
        views += [
            (f"{side}_tip_top.png", (140, sg * 470, 210), (140, sg * 470, 40), (1, 0, 0), 38, (900, 900), (0.25, 0.6, 0.75), f"{side} tip from above"),
            (f"{side}_tip_oblique.png", (20, sg * 560, 140), (135, sg * 465, 42), (0, 0, 1), 40, (1000, 760), (0.3, 0.55, 0.78), f"{side} tip, oblique from ahead, outboard and above"),
            (f"{side}_tip_below.png", (230, sg * 560, -80), (140, sg * 465, 40), (0, 0, 1), 40, (1000, 760), (0.2, 0.7, 0.69), f"{side} tip, oblique from behind, outboard and below"),
            (f"{side}_te_corner.png", (205, sg * 452, 70), (188, sg * 444, 44), (0, 0, 1), 30, (800, 700), (0.3, 0.55, 0.78), f"{side} trailing-edge corner close-up"),
            (f"{side}_le_corner.png", (70, sg * 452, 55), (92, sg * 443, 38.8), (0, 0, 1), 30, (800, 700), (0.3, 0.55, 0.78), f"{side} leading-edge corner close-up"),
        ]
    for fname, eye, target, up, fov, size, d, title in views:
        print(render(data, eye, target, up, fov, size, d, fname, title), flush=True)


if __name__ == "__main__":
    main()
