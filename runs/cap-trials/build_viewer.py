"""Generate the offline four-model review from the actual exported STEP solids."""

import base64
import gzip
import hashlib
import json
from pathlib import Path

import cadquery as cq
import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve2d, BRepAdaptor_Surface
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepTools import BRepTools
from OCP.GeomAbs import GeomAbs_Plane
from OCP.TopAbs import TopAbs_REVERSED
from OCP.gp import gp_Pnt, gp_Vec

ROOT = Path(__file__).resolve().parent


def packed(array, dtype):
    raw = np.asarray(array, dtype=dtype).tobytes()
    return base64.b64encode(gzip.compress(raw, mtime=0)).decode("ascii")


def rectangular_face(face, bounds):
    """Use a UV grid only when the trim lies on the rectangle's outer sides."""
    u0, u1, v0, v1 = bounds
    for edge in face.Edges():
        curve = BRepAdaptor_Curve2d(edge.wrapped, face.wrapped)
        values = [curve.Value(t).Coord() for t in np.linspace(
            curve.FirstParameter(), curve.LastParameter(), 9)]
        uv = (np.asarray(values) - [u0, v0]) / [u1-u0, v1-v0]
        if not np.any(np.max(np.abs(np.column_stack((uv, uv-1))), axis=0) < 1e-6):
            return False
    return True


def sample_face(face):
    surface = BRepAdaptor_Surface(face.wrapped)
    bounds = BRepTools.UVBounds_s(face.wrapped)
    if surface.GetType() == GeomAbs_Plane or not rectangular_face(face, bounds):
        positions, indices = face.tessellate(0.015, 0.06)
        normals = [face.normalAt(p).toTuple() for p in positions]
        return [p.toTuple() for p in positions], normals, {"triangles": packed(indices, "<u4")}

    u0, u1, v0, v1 = bounds
    # Sample in both surface directions, retaining the CAD normal on each side
    # of a shared edge. Cluster a few extra samples near collapsed boundaries.
    us = np.linspace(0, 1, 193)
    vs = np.linspace(0, 1, 97)
    us = u0 + (u1-u0) * (us - 0.55*np.sin(2*np.pi*us)/(2*np.pi))
    vs = v0 + (v1-v0) * (vs - 0.55*np.sin(2*np.pi*vs)/(2*np.pi))
    positions, normals = [], []
    p, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
    sign = -1 if face.wrapped.Orientation() == TopAbs_REVERSED else 1
    for u in us:
        for v in vs:
            surface.D1(u, v, p, du, dv)
            positions.append(p.Coord())
            normal = du.Crossed(dv)
            if normal.SquareMagnitude() < 1e-16:
                # The point is exact. Use the one-sided normal at a singular vertex.
                surface.D1(np.clip(u, u0+1e-6*(u1-u0), u1-1e-6*(u1-u0)),
                           np.clip(v, v0+1e-6*(v1-v0), v1-1e-6*(v1-v0)), p, du, dv)
                normal = du.Crossed(dv)
            normal.Normalize()
            normals.append(normal.Multiplied(sign).Coord())
    return positions, normals, {"grid": [len(us), len(vs)]}


def model_payload(config):
    path = ROOT / config["step"]
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    report = json.loads((ROOT / config["report"]).read_text(encoding="utf-8"))
    if report.get("readback", report)["step_sha256"] != digest:
        raise ValueError(f"Measurement report belongs to a different STEP revision: {path}")
    solids = cq.importers.importStep(str(path)).val().Solids()
    if len(solids) != 1 or len(solids[0].Shells()) != 1 or not BRep_Tool.IsClosed_s(solids[0].Shells()[0].wrapped):
        raise ValueError(f"Expected one closed solid in {path}")
    if not BRepCheck_Analyzer(solids[0].wrapped, True, False, True).IsValid():
        raise ValueError(f"Invalid reimported solid in {path}")
    faces = []
    points, normals = [], []
    for face in solids[0].Faces():
        pos, norm, topology = sample_face(face)
        faces.append(dict(count=len(pos), **topology))
        points.extend(pos)
        normals.extend(norm)
    # The display uses right, up, forward axes. The signed map is recorded per STEP.
    axes = np.asarray(config["axes"])
    def transform(xyz):
        return np.asarray(xyz)[:, abs(axes)-1] * np.sign(axes)
    xyz = transform(points)
    low, high = xyz.min(axis=0), xyz.max(axis=0)
    extent = high-low
    tip_axis = 1 if config["fin"] else 0
    tip_region = xyz[xyz[:, tip_axis] >= config["cut"] - (20 if config["fin"] else 24)]
    xyz16 = np.rint((xyz-low) / extent * 65535)
    n16 = np.rint(transform(normals) * 32767)
    for name, arr in (("positions", xyz16), ("normals", n16)):
        if not np.isfinite(arr).all():
            raise ValueError(f"Nonfinite {name} in {path}")
    return dict(**config, sha256=digest,
                faces=faces, minimum=low.tolist(), extent=extent.tolist(),
                wholeBounds=[low.tolist(), high.tolist()],
                tipBounds=[tip_region.min(axis=0).tolist(), tip_region.max(axis=0).tolist()],
                positions=packed(xyz16, "<u2"), normals=packed(n16, "<i2"))


if __name__ == "__main__":
    configs = json.loads((ROOT / "cases.json").read_text(encoding="utf-8"))
    models = [model_payload(case) for case in configs]
    template = (ROOT / "viewer.template.html").read_text(encoding="utf-8")
    (ROOT / "review.html").write_text(template.replace("/* MODELS */ []", json.dumps(models)),
                                     encoding="utf-8")
