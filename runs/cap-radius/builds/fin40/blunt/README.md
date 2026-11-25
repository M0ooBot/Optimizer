# fin40 / blunt — rounded tip cap

\`wing.step\` is the locked parent fin (Z 0–240 mm, 160 mm root chord, 85 mm chord at Z = 240, 12 % elliptical sections, flat 0.4 mm trailing-edge face) with a rounded cap from Z = 240 to the 280 mm crown, as one closed solid. Units are mm, with X aft, Y across the thickness and Z up.

## What was built

| Face | Source |
|---|---|
| root plane, ±Y parent skins | \`parent.step\`, unchanged |
| trailing-edge plane | parent face (a flat B-spline in `parent.step`, 1.6e-12 mm from `parent.py`), merged with its continuation into the cap into one plane face that closes in a pointed arch 1.289 mm above the attachment |
| ±Y cap skins | new B-spline patches (degree 5 × 5, 59 × 74 poles), mirror images that meet along the lip |

**Shape** (\`cap_def.py\`). The parent is a family of horizontal 12 % ellipses, Y² = 0.0144 (x − LE)(E − x), cut by the trailing-edge plane x = 160 − 10 z/240. Its virtual trailing end E sits 0.0327 mm aft of that plane at Z = 240. The cap continues the same family. Each horizontal section is Y² = G(x, z), with zeros on the locked outline:

    G = 0.0144 [ (x − xL)(xT − x) + d(w) (x − xL) (1 − b(w) (l/(r + l))³) ],   r = xT − x,  w = z − 240

- xL and xT are the outline's LE and TE branches.
- d is the parent's scale excess E − TE, faded out over 20 mm.
- b (a quintic smoothstep over 1.5 mm) and the weight (l/(r + l))³ with l = 3 mm pull each section's trailing end through the trailing-edge plane. The weight keeps this near the trailing edge.

The outline ends differ from the parent's straight LE and TE lines by O(h³), and the blends are 1 − O(w³) or O(w³). The cap therefore matches the parent skin in position, tangent plane and curvature (G2) at every chord point. For w ≥ 1.5 mm the bracket is zero exactly on the outline, so the lip is the locked outline.

Below the closure height the section's trailing end lies aft of the trailing-edge plane. The plane cuts it there, as in the parent: the flat face continues upward and narrows to a point at Z = 241.289, X = 149.946. The parent's skin-to-face corner (18.09°) continues around that arch and fades to 0.63° at its tip.

**Construction** (\`build_cap.py\`, CadQuery/OCP 7.9.3):
1. Fit a least-squares tensor B-spline to the analytic skin. Parameters are t ∈ [0, π] around the section and v ∈ [0, 1] up the span. The fit holds these exactly: z depends on v only, Y = 0 on the lip columns, x has zero chordwise slope at the lip, and the v = 1 row collapses to the crown with a horizontal tangent plane. Maximum fit error is 2.2e-8 mm, and 4e-10 mm on the attachment row.
2. Sew the ±Y patches and a flat Z = 240 face into a cap solid. Subtract the half-space aft of the trailing-edge plane.
3. Fuse with \`parent.step\`, then unify the coplanar trailing-edge faces.
4. Rebuild the two attachment edges and the two arch edges with exact same-parameter curves (deviation 4e-9 mm and 3e-13 mm). The model then passes the exact BRepCheck as well as the standard one.
5. Write AP214 STEP.

## Measured on the reimported STEP (\`readback.json\`)

| Check | Result |
|---|---|
| Solid | 1 solid, 1 closed shell, 6 faces, no free edges; BRepCheck valid (standard and exact); max edge tolerance 1e-7 mm |
| Volume | 368 288.016 mm³ (parent 350 160.431 + cap 18 127.585) |
| Parent skin | max 8.9e-9 mm from \`parent.py\`, 1.4e-13 mm from \`parent.step\`'s surfaces; TE face width 0.400 mm at Z = 0 to 240 |
| Height, chords | crown (112.083, 0, 280.000); root chord 160, attachment chord 85 |
| Outline (Y = 0 section vs locked outline) | max 7.1e-4 mm, at the trailing-edge face tip. Below the tip the silhouette is the flat face, which lies on the parent's TE line, and the outline leaves that line by a·h³/2. Where the lip forms the silhouette: max 4.9e-4 mm (just above the tip) |
| Attachment joins (±Y), 2397 points each | separation ≤ 4.2e-9 mm, normal angle ≤ 1.5e-6°, curvature difference ≤ 8.8e-7 /mm (G2) |
| Lip between the two cap patches | normal angle ≤ 1.7e-6°, curvature difference 0 (mirror) |
| **Smallest principal radius on the cap** | **0.5765 mm** at (146.590, 0, 260.079), the lip on the TE branch. Attachment thickness at that chord position is 4.02 mm, so 0.4 mm is required: ratio 1.44. Map: 600 × 1217 grid per patch plus 12 refined minima; the collapsed crown point is excluded within 6.3e-4 mm |
| Thin chord positions (t < 0.8 mm: within 0.13 mm of the LE and 0.098 mm of the TE) | lip radius 0.62–0.71 mm, larger than t/2 |
| Cuts normal to the outline (31) | min 0.5765 mm; in every cut the tightest point is the lip itself; no curvature reversal |
| Constant-X cuts (22) | normal-section radius ≥ 0.577 mm. The plane-section radius falls to 0.029 mm at chord fraction 0.999, because the cut grazes the lip (incidence 0.045); it is not a surface radius. Curvature reversals only at chord fraction ≥ 0.95, 0.8–1.5 mm above the attachment: the trailing-edge face closure, concave radius ≥ 11 mm |
| Zebra, both joins (25 views × 6 stripe axes, 36 stripes per turn) | phase mismatch ≤ 1.2e-6 stripe; stripe direction change at 6570 light/dark crossings ≤ 0.001°; gradient mismatch ≤ 2.9e-4 stripes/mm |

The renders (\`zebra_*.png\`) use CAD normals with no edge lines: the whole fin, the tip, the join at mid-chord and near the LE, and the trailing-edge face closure.

**Departure from the skill (see \`skill-notes.md\`).** The flat trailing-edge face ends in a pointed arch, not a semicircle of radius 0.2 mm. A semicircular end would force the lip to leave the locked outline at about 18° and would break G2 at the trailing-edge corners. Every turn around the closure is 0.6 mm or more, rounder than the 0.2 mm half-round.

## Files and how to run (from the project root)

    .venv/bin/python runs/cap-radius/builds/fin40/blunt/build_cap.py     # ~1 min  -> wing.step, build.json
    .venv/bin/python runs/cap-radius/builds/fin40/blunt/verify_cap.py    # ~10 min -> readback.json
    .venv/bin/python runs/cap-radius/builds/fin40/blunt/render_zebra.py  # ~5 min  -> zebra_*.png

- \`cap_def.py\`: the analytic cap, with exact derivatives from the parent's jets.
- \`occ_util.py\`: shared OCP helpers.
- \`build.json\`: build-side numbers: fit error, edge fixes, volumes, SHA-256.
