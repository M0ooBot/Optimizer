# fin40 / sharp: rounded tip cap

`wing.step` is the whole fin: the parent (Z = 0 to 240 mm, 12 % elliptical sections, 160 mm root chord, 85 mm chord at Z = 240) plus one rounded cap from Z = 240 to the crown at Z = 280 mm. It is one closed, valid solid of 4 faces: the two parent skins and the parent root plane, copied pole-for-pole from `parent.step`, plus one B-spline cap face. The parent's flat top face is gone.

Axes and units follow `spec.json`: X aft, Y thickness, Z up, millimetres.

## Shape

The cap uses the parent's own thickness law with the chord ends moved onto the locked outline:

    Y^2 = 0.0144 [ a(Z)^2 (1 - h^3) - (X - m(Z))^2 ],   h = (Z - 240)/40,
    m = 80 + 27.5 Z/240,   a = 80 - 37.5 Z/240

At h = 0 this is the parent, `Y^2 = 0.0144 (X - X_le)(X_te - X)`. So every horizontal section of the cap is still a 12 % ellipse, with its ends on the locked outline. The law differs from the parent's by `a^2 h^3`, which makes the cap match the parent in position, tangent plane and curvature at Z = 240. The first jump is in the third derivative. On the plane Y = 0 the surface follows the locked outline exactly, and it is smooth everywhere, including the crown.

The lip (the outline itself) has a normal-section radius between 0.577 mm and 0.717 mm. It is 0.634 mm at the LE end and 0.612 mm at the TE end, where it continues the parent's elliptical nose. It is 0.71 mm at the crown. The rest of the cap is flatter than the lip.

## Construction (`build_cap.py`)

1. **Profiles.** The two cut edges are the last pole rows (Z = 240) of the parent's +Y and -Y skins. They are refined exactly by knot insertion, giving 61 poles of degree 5. The middle profile is the locked outline at Y = 0, from `parent.tip_outline`. It is tangent to the parent LE and TE and is written in a parameter w that stays regular through the crown.
2. **Wrap guides.** There is one planar guide per chord station u, and u is the parent's chord parameter. In side view the guide runs straight from the cut point A(u) to the perimeter point L(u), where w(u) = -0.98 cos(pi (u - 0.1 sin 2 pi u)) + 0.02 (2u - 1). Along that line it follows the thickness law above, out to Y = 0 and back to the other skin. With this mapping each guide plane meets the lip within 1.15x of its normal plane, so a guide's own lip radius stays close to the true lip radius.
3. **Boundary surface.** The cap is one B-spline face of degree 5 x 5 with 61 x 29 poles. u runs along the chord and v from the lower cut edge (-1) through the lip (0) to the upper cut edge (+1). The u = 0 and u = 1 edges collapse to the LE and TE points. The boundary pole rows are the parent's cut-edge poles. Every interior pole comes from a least-squares fit to 721 x 601 samples of the analytic guide family; the maximum fit error is 1.5e-6 mm.
4. **Solid.** The parent skins, the parent root plane and the cap are sewn at 1e-6 mm, giving no free or multiple edges. The result is made into a solid and written as AP214 STEP in mm.

## Measured on the reimported STEP (`verify_cap.py` -> `readback.json`)

SHA-256 `66ad12ce7bb499980baf79ec0e81c4c3d605e27002abf9d5160aff09e8cc5f79`, 212 042 bytes.

| Check | Result |
|---|---|
| Topology | 1 solid, 1 closed shell, 4 faces, BRepCheck valid, 0 free-boundary wires |
| Volume | 368 159.32 mm^3 (parent alone 350 036.25) |
| Bounding box | X 0 to 160, Y ±9.6, Z 0 to 280.000001 |
| Parent skins | poles and knots identical to `parent.step` (difference 0); max deviation from `parent.py` 2.6e-8 mm |
| Outline | lip (v = 0) to outline max 9.0e-7 mm; outline to cap face max 9.0e-7 mm; no cap point outside the outline by more than 9.0e-7 mm; crown (112.083, 0, 280) |
| Joins (both skins, 1999 points each, u 0.001 to 0.999) | separation 2.3e-11 mm; tangent-plane angle at most 2.5e-5 deg; curvature in a common frame differs by at most 1.3e-3 /mm across the join, at 0.016 mm from the TE point (relative 7.8e-4), and by at most 4e-6 /mm elsewhere in the tensor |
| Smallest principal radius on the cap (1601 x 1201 grid plus refined maxima, to within 0.003 mm of the collapsed ends) | **0.5765 mm** at (146.59, 0, 260.08), on the lip at the TE side, chord fraction 0.75, attachment thickness 8.84 mm, required 0.4 mm: pass |
| Thin zone (attachment thickness < 0.8 mm) | smallest radius 0.612 mm (the parent's TE nose radius), required at most t/2 = 0.40 mm: pass |
| Grid points breaking the printable rule | 0 |
| Bending reversals | none in 77 section traces; principal k2 >= -1.2e-6 /mm away from the corners (the parent is a developable cone with k2 = 0 at the join) |
| Sections | lip-normal planes: 0.577 to 0.717 mm; guide planes and planes between guides: section radius >= 0.508 mm (thick zone), >= 1.36 t/2 (thin zone); smallest normal radius along any trace 0.5765 mm |
| Zebra across both joins (10 views x 8 stripe angles, 5 deg stripes, 2747 light/dark crossings) | phase mismatch at most 1.0e-5 stripe; stripe direction change at most 0.11 deg; gradient mismatch at most 0.0032 stripes/mm (4.2 % relative, 0.075 mm from the TE point) |
| STEP normals vs the design law | at most 7.1e-5 deg |

Constant-X and horizontal planes cut the lip obliquely, so their section radii are smaller than the surface radius: down to 0.036 mm at X = 149.8 and 0.10 mm at Z = 279.5. By Meusnier's theorem the radius of a section is the normal radius times cos(phi). The normal radius along those same traces is never below 0.5765 mm. `readback.json` lists every section with its section radius, normal radius, obliquity factor and lip crossing.

`review/` holds CAD-normal zebra renders of the whole fin and close-ups of the cap and both corners, with no edge lines, plus `radius_map_side.png`. The short stair-steps on some stripe edges in `zebra_cap_oblique_front.png` come from the renderer's per-triangle colour interpolation. A finer CAD-normal recheck along those lines finds a smooth phase (largest second difference 8.9e-7 stripe at 20 001 samples), and the normals match the design law to within 7.1e-5 deg.

## Files and how to run (from the project root)

    .venv/bin/python runs/cap-radius/builds/fin40/sharp/build_cap.py      # writes wing.step, build.json (about 6 s)
    .venv/bin/python runs/cap-radius/builds/fin40/sharp/verify_cap.py     # writes readback.json (about 1 min)
    .venv/bin/python runs/cap-radius/builds/fin40/sharp/review_render.py  # writes review/*.png (about 20 s)

`surfkit.py` reads the STEP and evaluates its B-spline faces with NumPy: points, outward normals, fundamental forms and zebra phase. `verify_cap.py` checks that evaluator against Open CASCADE's GeomLProp: points agree to 1e-13 mm, normals to 1.5e-6 deg and curvature to a relative 1e-11. A solid classifier confirms that the normals point outward. A rebuild writes a new STEP timestamp, so the SHA-256 changes and `verify_cap.py` must be rerun.
