# fin5 / sharp: rounded tip cap

The whole fin (parent from `runs/cap-radius/parents/fin5/sharp/`) closed at Z = 240 mm by a rounded cap reaching Z = 245 mm. Built with OCP / Open CASCADE 7.9.3. Units are mm, X runs aft, Y is the thickness direction and Z points up the span.

## Files

| file | what |
|---|---|
| `wing.step` | one closed solid: parent root face and both parent skins, unchanged, plus two cap patches |
| `cap_design.py` | the cap definition: wrap profile, spanwise functions and their B-spline fit (shared by both scripts) |
| `build_cap.py` | builds `wing.step` and writes `build.json` |
| `verify_cap.py` | measures the reimported STEP and writes `readback.json` and the `zebra_*.png` pictures |
| `build.json`, `readback.json` | build parameters and fit errors; all measurements, including the STEP's SHA-256 |

To rebuild and re-measure, run these from the project root:

```
.venv/bin/python runs/cap-radius/builds/fin5/sharp/build_cap.py
.venv/bin/python runs/cap-radius/builds/fin5/sharp/verify_cap.py     # about 2 minutes
```

Because the STEP header carries a timestamp, a rebuild changes the SHA-256. Run the verification again after any rebuild.

## Construction

The parent section at height Z is an ellipse with centre m(Z) = 80 + 27.5 Z/240, semi-chord a(Z) = 80 - 37.5 Z/240 and half thickness 0.12 a(Z), parametrised by the angle theta. The locked outline is X = m ± a sqrt(1 - h^3) with h = (Z - 240)/5. This is the parent's own chord, shrunk by sqrt(1 - h^3). The cap uses that directly. Every horizontal section of the cap is an ellipse:

    X = m(Z) - a(Z) sqrt(1 - h^3) cos(theta),   Y = ± Yp(h) sin(theta),   240 <= Z <= 245

- **Perimeter.** theta = 0 and theta = pi lie in Y = 0 and are the locked outline exactly. The outline leaves the leading and trailing edges tangent to them in 3D and with zero curvature.
- **Wrap profile Yp.** This is the mid-chord half thickness. It is designed as a plane curve from its curvature. It leaves the skin with the skin's inward lean (1.07 deg) and zero curvature, and ramps its curvature up with a smoothstep over the first 12 % of its 8.02 mm arc. It peaks at a radius of 4.27 mm, about 1 mm above the attachment. Then it relaxes smoothly to a 5.77 mm crown radius and arrives horizontal at (Y = 0, Z = 245). At other chord positions the wrap follows from the section rule.
- **Patches.** There are two B-spline patches, +Y and -Y, mirror images of each other. In u each uses the parent skin's own u parameter, knots and degree (5, 33 poles). In v each is degree 5 with 39 poles and h = 1 - (1 - v)^2. Each pole is (Mx_j + Ax_j cx_i, Ay_j sy_i, Z_j). Here cx_i and sy_i are read from the parent's top-edge poles, and Mx, Ax, Ay, Z are one-dimensional fits in v.
- **G2 by construction.** The cap's first pole row is the parent's top-edge poles, copied. The v-fits are constrained to the analytic value and first and second derivatives at v = 0, so tangency and curvature match the ruled parent exactly. The crown h = 1 is a collapsed pole. The two patches meet on the perimeter perpendicular to Y = 0, so they share a tangent plane and, by mirror symmetry, their curvature.
- **Assembly.** The parent's flat face at Z = 240 is dropped. Its two skins and root face are sewn with the two cap faces at 1e-7 mm.

## Measured on the reimported STEP (`readback.json`)

| check | result |
|---|---|
| solid | 1 solid, 1 shell, 5 faces, BRepCheck valid, 0 free edges; volume 352,471.23 mm³ (parent 350,036.25 mm³) |
| SHA-256 | `1a1fd864578361a11637fc53015c09a11635f47662e0b848b330b5fcd0e4c2a2` (288,894 bytes) |
| parent skin | poles and knots identical to `parent.step` (0 deviation); 2.7e-8 mm from `parent.py` |
| locked dimensions | height 245.000 mm; attachment edges at Z = 240 to 1.4e-13 mm; outline matched to 4.8e-9 mm (spec samples 3.7e-9 mm, crown 6.7e-12 mm) |
| skin-to-cap joins (2399 points each) | gap ≤ 6.4e-14 mm, tangent-plane angle ≤ 1.5e-6 deg, second fundamental form difference ≤ 9.6e-10 /mm (3.7e-8 relative) |
| perimeter between cap patches | angle ≤ 1.7e-6 deg, curvature difference ≤ 4.7e-6 relative (sampled to 7 µm from the crown) |
| smallest principal radius on the cap | **0.612 mm** at (149.997, 0, 240.080): the trailing-edge end of the perimeter, 0.08 mm above the attachment. It continues the parent's own trailing-edge nose radius (0.6125 mm). The leading-edge end is 0.634 mm. |
| printable rule (0.4 mm where t ≥ 0.8 mm, t/2 below) | met at every one of 2 × 652,740 grid points and at the refined minima. Smallest margin R/required = 1.68 (R = 0.670 mm where t = 0.811 mm) |
| crown radii | 5.76 mm across, 522 mm along the chord |
| sections (80 traces) | no curvature reversals. Surface radius along every trace ≥ 0.612 mm. The lip radius in planes normal to the perimeter is 0.635 → 5.69 mm (LE branch) and 0.612 → 5.67 mm (TE branch) from h = 0.02 to 0.99. At mid-chord (X = 107.5) the tightest point is 4.27 mm, 0.96 mm above the attachment, and the lip is 5.77 mm |
| zebra (40 views × 6 stripe axes) | phase jump ≤ 1.4e-9 periods and stripe-direction change ≤ 1.4e-5 deg on the skin joins. On the perimeter, 9.2e-8 periods and 5.6e-5 deg |

## Things to know

- **The wrap is not monotone.** At every chord station the tightest point lies below the lip, and the turn opens slightly toward the lip: at mid-chord from 4.27 mm to 5.77 mm, and at 10 % chord from 2.81 mm to 3.22 mm. There are no reversals. The locked 5 mm reach is less than the 5.1 mm half thickness (5.006 mm after the skin's lean), which leaves room for only a near-circular quadrant. Any G2 lead-in must turn tighter than the half thickness somewhere and then relax. The ramp length trades this: 10 % of the arc gives 4.39/5.63 mm, 15 % gives 4.10/6.02 mm, 20 % gives 3.84/6.54 mm.
- **Grazing section cuts read tight.** Planar cuts at constant X, or along the parent's rulings, near the leading and trailing edges meet the surface at a grazing angle at the lip. Their in-plane radius is therefore as small as 0.10 mm, and 0.25–0.39 mm where t = 0.91–1.12 mm. The surface's own radius at those points is 0.62–0.82 mm. The printability judgement uses principal curvature (see `skill-notes.md`).
- **The crown is a collapsed pole.** Derivative-based measures stop 7 µm (curvature map: 0.8 µm) from it.
- **Zebra pictures.** `zebra_*.png` are orthographic renders from exact CAD normals, with no edge lines and the camera fixed in `verify_cap.py`.
