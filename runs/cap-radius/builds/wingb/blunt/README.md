# Wing B (blunt) — rounded tip cap

`wing.step`: the Wing B half-wing (bowed SD7037, root 200 mm, blunt 0.4 mm trailing edge) with a rounded
tip cap from the attachment (X = -540) to the locked outline (600 mm projected half-span). One closed, valid
solid. Axes and units as the parent: mm, X = -span, Y = up, Z = forward.

## Files

| file | purpose |
|---|---|
| `capgeom.py` | geometry definition of the cap (base continuation, lip schedule, patch parametrisation) |
| `build_cap.py` | fits the cap patches, sews them to the parent, writes `wing.step` (~1.5 min) |
| `verify_cap.py` | reads `wing.step` back and measures everything below; writes `readback.json` and the `zebra_*.png` renders (~7.5 min) |
| `readback.json` | all measurements on the reimported STEP, including its SHA-256 |

Run from the project root with `.venv/bin/python runs/cap-radius/builds/wingb/blunt/build_cap.py`, then
`verify_cap.py`. A rebuild gives the same geometry but a new SHA-256, because the STEP header carries a timestamp.

## Construction

- **Parent untouched.** The parent's skin, trailing-edge (TE) face and root face are taken from `parent.step`
  as they are. Their poles, knots and weights read back identical (difference 0). Only the flat attachment face is removed.
- **Base surface.** Past the attachment, the parent skin is continued by its exact second-order Taylor expansion
  in span at fixed native u. Near the attachment the cap is exactly this continuation, so it matches the parent
  in position, tangent and curvature. The mean line (camber, 10° dihedral rise, sweep) comes from this continuation.
- **Closure.** Over the plan region inside the locked outline, the cap is `Y = m ± h`. Here `m` is the base mean line
  and `h` a half-thickness that reaches zero on the outline as a G2 round lip: `h² = H² − H_out²·F(δ/ℓ)`, where
  δ is the plan distance to the outline. F is a profile whose curvature ramps up from 0 over the first 29° of turning
  and then stays constant to the apex, so the lip is G2 everywhere and its smallest radius is 0.989 × the half-thickness.
  The lip width ℓ(p) is scheduled along the outline. On the leading-edge side and at the tip, the roll is spread
  over 0.6 × the medial distance (ℓ ≈ 17–25 mm, lip radius 0.6–1.0 mm). On the TE side, the lip is the roundest
  the thickness allows: a full half-round where the section is thinner than 0.8 mm.
- **Thin TE region.** Where the chord position is above about 0.9, the base thickness is blended (C2-flat at the
  attachment, over about 6 mm of span) to the attachment thickness at the same chord position plus 8 %. Across the
  lip there, the mean line is levelled (camber slope removed within 1.2ℓ and restored by 4ℓ) so the thin
  half-round is symmetric.
- **TE corner.** The parent's flat 0.4 mm TE face continues as a flat strip face (0.40 mm², ruled, tangent to the
  TE face and to both cap patches). The strip tapers to a point 2 mm outboard, where the lip has become a full
  half-round of radius about 0.2 mm.
- **Topology.** Six faces in all: the three parent faces, plus the upper and lower cap patches (bicubic B-splines,
  613 × 37 poles, collapsed to points at the nose and at the TE corners) and the strip. The patch rows at the cut
  are the parent's own edge (error 8e-13 mm), and rows 1–2 come from the exact first and second derivatives
  (G1/G2 by construction). The two patches share the perimeter row exactly and have vertical tangents there.
  Sewn at 1e-6 mm.

## Measured on the reimported STEP (`readback.json`)

**Solid:** 1 solid, 1 shell, 6 faces, no open edges. BRepCheck valid, maximum edge tolerance 1e-7 mm.
Volume 769 329.9 mm³ (parent 747 480.8). Size 7.2 MB.
SHA-256 `2312ee5e4e4c252f58fcf3b3163377efea1a1224062c1d359223d64a444891b4`.

**Locked dimensions:** projected half-span 600.0000001 mm. All join edges lie at X = -540 to within 1e-10 mm.
No cap point lies outside the outline by more than 1.9e-6 mm, and the perimeter edges lie on the outline to
within 1.9e-6 mm over p = 0…1. Perimeter elevation: 48.63 mm at the nose, 62.34 mm at the tip, 48.74 mm at the TE end.

**Joins** (600 samples per edge, excluding 0.05 mm around the collapsed corners):

| join | gap mm | tangent angle | cross-edge curvature diff (median / max, 1/mm) | zebra: phase / direction change |
|---|---|---|---|---|
| parent skin – upper cap | 1.7e-10 | 0.00019° | 8e-8 / 0.0022 (at nose) | 4e-5 stripe / 0.23° |
| parent skin – lower cap | 1.3e-10 | 0.00054° | 2e-8 / 0.0021 (at nose) | 6e-5 stripe / 0.18° |
| parent TE face – strip | 2e-12 | 0.0000065° | 5.5e-5 | 1e-6 stripe / – |
| upper – lower cap (perimeter) | 4.5e-8 | 0.0000015° | 4e-5 / 22.6 (0.09 mm from the nose point) | 5e-8 stripe / 0.05° |
| strip – upper / lower cap | 6e-11 | 0.0000017° | 17–20 / 2e8 (corner) | G1 only, by design (corner transition) |

The zebra check covers 6 views × 4 stripe axes, 8 stripe pairs per half turn; the largest gradient mismatch on
the parent joins is 3 %.

**Printable rule.** The principal curvature map uses 2798 × 1256 points on each patch, then refines the
60 worst local minima with Nelder-Mead. Each point is judged against the attachment thickness at its chord
position (spec.json table):

- Outside the two zones below, every radius meets the rule. The worst ratio is R/rule = 1.017: R = 0.2426 mm
  against a rule of 0.2385 mm (attachment thickness 0.477 mm), at (-555.635, 52.073, -154.770) on the TE lip,
  15.6 mm outboard, chord position 0.994. The smallest radius anywhere outside the zones is 0.2065 mm at
  (-542.02, 49.31, -155.42), on the TE half-round (rule 0.200).
- **TE corner zone** (within 2 mm of span from the attachment, chord position > 0.98): the parent's sharp TE corners
  (radius 0) fade into the half-round. The radius is below the rule for up to 1.71 mm (upper) and 1.78 mm (lower)
  outboard of the attachment.
- **Nose corner** (within 1 mm of span from the attachment, chord position < 0.01): below the rule only within
  0.028 mm of the nose point. The parent's own nose has a 0.064–0.08 mm radius spot within a few µm of its nose
  line (B-spline micro-structure), which the cap carries for its first 1 mm of span. The patch corner collapses at
  the nose point; the grid reaches to 1.3e-6 of the parameter range from it.

**Sections** (planes perpendicular to the outline at 30 values of p; also span stations X = const and chord stations
Z = const; radius = surface normal-curvature radius along the trace):

- Perpendicular sections, outside the TE corner zone: min R/rule ≥ 1.027.
- Lip radius: 0.80 mm at p = 0.1, 1.02 at 0.3, 0.83 at the tip (p = 0.5), 0.67 at 0.7, 0.357 at 0.9 (attachment
  thickness 0.61), 0.214 at 0.99 (thickness 0.40).
- X and Z sections: ratio ≥ 1.02 outside the zones. Their in-plane radii are smaller where the plane crosses the lip
  obliquely, down to 0.09 mm in the Z = const sections near the TE; these are not surface radii.

**Bending reversals.** The sections on the TE side (p ≥ 0.8) are concave between 1.2ℓ and 4ℓ inside the lip, at
radius ≥ 1.03 mm. This is the deliberate mean-line levelling S-bend. Gentler concavity (radius ≥ 20 mm) continues
the parent's own concave aft surfaces: the SD7037 lower surface is concave from 0.54c to the TE, radius ≥ 97 mm.
Elsewhere there are no reversals with radius below 150 mm.

**Renders:** `zebra_full_top.png`, `zebra_tip_*.png`, `zebra_te_corner.png` and `zebra_le_corner.png` show the
reimported surfaces with stripes and no edge lines. Speckle in flat regions is undersampling in the point renderer,
not the surface. One stripe jog remains near the TE, 2–6 mm outboard, where the 8 % TE thickening blends in.

## Deviations to note

- The thin TE region of the cap is 8 % thicker than the parent at the same chord position: about 0.43 mm instead of
  0.40 mm at the TE. Without it, a G2 lip cannot reach radius = half the thickness.
- In the TE corner zone, the rule is not met by construction (inherited sharp corners), as described above.
