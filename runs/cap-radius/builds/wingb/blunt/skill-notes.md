# Skill notes: close-rounded-wingtip on Wing B (blunt)

This is an honest account for the maintainer of `.agents/skills/close-rounded-wingtip/SKILL.md`, from one run on
the blunt Wing B parent with CadQuery/OCP only (no SolidWorks).

## How each step was followed

**Build 1, attachment curves.** Followed. The cap rows at the cut are the parent's own iso-curve at v_attachment,
fitted with knots that include all of the parent's u-knots (error 8e-13 mm). One problem: the parent skin is a
single face that wraps around the nose, so there is no separate "upper cut edge" and "lower cut edge". I built two
patches and let sewing split the parent's edge at the nose. That worked; the skill does not mention it.

**Build 2, 3D perimeter.** Plan shape: locked. Elevation: the continued parent mean line at the outline, so the
perimeter rises with the 10° dihedral bow (48.6 mm → 62.3 mm at the tip). Tangency to the LE and TE is inherited
from the locked outline. No problems.

**Build 3, wrap guides.** Partly followed. The skill pictures guides at fixed chord locations running from the
upper cut, through the perimeter, to the lower cut. Near the LE and TE ends this fails: the locked outline continues
the LE and TE lines, so a guide at a fixed chord station meets the perimeter almost tangentially (about 1° at the TE
end). A half-round drawn in that guide's plane becomes a knife edge in the section perpendicular to the lip. I worked
this out before building and avoided it, but the skill does not warn about it.
- My guides are plan Hermite paths that leave the cut spanwise and reach the perimeter along its normal.
- The lip shape is defined against the distance to the perimeter, not along the guide.
- **Suggested wording:** "guides must cross the perimeter at right angles; shape the lip radius in sections
  perpendicular to the perimeter."

**Build 4, boundary surface.** Improvised. OCP has no boundary-surface tool with guide curves and G2 conditions.
- I defined the cap as a smooth field (`capgeom.py`) and least-squares fitted bicubic B-spline patches to it.
- Exact constraints: cut row = parent edge; rows 1–2 = exact first and second cross derivatives, giving G1/G2;
  shared perimeter row; zero plan derivative at the perimeter, giving vertical, common tangent planes.
- "Match each skin's curvature" was achieved by continuing the parent skin with its own second-order Taylor
  expansion in span. Measured: tangent angle ≤ 0.0005°, cross-edge curvature difference ≤ 0.0022 /mm.

**Build 5, printable rule.** Met everywhere except two corner zones (numbers below). I had to decide four things
the skill does not say.
1. *What "chord position" means for a cap point 40 mm outboard of the attachment.* I used the section-frame chord
   fraction at the point's own span (the spec's definition carried outboard). In the build I used the continued
   section's own nose-to-TE fraction. The two differ by about 1e-4.
2. *"A full half-round whose radius is half the local thickness" cannot be G2.* Only an exact circle turns 90°
   within a rise equal to its radius, and a circle meets the flat skin with a curvature jump. With a G2 lip, the
   smallest radius is at best about 0.99 × the half-thickness (my profile: curvature ramp over the first 29° of
   turning). Also, rule points sit slightly inboard of the apex, where the parent is thicker.
   - Result: the rule as literally stated could not be met with a G2 lip. I thickened the thin TE region of the cap
     by 8 % over the parent's thickness at the same chord position.
   - Needs a ruling: either allow G1 (circular) lips in thin sections, or allow a small thickening.
3. *Camber and dihedral tilt the thin lip.* With a vertical apex, one side of the lip must turn 90° plus the slope
   angle inside the same thickness, which tightened that side by about 15 %. I levelled the mean line across the
   thin lip. That is a deliberate S-bend with concave radius ≥ 1.03 mm, between 1.2ℓ and 4ℓ inside the lip, and it is
   documented in the README. The skill could mention orienting the thin lip along the camber line.
4. *The blunt TE end.* The sentence "the half-round closes the end of the flat trailing-edge face, so its radius
   there is half the trailing-edge thickness" does not work as written here.
   - The locked outline ends at the midpoint of the TE face at X = -540, so the half-round cannot use the face as its
     diameter: its apex would have to sit 0.2 mm aft of, or outboard of, the outline's end.
   - The parent's TE corners are sharp (about 90°), and no surface that is tangent to both the skin and the TE face can
     have a finite radius at those corners.
   - What I did: the flat face continues into the cap as a tapering strip that ends 2 mm outboard, where the lip is a
     full half-round (about 0.2 mm). The cost is a radius below the rule within 1.78 mm of the attachment near the
     TE corners, falling to 0 at the corners themselves.
   - The skill should state what is expected here, and whether this zone is exempt.

**Build 6, progressive turn.** Followed. The lip profile's curvature rises steadily to the apex. The perpendicular
sections show no tight bend before the perimeter, and the minimum radius is at or next to the apex. "Spread the
roll" was implemented as a lip width of 0.6 × the medial distance on the LE side and at the tip. The profile is
capped by thickness on the TE side.

**Verify 1–5.** All done on the reimported STEP (`verify_cap.py`, `readback.json`): closed and valid solid, parent
poles identical, locked outline within 1.9e-6 mm, joins (gap, angle, curvature in a common frame), sections,
principal-curvature map with refinement, numerical zebra (6 views × 4 stripe axes) and renders. Notes:
- **Fixed-plane sections.** X = const and Z = const planes cut the lip obliquely, and their in-plane radii are
  smaller than the surface radius (down to 0.09 mm against a surface radius of 0.36 mm). I report the
  normal-curvature radius along the trace and judge with it. The skill should say which radius to judge, or
  require perpendicular sections only.
- **"Exclude collapsed endpoints".** This needed a definition. I excluded the degenerate patch corners and reported
  two zones separately: the TE corner zone above, and the nose corner, where the parent's own nose has a
  0.064–0.08 mm radius spot within a few µm of its nose line (B-spline micro-structure from the source points).
  That spot is a parent property that the cap inherits near the attachment. Points closer than 0.028 mm to the nose
  point fall below the rule.
- **Zebra "direction change at crossings".** This is ill-conditioned where the stripes are very wide. I count
  crossings only where the stripes are narrower than 20 mm, and measure the gradient mismatch against that floor.
- **Review in context.** I had no CAD viewer, so I used my own point-splat renders with edge lines off. They are
  coarse, but they did reveal two stripe jogs that the numbers alone did not make obvious (see iterations).

## Iterations (about 15 build/verify cycles)

1. A multiplicative lip function: fit failures and knife-like micro-features at the LE corner, where the outline runs
   on the parent's nose. I switched to subtracting in h² so the nose radius is preserved.
2. Micro-structure in the parent's nose spoiled the lip apex near the LE corner (graph pairing at equal Z amplifies
   it). I tried pinning the base nose to the outline, then replaced it with a 0.05 mm apex regularisation faded in
   over 1 mm of span. The pinning was later removed.
3. The thin TE lip was asymmetric (lower side 15 % tight) because of camber slope. I added mean levelling. The first
   zone was too wide and crossed the medial axis at the tip, giving a 1 mm data error, so I capped it.
4. The shoulder of the old profile was tighter than its apex. I replaced it with the curvature-ramp G2 profile plus
   TE thickening (4 %, then 8 %).
5. The thickening was first applied near the nose too, which created a 0.27 mm spot. I restricted it to chord
   positions above 0.85.
6. G1 error of 0.1° at the cut from the least-squares fit. Fixed by setting rows 1–2 exactly from the derivatives.
   I added t-knots near the cut.
7. The zebra renders showed stripe jogs about 4 mm outboard all along the cap, caused by undoing the parent's
   thickness taper. I now keep the taper except in the thin TE region.
8. The renders also showed a stripe step near p ≈ 0.75 from the levelling weight, which I now apply only where the
   lip is thin.
- **Fitting issues along the way:** the locked outline is only C2 at p = 0.18 and 0.82 (fairing), so knots are
  needed there. Interleaved knot sets needed care.

## Requirements: met / partly / not met

- **Met:** single closed valid solid; parent skin and faces preserved exactly; span, attachment and outline locked;
  tangent and curvature match to both skins; common tangent plane and smooth turn on the shared perimeter
  (0.0000015°); finite lip radius everywhere outside the TE corner zone; rule outside the two zones, with minimum
  R/rule = 1.017 (R = 0.2426 mm vs 0.2385 mm at chord 0.994, 15.6 mm outboard) and smallest radius 0.2065 mm on the
  TE half-round (rule 0.200).
- **Partly:**
  - TE corner zone: rule not met within 1.71 mm (upper) and 1.78 mm (lower) of the attachment; radius goes to 0 at
    the parent's sharp corners.
  - Nose point: rule not met within 0.028 mm of it.
  - The thin TE region is 8 % thicker than the parent at the same chord position.
  - The strip-to-patch edges in the TE corner are G1, not G2.
- **Not done:** a SolidWorks / real CAD visual review (no SolidWorks available); the renders are my own.
