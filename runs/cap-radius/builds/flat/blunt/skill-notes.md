# Skill notes: close-rounded-wingtip on the flat / blunt parent

## What I followed, and how

**Build 1: attachment curves.** I followed this step. The cap's two cut-edge rows of poles are the parent skin's own edge curve at |Y| = 350, read from `parent.step`, with extra knots inserted. The cap and the parent therefore share the edge exactly.

**Build 2: perimeter.** I followed this step. The perimeter is the locked outline at Z = 0 and is tangent to the leading edge at its start. There is no twist or rise, so the elevation check is trivial.

**Build 3 and 4: wrap guides and boundary surface.** I followed these steps in substance, not with the same tool.
- CadQuery/OCP has no boundary-surface-with-guides feature, so I wrote the surface analytically as an exact tensor-product B-spline.
- Its iso-curves at w = +1, 0 and -1 are exactly the upper cut edge, the perimeter and the lower cut edge. Its iso-curves in the other direction are the wrap guides: one in every constant-X plane, not a handful.
- Upper and lower form one face, so there is no shared-perimeter seam to check.
- The match is G2 to the parent at the cut edges.

**Build 5: printable radius.** I treated it as a design requirement from the start, using an analytic principal-curvature map before any CAD was built.

**Build 6: progressive turn.** I checked it. Over the thick part of the chord, the normal curvature along every constant-X wrap rises monotonically to the lip.

**Verify 1 to 5.** I did all five on the reimported STEP. B-spline faces were evaluated exactly from their STEP poles and knots, and cross-checked against Open CASCADE's D2 to within 2e-13 (relative). The numbers are in `README.md` and `readback.json`.

## Where the skill was unclear, ambiguous or wrong for this parent

1. **Blunt trailing edge with a locked outline.** Build 5 says the half-round closes the end of the flat trailing-edge face with radius = half the trailing-edge thickness. But this locked outline reaches the trailing edge with zero reach, tangent to the trailing edge: it runs 1200 (1-p)^3 mm forward of X = 120. Any closure of finite size therefore has to leave the outline.
   - The skill does not say whether the outline may move near the trailing edge, or by how much.
   - I extended the flat face with a planar rounded closure, apex radius 0.200 mm, reaching 0.352 mm. I shifted the lip aft by at most 2.7 µm over the last 6.8 mm of reach.
   - The closure is a rounded end, not an exact semicircle: its outline has radius 0.200 mm at the apex and 0.1725 mm at its shoulders. Making it an exact circle would conflict with the G2 wrap shape at the corners.

2. **Which chord position judges a point on the cap.** The skill says "at the same chord position" but does not define it for a point on the cap. I used the point's X coordinate, which is the chord station of the wrap it lies on.

3. **"Full half-round" in the thin zone.** Here the thin zone is X > 119.977 mm. The wraps there are 0.35 to 5.3 mm long but only 0.4 to 0.8 mm tall, and the lip runs nearly parallel to the trailing edge. The true lip radius there is large: the cap continues the parent's blunt, near-vertical trailing-edge face outboard. So the lip is not a half-round of radius t/2.
   - I judged the zone by "no radius below t/2". That holds: the smallest radius is 0.356 mm against 0.2 mm required.
   - The literal "full half-round of radius t/2" is **not** met in the thin zone.
   - A literal half-round would make the cap thinner than the parent right at the cut, which tangency forbids.

4. **The parent's crease.** The parent has a 2.85° crease between skin and trailing-edge face, and the skill does not mention it.
   - My first design closed the cap tangent to the trailing-edge plane. Its radius collapsed to 0.001 mm next to the corner, because the crease had to be smoothed out within a vanishing distance.
   - The fix was to let the crease continue along the closure, decaying from 2.85° at the corners to 0.10° at the apex. The skill should say that a creased parent edge must either be carried into the cap or will force unbounded curvature at its end.

5. **Section checks in oblique planes.** On this wing the lip runs almost parallel to X near the leading and trailing edges, so constant-X sections (the guide planes) cut the lip very obliquely.
   - The in-plane radius of such a trace falls to 0.05 mm while the surface's own radius there is 2 to 55 mm (Meusnier).
   - The skill asks for "signed curvature along each trace" and mentions sweep-following cuts only for swept wings. Read literally, it would fail a good cap.
   - I reported trace radius, surface normal-curvature radius and principal radius on each section, and added cuts normal to the lip. The map is what I judged by.

6. **Missing references and tolerances.**
   - The skill links the Wing B study for guidance; I was not allowed to read it.
   - The zebra step says "judge against the design's tolerances", but none are given. I reported raw numbers: phase 7e-8 stripe, gradient 6e-6.

## Decisions and improvisation the skill did not cover

- **Wrap shape.** E = 1 - w^2 and S = the quintic odd smoothstep. I scanned E = (1 - w^2)(1 + A w^2) for A = 0, 0.3 and 0.6 and kept A = 0. A rounder lip pushed the tightest bend off the lip into the wrap, which is the shoulder the skill warns about.
- **Trailing-edge reach law.** At the trailing edge the reach is `R_lock(eps + 2e-3·exp(-eps/6e-3)) + (0.352 - R_lock(2e-3))·exp(-eps/6e-3)`, with eps = 120 - X. It was chosen from a scan of the shift and decay lengths, trading radius margin against outline deviation.
- **TE closure face.** A separate planar face, coplanar with the parent's trailing-edge face. I did not merge the two because the parent's face is a B-spline.
- **Rendering.** Zebra images use exact CAD normals splatted point by point, with no mesh normals, so mesh triangles cannot produce false marks.

## Iterations

There were 3 design iterations on the analytic surface, then 2 CAD-fit iterations.
1. Locked outline exactly, plus a smooth-max cone closure at the trailing edge. **Failed**: radius 0.001 mm at the corner of the trailing-edge face.
2. Locked outline shifted aft globally. Fixed the trailing edge but moved the leading-edge reach. **Rejected**.
3. Local shift, with (shift, decay) scanned:

   | Shift, decay (mm) | Minimum ratio | Lip deviation | Result |
   |---|---|---|---|
   | 5e-4, 1.5e-3 | 0.77 | | Fail |
   | 1e-3, 3e-3 | 1.17 | 1.4 µm | |
   | 1e-3, 1e-2 | 1.46 | 3.4 µm | |
   | 2e-3, 6e-3 | 1.78 | 2.7 µm | Chosen |

CAD fit:
- **Knot refinement.** Spacing ratio 1.45 gave a reach fit error of 1.9e-4 mm; ratio 1.15 gave 1.1e-6 mm.
- **Leading-edge precision bug.** Computing p from 120 - X cancelled near the leading edge. I fixed it by computing p from X there.

## Requirements status

| Requirement | Status |
|---|---|
| Closed valid solid, parent unchanged, span 800, attachment at ±350 | Met. Skin and trailing-edge face poles identical to parent.step |
| Locked outline | Met to 1e-5 mm except near the trailing edge, where it is off by up to **0.0027 mm** (X > 119.95). Needed for the rounded trailing-edge closure |
| Tangent (and curvature) match to both skins | Met. 1.5e-6°, 2.2e-6 /mm |
| Printable radius ≥ 0.4 mm where t ≥ 0.8 mm | Met. Smallest **0.743 mm** at (111.33, ±384.04, 0) on the lip, t = 6.25 mm |
| Thin zone: radius t/2 | Met as "no radius below t/2": smallest 0.356 mm against 0.2 mm required, at (120, ±350.104, ±0.198) |
| Thin zone: lip a full half-round | **Not met literally**: the lip there is a flatter continuation of the blunt trailing-edge face |
| Trailing-edge face closed by a half-round of radius 0.2 mm | **Partly met**. Apex radius 0.200 mm, but not an exact semicircle (reach 0.352 mm, shoulder radius 0.1725 mm in the face plane) |
| No reversals, progressive turn | Met over the thick zone. At X > 119.9999 the wrap bends before the flat aft face, as described above. A tiny saddle region lies at X > 119.9993 |
| Zebra continuity across the joins | Met (G2) |
| Review in context | Done: images/ |
