# Skill notes: close-rounded-wingtip on fin40 / blunt

These notes are for the skill's maintainer. I read only the skill, my parent folder and the project's \`src/\` and \`tools/\` listings, and built with CadQuery/OCP. There was no SolidWorks.

## Build steps: what I did

**1. Attachment curves.** I followed this. The cap's bottom boundary is the parent's own cut edges: the ±Y skin edges at Z = 240 in \`parent.step\`, which the solid shares after the fuse. I fitted the cap's v = 0 row to the analytic attachment section (4e-10 mm) and gave the shared edges exact pcurves (4e-9 mm). The skill does not say what happens to the **trailing-edge segment** of a blunt cut. The cut has three pieces (upper edge, the 0.4 mm segment, lower edge), but the boundary surface has only three profiles (upper, perimeter, lower). I had to add a face for the segment; see step 5.

**2. Perimeter curve.** I used the locked outline in the plane Y = 0 as it stands. It is already tangent to the LE and TE lines, with zero curvature at both ends (its h³ term). The fin is symmetric with no twist, so the elevation check is trivial. One thing is unclear: on a blunt parent the outline ends at the *midpoint* of the 0.4 mm trailing-edge segment, while the skill speaks of arriving "at the trailing-edge end of the cut".

**3. Wrap guides.** I improvised here. I drew no discrete guide curves. OCP has no boundary-surface-with-guides tool, and its filling and plate tools give neither exact boundaries nor G2. Instead I defined the whole cap analytically: horizontal sections continue the parent's 12 % ellipse family, with ends on the outline (formula in README). That is a continuum of wraps. The lip radius came out of the section shape and the outline (0.58–0.72 mm); it was not a free parameter, so "give the turn a finite lip radius" needed no separate decision.

**4. Boundary surface.** I followed this in substance, not in method. The cap has two patches (±Y): a least-squares B-spline fit of the analytic surface, accurate to 2.2e-8 mm. Tangency and curvature match both skins by construction: the difference from the parent is O(w³). The two patches share the lip with a common tangent plane (1.7e-6°) and equal curvature across it (mirror images).

**5. Printable radius.** Met for the surface, but the trailing edge was not built as written.
- Smallest principal radius on the cap: 0.5765 mm at the TE-branch lip, Z = 260.08 mm, X = 146.59 mm. Attachment thickness there is 4.02 mm, so the requirement is 0.4 mm (ratio 1.44).
- Thin chord positions (t < 0.8 mm: chord fraction < 0.00154 or > 0.99884) have lip radii of 0.62–0.71 mm, more than t/2.
- **I did not build the "half-round that closes the end of the flat trailing-edge face" with radius 0.2 mm.** Read literally, the face end is a semicircle of radius 0.2 mm in the face's plane. That conflicts with the locked outline. Where the face end meets the lip, the lip leaves the TE line at angle α with tan α = r_end / r_nose. The nose radius is fixed at about 0.61 mm by the parent's 12 % section and G1 at the attachment, so a 0.2 mm end gives α ≈ 18°. The lip would then have to leave the locked outline, or kink back onto it. A semicircular end also breaks G2 at the trailing-edge corners: spanwise curvature 5 /mm against about 0 on the parent.
- What I built instead: each section's virtual trailing end retreats through the trailing-edge plane over about 1.3 mm. The flat face narrows to a pointed arch 1.289 mm above the attachment. The skin-to-face corner fades from the parent's 18.09° to 0.63° at the tip. No radius anywhere in that closure is below 0.6 mm.
- The outline cost is 7.1e-4 mm at the tip. That is the outline's own departure from the TE line there; the flat face holds the silhouette below the tip.
- In intent this is "as round as the section allows"; literally, the face end is not a half-round.

**6. Progressive turn.** Met.
- In all 31 cuts normal to the outline, curvature rises steadily into the lip. The lip is the tightest point and nothing tightens earlier.
- Constant-X cuts show curvature reversals only at chord fraction ≥ 0.95, 0.8–1.5 mm above the attachment. This bend is deliberate: the face closure has to slow down to meet the outline tangentially. Its concave radius is at least 11 mm (at chord fraction 0.998) and about 1.2 m at 0.95.

## Verification steps

1. **Reimport checks.** Done. The model is valid under both the standard and the exact BRepCheck. Getting the exact check to pass took extra work: OCC's boolean leaves pcurves and intersection curves that the exact check rejects, so I rebuilt those edges. The skill could say which validity level is expected.
2. **Sections.** Done: 22 constant-X cuts, 31 cuts normal to the outline (both branches) and a crown cut, with signed curvature, reversals, smallest radius, lip radius, attachment thickness and verdict for each.
   - **The skill does not say which radius to judge.** Near the LE and TE, a constant-X plane crosses the lip at a grazing angle (incidence down to 0.045). The plane-section curve then has radii down to 0.029 mm, although the surface's smallest radius there is about 0.6 mm. I report the plane-section radius, the normal-section radius along the trace and the incidence, and judge by the normal section and the principal-curvature map. The skill should say that the printable rule applies to surface curvature, or to sections normal to the lip, and that grazing cuts mislead.
   - "Chord position" above the attachment is ambiguous on a swept planform. I used the attachment frame (X − 65)/85, as \`spec.json\` tabulates it.
   - The crown cut passed exactly through the collapsed pole and returned a sign-flipped curvature there. I excluded points within 1e-5 in v of the pole.
3. **Curvature map.** Done: 600 × 1217 points per patch plus Nelder–Mead refinement of 12 seeds. The crown pole is excluded within 6.3e-4 mm. Smallest radius 0.5765 mm, 1.44 × the requirement.
4. **Numerical zebra.** Done, on both joins.
   - Set-up: 25 view directions × 6 stripe axes, 36 stripes per turn, CAD normals with Weingarten derivatives.
   - Results: phase mismatch ≤ 1.2e-6 stripe; stripe direction change at 6570 crossings ≤ 0.001°; gradient mismatch ≤ 2.9e-4 stripes/mm. Relative mismatch is 0.5% where the gradient is at least 10% of that view's median. The skill gives no tolerance, so I judged these as continuous.
   - A relative gradient-angle metric spikes (0.9°) where the stripe gradient is almost zero (9e-6). The skill could warn about this.
5. **Review in context.** I rendered zebra images with my own depth-buffered renderer, using CAD normals and no edges, at fixed cameras. The views are the whole fin, the tip, the join at mid-chord and near the LE, and the trailing-edge closure. Iteration 1's concave band showed as a visible notch in the stripes, which confirmed the numerical finding. The final renders show no jog at the join.

## Iterations: 3 design iterations

1. **Built.** Pure ellipse sections: one trailing end per height, retreating the parent's 0.04 % chord excess over 1.5 mm. Radius rule and G2 were met, but because the retreat acted on the whole chord it put a concave band (−0.0047 /mm at mid-chord, −0.1 /mm near the TE) 0.7–1.4 mm above the attachment across the whole span of chord. It showed as a zebra notch. That version also failed the exact BRepCheck, which I fixed in the build.
2. **Analytic only.** I confined the closure to the trailing edge with a weight l/(r + l), l = 1 mm. That was still concave back to X ≈ 133. The same term had a square-root singularity at the crown, so I added the 20 mm fade of the excess.
3. **Final.** Weight (l/(r + l))³, l = 3 mm, a general B-spline fit, and refined span knots (fit error from 1.2e-7 to 2.2e-8 mm). The concave band is now limited to the last 5 % of chord.

## Requirements: met, partly met, not met

- **Met:**
  - one closed valid solid;
  - parent skin, attachment, span (280.000) and chords preserved;
  - G1 and G2 at both joins (4e-9 mm, 1.5e-6°, 8.8e-7 /mm);
  - shared perimeter;
  - printable radius ≥ 0.4 mm wherever t ≥ 0.8 (min 0.5765 mm);
  - progressive turn;
  - all verification steps.
- **Partly met:**
  - Thin-section rule: radii exceed t/2, but they are not "a half-round of radius t/2".
  - Locked outline: within 7.1e-4 mm measured. By design it is exact from 1.5 mm above the attachment; the deviation comes from the trailing-edge face closure below that.
- **Not met as written:** the 0.2 mm half-round closing the flat trailing-edge face, for the reason under build step 5.
