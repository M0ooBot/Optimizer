# Skill notes: close-rounded-wingtip on the flat wing with a sharp TE

This is for whoever maintains the skill. Work and numbers are in `README.md` and `readback.json`.

## What I followed, and how

- **Build 1 (attachment curves).** I used the parent's actual cut edges: the Y = ±350 pole columns of
  the skin B-splines in `parent.step`. The cap's first and last pole rows are those poles, copied
  bit for bit. Span, attachment and outline are unchanged.
- **Build 2 (3D perimeter).** The perimeter is the locked outline. It is already tangent to the LE and TE
  in plan and lies in Z = 0. On a flat, untwisted wing the elevation check is trivial.
- **Build 3 (wrap guides).** Every iso-u curve of the cap is a guide in the plane X = X(u), and it crosses
  all three profiles in order. I did not place a few discrete guides; there is effectively one per chord
  station.
- **Build 4 (boundary surface).** OCCT has no boundary-surface-with-guides feature that matches the
  skill's description. I wrote the B-spline surface directly: profile rows are upper edge, outline,
  lower edge, and the guides are its iso-curves. The surface is one patch through the lip, so there are
  no separate upper and lower patches to match. G2 to both skins is built into the poles: the first three
  wrap rows differ only in Y.
- **Build 5 (printable radius).** I designed for it from the start, but the first design still failed (see
  Iterations).
- **Build 6 (progressive turn).** In every guide the in-plane curvature rises monotonically from the
  attachment into the lip, and the lip is the tightest point of each guide.
- **Verify 1–5.** All checks run on the reimported STEP:
  - topology, volume and extents;
  - skins and outline compared with `parent.py`;
  - separation, tangent angle and curvature in a common frame along all four joins;
  - X = const sections, plus sections normal to the lip;
  - a principal-curvature map with refinement in OCCT;
  - numerical zebra over 7 views × 5 stripe axes;
  - rendered zebra images (centred and mirrored cameras, no edges).

## Where the skill was unclear, ambiguous or wrong

1. **Which radius the printable rule judges.** Build 5 and verify 2 talk about "the turn through the
   lip" in a section "at the same chord position". That assumes the section crosses the lip roughly at
   right angles. Here the outline is tangent to the LE and TE, which build 2 requires. So near both corners
   every chord-normal plane is almost tangent to the surface at the lip: the cosine between the plane and
   the surface normal falls to 0.002. The in-plane radius of such a section measures the grazing angle,
   not the bend.
   - Iteration 1 met the rule in every X = const section (lip radius 0.958 h) but had a 0.122 mm principal
     radius across the lip.
   - Iteration 2 is the reverse: the principal map and the across-lip normal sections pass with ≥ 0.905 mm,
     while X = const in-plane radii fall to 0.025 mm near the corners.
   - I treated the principal map (verify 3) and normal sections across the lip as authoritative, and report
     the X = const in-plane numbers separately. The skill should say which one counts. It should also say
     to cut normal to the outline wherever the outline is oblique to the guide planes. The existing "on a
     swept wing, cut along the sweep" hint is the same idea, but it does not cover outlines that become
     tangent to the LE or TE.
2. **The 0.8 mm threshold cannot be met literally in a G2 chord-normal section.** For a turn that drops
   h = t/2 while every radius stays ≥ 0.4 mm = h, the turn must be an exact quarter circle starting at the
   attachment. That is only G1. Iteration 1 measured 0.383 mm at t = 0.801 mm. The rule needs either a
   tolerance band or a statement that it applies to principal radii.
3. **"Half-round whose radius is half the local thickness" versus "end the cap the way the parent's TE
   ends".**
   - This parent's sharp TE (and its LE) closes with a vertical tangent and zero curvature: X ∝ Z³. It is
     not a half-round.
   - Here the thin zone (t < 0.8 mm) is only the first 0.002 mm and the last 0.026 mm of the chord.
   - I followed "end it the same way": near the corners the cap reproduces the parent's own closure along
     the outline.
   - Principal radii in the thin zone are at least 2.26 × t/2 (the global minimum ratio on the map), and
     6.5–13.8 mm at the sections I tested there. But the across-lip lip there is not a half-round of
     radius t/2; it is flatter.
   - I count this as met in intent but not to the letter.
4. **Step 6's failure mode is described only for the guide.** Iteration 1 had no such shoulder in any guide.
   Its tight spot was across the lip, near the corners. A warning is needed: "when the outline approaches
   the LE/TE tangentially, a roll packed into a short band of the guide is compressed across the lip by the
   sine of the outline's angle to the guide plane."
5. **Zebra crossings (verify 4).** Read literally, the step compares stripe direction at every crossing.
   That gave false 175–180° "direction changes" on back-facing samples and where the stripe gradient
   vanishes (lighting singularities). I restricted the check to front-facing samples, with a gradient
   floor of 1e-3 stripe/mm. The skill should say so.
6. Verify 1 asks for curvature "in a common tangent frame" without defining one. I used e1 = the cut edge
   tangent and e2 = n × e1, and compared normal curvatures k(e1), k(e2) and the twist.

## Decisions the skill did not cover

- Guide family. The skill gives no family. I chose η⁵ + ζ² = 1 in units of (reach, half thickness):
  - **G2 start.** ζ = 1 − η⁵/2 at the attachment.
  - **Round lip.** In-plane lip radius 5h²/(2r).
  - **Corners.** Near them the cap reduces to the parent's LE/TE closure moved along the outline.
  - **Exponent.** Swept on the exact surface: P = 3 gives 0.63 mm minimum radius, 4 gives 0.85, 5 gives
    0.905, 6 gives 0.85. I chose 5.
- Guides are planar at X = const. That makes the chord position of every cap point well defined for the
  thickness rule.
- One cap patch through the lip rather than two.
- A sewing tolerance of 1e-6 mm, and degenerate edges at the collapsed LE/TE corners.

## Iterations

Two built designs; the first was rejected by the verification.

1. **Flat crown plus a curvature-ramped, near-constant-radius roll in the outer 1.7 h of each guide.**
   - Before this I spent about four design-only attempts at fitting it as a separable B-spline; the
     parameterisation problems were mine, not the skill's.
   - The build was valid and exactly G2.
   - Principal map: 0.122 mm at X = 119.89, |Y| = 357.96, attachment thickness 1.31 mm. **Fail.**
2. **η⁵ + ζ² = 1** (delivered). The smallest principal radius is 0.905 mm.

## Requirements: met, partly met, not met

| requirement | status |
|---|---|
| closed valid single solid, parent skin preserved, span/attachment/outline locked | met: BRepCheck valid, 0 free edges; skins bit-identical; outline within 7e-11 mm |
| tangent to both skins at the attachment | met: ≤ 1.5e-6° |
| curvature continuity at the attachment (smooth reflection) | met: Δk ≤ 1.2e-7 /mm; zebra phase ≤ 7e-11 stripe; direction change ≤ 7e-9° |
| perimeter tangent to LE and TE, upper and lower share one tangent plane at the lip | met: one patch; lip normal horizontal to 5e-16 |
| finite lip radius, progressive turn, no reversal | met in every guide and across-lip section. The lip radius is finite everywhere except at the collapsed LE/TE points. Across-lip lip radius is ≥ 1.06 mm, and grows very large towards the corners. |
| printable radius ≥ 0.4 mm where t ≥ 0.8 mm | **met by principal curvature** (min 0.905 mm, ratio 2.26) and by across-lip normal sections (min 0.952 mm). **Not met in X = const in-plane terms near the corners** (down to 0.077 mm at chord fraction 0.99975, where the plane meets the surface at cos 0.016). |
| thin zone: lip a half-round of radius t/2 | **partly met**: every principal radius is ≥ 2.26 × t/2 (6.5–13.8 mm at the tested sections), but the shape is the parent's flat cubic closure rather than a half-round |
| zebra review in context | done: `zebra_*.png`, nothing flagged |
