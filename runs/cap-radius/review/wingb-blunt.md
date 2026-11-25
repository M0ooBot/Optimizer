# Review: wingb / blunt

**Verdict: accept with caveats.** The cap is strictly valid and G1 at every join. The thick-zone rule is met with a worst value of 0.427 mm. Three things count against it: the rule fails over about 1.7 mm at each TE corner, the TE band is thickened by 8 %, and the mean line is levelled with an S-bend. The audit's "TE closure: closes False" is wrong for this model.

Governing numbers:
- **Thick-zone worst:** 0.427 mm (margin +0.027 mm) at x/c 0.9723.
- **TE corner (my measurement):** the radius is below 0.9 × t/2 for 1.6–1.75 mm outboard. In the first 0.5 mm it is at most 5e-5 mm, which is effectively a crease.
- **TE band thickening:** +7.1 to +8.1 % over x/c 0.92 and aft.
- **Joins:** 0.0067°.

## 1. Deliverables and claims
- All files are present. The SHA-256 `2312ee5e…891b4` is the same in the README, `readback.json` and the audit. The volume is 769,329.87 mm³.
- The claims are consistent with each other and with the STEP. One thing is not reported in the README: the build log shows the patch fit misses the builder's own analytic design by up to 2.3e-3 mm in Y near the nose end of the perimeter (s = 0.009). This is not a requirement failure.

## 2. Builder vs audit
| Item | Audit | Finding |
|---|---|---|
| Validity | PASS (strict) | Agree. |
| Joins and cap seams | PASS (0.0067°; seams 1.1e-6°) | Agree. The seam curvature mismatch of 1.9e8 /mm is at the TE-corner collapse. |
| Map | FAIL, 4.0e-9 mm at 0.05 mm from the TE corner | **Correct, and the README says the same.** My radius against distance outboard of the corner: 3e-13 mm (0–0.05 mm), 5e-5 (to 0.5 mm), 0.0023 (0.5–0.75), 0.017 (0.75–1.0), 0.056 (1.0–1.25), 0.115 (1.25–1.5), 0.173 (1.5–1.75), then at least 0.202 mm. The README gives 1.71 mm (upper) and 1.78 mm (lower). |
| TE closure | FAIL, "closes False" | **Audit tool error.** See below. Separately, the approach cuts fail only because of obliquity (79–89°); their lip-normal radii are 0.234–0.33 mm, which is at least t/2. |
| Section rule | FAIL | **Audit artefact (obliquity).** Line 4b gives at least 0.477 mm in the thick zone and at least 1.181 × t/2 in the thin zone. |
| Skin, outline, extent | PASS | Agree. |

**What the audit tool gets wrong.** `audit.py` recognises the parent TE face only if it is planar: the face must be of type `GeomAbs_Plane` or pass `_flat()`. On Wing B the blunt TE face is a ruled B-spline that follows the bowed TE. Its normals spread over 3.8° and it is up to 0.12 mm off its best-fit plane along its 540 mm.
- The TE-face-to-strip join is therefore filed as "parent skin to cap". This is audit edge faces [1, 5], 0.39999 mm long, with a gap of 2.0e-12 mm and an angle of 6.5e-6°.
- `te_face_joins` comes out empty, so `closes = False`.
- The builder's 2 mm strip is planar to 0.013°, but it also fails `_flat()`. So "seams with planar cap faces" is empty and the strip is counted as a curved cap face.

In fact the TE face is closed by a G1 strip.

## 3. Classification
- **TE-corner sub-rule zone.** Partly forced: the parent's 0.4 mm TE face has sharp corners, and a G1 cap must dissolve them. The 1.7 mm length is the builder's choice: swept55 blunt brings the radius above t/2 within about 0.5 mm, keeping an inherited crease. The skill does not say how to handle a creased TE-face corner (skill gap).
- **8 % TE thickening: builder choice, driven by a skill gap.** The builder's own lip profile reaches 0.989 × h, which already meets the audit's 0.9 gate without thickening. The skill states no tolerance and does not say whether the section may be thickened. Swept55 blunt shows no thickening.
- **Approach cuts and section rule:** audit artefacts.

## 4. Design quality
- **Thickening (measured).** I compared the thickness along constant-chord-fraction lines with the attachment thickness:

  | Chord fraction | Thicker by | Where |
  |---|---|---|
  | 0.92 | 7.1 % | peak 11.5 mm outboard |
  | 0.94–0.99 | 7.9–8.1 % | peaks 1.6–2.6 mm before the lip |
  | 0.88 and below | none | |

  The cap matches the parent exactly at xi 0.7 and below. At xi 0.9–0.95 the spanwise curvature swings between -2.7e-3 and +3.4e-3 /mm over the first 6 mm, where the thickening blends in.
- **Mean-line levelling S-bend.** The README gives concave R of 1.03 mm or more, between 1.2ℓ and 4ℓ inside the lip, for p ≥ 0.8. The audit counts 15–47 reversals on the TE-side cuts.
- **Zebra.** There are stripe jogs near the TE, 2–6 mm outboard (`zebra_tip_top`, `zebra_tip_bottom`, `zebra_te_corner`), and a small stripe step at the nose end of the lip (`zebra_tip_oblique`).
- **TE end.** The parent's flat TE face carries on as a ruled strip that tapers to a point 2 mm outboard. This is one reading of "the half-round closes the end of the flat TE face".

## 5. Reproducibility
Rebuilt in a scratch copy (/tmp/claude-1000/review-A/wingb-blunt/, parents copied beside it). The build took 78 s. On the same (u, v) grid, every face of the rebuilt STEP matches the delivered one to 0 mm, and the volume difference is 0.000 mm³. Only the SHA-256 changes, because of the header timestamp.

## 6. Sources
None found.
