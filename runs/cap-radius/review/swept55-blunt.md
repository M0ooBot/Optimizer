# Review: swept55 / blunt

**Verdict: accept with caveats.** Away from the TE corners the cap meets the rule with margin. Three caveats:
- An 82° crease runs round the 0.31 mm half-disc, and below t/2 there is a sub-rule zone at each TE corner, about 0.5 mm long. Both come from the parent's sharp TE-face corners.
- The taper let-out over the first 3 mm is short and visible in the zebra.
- The STEP fails strict BRepCheck, on tolerance bookkeeping only.

Governing numbers:
- **Thick-zone minimum:** 0.412 mm at x/c 0.980 (builder and the audit map).
- **TE-corner minimum:** 0.0665 mm, which is 0.333 × t/2, 0.04 mm from the corner.
- **Crease:** 81.9° at the corner, median 33.9° along the half-disc edge.
- **Strict BRepCheck:** false. The pcurve deviates by 4.5e-6 mm against an edge tolerance of 1.7e-6 mm.

## 1. Deliverables and claims
- All files are present. The SHA-256 `7f7e6fb6…7b17b7a` is the same in the README, `readback.json` and the audit. The volume is 1,239,793.36 mm³.
- The claims are consistent, with two exceptions:
  - "BRepCheck valid" holds for the standard check only.
  - The README describes the taper let-out bend as concave R ≈ 850–910 mm. I measure R ≈ 300 mm at xi 0.3 and 500 mm at xi 0.7, peaking 2.2 mm outboard of the join.

## 2. Builder vs audit
| Item | Audit | Finding |
|---|---|---|
| Strict validity | FAIL | **Audit right.** On the four roll-patch join edges the cap-side pcurve is up to 4.5e-6 mm from the 3D edge, against an edge tolerance of 1.7e-6 mm. The exact check flags InvalidCurveOnSurface, and BRepCheck marks the faces UnorientableShape. The physical gap is 3.6e-6 mm. |
| Joins | PASS: 3.6e-6 mm, 3.6e-5°; curvature mismatch 10.3 /mm | Both are right. The audit's 10.3 /mm sits 0.0006 mm from the TE corner (median 1e-7 /mm). The README's 0.18 /mm excludes the corner. |
| Map | FAIL, 0.0665 mm at x/c 1 | **Correct; the README says 0.067 mm.** My radius against distance outboard of the corner: 0.067 (0–0.05 mm), 0.078, 0.083 (to 0.25), 0.164 (0.25–0.5), 0.200 (0.5–0.75), then at least 0.213 mm. |
| TE closure | FAIL | The cap closes: the TE face meets the coplanar half-disc at 0°. The failure comes from the approach cuts at 86–89.5° obliquity, an artefact; their lip-normal radii are 0.227–0.33 mm, at least t/2. The 81.9° crease is real; see section 3. |
| Section rule | FAIL | **Audit artefact (obliquity).** Line 4b gives at least 0.430 mm in the thick zone and at least 1.075 × t/2 in the thin zone. |
| Skin, outline, extent | PASS | Agree: 1.2e-7 mm and 1.5e-9 mm. |

## 3. Classification
- **TE-corner crease and the sub-rule zone of about 0.5 mm.** Forced by the locked parent: its TE face meets the skins at an 84° crease. The handling is a skill gap, because the skill never mentions the TE-face corners. "Half-round closes the end of the flat TE face" is here read as an in-plane half-disc of radius 0.2 mm. Wing B blunt reads the same sentence as a tapering strip.
- **Strict invalidity: builder error.** The faces were sewn at 1e-5 mm without restoring exact edge tolerances. It has no geometric effect.
- **Short taper let-out: builder choice.** The sharp sibling eases its taper with no concave bend at xi 0.3.

## 4. Design quality
- **Taper let-out.** A concave band reaching R ≈ 300 mm lies in the first 3 mm outboard of the join. It shows as a stripe bend in `tip_top_zebra`. The cap is not thicker than the attachment anywhere at x/c 0.3–0.99.
- **TE corner.** Within about 1 mm of the corner the stripes are wavy and break into islands (`te_corner_top_zebra`, `te_corner_oblique_zebra`). The station cut at Y = 390.5 shows a concave R of 2.5 mm 0.5 mm from the corner (README). This is a local curvature ripple where the corner fades out.
- **Lip.** In lip-normal cuts the roll is otherwise convex. The lip radius is 0.60 mm over the thick tip and 0.23–0.24 mm near the TE end.

## 5. Reproducibility
Rebuilt in a scratch copy (/tmp/claude-1000/review-A/swept55-blunt/, parents copied beside it). The build took 23 s. On the same (u, v) grid, every face of the rebuilt STEP matches the delivered one to 0 mm, and the volume difference is 0.000 mm³. Only the SHA-256 changes, because of the header timestamp.

## 6. Sources
None found.
