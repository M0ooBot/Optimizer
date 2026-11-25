# Review: flat / blunt

**Verdict: accept with caveats.** The radius rule is met everywhere. The one judged failure is the locked outline: near the TE the lip sits up to 2.7 µm off it, against a 1 µm tolerance. The builder chose this as a trade.

Governing numbers:
- **Thick-zone minimum:** 0.743 mm at x/c 0.928.
- **Tightest point:** 0.356 mm at the edge of the TE closure, where t = 0.4 mm (1.78 × t/2).
- **Outline-to-surface:** 2.68e-3 mm (audit), 2.7e-3 mm (builder), both at X above 119.95 mm.
- **Strict BRepCheck:** passes.

## 1. Deliverables and claims
- All files are present. The SHA-256 `fdf987b5…6a0fd5` is the same in the README, `readback.json` and the audit. The volume is 958,092.90 mm³.
- The claims are consistent, and the README states the outline deviation openly.

## 2. Builder vs audit
| Item | Audit | Finding |
|---|---|---|
| Validity, extent, skin, joins | PASS | Agree. The joins are G2: 1.1e-11 mm and 1.3e-7°. |
| Locked outline | FAIL, outline-to-surface 2.68e-3 mm | **Audit right; it is a real deviation.** The constant-X cuts near the TE show the lip up to 0.99 mm off in Y. That is only because the outline and the lip run almost parallel to Y there. The true distance between them is 2.7 µm. |
| Map | PASS 0.743 / 0.356 | Agree. |
| TE closure | FAIL | It closes: the TE face joins the planar closure with a gap of 0 mm and an angle of 0°. The failure comes from the approach cuts at 80–89° obliquity, an artefact; their lip-normal radius is at least 1.1 mm. The closure–cap seam has a crease of 2.85°, the same as the parent's own skin-to-TE-face crease, fading to 0.10° at the apex. |
| Section rule | FAIL | **Audit artefact.** Line 4b is at least 0.756 mm in the thick zone and at least 18.7 × t/2 in the thin zone. |

## 3. Classification
- **Outline deviation: builder choice, made under a skill gap.**
  - The outline approaches the TE as 1200(1 − p)³ mm in front of X = 120, so it hugs the TE plane.
  - In the builder's own scan, a 1e-3 mm shift gives a 1.4 µm deviation at 1.17 × the rule, and a 5e-4 mm shift fails at 0.77. A setting that stays within 1 µm and still passes 0.9 very likely lies between these. The builder took 2.7 µm to get a margin of 1.78.
  - The skill gives no outline tolerance where a blunt TE meets a TE-tangent outline.
  - The geometric effect is nil.
- **Thin-zone "full half-round" not literal: skill wording.** The thin zone is the last 0.023 mm of chord. There the cap continues the parent's almost vertical TE face, so the lip is flatter than t/2, which is a larger radius. The closure is a rounded planar end: apex radius 0.200 mm, shoulder radius 0.1725 mm, reach 0.352 mm. It is not an exact semicircle.
- **2.85° crease: forced.** It is inherited from the parent.

## 4. Design quality
The zebra renders are clean, with no reversals and no bulges. The TE close-up shows only the inherited crease. The cap is not thicker than the parent anywhere.

## 5. Reproducibility
Rebuilt in a scratch copy (/tmp/claude-1000/review-A/flat-blunt/, parents copied beside it). The build took 2 s. On the same (u, v) grid, every face of the rebuilt STEP matches the delivered one to 0 mm, and the volume difference is 0.000 mm³. Only the SHA-256 changes, because of the header timestamp.

## 6. Sources
None found. The skill notes say the builder did not read the linked Wing B study.
