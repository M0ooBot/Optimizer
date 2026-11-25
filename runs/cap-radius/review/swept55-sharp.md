# Review: swept55 / sharp

**Verdict: accept.** Every requirement the audit judges in an angle-independent way passes. The section-rule and TE-closure failures are audit artefacts.

Governing numbers:
- **Map:** worst 0.4133 mm (margin +0.013 mm) at x/c 0.9604, t = 0.800. The builder's value is 0.4148 at ξ 0.960.
- **Joins:** 8e-12 mm and 7.6e-8°.
- **Strict BRepCheck:** passes.

## 1. Deliverables and claims
- All files are present. The SHA-256 `5d511033…6f2420` is the same in the README, `readback.json` and the audit. The volume is 1,208,744.50 mm³.
- The README matches the readback and my checks.

## 2. Builder vs audit
| Item | Audit | Finding |
|---|---|---|
| Validity, extent, skin, outline | PASS | Agree: strict valid; 1.9e-7 mm (the parent's own fit); 1.2e-6 mm. |
| Joins | PASS, curvature mismatch 1e-7 /mm | Agree. |
| Map | PASS 0.4133 | Agrees with 0.4148. |
| Section rule | FAIL | **Audit artefact (obliquity).** Line 4b is at least 0.4695 mm in the thick zone and at least 1.032 × t/2 in the thin zone. The builder's 64 sweep-following and 34 lip-normal cuts judge the surface radius and pass. |
| TE closure | FAIL | **Audit artefact.** The cap closes at the TE point, and the principal radius within 0.1 mm of it is at least 0.91 mm. Only the approach cuts fail, at 86–89° obliquity. |

## 3. Classification
Nothing unmet except the audit artefacts. The thick-zone margin is thin (+0.013 mm at the 0.8 mm crossing) but positive.

## 4. Design quality
- **Aft lip.** The aft half of the cap holds the parent's full thickness (taper eased from q = 0.6) and rolls over only the last ~7 % of the reach. The zebra shows this as a hook near the aft lip (`tip_top_zebra`, `tip_aft_lip_top_zebra`, `tip_rear_oblique_zebra`). It is a consequence of the rule, not an error: where t is about 0.8 mm, a 0.4 mm radius needs a full-thickness half-round at the lip. Thickness along constant-chord-fraction lines never exceeds the attachment value (no bead). There is no concave bend at xi 0.3 in the first 12 mm.
- **Reversals.** They are slight: radius at least 1.35 m in the flat run. The roll-start fit ripple is at most 1.1e-3 × the local lip curvature. The saddle shape (31 % of cap points with K < 0) is inherited from the parent's ruled skin.
- **Compared with the accepted swept55 cap.** The accepted cap's worst radius is 0.0385 mm at the same x/c, so it fails the new rule. Join quality is equal.

## 5. Reproducibility
Rebuilt in a scratch copy (/tmp/claude-1000/review-A/swept55-sharp/, parents copied beside it). The build took 72 s. On the same (u, v) grid, every face of the rebuilt STEP matches the delivered one to 0 mm, and the volume difference is 0.000 mm³. I also ran design.py (45 s) first; it regenerated design_table.json byte-identical. Only the SHA-256 changes, because of the header timestamp.

## 6. Sources
None found.
