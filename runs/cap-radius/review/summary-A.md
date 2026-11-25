# Review summary A: wingb, swept55, flat (sharp and blunt)

Written by the coordinator from reviewer A's report; the harness blocked the reviewer from writing this file. Per-build detail is in the S-V.md files beside it.

All six STEP hashes match their README, readback.json and the audit. All six rebuilt bit-identically from their scripts in scratch copies (point deviation 0, volume difference 0, 2–78 s each). No build references earlier trials, docs/studies, git history or other worktrees.

| Build | Verdict | Whole-cap worst radius | Unmet and why |
|---|---|---|---|
| wingb sharp | accept with caveats | 0.438 mm (needs 0.4) | Join not tangent near the nose: up to 0.33° within 0.06–0.25 mm of it (builder error). Fails strict BRepCheck (builder error). Lip near the TE thickened 6–15 % (builder choice prompted by a skill gap). |
| wingb blunt | accept with caveats | 0.427 mm; 4e-9 mm at the TE corner | Below the rule for 1.7 mm at each TE corner (partly forced; the length is a builder choice). TE band thickened 8 % with an S-bend (builder choice, skill gap). Audit "closes False" is a tool error. |
| swept55 sharp | accept | 0.413 mm | Audit artefacts only. |
| swept55 blunt | accept with caveats | 0.412 mm; 0.067 mm at the TE corner | 82° crease round a 0.31 mm half-disc, ~0.5 mm below the rule at each corner (forced by the parent's corners; skill gap). Fails strict BRepCheck with no geometric effect (builder). Taper let-out squeezed into 3 mm (builder choice). |
| flat sharp | accept | 0.905 mm | Audit artefacts only. |
| flat blunt | accept with caveats | 0.743 mm; 0.356 mm at the closure | Lip 2.7 µm off the locked outline near the TE against a 1 µm tolerance (builder trade; skill gap). 2.85° crease inherited from the parent. |

## Builder claims against the audit

- Section-rule and TE-approach-cut failures are audit artefacts on all six: those cuts cross the lip 78–90° from square. The radius across the lip, perpendicular to the outline, passes everywhere.
- Strict validity: the audit is right on wingb sharp (cap-side edge curve 4.4e-4 mm from the true edge against a 4e-6 mm tolerance) and swept55 blunt (4.5e-6 against 1.7e-6). Both builders reported only the standard check.
- wingb sharp join: the audit (0.222°) is right; independent projection gives 0.33°. The builder's 0.032° sampled from about 0.13 mm from the nose and missed the worst zone.
- Audit tool error on wingb blunt: audit.py recognises the parent TE face only if it is planar or passes its flatness test. Wing B's blunt TE face is a curved ruled surface (normals spread 3.8°, 0.12 mm off a plane over 540 mm), so the tool files that join as a skin join and reports "closes False". The face is in fact closed and tangent.
- The audit judges radius against the parent's attachment thickness and never checks the cap's own thickness, so the two Wing B thickenings pass unnoticed.

## Recurring skill problems

1. The thin-zone rule cannot be met as written and has no tolerance. A lip of radius exactly t/2 is a circle, which can meet the skin only tangentially, not with matched curvature as build step 4 asks. Both Wing B builders thickened the section to exceed 1.0 × t/2; swept55 sharp reached 1.018 × t/2 without thickening.
2. The sharp corners of a blunt TE are never mentioned. Three blunt builds, three treatments: swept55 carried the crease round a half-disc (~0.5 mm below the rule); wingb smoothed it and tapered the face over 2 mm (1.7 mm below); flat kept the parent's 2.85° crease (none below). "The half-round closes the end of the flat TE face" was read both as a flat half-disc and as a tapering strip ending in a rounded edge.
3. The section checks assume the cut crosses the lip square. The prescribed outline runs tangent to the LE and TE, so literal sections fail good caps. Name the whole-cap curvature map and cuts across the lip as the measure.
4. No tolerances for join angle, zebra or validity, and nothing requires the cap edge to be the parent's own edge — how wingb sharp shipped a 0.33° join and two builds failed the strict check.
5. No outline tolerance where a blunt TE meets an outline tangent to it (flat blunt's 2.7 µm trade).
6. Smaller: "chord position" undefined on the cap; "full half-round" conflicts with "end as the parent ends" on flat's vertical LE/TE; the skill links a study builders were not allowed to read.

## Design quality beyond the radius numbers

- Clean: both flat caps and swept55 sharp (which holds full thickness at the aft lip and rolls late, a zebra hook the 0.4 mm rule makes necessary).
- Hidden by the radius numbers: wingb sharp has a bead at the lip near the TE with a concave band behind it; wingb blunt has the 8 % band and S-bend (stripe jogs near the TE); swept55 blunt's 3 mm taper let-out is concave to about 300 mm radius (its README says 850–910 mm); the blunt corner zones show wavy stripes.
- Against the previously accepted caps: swept55 sharp now passes the rule (accepted 0.0385 mm); flat sharp is on a par (0.905 vs 0.845 mm); wingb sharp passes (accepted 0.199 mm) but its join is worse than the accepted cap's 2.2e-4°.
