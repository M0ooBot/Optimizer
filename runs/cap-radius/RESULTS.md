# Cap radius trial: results so far (stopped early)

Test of the 0.4 mm printable-radius rule added to .agents/skills/close-rounded-wingtip/SKILL.md (commit d88407d6), on six shapes × two trailing edges (sharp, and blunt 0.4 mm via the Airfoil Converter's keep-chord cut). Each cap was built by a fresh agent given only the skill, a locked parent (parents/) and a short brief. Independent measurements: audit/ (tool, calibration, results/batch1.md and batch2.md). Reviews: review/ (reviewer A covered wingb, swept55 and flat; reviewer B for twisted and the fins was not run).

The trial was stopped at the user's request before the twisted blunt build finished; builds/twisted/blunt/ is incomplete and unreviewed.

## Whole-cap minimum radius (audit, angle-independent map)

| Shape | Sharp TE | Blunt TE (0.4 mm) |
|---|---|---|
| Wing B | 0.438 mm pass | 0.427 mm, but below the rule for 1.7 mm at each TE corner |
| 55° swept | 0.413 mm pass | 0.412 mm, but 0.067 mm within ~0.5 mm of each TE corner |
| Straight | 0.905 mm pass | 0.743 mm pass (0.356 mm at the TE closure, rule there 0.2) |
| Twisted | 0.412 mm pass | incomplete |
| Fin 40 mm | 0.577 mm pass | 0.577 mm pass (audit misreads the extended TE face; unreviewed) |
| Fin 5 mm | 0.668 mm pass | below the rule within 0.2 mm of each TE corner |

Previously accepted caps for comparison: swept55 0.039 mm, twisted 0.048 mm, Wing B 0.199 mm (all fail the rule); flat 0.845, fin40 0.577, fin5 0.668 (pass).

## Findings

1. Sharp trailing edges: the skill now produces caps that meet 0.4 mm on every shape tried (5/5 finished), including the two where the earlier accepted caps were far below.
2. Blunt trailing edges: the cap meets the rule everywhere except at the corners where the flat TE face meets the skins. The parent's crease there (18–84°) cannot be rounded to 0.2–0.4 mm within the cap, and the skill says nothing about it; every builder handled it differently.
3. Skill wording to fix:
   - Judge the rule on the surface (whole-cap principal curvature, and cuts square to the lip). Constant-chord or sweep-following cuts graze the lip near the LE/TE and read far too tight; every builder and the auditor hit this.
   - The thin-zone "full half-round of radius t/2" cannot coexist with curvature matching at the join; give a tolerance (e.g. ≥ 0.9 × t/2) and forbid thickening the section to reach it (both Wing B builders thickened the lip 6–15 %).
   - Specify the blunt-TE corner treatment and the outline tolerance where a blunt TE meets an outline tangent to it.
   - Add tolerances for join angle and strict BRepCheck validity (Wing B sharp shipped a 0.33° join near the nose; two builds fail the strict check).
   - "Spreading the roll over more of the reach" was wrong advice for the aft lip: meeting 0.4 mm there needed full thickness to 80–93 % of the reach and a late roll.
4. Audit tool limits found in review: it does not recognise a curved (ruled) blunt TE face (Wing B blunt), misreads a TE face extended into the cap (fin40 blunt), and does not check the cap's own thickness.
