# Review: wingb / sharp

**Verdict: accept with caveats.** The printable-radius rule is met. Two construction errors by the builder remain: the join is not G1 within about 0.25 mm of the nose, and the STEP fails strict BRepCheck. The builder also added a deliberate thickening bead of 6–15 % at the TE lip, which a better design could probably have avoided.

Governing numbers:
- **Whole-cap map:** worst 0.438 mm at x/c 0.9337, t = 0.800 mm (audit). The builder's refined value at the same spot is 0.4242 mm.
- **Join tangent angle** (my measurement): 0.33° at 0.059 mm from the nose on the upper side and 0.10° at 0.077 mm on the lower side. It is 0.03° from 0.25 to 1 mm, and still 0.014° on the lower side 3.1 mm from the nose. The tolerance is 0.01°.
- **Strict BRepCheck:** false.

## 1. Deliverables and claims
- All files are present. The SHA-256 `18881aab…836a20` is the same in the README, `readback.json`, the audit and the file itself. The volume is 731,173.83 mm³ in all of them.
- Two README claims go further than the STEP supports:
  - "BRepCheck valid" holds for the standard check only.
  - "Max tangent-plane angle 0.032° / 0.025°": the builder's join samples start at xi 3.55e-4, about 0.13 mm of arc from the nose, so they miss the worst zone.

## 2. Builder vs audit
| Item | Audit | Finding |
|---|---|---|
| Strict validity | FAIL | **Audit right.** On both cap–skin join edges, the cap-side pcurve is up to 7.8e-5 mm (upper) and 4.4e-4 mm (lower, next to the nose) from the 3D edge. The edge tolerances are 8.8e-6 and 4.1e-6 mm. The exact check therefore reports InvalidCurveOnSurface and InvalidSameParameterFlag, and BRepCheck marks the cap face UnorientableShape. The parent-side pcurves are exact (3e-13 mm). |
| Joins G1 | FAIL, 0.222° at x/c 1.4e-4 | **Audit right.** I projected 4,000 parent-edge points onto the cap surface and compared normals at the true foot points. This gives 0.33° (upper) and 0.10° (lower) beyond 0.05 mm from the nose, with a gap of at most 7.8e-5 mm. The accepted Wing B cap measures 2.2e-4°. |
| Section rule | FAIL | **Audit artefact.** 19 of the 24 failing cuts cross the lip more than 60° from square. The lip-normal radius (4b) is at least 0.459 mm in the thick zone, and at least 1.145 × t/2 in the thin zone. |
| Map | PASS 0.438 | Agrees with the builder's 0.4242. |
| TE closure | FAIL | **Audit artefact.** The cap closes at the TE point. Only the approach cuts fail, at 79–89° obliquity; their lip-normal radii are 1.27–1.83 × t/2. The principal radius within 0.1 mm of the TE point (1e-7 mm) comes from the collapsed sharp TE, where t is below 1.5 µm. That is forced. |
| Skin, outline, extent | PASS | Agree: 2.4e-13 mm, 1.2e-5 mm and 1e-7 mm. |

## 3. Classification
- **Join angle and strict invalidity: builder error.** The cap was interpolated through a 241 × 201 guide net and sewn at 1e-4 mm. Its boundary rows are not the parent's edge. Wing B blunt, swept55 sharp and flat sharp use the parent edge poles directly, and all three pass strict BRepCheck with join angles of 1e-6° or less.
- **TE-lip bead: builder error, prompted by a skill gap.** The builder argued that a G2 lip cannot reach R = t/2. That is true only for parallel sides, and the skill gives no tolerance. Swept55 sharp has the same kind of TE-tangent outline and gets at least 1.018 × t/2 with no thickening.
- **Thin zone within 2.2 mm of the TE point: forced** by the sharp TE.

## 4. Design quality
- **Bead (measured).** I compared the thickness along constant-chord-fraction lines with the attachment thickness:

  | Chord fraction | Thicker than the attachment by | Where |
  |---|---|---|
  | 0.92 | 5.9 % | 1.46 mm before the lip |
  | 0.94 | 8.1 % | 0.99 mm before the lip |
  | 0.96 | 8.8 % | 0.64 mm before the lip |
  | 0.98 | 14.6 % | 0.36 mm before the lip |
  | 0.88 and below | none | |

  Behind the bead is a concave band; the README gives its radius as 14 mm at p = 0.85 and 1.2 mm at p = 0.92. It shows as stripe jogs at the TE end of the lip in `zebra_top.png`, `zebra_bottom.png` and `zebra_tip-on.png`.
- **Thickness taper levelling (tanh).** At xi 0.7 the spanwise curvature rises from 3e-4 to 1.0e-3 /mm about 4 mm outboard of the join, on the upper side. That is smooth and mild.
- **Compared with the accepted Wing B cap.** The accepted cap re-scales the section into the outline: at xi 0.8 it is 1.41 × the attachment thickness. This cap keeps the parent section and rolls it at the outline, which departs far less from the parent.

## 5. Reproducibility
Rebuilt in a scratch copy (/tmp/claude-1000/review-A/wingb-sharp/, parents copied beside it). The build took 6 s. On the same (u, v) grid, every face of the rebuilt STEP matches the delivered one to 0 mm, and the volume difference is 0.000 mm³. Only the SHA-256 changes, because of the header timestamp.

## 6. Sources
No references to earlier trials, docs/studies, git history or other worktrees in the scripts, README or notes.
