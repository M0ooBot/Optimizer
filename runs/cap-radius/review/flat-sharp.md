# Review: flat / sharp

**Verdict: accept.**

Governing numbers:
- **Map:** worst 0.905 mm at x/c 0.999, t = 1.355 mm. The margin is +0.505 mm, and the builder's figure agrees.
- **Joins:** exactly G2. The gap is 0 mm, the angle 2e-14° and the curvature mismatch 3e-16 /mm.
- **Strict BRepCheck:** passes.

## 1. Deliverables and claims
- All files are present. The SHA-256 `d32211f9…5671a6` is the same in the README, `readback.json` and the audit. The volume is 962,251.54 mm³.
- The claims are consistent.

## 2. Builder vs audit
| Item | Audit | Finding |
|---|---|---|
| Validity, extent, skin, outline, joins | PASS | Agree. |
| Map | PASS 0.905 | Agree. |
| Section rule | FAIL (constant-X cuts near the corners) | **Audit artefact.** Cuts there run 78–90° from square, so the planes are almost tangent to the surface. Line 4b is at least 1.07 mm in the thick zone and at least 12.65 × t/2 in the thin zone. The builder reports the same effect. |
| TE closure | FAIL | **Audit artefact.** The cap closes at the TE point, and the principal radius within 0.1 mm of it is at least 4.56 mm. The approach cuts fail at 80–89° obliquity; their lip-normal radii are 1.6–5.5 mm, at least 1.1 × t/2. |

## 3. Classification
Nothing unmet beyond the audit artefacts. The literal "half-round of radius t/2" in the thin zones (the first 0.002 mm and last 0.026 mm of chord) is not built. There the cap continues the parent's own flat cubic closure, and its radii are at least 2.26 × t/2. This is a skill wording problem: the rule assumes a round-nosed thin end, and this parent's LE and TE are vertical.

## 4. Design quality
- **Wraps.** Each wrap follows η⁵ + ζ² = 1, with a flat crown and a progressive roll. The lip radius is 1.8 mm at mid-chord. There are no reversals.
- **Zebra.** The renders are clean, with no jogs.
- **Compared with the accepted flat cap.** The accepted cap's worst radius is 0.845 mm; this one's is 0.905 mm. The look is similar.
- No local thickening. No departure from the parent's section.

## 5. Reproducibility
Rebuilt in a scratch copy (/tmp/claude-1000/review-A/flat-sharp/, parents copied beside it). The build took 2 s. On the same (u, v) grid, every face of the rebuilt STEP matches the delivered one to 0 mm, and the volume difference is 0.000 mm³. Only the SHA-256 changes, because of the header timestamp.

## 6. Sources
None found.
