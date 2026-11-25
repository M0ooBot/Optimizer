# swept55 / blunt — rounded tip caps

**Model:** `wing.step` is one closed, valid solid. It is the full symmetric wing: the parent, untouched, from Y = -390 to +390, plus a rounded cap at each tip. The span is 900 mm.
AP214, mm, 3.40 MB, SHA-256 `7f7e6fb63931f08bb8c52ca41814bcfd0e18149c1a56afec1cc75172e7b17b7a`.

The skill used is `.agents/skills/close-rounded-wingtip/SKILL.md`. The parent is `runs/cap-radius/parents/swept55/blunt/`.

## What was built

- **12 faces.**
  - The 6 parent B-spline faces (the upper and lower skins and the flat TE face, each split at Y = 0) are copied from `parent.step` unchanged; their poles are identical. The two flat end faces at Y = ±390 are removed.
  - Each tip has three faces: an upper roll patch, a lower roll patch and a planar half-round face.
  - The roll patches are quintic B-splines with 200 × 40 poles. The lower patch and the left tip are exact mirrors of the right upper patch.
- **Profiles.** The boundary-surface profiles are:
  - the parent's upper cut edge at Y = 390,
  - the locked outline (the lip, at Z = 0),
  - the lower cut edge.

  The upper and lower patches share the lip edge. All lip poles lie at Z = 0, and the next pole row sits straight above them, so the two patches share a vertical tangent plane there: G1 by construction and G2 by symmetry.
- **Guides.** Spanwise guides follow the sweep: they run along lines of constant chord fraction of the continued planform. Near the trailing edge they fan forward, so that the lip, which runs almost along the TE line there, gets a real half-round.
- **Trailing-edge end.** The flat TE face (0.4 mm) is closed by a coplanar half-disc of radius 0.2 mm, which continues the TE plane outboard. The cap patches meet its semicircular edge. The parent's 84° TE-corner crease continues around that edge and fades to 0° at the apex, 0.2 mm outboard along the TE line, where the lip begins.
- **Height function.** In the cap,

  `z = Zp(P) − H_lip · Δ(τ)`

  - `Zp` is the parent's own half-thickness, continued outboard. The chord taper is let out over 3 mm.
  - `τ` is a roll coordinate: 0 at the lip, ≥ 1 in a flat run next to the cut.
  - `Δ` is a stretched circular arc. It is vertical at the lip and tangent to the parent's slope at the roll top.

  As a result the cap is exactly the parent surface next to the cut (G2), the roll is convex, and the curvature rises progressively to the lip.
- **Lip-radius targets.**
  - 0.6 mm where the section is thick.
  - 1.08 × the local half-thickness where it is thinner. Near the TE this is multiplied by (1 + μ), using the forward-thickening wedge, which gives about 1.17 × t/2.
  - Near the leading edge, the lip continues the parent's nose, with a radius of about 1.9 mm.
- **Construction method.** The analytic surface (`cap_model.py`) is sampled and fitted with constrained least squares (`build_cap.py`), which also enforces the parent's first and second cross-derivatives at the cut. The faces are sewn with a 1e-5 mm tolerance and exported. Fit error to the analytic surface: max 0.97 µm, rms 0.04 µm.

## Results (measured on the reimported STEP, `verify.py` → `readback.json`)

| Check | Result |
|---|---|
| Solid / shells / faces | 1 / 1 / 12. BRepCheck valid. Every edge has 2 faces (4 degenerate edges at the collapsed LE corners). Volume 1 239 793.36 mm³ (parent 1 191 108.33). |
| Parent preserved | 6/6 faces have poles identical to `parent.step`. Max distance to `parent.py`: 1.25e-7 mm. |
| Locked dimensions | Span 900.0000 mm. Attachments at Y = ±390.0000. Lip Z = 0. Lip vs locked outline ≤ 0.16 µm in plan, all of it in the last 0.5 mm before the TE end, where the lip is nudged onto the TE plane. |
| Cap/parent joins (4) | Gap ≤ 1.0e-7 mm. Tangent-plane angle ≤ 3.6e-5° (≤ 5e-4° at the collapsed LE point). Curvature difference ≤ 0.18 mm⁻¹, worst 0.21 mm from the TE corner; far from the corners it is at fit-noise level (see zebra). |
| Lip seam (upper/lower) | Angle ≤ 3e-4°. Zebra phase continuous. Cross-gradient jump p95 0.10 °/mm against a median gradient of 59 °/mm. |
| Half-disc / parent TE face | Coplanar (0.000°). |
| Crease along the half-round (designed) | 80.9° at the TE corner, falling to 1.2° next to the apex. |

### Printable radius rule (larger principal curvature, 900 × 260 grid plus Nelder–Mead refinement)

| Zone (attachment thickness at the point's chord fraction) | Requirement | Measured |
|---|---|---|
| t ≥ 0.8 mm | every radius ≥ 0.4 | **min 0.412 mm** at (671.04, 408.40, 0.106), chord fraction 0.97996, t = 0.800. Lip-normal sections give lip radii of 0.60–0.61 over the tip and 1.8–1.9 near the LE. |
| LE, t < 0.8 (chord fraction < 0.001) | lip half-round with R = t/2 | lip R 1.92 mm, which is ≥ 4.8 × t/2 (the parent's nose). Interior ≥ 1.34 mm beyond 0.05 mm of the LE point, where the patch is collapsed (excluded; the grid reaches 1e-12 mm of it). |
| TE, t < 0.8 (chord fraction > 0.98) | lip half-round with R = t/2 | lip R / (t/2) = 1.07–1.17 more than 0.1 mm beyond the apex, and 1.3–2.0 at 0.05–0.1 mm. Within 0.05 mm of the apex the lip is flat, where the TE face fades into it. |
| TE, t < 0.8, interior | (not stated; R = t/2 = 0.2 taken as the reference) | **min 0.067 mm**, on the half-round edge 0.04 mm from the TE corner. 0.083 beyond 0.1 mm, 0.107 beyond 0.2 mm, 0.153 beyond 0.3 mm, 0.186 beyond 0.5 mm, and ≥ 0.210 (≥ t/2) beyond 1.0 mm from the TE corner. |

**Below t/2:** the only radii below t/2 are in a region about 1 mm across at each TE corner. There the cap has to dissolve the parent's 84° TE-corner crease along the 0.31 mm half-round. See `skill-notes.md`.

### Sections (`readback.json` → `sections`)

The cut sets are:

- 20 lip-normal cuts (p = 0.01 … 0.998),
- 11 station cuts (Y = 390.5 … 449.5),
- 17 cuts along constant chord fraction (sweep-following).

In the lip-normal cuts the roll is convex from the parent to the lip; the only "reversals" are at R ≥ 29 mm, which is numerical level. Lip-normal lip radii, as p runs from 0.01 to 0.998 (LE end to TE end):

1.93, 1.91, 1.88, 1.82, 1.54, 1.14, 0.76, 0.61, 0.60, 0.60, 0.60, 0.60, 0.51, 0.37, 0.29, 0.24, 0.23, 0.23, 0.24, 0.24 mm

The deliberate taper let-out shows as a very gentle concave bend along the sweep in the first 3 mm of the cap (R ≈ 850–910 mm). The station cut at Y = 390.5 crosses the TE-corner transition: concave R = 2.5 mm, 0.5 mm from the corner.

The sweep-following and station cuts meet the lip obliquely near the LE and TE corners, so their in-plane lip radii (down to 0.02 mm) are not surface radii; the lip-normal cuts and the principal-curvature map are the measures used.

### Numerical zebra (CAD normals of the reimported faces)

Setup: 5 view directions × 4 stripe axes, 5° stripes, `r = d − 2(d·n)n`, phase `asin(r·e)`.

- **Cap/parent joins.**
  - Phase jump ≤ 6.5e-5° everywhere.
  - Cross-join gradient jump, against a median gradient of 0.28 °/mm:
    - ≤ 3.5e-4 °/mm more than 3 mm from the corners,
    - ≤ 1.5 °/mm at 1–3 mm and ≤ 2.9 °/mm at 0.5–1 mm from the TE corner,
    - 24.5 °/mm maximum, 0.26 mm from the TE corner.
  - Stripe direction change at stripe boundaries on the join: ≤ 0.02° more than 3 mm from the corners, up to 138° within 1.5 mm of the TE corner, where the stripes kink.
- **Lip seam.** Max gradient jump 1.1 °/mm against a 59 °/mm median, and direction change ≤ 0.014°.

### Renders (`figures/`, from the STEP poles; review in context)

- `wing_top_zebra*.png`: the whole wing, centred camera, mirrored stripe lighting. The two tips match.
- `tip_top_zebra*.png`, `tip_bottom_zebra.png` and `tip_oblique_aft_zebra.png`: the right tip.
  - Stripes cross the join without a jump.
  - Within about 3 mm outboard of the join they bend. This is the deliberate taper let-out: the chord taper stops over 3 mm, so the section keeps its thickness for the lip. The join itself is G2, but the let-out is short, so the easing is visible in the stripes.
  - The parent area in these close-ups is sampled coarsely and looks hatched. That is a render artifact.
- `te_corner_*_zebra.png`: a 2 mm field at the TE corner, with the half-round face (flat, white). Stripes converge on the half-round edge (the crease) and are wavy within about 1 mm of the corner; this is the zone where the radius map and the numerical zebra flag the corner. Fine ring patterns in the dark bands are aliasing.

## Files

- **`cap_model.py`**: the analytic cap (guides, TE-corner field, height function, design parameters).
- **`build_cap.py`**: fits, sews and exports `wing.step`. Also writes `build_report.json` and `cap_patch.npz` (the poles).
- **`bspline_fit.py`**: the B-spline basis and the constrained least-squares fit.
- **`verify.py`**: every check above, on the reimported STEP. Writes `readback.json` (about 11 min) and `radius_map_R_upper.npz`.
- **`radius_map.py`**: the principal-radius grid and the rule judgement (used by `verify.py`).
- **`render.py`**: `figures/*.png`, zebra renders of the STEP. Orthographic, no edges, no overlays, mirrored stripe lighting.
- **`design_scan.py`**: the design-loop radius scan on the analytic surface (finite differences).

Run from the project root with `.venv/bin/python <file>`, in this order: `build_cap.py`, `verify.py`, `render.py`.
