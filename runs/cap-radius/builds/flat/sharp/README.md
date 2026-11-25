# Flat wing, sharp trailing edge: rounded tip caps

Straight 800 mm wing, 120 mm chord, symmetric polynomial section, attachments at Y = ±350 mm,
locked plan outline reaching 50 mm. `wing.step` is the whole wing (parent skins plus a cap at each
tip) as one closed solid.

Files (run from the project root with `.venv/bin/python`):

| file | what |
|---|---|
| `wing.step` | the model, 1 solid, 4 B-spline faces, 232 522 bytes, SHA-256 `d32211f9c4f4d57c49638fcf189fb86ae2ed59961b16527d6657d1cd445671a6` |
| `wrap_design.py` | 1D design of the wrap guide (prints its own checks) |
| `build_wing.py` | builds `wing.step` (and `build.json`) from `parents/flat/sharp/parent.step` |
| `verify_wing.py` | all checks on the reimported STEP; writes `readback.json`, `zebra_join_samples.json`, the PNGs (about 13 min) |
| `readback.json` | the measurements quoted below |
| `zebra_*.png` | CAD-normal zebra renders of the reimported STEP (no edge lines, fixed cameras, mirrored for the left tip) |

## Construction

- **Parent skins kept as they are.** The upper and lower skin faces are copied from `parent.step`
  (same poles and knots, bit for bit). The two flat end faces are dropped.
- **One B-spline face per cap**, wrapping upper cut edge → lip → lower cut edge. In the chord direction it
  uses the parent's own basis (native chord parameter u, degree 5, 33 poles), so its first and last pole rows
  *are* the parent's cut-edge poles at Y = ±350 and the joins are exact.
- **Wrap guides** are the cap's iso-u curves. Each lies in the plane X = X(u) and is the curve
  η⁵ + ζ² = 1, with η = (|Y| − 350)/r(u), ζ = Z/h(u), r(u) = 200u(1 − u) the locked outline reach and h(u)
  the parent half thickness. The cap poles are P = (X_j, ±(350 + r_j η_i), h_j ζ_i): X_j and h_j are the parent's
  edge poles, r_j the exact coefficients of r(u), and η_i, ζ_i a 33-pole degree-5 fit of the guide curve
  (max deviation 3.5e-8 in units of h).
- What the pole layout guarantees:
  - The first three wrap pole rows differ only in Y. That makes the cap tangent and curvature-continuous
    (G2) with the parent skin, whose spanwise curvature is zero.
  - The wrap is mirror-symmetric, so the lip lies in Z = 0 with a horizontal normal.
  - η = 1 at the lip, so the perimeter is the locked outline exactly.
- **Shape.** The crown stays almost flat (ζ ≈ 1 − η⁵/2), then turns progressively into a round lip. In the
  guide plane the lip radius is 5h²/(2r), and the curvature rises steadily from the attachment to the lip.
  Near the LE and TE corners the outline runs almost spanwise. There the cap reproduces the parent's own
  LE/TE closure (Z ∝ distance^(1/3)), moved along the outline, so the surface ends the same way the
  parent does.
- The four faces are sewn into one shell and closed into one solid (OCCT, sewing tolerance 1e-6 mm).
  Each cap has two degenerate edges, at the LE and TE points.

## Measured on the reimported STEP (`readback.json`)

| check | result |
|---|---|
| topology | 1 solid, 1 shell, 4 faces, BRepCheck valid, 0 free edges, 4 degenerate edges (cap corners) |
| volume | 962 251.5 mm³ (parent 864 000 + 2 × 49 125.8) |
| extents | X 0 … 120, Y −400.000 … +400.000, Z ±6.2413 mm; joins at |Y| = 350.000 |
| parent skins | poles and knots identical to `parent.step`; max 2.0e-11 mm from `parent.py` |
| outline | lip within 7.2e-11 mm of the locked outline; lip normal ⟂ outline and horizontal to 4e-12; no cap point outside the outline |
| joins (4) | separation ≤ 7.6e-12 mm; tangent-plane angle ≤ 1.5e-6°; chordwise curvature difference ≤ 7.4e-12 /mm; spanwise ≤ 1.2e-7 /mm; twist ≤ 1e-16 /mm (u from 1e-4 to 1 − 1e-4) |
| principal-curvature map | 1374 × 1401 grid per cap, local minima refined in OCCT. Smallest radius **0.905 mm** at X = 119.880, |Y| = 355.98, Z = −0.634 (attachment thickness there 1.355 mm, limit 0.4 mm, ratio 2.26). The grid reaches within 1.2e-8 mm of the collapsed LE/TE points. |
| sections across the lip | vertical planes normal to the outline at 19 outline points per cap. Worst 0.952 mm at p = 0.98 (limit 0.4). No bending reversals. Lip radius 1.80 mm at mid-chord, smallest 1.06 mm at p = 0.8, and much larger towards both corners (14 mm at p = 0.05 and at p = 0.99). |
| chord-normal sections (X = const) | 32 per cap. No bending reversals. Near the corners these planes are almost tangent to the surface at the lip (cos ≥ 0.002), so their in-plane radius drops below 0.4 mm where the outline runs almost spanwise. Of the stations tested, X ≤ 0.012 mm and X ≥ 118.8 mm are below; X = 0.12 mm and X = 117.6 mm pass. The surface's principal radius at those same points is 1.5 to 65 mm. |
| zebra on the joins | 7 views × 5 stripe axes, front-facing samples only. Phase mismatch ≤ 6.7e-11 stripe; gradient mismatch ≤ 2.3e-5 stripe/mm; stripe-direction change at 1766 light/dark crossings ≤ 7e-9°. |

## Iterations

1. **Flat crown plus constant-radius roll in the guide plane** (lip radius 0.958 h in every X = const plane).
   It was G2 and valid, but the principal map found a **0.122 mm** radius at X = 119.89, |Y| = 357.96
   (thickness 1.31 mm). Near the corners the outline runs almost spanwise, so a roll packed into the last
   1.7 h of an X = const guide is squeezed across the lip. Rejected.
2. **η⁵ + ζ² = 1 wrap** (this build). The turn is spread over the whole reach in a way that matches the
   parent's own LE/TE near the corners and gives a round lip in between.

See `skill-notes.md` for how the skill's printable-radius rule was interpreted.
