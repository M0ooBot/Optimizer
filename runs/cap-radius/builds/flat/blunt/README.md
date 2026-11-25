# Rounded tip caps: flat / blunt parent

`wing.step` is the 800 mm wing with a rounded cap on each tip. The parent is the 700 mm straight wing from `runs/cap-radius/parents/flat/blunt/`: 120 mm chord, symmetric polynomial section, and a flat trailing-edge face 0.4 mm thick. The parent's two skins and its trailing-edge face are copied unchanged from `parent.step`. Each cap adds one B-spline face plus a small planar face that rounds off the end of the trailing-edge face. The result is one closed, valid solid with 7 faces.

## Construction

- **Wraps (guides).** Each wrap lies in a plane of constant X, at the chord station of parent chord parameter `s`. It runs from the upper cut edge at |Y| = 350 out to the lip and back to the lower cut edge:
  `|Y| = 350 + R(s) E(w)` and `Z = h(s) S(w)`, for `w` from +1 to -1.
  - `h(s)` is the parent's half thickness.
  - `E = 1 - w^2` and `S = (15w - 10w^3 + 3w^5)/8`.
  - Because `S'(±1) = S''(±1) = 0`, the cap meets the parent skin with the same tangent plane and the same curvature (G2).
  - `E'(0) = 0`, so the wrap crosses the lip vertically. The upper and lower halves are one smooth face.
- **Lip.** The lip is the `w = 0` curve, at Z = 0. `R(s)` is the locked outline's reach at the same X, so the lip is the locked outline.
- **Trailing edge.** Here the locked outline arrives tangent to the trailing edge with zero reach, so a rounded closure cannot sit on it exactly. In the last few millimetres of reach, `R` is shifted aft by at most 2e-3 mm. The shift fades out as `exp(-(120 - X)/6e-3)`, and at X = 120 the reach is 0.352 mm.
  - The flat trailing-edge face is extended at each end by a planar rounded closure: 0.4 mm high, reaching 0.352 mm, with apex radius 0.200 mm (half the trailing-edge thickness).
  - The cap meets this closure at 2.85° at the corners, the same angle as the parent's crease between skin and trailing-edge face. The angle falls to 0.10° at the apex.
- **B-spline form.** The cap is an exact tensor-product B-spline: degree 5 × 5, with 82 × 6 poles.
  - Its cut-edge rows of poles are the poles of the parent's own edge curve (after knot insertion), so the cap and the parent share that edge exactly.
  - `R(s)` is interpolated to within 1.1e-6 mm.
  - The leading-edge row collapses to the leading-edge point.
  - The left cap mirrors the right one.

## Results, measured on the reimported STEP (`readback.json`)

| Check | Result |
|---|---|
| STEP | SHA-256 `fdf987b5931fdd4719e184dc9dc414b82b53cd2f9e3445a92d15ce32ec6a0fd5`, 140 110 bytes |
| Solid | 1 solid, 1 closed shell, 7 faces, BRepCheck valid, 0 free edges, maximum tolerance 1e-7 mm, volume 958 092.90 mm³ |
| Extents | X 0 to 120, Y ±400.000000 (span 800.0000002), Z ±6.2415 |
| Parent preserved | Skin and trailing-edge face poles identical to `parent.step` (0.0 difference); within 2.1e-11 mm of `parent.py` |
| Cut edges | At |Y| = 350 exactly |
| Joins (4) | Separation ≤ 1.1e-11 mm; tangent-plane angle ≤ 1.5e-6°; curvature difference ≤ 2.2e-6 /mm chordwise, ~1e-18 /mm spanwise (G2) |
| Lip and outline | Lip at Z = 0; tip at Y = 399.9999999972, X = 60.0008. Lip within 1e-5 mm of the locked outline everywhere except X > 119.95, where the trailing-edge shift puts it up to 0.0027 mm aft (at X 119.9938, Y ±352.83) |
| Smallest radius where t ≥ 0.8 mm | **0.743 mm** on the lip at (111.33, ±384.04, 0), where t = 6.25 mm (required 0.4 mm) |
| Smallest radius anywhere on the cap | **0.356 mm** at (120.000, ±350.104, ±0.198), on the edge of the closure, where t = 0.40 mm (required t/2 = 0.20 mm) |
| Map | 1580 × 919 grid per cap plus local refinement; 0 points below the rule; no concave points. A small saddle region lies only at X > 119.9993. The grid comes within 0.002 mm of the collapsed leading-edge point |
| Sections | 26 at constant X, 14 normal to the lip, 12 at constant Y: no curvature reversals. In normal-to-lip cuts the lip radius ranges from 0.748 mm (p = 0.8) through 1.27 mm (p = 0.5) to 27 mm near the ends |
| Zebra | 50 views × 8 stripe angles over the 4 joins, 34 032 light/dark crossings. Maximum phase mismatch 7e-8 stripe, stripe-gradient mismatch 6e-6 (relative), stripe-direction mismatch at crossings 6e-13° |

For context, the parent's own skin has its tightest radius, 0.991 mm, at X = 119.92. The cap's tightest point near the trailing edge is 0.974 mm at (119.90, 354.5, ±0.60), the same feature carried outboard. Images are in `images/`:
- zebra stripes on the whole wing and on both tips (centred camera, mirrored lighting);
- oblique views;
- a trailing-edge close-up showing the rounded end of the flat face;
- `radius_ratio_map.png`.

## Files and how to run them (from the project root)

```
.venv/bin/python runs/cap-radius/builds/flat/blunt/build_cap.py         # writes wing.step and build.json
.venv/bin/python runs/cap-radius/builds/flat/blunt/verify_cap.py        # writes readback.json (about 1 min)
.venv/bin/python runs/cap-radius/builds/flat/blunt/render_zebra.py      # images/zebra_*.png
.venv/bin/python runs/cap-radius/builds/flat/blunt/render_te_corner.py  # images/*te_corner*.png
.venv/bin/python runs/cap-radius/builds/flat/blunt/plot_radius_map.py   # images/radius_ratio_map.png
.venv/bin/python runs/cap-radius/builds/flat/blunt/explore_te_law.py 1e-3,3e-3 2e-3,6e-3   # design scan
```

`design.py` holds the analytic cap definition and its exact derivatives. Rebuilding changes the STEP header timestamp, and with it the SHA-256, so rerun `verify_cap.py` after any rebuild.
