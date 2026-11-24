# Wing B rounded wingtip

[Open the paired zebra viewer](../../runs/wing-b-rounded-tip-g2/zebra-viewer.html) by downloading the HTML and opening it locally. It starts with the whole 1,200 mm wing. Switch between **Previous** and **Revised** at the same camera and stripe angle. The [section viewer](../../runs/wing-b-rounded-tip-g2/review.html) also includes fixed-aft cuts through the lip. The [revised STEP](../../runs/wing-b-rounded-tip-g2/wing-b-rounded.step) contains one 600 mm half-wing. Both viewers mirror it at the root.

## What the designer asked for

The discussion started with three possibilities for a wingtip device: an airfoil-sectioned winglet, a plate, and a sharp terminal cap. The Spitfire reference clarified the desired first example. The wing itself should finish with a broad oval in plan, follow its existing bow, and have a little roundness through its thickness. The designer proposed a boundary surface that wraps from the upper skin around a perimeter curve to the lower skin. That became the construction brief.

The parent was the [untwisted, bowed SD7037 Wing B](../../runs/paired-wing-trial/a/measured/01-wing-occt.step) from the paired silhouette trial, with a 75 mm tip chord. The later twisted Wing B spar study has a separate parent. The cap starts at 540 mm, beyond the 510 mm end of the hinge, and closes at the original 600 mm projected half-span.

## What changed between the two builds

The [first build](../../runs/wing-b-rounded-tip/build.py) took the upper and lower cut edges directly from the parent. A curved perimeter set the planform. Wrap guides set the thickness and turned around the rim, tightening near the trailing edge. The resulting [first rounded STEP](../../runs/wing-b-rounded-tip/wing-b-rounded.step) had closely matched tangent planes at the join. A full-wing zebra review drew attention to the abrupt change in how quickly the skin bent immediately outboard.

The [revised build](../../runs/wing-b-rounded-tip-g2/build.py) samples the original skin at the cap's projected coordinates and blends those heights into the rolled guide. It matches first and second derivatives at the cut and redistributes aft guides as thickness falls. Fairing the leading and trailing outlines moved the first perimeter by up to 0.84 mm in plan. The exporter reconciles the shared 3D cut edge with both faces' UV curves, verifies the reimported solid, then publishes the STEP.

## What the checks showed

| Check | Result |
| --- | --- |
| Reimport | One closed four-face solid; standard and strict Open CASCADE checks passed. The retained skin's B-spline poles, knots and weights matched the original. |
| Upper join, 150 mm aft | The mismatch in the local curvature tensor fell from 0.03350 to 0.00000173 per mm. Sections beyond the join showed the turn building over several millimetres. |
| Full-wing zebra comparison | Across 1,001 locations per skin and 364 combinations of camera, stripe angle and displayed wing half, the 95th-percentile stripe-gradient mismatch fell from 0.0613 to 0.0000400 cycles per mm of surface distance. |
| Detected stripe edges | At 11 upper-skin light/dark crossings in the default centered view with stripes at 51°, the largest change of stripe direction fell from 26.30° to 0.00693°. |

The zebra analyzer computes the paired viewer's continuous reflection phase from CAD normals, then its gradient from surface derivatives. A local stripe direction can be evaluated between actual bands. The 11-crossing result uses solved light/dark transitions. Fine CAD-normal sampling removed false triangular islands seen in a coarse mesh. Switching revisions at a fixed pose kept real differences visible.

The numerical zebra sweep covers the middle 99.8% of the common fore-to-aft interval along the join. The UV patches collapse at the terminal vertices, so ordinary two-direction surface derivatives are undefined there. Smoothly turning bands can persist in grazing views. These findings support the cap blend over the measured region. The mirrored root, printability, and aerodynamic performance are separate questions. [STEP readback](../../runs/wing-b-rounded-tip-g2/g2-readback.json) · [Zebra measurements](../../runs/wing-b-rounded-tip-g2/zebra-analysis.json) · [Fixed-aft trace](../../runs/wing-b-rounded-tip-g2/zebra-trace.png)

## Reproduce the checks

From the repository root, use Python 3.12 or newer with the project's CAD dependencies. The study used CadQuery 2.8 and OCCT 7.9.3. Run the analyzers on the committed STEP files first. Supply output paths outside the checkout to keep the recorded evidence.

```sh
python runs/wing-b-rounded-tip-g2/check_g2.py --output g2-check-local.json
python runs/wing-b-rounded-tip-g2/analyze_zebra.py --output zebra-local.json
```

To regenerate either CAD revision, run its `build.py` in a working copy. It writes `wing-b-rounded.step` in its own study directory. STEP headers include export timestamps, so rerun the checks and record the hash of the file measured. The archived numerical zebra results use the [paired viewer's shader](../../runs/wing-b-rounded-tip-g2/zebra-viewer.html); the separate section viewer has its own lighting model.

Use the [wingtip cap skill](../../.agents/skills/close-rounded-wingtip/SKILL.md) for future designs.

For subsequent winglet work, start with the [accepted whole-wing and fin models](../../runs/cap-trials/README.md).
