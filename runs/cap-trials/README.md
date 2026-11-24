# Rounded wingtip skill trials

## Starting point for winglet design

The designer accepted **Rebuilt** for all five geometries. The source of record for each is its builder and whole-model STEP under `retest/`. **Control** identifies the earlier export used for comparison. The next session can start from the accepted model and its parent-to-cap attachment geometry.

| Accepted whole model | Construction and measurements |
| --- | --- |
| [Straight wing](retest/flat/wing.step) | [Straight wing study](retest/flat/README.md) |
| [Twisted wing](retest/twisted/wing.step) | [Twisted wing study](retest/twisted/README.md) |
| [55° swept wing](retest/swept55/wing.step) | [Swept wing study](retest/swept55/README.md) |
| [Vertical fin, 40 mm cap](retest/fin/40/wing.step) | [40 mm fin study](retest/fin/40/README.md) |
| [Vertical fin, 5 mm cap](retest/fin/5/wing.step) | [5 mm fin study](retest/fin/5/README.md) |

The existing skill covers closing an exposed tip. Winglet geometry is the next design decision. The reusable pieces are the parent skin edges, 3D guides, and reimported STEP checks for attachment continuity and local radius. On the accepted 55° wing the reversals disappeared in the nine paired cuts, and the 97%-chord cut tightened to a 0.000758 mm lip radius. The designer reviewed and accepted that tradeoff. The rebuilt 5 mm fin is symmetric across its center plane.

Open `runs/cap-trials/review.html` locally from a checkout to orbit the models and compare revisions.

## Original trial controls

[Open the 3D review](review.html). It opens on the complete rebuilt 55° wing. Switch **Skill version** between **Control** and **Rebuilt** at the same camera position and scale. Select **Vertical fin** and change **Cap reach** to compare the two fin heights. Each selection links to the STEP and its own measurements. [Paired skill retest](retest/README.md)

Four fresh Sol agents received dimensional briefs, the [rounded wingtip skill](../../.agents/skills/close-rounded-wingtip/SKILL.md), and access to the project tools. Each agent chose its construction independently. Fresh reviewers checked the exports and the evidence. The 5 mm fin was a follow-up to the original four trials.

| Study | Geometry exercised | Measured result |
| --- | --- | --- |
| [Straight wing](flat/README.md) | 800 mm span, constant chord, original symmetric polynomial section | The checker found 3,924 stripe crossings across the four skin-to-cap joins and measured no direction change at those crossings. |
| [Twisted wing](twisted/README.md) | 1,000 mm span, 5° dihedral, 3° washout, modified NACA 2412 | The largest measured stripe turn at an attachment was 0.00188°. At 99% chord, the guide radius at the perimeter was about 0.0052 mm. |
| [55° swept wing](swept55/README.md) | 900 mm span, level, closed-trailing-edge NACA 0009 | The attachment stripe change peaked at about 0.000000103°. At the 97% chord section, the turn tightens inside the cap before reaching the lip. |
| [Vertical fin](fin/README.md) | 280 mm height, 12% elliptical section, 40 mm cap | The sampled joins match in position, tangent plane and curvature. The ten measured section radii at the lip ranged from 0.262 to 0.702 mm. |
| [Short-cap fin](fin/cap-5.0mm/README.md) | Same parent and attachment, 245 mm height, 5 mm cap | The tightest measured radius along the ten wraps was 1.88 mm, about 2.17 mm above the attachment. |

All five STEP files reimport as single, closed, strictly valid solids. Each builder wrote a geometry-specific surface construction and measurement adapter on top of CadQuery and Open CASCADE. These runs demonstrate that agents can use the skill to design and inspect new analytic parents. Applying the process to an unfamiliar supplied STEP, including identifying corresponding edges and positions on differently parameterized faces, remains untested. Stripe frequencies and sample coverage vary by trial. Ordinary two-direction derivatives also break down at collapsed endpoints. Read the individual reports for those limits.

To generate another fin with the parent held fixed, run this command from the repository root and substitute the desired reach in millimetres:

```powershell
python runs/cap-trials/fin/build_variant.py --cap-reach-mm 5
```

The command exports and measures the resulting fin in its own reach-specific folder. To regenerate the review from the models listed in `cases.json`, run `python runs/cap-trials/build_viewer.py`. The generator verifies each STEP against the hash in its measurement report. The viewer displays sampled CAD face normals and applies mirrored diagnostic lighting to the wings. The numerical reports evaluate surface derivatives directly on the STEP.
