# CAD-based wing and tail optimization

**Working design spec · M0ooBot fork · 2026-11-23**

## Purpose

Start with an aircraft the designer has drawn in SolidWorks. Answer a specific question about changing its wing or tail, such as “At my cruise condition, can a different wing planform or dihedral reduce the power this airplane needs?” Return candidate geometries the designer can inspect and a rough whole-aircraft performance comparison. The designer chooses which candidate, if any, to carry into CAD.

The current optimizer searches a parametric aircraft through a series of studies defined in Python. Its STEP importer measures a fuselage. This draft describes the intended overhaul. The proposed Airfoil Converter MCP interface, lifting-surface import and candidate verification are design work. [Current solver](../src/planeopt/solve.py) · [Current STEP importer](../src/planeopt/cadimport.py)

**Flow deflectors belong to CFD.** When the design question is whether a vane redirects airflow onto another surface or around a junction, the coordinator takes the placement and shape study to the CFD tool. For example, a vane intended to improve flow onto the tail needs evidence of what changed at the tail. A canard being studied for its contribution to lift and pitch trim falls within the optimizer's lifting-surface scope.

## What the designer authorizes

A study records a CAD revision, the flight conditions and objective, which regions or dimensions may change, and the limits that must hold for every candidate. The original aircraft is the zero-change candidate. Show the reconstructed baseline overlaid on the CAD geometry and quantify the discrepancy before using it to score changes.

For example, the designer might free outer-wing sweep, dihedral and twist while locking the projected span, root attachment, airfoil profiles and aileron hinge axis. A fixed tail can still use elevator deflection to trim each candidate. Changing elevator chord, tail incidence or tail shape is a separate design permission. The same study can be run with different permissions, recorded as a new comparison.

The coordinator can generate and evaluate candidates within an approved study. A wider envelope or an additional editable feature goes back to the designer. This revises the harness charter's earlier assumption that the designer reshapes the plane at every iteration. The existing design-loop handoff provides the context for that change. [Harness charter](https://github.com/M0ooBot/Engineering-Harness/blob/HEAD/CHARTER.md) · [Optimizer handoff](https://github.com/M0ooBot/Optimizer/blob/handoff-loop/docs/HANDOFF-LOOP.md)

## Geometry shared by the Converter and optimizer

The proposed Airfoil Converter MCP server exposes a lifting surface through its sections, placements, leading and trailing edges, and surface guides. The agent can propose a particular loft or a family with bounded parameters. The numerical evaluator samples that same definition. Shortlisted candidates are rebuilt in a study copy through the harness's SolidWorks bridge and checked against the geometry that earned the score. The Converter already supports section placement, pitch, edge curves and derived surface guides. [Airfoil Converter wing workflow](https://github.com/Trying-to-Fly/Airfoil-Converter#shape-between-ribs)

- **Driving geometry.** Record which inputs drive the loft and which quantities follow from them. If two prescribed edges determine a section's chord and pitch, flag a conflicting pitch request at that station. A guide on the skin and an internal hinge axis constrain different parts of the geometry.
- **Section orientation.** At each station the section plane is perpendicular to the local tangent of one designated spanwise guide. The quarter-chord curve is the default. Twist rotates the profile within that plane. Root and tip planes can differ or be parallel, depending on their local tangents. Report tip rise and tip tilt separately, with the guide and reference frame used. A raised tip that levels off can have zero tilt relative to a horizontal root.
- **Between-station smoothness.** Check the leading and trailing edges, spanwise bend, twist and resulting skin between defining ribs. For a single smooth upsweep, the local dihedral angle should increase monotonically through the bend, subject to the study's minimum bend radius and twist-rate limit. A gull wing or deliberate polyhedral break requires a study definition that allows it. Preserve declared sharp features such as the trailing edge and hinge gap.
- **Proposed transition rule.** Where an editable region leaves a fixed flat panel in a smooth bend, match position, tangent and curvature at the join. Express the initial search with a few controls such as onset, rise and bend distribution. Validate intermediate geometry for ripples and local curvature peaks even when all control stations pass.
- **CAD readback.** Inspect the actual loft between its sections and guides. Report the worst sampled surface deviation and its location. Compare rebuilt section-plane normals against the designated guide tangents and report the worst alignment error. Record the sampling and tolerance used so the designer can interpret a pass.

## Straight spar and printable wing sections

Treat spar fit as a hard geometric constraint on each half-wing. The working interpretation of a printable section is a physical print piece with a maximum **250 mm build height** in its print orientation. Use the piece boundaries and cut planes from the study.

| Requirement | Rule |
| --- | --- |
| Continuous path | One straight spar axis runs from the root through every intervening print joint into the final tip piece. The coordinator or solver may choose its placement and direction. |
| Tip engagement | The full spar extends at least **20 mm** into that final piece, measured along its axis from the entry cut face. Account for the full cross-section at an oblique joint. |
| Stock dimensions | Allow round or square stock in whole-millimetre increments, provisionally **6 mm minimum**. Round sizes use outside diameter. Square sizes use side length, such as 6×6 or 7×7 mm. |
| Envelope | Check the full circular or square envelope along its path with the study's required clearance. A square spar has one rotation about its axis along the full length. Include its corners. |
| Evidence | Show the axis and cut planes on the wing. Report chosen stock, minimum clearance and its location, and actual tip engagement. A fit check establishes geometry. Strength and stiffness require separate evidence about material, wall thickness and loading. |

A wing may keep bending beyond the end of the required spar engagement. The proposed first study fixes the print split so every candidate faces the same tip joint. Whether the optimizer may later move joints is open.

## Who chooses the next shape

1. **Evaluate the baseline.** Record trim and the performance measures relevant to the question. Apply common flight conditions, mass assumptions and trim rules to the comparison. If a candidate changes estimated mass or CG, record the change and its source.
2. **Let the coordinator choose an experiment.** It can submit a single candidate or ask the numerical solver to tune a family. For example, a first family frees only dihedral. A second can free taper if the earlier results suggest the spar limits any further useful rise.
3. **Use code for numerical search and enforcement.** The solver tunes continuous variables, can enumerate small discrete choices such as stock sizes, and records infeasibility or failed convergence. The Converter's geometry checks apply to candidates from every search method. Save enough inputs, settings and versions to repeat the evaluation.
4. **Spend the next run on a question the results leave open.** The coordinator can change family within the authorized freedoms, hold one variable fixed to test an explanation, or stop when remaining differences are small relative to model sensitivity. A finer aerodynamic evaluation or a targeted CFD comparison can check a shortlisted result. The CFD app supplies aerodynamic evidence at recorded poses and conditions. The coordinator interprets that evidence alongside whole-aircraft performance.

## What comes back after a run

Every claim in the report should say **how the tool knows it**. The coordinator gets a compact comparison first and can request spanwise curves, geometry, solver records or CAD evidence as needed.

| Question | Evidence returned |
| --- | --- |
| What was tested? | Baseline and candidate revisions, editable parameters and locks, study objective, flight cases, input assumptions, solver status and check coverage. |
| Where did the wing change? | Front and plan overlays, changes in chord, twist and dihedral along the span, and the reference guide. Identify whether alignment was checked on generated geometry or measured after the CAD rebuild. |
| What did the airplane gain or lose? | Baseline and candidate values with units and deltas at matched conditions. Include the relevant measures, such as trimmed power, lift-to-drag ratio, control deflection and stability margin. Identify any condition where trim failed. |
| What limited the result? | Constraint margins, locations and the tested settings. A spar placement that fails establishes a failure for that placement. Ruling out every placement requires a separate feasibility result. |
| What else was tried? | For a family search, show nearby feasible and rejected candidates, including how much performance and clearance changed. Preserve useful simpler alternatives when their scores are close. |
| How reliable is the comparison? | Identify the model, geometry and resolution behind the prediction. Give the result of any finer-resolution, second-condition, CFD or physical check and list outstanding checks. Label causal explanations as hypotheses until a controlled comparison supports them. |

## First end-to-end study

The proposed first slice uses one designer-supplied wing, a fixed tail geometry and a declared cruise condition. Evaluate the baseline, tune one authorized wing family, rebuild a shortlisted loft, and deliver the baseline-versus-candidate report. Tail tuning uses the same design permissions in a subsequent study.

- With all freedoms locked, reproduce the baseline within a stated tolerance.
- Exercise a flat wing, a smooth upsweep, and a wing with a ripple between ribs. Locate the ripple. Check both a tilted tip and a raised tip whose guide levels off.
- Exercise the spar engagement boundary at 20 mm and just below it. Include an oblique cut and a square spar whose centerline clears but whose corner intersects the skin.
- Compare an agent-guided study with a fixed numerical search using the same freedoms, constraints and evaluation budget. Judge feasible performance, persistence of the gain under cross-checks, and which evidence helped the designer decide.

## Decisions still open

- How the CAD revision, section profiles, guide curves and locks travel between the Converter and optimizer, and how custom intermediate guides change the section definition used for aerodynamics.
- Numerical tolerances for baseline reconstruction, orientation, curvature, surface deviation and spar clearance. The bend-radius and twist-rate limits need a physical basis for each design.
- Source of print orientation and cut planes, whether studies can vary them, and which available stock sizes, wall thicknesses and materials to supply.
- Which aircraft revision and flight case to use first, which performance gain matters to the designer, and when the remaining uncertainty justifies commissioning CFD.
