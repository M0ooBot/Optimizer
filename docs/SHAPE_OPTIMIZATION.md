# CAD-based wing and tail optimization

**Working design spec · M0ooBot fork · 2026-11-23**

## Purpose

Start with an aircraft the designer has drawn in SolidWorks. Answer a specific question about changing its wing or tail, such as “At my cruise condition, can a different wing planform or dihedral reduce the power this airplane needs?” Return candidate geometries the designer can inspect and a rough whole-aircraft performance comparison. The designer chooses which candidate, if any, to carry into CAD.

The current optimizer searches a parametric aircraft through a series of studies defined in Python. Its STEP importer measures a fuselage. This draft describes the intended overhaul. The proposed Airfoil Converter MCP interface, lifting-surface import and candidate verification are design work. [Current solver](../src/planeopt/solve.py) · [Current STEP importer](../src/planeopt/cadimport.py)

**Flow deflectors belong to CFD.** When the design question is whether a vane redirects airflow onto another surface or around a junction, the coordinator takes the placement and shape study to the CFD tool. For example, a vane intended to improve flow onto the tail needs evidence of what changed at the tail. A canard being studied for its contribution to lift and pitch trim falls within the optimizer's lifting-surface scope.

The optimizer's design variables focus on lifting surfaces and their controls. Carry the designer's fuselage geometry, mounting interfaces, mass and estimated aerodynamic contribution into the whole-aircraft calculation. Route questions about local fuselage shape, flow separation and wing-root fairings to CFD. Record when a change in wing placement makes an assumed wing–body interaction uncertain.

## Coordinator and agent handoffs

The coordinator is an AI agent. Starting with the rough mission, it works with the designer to record requirements, identify assumptions and choose the next engineering question. Different agents may handle CAD, optimization and CFD, and another agent may take over coordination later.

Each agent updates the aircraft's shared design record when a decision is made or work passes to another agent. Keep the handoff short: identify the aircraft revision, the question and constraints, what the designer approved, what was tried and why, and the next unresolved question. Link the CAD revision or run that supports each finding and identify assumptions still awaiting a check. A new agent should be able to resume from that record without access to an earlier chat.

Save the input revisions actually used by each run. A change to a dependency such as measured mass, available elevator travel or airfoil evidence marks the affected verdicts as needing re-evaluation. Keep the original run available for comparison.

**Aircraft workbook (accepted design).** Keep one master XLSX per aircraft with separate sheets for mission and flight cases, geometry and mounts, controls, component masses and loading-case CGs, and results and decisions. Identify units, frame or sign convention, source and revision. Distinguish editable design intent from values read from CAD, Systems or the slicer. Keep candidate values separate from the accepted design. Store predicted control deflections by candidate, flight case and command, with a link to the originating run. The coordinator updates the workbook from approved edits and revisioned results. Workbook editing and synchronization details are open.

## Optimizer MCP server

The proposed optimizer exposes numerical aircraft studies through an MCP server. The designer asks the coordinator an engineering question in ordinary language. The coordinator uses the aircraft's intent, revision and approved design freedoms to formulate an experiment, interpret the results and explain the trade-offs. The server enforces study bounds, computes comparable results and preserves the inputs and evidence behind them.

The coordinator needs a small set of operations it can combine across questions:

- Inspect the baseline, flight cases, available analyses, adjustable parameters and their limits.
- Find earlier runs for the same geometry and conditions, or evaluate a specified candidate.
- Sweep a parameter or optimize a bounded family against a named objective and constraints.
- Compare saved candidates at matched conditions and retrieve the margins, spanwise data, model limits and run references behind a result.

**Example study.** A designer asks whether changing wing dihedral helps. For a cruise-power objective, compare the baseline with several permitted dihedral angles while holding other design choices fixed and re-trimming each candidate. Record whether span is fixed in projection or along the wing. Check spar fit and the available static stability measures, then refine promising ranges if the initial differences warrant it. Report the best feasible angle tested for that objective, or a range with indistinguishable predicted performance. If the designer means handling in a gust, identify the required dynamic evidence and the extent of the available model before drawing a conclusion.

### From mission requirements to flight cases

The coordinator translates the designer's brief into numerical requirements before submitting a study. For a camera mission, record the payload, outbound and return distances, time and task on station, applicable speed and wind assumptions, and minimum reserve on landing. Keep assumptions distinguishable from designer requirements. Record the units and link each requirement to the flight cases or mission segments that test it. The designer can bring a rough CAD layout directly to this step. Preliminary sizing from the mission alone is a separate capability whose scope remains open.

For each requirement the optimizer returns its target, predicted value, margin, and the case or segment that determines the verdict. Aggregate energy across the mission before checking landing reserve. A segment outside the model's coverage leaves the dependent mission verdict unresolved. For acceptance, exercise a mission whose individual legs trim successfully but whose combined energy use misses the requested landing reserve.

## What the designer authorizes

For the first import trial, the designer will supply a SolidWorks aircraft with the wings, tail surfaces and fuselage as separate bodies. Use that real design to decide how bodies are identified, how their placement is read, and how they are matched to Airfoil Converter sections. Body naming, part or assembly organization, and the exchange format remain open until that trial.

A study records a CAD revision, the flight conditions and objective, which regions or dimensions may change, and the limits that must hold for every candidate. The original aircraft is the zero-change candidate. Show the reconstructed baseline overlaid on the CAD geometry and quantify the discrepancy before using it to score changes.

For example, the designer might free outer-wing sweep, dihedral and twist while locking the projected span, root attachment and airfoil profiles. If the study also defines an aileron hinge axis, the designer can lock it. A fixed tail with a modeled elevator can still use elevator deflection to trim each candidate. Changing elevator chord, tail incidence or tail shape is a separate design permission.

The coordinator can generate and evaluate candidates within an approved study. A wider envelope or an additional editable feature goes back to the designer. This revises the harness charter's earlier assumption that the designer reshapes the plane at every iteration. The existing design-loop handoff provides the context for that change. [Harness charter](https://github.com/M0ooBot/Engineering-Harness/blob/HEAD/CHARTER.md) · [Optimizer handoff](https://github.com/M0ooBot/Optimizer/blob/handoff-loop/docs/HANDOFF-LOOP.md)

## Controls on continuous CAD surfaces

The expected input can have a continuous wing, horizontal tail and vertical fin. Import accepts these surfaces without modeled ailerons, elevators, rudders, hinge gaps or moving parts.

For a study that needs a control, define it on its parent surface for analysis. Record its role, inboard and outboard limits, hinge line or chord fraction, deflection direction and available travel. A designer can supply the layout or authorize the coordinator to propose one. Show the outline on the imported surface and label its source. Placement and sizing become search variables only when the study frees them. Deflection during a flight case is a separate quantity.

Use the same declared control layout for baseline and candidate when comparing wing shapes, and record how its span and chord coordinates map onto each geometry. When comparing control layouts, identify the layout on each candidate. Tie predicted trim deflection and control margins to the layout and travel used. If the study has no declared pitch control, it can report the pitching moment at a specified condition and flag elevator trim and remaining travel as unevaluated. A proposed layout can be handed back to the designer as an annotated outline for detailed CAD work.

### First aileron and rudder studies

Start with the accepted wing and tail placement and a trimmed baseline for each declared speed and loading case. The designer supplies a required authority or asks for a comparison of available authority. With no target, report the trade-off curve and the limiting condition. Keep static stability, control effectiveness and time response as separate reported quantities. The current sideslip check measures static stability derivatives. Control-deflection sweeps need their own implementation and verification. [Current model limits](MODEL_DETAILS.md#34-stability-quantities)

- **Ailerons.** Vary inboard and outboard stations and hinge position within the permitted wing region. Declare independent up/down travel and any differential. Sweep both command directions at each case. Report the change in rolling and yawing moment about the CG, individual left/right surface deflections required for a target rolling moment, remaining travel and proximity to the model's section-lift limits. Check a physical hinge axis and clearances against the curved wing. Report reduced confidence where the deflected outboard section approaches stall.
- **Rudder.** For a fixed fin and tail placement, vary the movable span and hinge position. Specify a sideslip range or an external yawing moment to oppose. Sweep rudder both ways, find the deflection that balances yawing moment where attainable, and report remaining travel, side force and accompanying roll moment. Identify fin or fuselage wake uncertainty. For mixed controls such as ruddervators, enforce physical travel on each surface under simultaneous commands, or mark the configuration unsupported.
- **How to choose.** Compare feasible outlines by the designer's priority, for example least movable area that achieves the specified moment and reserve across all cases. Treat predicted roll rate, bank response time and gust recovery as separate dynamic questions requiring suitable inertia, damping and validation. A candidate needing uncertain local wake or hinge-load predictions can be referred for focused CFD evidence.

## Geometry shared by the Converter and optimizer

The proposed Airfoil Converter MCP server exposes a lifting surface through its sections, placements, leading and trailing edges, and surface guides. The agent can propose a particular loft or a family with bounded parameters. The numerical evaluator samples that same definition. Shortlisted candidates are rebuilt in a study copy through the harness's SolidWorks bridge and checked against the geometry that earned the score. The Converter already supports section placement, pitch, edge curves and derived surface guides. [Airfoil Converter wing workflow](https://github.com/Trying-to-Fly/Airfoil-Converter#shape-between-ribs)

- **Driving geometry.** Record which inputs drive the loft and which quantities follow from them. If two prescribed edges determine a section's chord and pitch, flag a conflicting pitch request at that station. A guide on the skin and an internal hinge axis constrain different parts of the geometry.
- **Section orientation.** Preserve the as-drawn section planes when evaluating an imported baseline or swapping profiles on it. Report any mismatch with the guide-normal convention. If the evaluator cannot represent those planes faithfully, flag the affected analysis as unsupported. A refit requires designer authorization and its own candidate revision. For generated bends, each section plane is perpendicular to the local tangent of the designated spanwise guide (quarter-chord by default). Twist rotates the profile within that plane. Report tip rise and tip tilt separately with the guide and reference frame. A raised tip that levels off can have zero tilt relative to a horizontal root.
- **Between-station smoothness.** Check the leading and trailing edges, spanwise bend, twist and resulting skin between defining ribs. For a single smooth upsweep, the local dihedral angle should increase monotonically through the bend, subject to the study's minimum bend radius and twist-rate limit. A gull wing or deliberate polyhedral break requires a study definition that allows it. Preserve declared sharp features such as the trailing edge and hinge gap.
- **Proposed transition rule.** Where an editable region leaves a fixed flat panel in a smooth bend, match position, tangent and curvature at the join. Express the initial search with a few controls such as onset, rise and bend distribution. Validate intermediate geometry for ripples and local curvature peaks even when all control stations pass.
- **CAD readback.** Inspect the actual loft between its sections and guides. Report the worst sampled surface deviation and its location. Check rebuilt section planes against the source planes for an as-drawn baseline or profile swap, and against the guide tangents for a generated bend. Record the sampling and tolerance used so the designer can interpret a pass.

On each shortlisted CAD rebuild, recheck the selected spar's full envelope, locked interfaces, minimum clearance and 20 mm tip engagement against the rebuilt geometry and its actual cut planes. Record the CAD revision and limiting location. Until that check runs, label a passing fit as verified on generated geometry only. A clearance failure on the rebuild fails CAD fit even when the loft meets its overall surface-deviation tolerance.

## Airfoil sources and wind-tunnel evidence

The imported aircraft's baseline uses the section coordinates that shaped its CAD wing. When the designer opens airfoil selection to the coordinator, the default source for new candidate profiles is the UIUC Airfoil Coordinates Database. Save the exact coordinates used and their source with the study. [UIUC coordinate database](https://m-selig.web.engr.illinois.edu/ads/coord_database.html)

For a profile-only swap, lock the leading-edge and trailing-edge guides, section stations and planes, chord lengths, and installed incidence and twist. Use the same trailing-edge chord reference for finite-thickness profiles. Rebuild derived upper and lower skin guides from the new profile. Check any independently locked guide or hinge for compatibility and report a conflict. Recheck spar clearance and re-trim the airplane at matched flight conditions. An experiment that also varies planform or twist is a separate study.

For each candidate, look for UIUC Low-Speed Airfoil Tests covering the relevant section operating range. Link the physical model and test run, and report the Reynolds number, incidence or section lift coefficient, and surface condition used in the comparison. Coverage varies by profile and condition. Where a suitable measurement is unavailable, record the gap and label any NeuralFoil or CFD result as a prediction. [UIUC low-speed airfoil tests](https://m-selig.web.engr.illinois.edu/uiuc_lsat.html) · [UIUC wind-tunnel data](https://m-selig.web.engr.illinois.edu/pd.html)

## Straight spar and printable wing sections

Treat spar fit as a hard geometric constraint on each half-wing. The working interpretation of a printable section is a physical print piece with a maximum **250 mm build height** in its print orientation. Use the piece boundaries and cut planes from the study.

| Requirement | Rule |
| --- | --- |
| Continuous path | One straight spar axis runs from the root through every intervening print joint into the final tip piece. The coordinator or solver may choose its placement and direction. |
| Tip engagement | The full spar extends at least **20 mm** into that final piece, measured along its axis from the entry cut face. Account for the full cross-section at an oblique joint. |
| Stock dimensions | Allow round or square stock in whole-millimetre increments, provisionally **6 mm minimum**. Round sizes use outside diameter. Square sizes use side length, such as 6×6 or 7×7 mm. |
| Outer skin | Leave at least **1.0 mm** from the finished spar channel to the exterior wing skin everywhere along the path. Measure the shortest 3D distance, including square-channel corners; a square spar keeps one rotation about its axis. The channel envelope includes any installation allowance around the stock. Intentional channel openings at root and print-joint cut faces are excluded from the skin measurement. |
| Evidence | Show the axis and cut planes on the wing. Report chosen stock, minimum clearance and its location, and actual tip engagement. A fit check establishes geometry. Strength and stiffness require separate evidence about material, wall thickness and loading. |

A wing may keep bending beyond the end of the required spar engagement. The proposed first study fixes the print split so every candidate faces the same tip joint. Whether the optimizer may later move joints is open.

The 1.0 mm floor is a geometric wall allowance intended to accommodate two adjacent perimeter paths. If the selected print profile requires a thicker wall for two paths, use that larger minimum. Record minimum wall thickness and location separately from any slicer verification of two paths in the critical region. Until the print orientation, settings and toolpaths are available, only the geometric clearance can be verified.

For a profile-only airfoil swap, spar placement is free by default. Passing requires a verified placement of the full straight spar through the candidate wing. Start with the baseline axis, then search permitted positions and directions if needed. Honor any explicitly fixed interface, such as a root socket. Apply the same stock assumptions to all profiles; resizing stock is a separate study freedom. Save a passing axis, its change from baseline, worst clearance and tip engagement. Carry the actual spar position and length into the candidate mass and CG calculation before re-trimming. If a bounded search finds no passing placement, report what was searched and leave feasibility unresolved unless a geometric check establishes that none can fit.

## Candidate mass, CG and endurance

For each loading case, the optimizer reads component masses and positions, payload and battery configuration, propulsion data, usable battery energy and the mission reserve. Record the source of each input. Designer estimates can supply the first study. Slicer estimates and weighed parts from the Systems mass ledger can refine later revisions. Record whether usable energy already includes the reserve so the calculation charges it once.

Estimate printed wing and tail mass from candidate geometry and a construction profile for the material and print settings. For a hollow skin-and-rib build, use skin mass per unit external surface area plus estimated ribs, joints and local reinforcement. Account separately for the selected spar and hardware and identify any items already included in a weighed assembly. Wing and tail may use different profiles. Start with provisional coefficients, then calibrate against representative sliced and weighed sections as the manufacturing method settles. The current per-surface mass model is a starting point. [Printed-surface calculation](../src/planeopt/massmodel.py)

Recompute each affected component's mass and mass location when its geometry or placement changes, then update aircraft mass, CG, trim, power and endurance for each loading case. A rigid translation of the same wing changes its mass location while preserving its own mass. Report which values were fixed, estimated, recomputed or measured. When plausible print-mass variation changes the preferred candidate, show that sensitivity alongside the ranking.

## Who chooses the next shape

1. **Evaluate the baseline.** Record trim and the performance measures relevant to the question. Apply common flight conditions, mass assumptions and trim rules to the comparison. If a candidate changes estimated mass or CG, record the change and its source.
2. **Let the coordinator choose an experiment.** It can submit a single candidate or ask the numerical solver to tune a family. For example, a first family frees only dihedral. A second can free taper if the earlier results suggest the spar limits any further useful rise.
3. **Use code for numerical search and enforcement.** The solver tunes continuous variables, can enumerate small discrete choices such as stock sizes, and records infeasibility or failed convergence. The Converter's geometry checks apply to candidates from every search method. Save enough inputs, settings and versions to repeat the evaluation.
4. **Spend the next run on a question the results leave open.** The coordinator can change family within the authorized freedoms, hold one variable fixed to test an explanation, or stop when remaining differences are small relative to model sensitivity. A finer aerodynamic evaluation or a targeted CFD comparison can check a shortlisted result. The CFD app supplies aerodynamic evidence at recorded poses and conditions. The coordinator interprets that evidence alongside whole-aircraft performance.

### When the rough layout misses a requirement

The first CAD layout is a starting estimate. During study setup the designer identifies which preliminary dimensions the coordinator may resize or move, the permitted ranges, and which features are fixed. If a required flight case misses a modeled limit, first search for layouts that meet the requirements. Rank those layouts on the mission objective afterward.

Identify the case, the failed check and its margin. For example, at minimum speed distinguish an estimated section lift limit from insufficient pitching moment at the allowed elevator travel. If the solver fails to converge, check the force and moment balance across the permitted angle and control range before assigning a cause. Record an unresolved case when the available model or numerical evidence is insufficient.

- **Wing example.** Compare additional wing area with a planform or twist change at fixed area. This tests whether lowering the overall lift demand per unit area or redistributing it along the span addresses the limiting section.
- **Trim example.** Compare increasing elevator chord or span with changing wing position or tail incidence. Record required deflection, remaining travel and static margin across the declared CG and speed cases. Check control effectiveness within the model's range.
- **Coupled changes.** Use single changes to understand the response and test combinations when one change cannot satisfy all cases. Recompute dependent mass, CG, trim and spar fit for each candidate and identify estimates. Report the least extensive tested change that passes and any useful alternatives. If the search finds none, report the ranges explored and propose a specific freedom to revisit with the designer.

## Using CFD evidence in an optimizer result

Use CFD as a comparison by default. A model update must name the contribution it replaces or corrects, its source run, geometry revision, configuration, operating range, and numerical uncertainty or outstanding checks. Keep the original model result alongside the updated result.

- For a 2D section, match the profile coordinates, Reynolds number, surface condition and incidence or section lift range. Apply supported sectional lift, drag or moment data only to the corresponding section-model terms over that range.
- For a whole aircraft, match the assembled geometry, control positions and flight condition. Reconcile units, force axes and moment reference. A CFD total replaces the corresponding total or supplies an explicitly defined difference against the original model at the same condition. Count each contribution once.
- Re-trim when a change affects force or moment balance, and recompute dependent power and mission results. If the new trim state leaves the supported range, report the corrected prediction as unresolved and retain the CFD measurement at its actual pose.

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
- Import a swept reference wing with parallel ribs that depart from the default guide normal. Preserve the planes during baseline evaluation and a locked profile swap, or report unsupported analysis. Check that any authorized refit has a separate candidate ID.
- Import a continuous tail, evaluate its pitching moment, then add a proposed elevator in the study and evaluate trim. The report must identify the elevator as an analysis assumption and show the travel used.
- Start a study with a known control-travel shortfall. Recover the limiting case and compare two authorized remedies across the full set of flight cases. Include a run that fails numerically and verify the report marks its outcome as unresolved.
- Exercise a flat wing, a smooth upsweep, and a wing with a ripple between ribs. Locate the ripple. Check both a tilted tip and a raised tip whose guide levels off.
- Exercise the spar engagement boundary at 20 mm and just below it. Include an oblique cut and a square spar whose centerline clears but whose corner intersects the skin.
- With the 1.0 mm minimum active, test a square channel whose closest corner leaves exactly 1.0 mm of wall and one just below it. Include the stock-to-channel allowance. Verify that a thicker two-perimeter requirement from a print profile controls the verdict when supplied.
- Swap to a profile that fails at the baseline spar axis but admits the same stock at another position. The candidate must pass with the verified new axis and retain its locked external geometry.
- Rebuild a shortlisted wing with a local deviation inside the permitted loft tolerance that consumes the spar clearance. Verify the report retains the generated fit and flags the CAD interference at its location.
- Run left and right aileron and rudder commands on a symmetric reference aircraft with independently specified expected surface motions and moment signs. For a mixed tail, include commands that fit separately but exceed one named surface's travel when combined. Check that the report identifies that surface and the exceeded limit, or declares mixed controls unsupported.
- Replay a CFD result in comparison-only mode and verify the original power estimate is preserved. Then replace whole-aircraft drag at a matched condition with a known value and verify total drag equals that value. Exercise a section correction outside its Reynolds range and a moment correction that requires re-trim.
- Compare an agent-guided study with a fixed numerical search using the same freedoms, constraints and evaluation budget. Judge feasible performance, persistence of the gain under cross-checks, and which evidence helped the designer decide.

## Decisions still open

- How the CAD revision, section profiles, guide curves and locks travel between the Converter and optimizer, and how custom intermediate guides change the section definition used for aerodynamics.
- Numerical tolerances for baseline reconstruction, orientation, curvature, surface deviation and spar clearance. The bend-radius and twist-rate limits need a physical basis for each design.
- Source of print orientation and cut planes, whether studies can vary them, and which available stock sizes, wall thicknesses and materials to supply.
- Which aircraft revision and flight case to use first, which performance gain matters to the designer, and when the remaining uncertainty justifies commissioning CFD.
