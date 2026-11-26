# Optimizer studies for agent-led aircraft design

**Working draft · 2026-11-26, revised the same day after review.** This is the proposed role of the Optimizer in the aircraft design harness. The designer and agent may begin with a brief and no geometry. Later the same tool must evaluate a particular wing, tail and human-designed fuselage as drawn. This document defines what questions the tool should answer, the inputs and evidence needed, and how to tell whether the answers are usable. [Harness design process](https://github.com/M0ooBot/Engineering-Harness/blob/main/AIRCRAFT_DESIGN_PROCESS.md) · [Decisions from the review](https://github.com/M0ooBot/Engineering-Harness/blob/main/CHARTER.md)

Each section separates two things. **The Optimizer provides** is what gets built into the app and its MCP server. **The agent** is instruction for whoever runs the study, and it is not built into the tool.

## 1. Contract and current starting point

The agent proposes an aircraft, chooses the next experiment and explains trade-offs. The designer sets delegations and adopts revisions. The Optimizer supplies deterministic sizing, performance calculations, comparison and bounded numerical search. It must support a plain evaluation that holds the entire geometry fixed. A designer asking whether their wing trims has authorized analysis of that wing.

The existing library already has mass and propulsion models, propeller performance tables, fixed-design speed sweeps, longitudinal trim, stability checks, airfoil studies and a gradient-based optimizer. The existing search starts from a Python aircraft definition and can run a scripted sequence of discrete studies. Those are useful calculations and search options to expose independently. The stages below describe the target interface. [Existing concept](OPTIMIZATION_CONCEPT.md) · [Model details](MODEL_DETAILS.md) · [Reviewed solver revision](https://github.com/M0ooBot/Optimizer/blob/1c6301660c08fbee26e771fab6d461602099fd58/src/planeopt/solve.py)

| Reuse or adapt | Capability the harness still needs |
|---|---|
| Fixed-geometry pitch trim, mass estimates and propulsion operating points | Candidate inputs from WingBuilder and the aircraft workbook, arbitrary installation poses, mass-case sweeps and an agent-callable evaluation |
| Airfoil comparisons and numerical optimization | Provisional profile screening from an estimated operating range, assessment of the as-built profile, and bounded studies on editable regions of a supported wing |
| Existing mission objectives and speed sweeps | Concept sizing from partial inputs, ground takeoff, climb and complete mission energy accounting |
| Static stability and pitch-control travel checks | Physical control definitions, mixed commands, roll and yaw authority studies, with separate evidence for time-response requirements |
| Propeller tables | A catalog of exact purchasable propellers, motors, ESCs and batteries |
| Run artifacts and diagnostics | Results linked to design, geometry, data and models across tools, comparison views, and condition-specific CFD updates |

The prior code review recorded concrete regressions to carry forward, including incompatible checkpoint reuse, disagreement between solver and final sweep limits, changed propeller data retaining a fingerprint, and regenerated geometry ignoring saved design variables. These matter directly to the new interface. [Optimizer issue tracker](https://github.com/M0ooBot/Optimizer/issues)

**Phased out or converted.** Some existing parts conflict with the decisions above. They are dropped or changed as the new interface is built.

| Existing part | What happens to it |
|---|---|
| The Optimizer's own equipment list | Converted: it reads the aircraft workbook, so there is one parts list |
| Spar fit, stress and stiffness constraints | Dropped: WingBuilder owns the spar rule, and strength waits for FEA |
| Parametric fuselage generator and design-brief generator | Dropped from design use: the designer owns the fuselage |
| T-tail option | Dropped from the first version |
| STEP import and fuselage drag check | Kept, and reused in section 7 |

### What the first version covers

- **Configurations.** A conventional tail and a V-tail. Each further configuration needs its geometry in WingBuilder, its trim and control model here, and its test cases, so the others are added one at a time.
- **Sizing** from a partial brief, shown as a sizing chart, and checks of complete power sets from the catalog.
- **Evaluation of a fixed candidate:** mass and CG for each mass case, pitch trim, static stability and the default control criteria in section 6.
- **Mission segments:** ground takeoff, climb, transit and return, time on station, whole-mission energy and approach. Hand launch gets screening criteria only.
- **Left open.** A requirement stated in time, such as a roll time, a gust recovery or a launch transition. CFD can improve the inputs to these questions. It does not answer them, and flight test is not part of the plan for now.

The aerodynamic method is not fixed by this document. It is settled when that code is built, and it decides how much of sections 4 and 6 the first version can answer. Until then, a question the model cannot represent returns "cannot assess" and never a number.

The structural envelope, wing loads for FEA and spar sizing are withdrawn from this tool until FEA development starts.

## 2. Design data crossing the tool boundary

A study identifies the aircraft revision, the requirements and primary measure from the intent file, the flight cases and mass cases, the comparison question, the delegation it runs under and the expected output. Every physical quantity crossing tool boundaries has units and a coordinate convention. Persist the transformation to the fixed aircraft datum used by the board and CAD.

**The Optimizer provides** a reader for each of these inputs and keeps the propulsion catalog.

- **Geometry.** Read the actual section profiles, placements and planes, planform and skin guides, incidence, twist and control definitions from the candidate. Link the geometry revision used to create the analysis mesh. An early layout may consist of lifting surfaces, a provisional body estimate and dimensioned equipment envelopes, all built in WingBuilder. A later evaluation may reference the designer's CAD assembly.
- **Parts and poses.** The aircraft workbook keeps the parts, masses, positions and mass cases. A candidate places instances in 3D with an orientation, shape envelope and local mass centre where known. Two identical motors use two instances of one component definition. Payload and battery positions may vary by mass case. When the designer places the real component in CAD, use that pose for the assembled revision. Until Systems is built, the Optimizer checks the workbook each time it reads it: units, each component once, and measured masses kept beside estimates.
- **Propulsion catalog.** The Optimizer keeps the catalog of propellers, motors, ESCs and batteries and evaluates them. Every entry is an exact part that can be bought. An exact part may take its performance data from a labelled stand-in, such as a folding propeller modelled by a fixed-blade table with a derate. Each entry carries the maker's part number, the source of its data, units, the test or fit range, any override revision and a status: owned, buyable or discontinued. The workbook keeps the identity, mass and dimensions of the parts chosen for a design and refers to their catalog entries. A flight-case query combines the propeller table with motor, ESC and battery data to return thrust or required RPM, current, voltage, electrical power and applicable margins. The Optimizer never evaluates a part that is not in the catalog. When a study needs one the catalog lacks, it reports the gap: what was needed and the nearest entries. The catalog starts small, with parts the designer owns or would buy. The designer approves every addition, and its source is recorded. The propeller tables that already ship with the library stay as data, and a propeller becomes selectable when it is approved into the catalog.
- **Installation.** Accept a motor axis, thrust application point and propeller disc at any candidate pose. WingBuilder checks swept-disc clearance and the board shows the result. Use the force direction and moment arm about the case's CG in aircraft equilibrium. A free-air propeller table can support a free-air prediction while a wake-interaction claim remains open. Representing a pose geometrically does not itself establish an installation correction.
- **Early mass estimates.** Use process-specific evidence where available: skin mass per unit surface area, rib and joint estimates, spars, hardware and local allowances. Count each item once. Take the spar's stock, position and length from the WingBuilder candidate that passed the spar rule. Recalculate mass and its location when shape changes; a rigid translation moves the mass centre while preserving component mass. Sliced estimates and measured masses can subsequently replace assumptions for new revisions.
- **Ranges.** Every estimated input carries a low, an expected and a high value with its source. A measured input carries its value and source. Section 5 says how the ranges reach a verdict.

**The agent** assembles the study from the intent file and the delegation. It keeps to the delegation's limits itself; the Optimizer does not enforce them. It says which inputs are estimates and gives their ranges and sources. Before the run it writes down the answer it expects and how it got there, and the Optimizer stores that prediction beside the result.

## 3. First call with no aircraft definition

**The Optimizer provides** sizing from a partial brief.

- **Least it needs.** The payload mass, the endurance or distance with a cruise speed, and the departure: the available ground roll and surface, or a hand launch. It also takes known hardware, a stall or approach speed limit, span or packed-size limits, and bounded estimates for missing inputs. A designer may bring dimensions of a payload box before knowing its final station.
- **Defaults for what the brief lacks.** A table of default values for the maximum lift coefficient, the drag level, the airframe's share of total mass for printed construction, the efficiency of the propulsion chain and the usable share of battery energy. The usable share is where the landing reserve is counted, so the two are one budget. Each default has a low, an expected and a high value and a source. The designer's own past printed aircraft are the preferred source. A textbook value is labelled as one.
- **The calculation.** For each assumed takeoff mass and maximum lift coefficient, estimate the minimum wing area for the stall requirement. Combine this with the applicable takeoff, climb and cruise power requirements. Estimate energy by segment, update battery and airframe mass assumptions, then repeat or expose the range over which the mass and energy budgets close. Use the propulsion model at each segment's own speed, because a propeller is far less efficient at takeoff than at cruise.
- **The answer.** A sizing chart with weight per wing area on one axis and power per weight on the other. Each requirement draws one line, and the allowed region is what the lines leave. The power axis is battery power, and the stall line is drawn as a band over the range of maximum lift. Beside the chart, show the mass and energy closure as its own output, because endurance usually sets the size of a battery aircraft and is not a line on the chart. With it come the ranges of flying mass, wing area, installed power and usable energy over which the budgets close, the requirement that sets the size, and the sensitivity to each uncertain input. If a segment cannot yet be evaluated, show its missing contribution and withhold any dependent mission verdict.
- **Power-set check.** Given complete sets from the catalog, each a motor, propeller, ESC and battery, evaluate every set at takeoff, climb and cruise across the sized mass range. Report thrust, current, voltage, the limit of each part that is closest, and the mission energy. Say which sets cover the sizing range and what limits each one.

**The agent** uses the chart to decide what to try next, such as requesting two WingBuilder concepts at comparable wing area with different mounting positions. The result must also support the outcome that several concepts remain credible, or that the declared mass, energy and operating limits conflict. The agent proposes two or three power sets as a shortlist, and the designer picks one after the configuration geometry exists. When the Optimizer reports a gap in the catalog, the agent looks for a real part and proposes adding it. The choice of layout stays with the agent and designer.

## 4. Profiles and curved lifting surfaces

**The Optimizer provides:**

- **Profile screening.** For an estimated lift requirement and ranges of speed and chord, return section data at the relevant root and tip Reynolds numbers, drag across the working lift range, pitching moment and available thickness at the spar location. Label each value as measured, predicted or extrapolated.
- **An analysis adapter that keeps the submitted shape.** Airfoil Converter prepares and places profiles. WingBuilder defines lifting surfaces and the SolidWorks bridge rebuilds shortlisted lofts. The adapter must preserve dihedral, sweep, varying twist, orientation and offsets along curved surfaces, including what happens between specified stations. Return the derived analysis mesh and its discrepancy from the source at a stated sampling resolution. A finer mesh or alternative method is appropriate when a design decision turns on a bend the current analysis poorly resolves. No refit to a fixed wing family may silently replace the geometry that was submitted. [WingBuilder study specification, kept as reference](https://github.com/M0ooBot/Wing-Builder/blob/main/docs/SHAPE_OPTIMIZATION.md)
- **Profile and configuration comparisons.** Re-trim shortlisted profiles on the aircraft, taking their effect on required trim, mass estimates, clearance and mission performance into account. Hold geometry fixed during a profile-only comparison and identify any change that forces a separate design revision. For a configuration trade, give each candidate its own reasonable dimensions and controls under common mission and hardware constraints.
- **Bounded search on an opened family.** A numerical study may free only the onset of a bend or a few guide coordinates. The existing greedy order of discrete studies is one optional search strategy.

**The spar rule belongs to WingBuilder.** WingBuilder checks that one straight spar runs from the root through every print joint into the tip piece, within the print build height, with the required tip engagement and wall around the channel. The Optimizer does not repeat the check. It treats a wing without a passing spar check as not eligible for acceptance, and it takes the spar's mass and position from the passing candidate. WingBuilder runs the check as soon as a wing candidate is made. Only the spar rule of the linked spec applies here. [Spar rule](https://github.com/M0ooBot/Wing-Builder/blob/main/docs/SHAPE_OPTIMIZATION.md#straight-spar-and-printable-wing-sections)

**The agent:**

- **Chooses a provisional profile.** Estimate required lift and the speed and chord ranges before detailed geometry exists. Research a documented reference profile for the surface's role. Tails deserve their own operating assumptions. Record the exact coordinates, data source, coverage, and the condition that would cause the selection to be reopened. Match print-related changes such as a blunt trailing edge to the profile actually passed to the analysis.
- **Changes shape inside its delegation.** The agent can create sections or change guides in an editable region. It can try two wings with identical projected span, area, and root and tip chords and still produce materially different sweep, section loading and curvature. Which freedoms a numerical study opens is a choice for that study.
- **Compares coherent alternatives directly** and reopens coupled choices when the results call for it.

## 5. Fixed-geometry studies and mission performance

**The Optimizer provides** the evaluation of a candidate as drawn. For a flight case and mass case, compute the CG, find the equilibrium within the supported model, and report its residuals and applicable limits. Existing longitudinal trim solves for angle of attack and pitch-control deflection at specified airspeed, weight and longitudinal CG. Its propulsion match supplies the thrust required to offset drag. The proposed installed-propulsion model must also account for thrust direction and pitching moment. A raised pusher and a nose puller need separate installation inputs even when their propeller tables match.

Every numerical path claiming a legal operating point must apply consistent voltage, current, placard, flight-envelope, stability and control-travel limits as applicable. Record a converged point that exceeds a limit as a violation.

**Verdicts.** Each requirement gets one of three verdicts.

| Verdict | Meaning |
|---|---|
| Passes across the range | Passes at the expected values and at the ends of every estimated input that matters to it |
| Passes at expected, fails at an extreme | The result names the input and the value at which it fails |
| Fails | Fails at the expected values. The result gives the governing case and margin |

A run can also end in one of two states that give no verdict. **Did not converge** means the solver did not reach an answer, which is not a failure of the design. **Cannot assess** means the model does not cover the question. An open verdict is the same thing.

Estimated inputs are varied one at a time by default. Inputs that move together, such as mass and energy, are varied together. Show the governing case and margin, including very small margins. Model error is separate from input ranges: each output carries a stated band, marked as assumed until the checks in section 9 supply it. A margin smaller than that band is marked as within model error.

**Musts, targets and preferences.** A must that fails rules the candidate out. A must that passes only at the expected values, or that cannot be assessed, is flagged for the designer, who approves or rejects the candidate case by case. A missed target is reported and does not rule a candidate out. A preference only informs the comparison.

**Comparisons.** One candidate beats another only when it does so across the range. Otherwise the comparison reports the two as indistinguishable and says which input reverses the ranking and where.

**Mission segments.** The pass rules below are defaults and their numbers are rules of thumb. The intent file can replace any of them.

| Mission segment | What is computed | Default pass rule |
|---|---|---|
| Ground takeoff | Ground roll to liftoff at 1.1 × stall speed, from thrust, the rolling resistance of the stated surface, wind and slope. The elevator needed to rotate at liftoff speed with the most forward CG. The airborne distance to a required obstacle height | The ground roll fits the available distance, the elevator can rotate the aircraft, and a required obstacle is cleared. Ground effect is not modelled and is listed as not checked |
| Climb | Best climb rate at full throttle from the departure state, with current, voltage, time and energy to the target altitude | Meets the climb rate in the intent file within motor, ESC and battery limits |
| Transit and return | Trimmed power and current at cruise speed for the segment's mass case, with wind and ground progress | Trims, stays within continuous limits and makes headway along the intended track |
| On station | Power at the loiter speed for the time on station | Time on station meets the requirement. It is reported apart from total endurance |
| Whole mission | Energy summed across the segments | The landing reserve is met. The reserve and the usable share of the battery are one budget: by default 80% of nominal capacity is usable and the other 20% is the reserve, counted once |
| Approach and landing | Stall margin and trim at 1.3 × stall speed with the most forward CG, in each mass case | Trims within the elevator limit in section 6. Landing distance is not computed and is listed as not checked |
| Hand launch, screening only | Stall speed at launch mass in the launch setting, and thrust at the release speed | Release speed is at least 1.2 × stall speed and thrust at release is at least half the weight, a rule-of-thumb placeholder. The verdict stays open and is flagged for the designer |

Honor the departure logic in the intent file: ground takeoff, hand launch, both, or either. For powered phases, preserve the flight case and installation assumptions. An energy-per-distance measure needs positive progress over the intended ground track. Check every must-requirement before ranking candidates on the primary measure.

**The agent** picks the flight cases and mass cases from the intent file, writes its prediction, runs the evaluation and reads the verdicts. It carries the list of what was not checked into what it tells the designer.

## 6. Settling surfaces and controls ahead of the fuselage

**The Optimizer provides** placement and shape studies against a fixed datum and explicit body estimates. Compare wing position and incidence, tail arm and placement, and surface shape with mass and CG updated for every revision. The Optimizer makes the provisional body drag and mass estimate from the body envelope and shows it with a range. Return viable ranges as well as points so the human knows what freedom remains for packaging.

A control definition belongs to a specific parent surface and revision. It includes its outline, hinge axis, inboard and outboard boundaries, positive direction, travel limits and mixing. Downstream outputs identify the physical left or right surface, commanded and achieved deflection, and remaining travel after trim. WingBuilder computes swept volumes on the editable geometry and again on the CAD rebuild. The analysis model assesses forces and moments where it has coverage.

| Study | Required comparison |
|---|---|
| Elevator | Trim deflection across speed and CG, and remaining nose-up and nose-down authority from each trimmed point |
| Ailerons | Rolling and yawing moments for both command directions, differential where relevant, and limiting outboard section conditions |
| Rudder | Yawing authority at the specified sideslip or command case and the condition that limits it |
| V-tail | Per-side travel under simultaneous commands, including baseline trim. Independent single-axis passes do not establish mixed-command feasibility. |
| Flaps, when present | Lift, drag, pitching-moment change, revised trim and remaining control authority for the applicable segment |

**Default static criteria.** A study needs a target to return a verdict. These apply unless the intent file sets others. Where the existing library already has a default, it is used. The others are textbook starting values that have not been checked against these aircraft. The slowest speed is the approach speed, 1.3 × stall speed. The fastest is the top speed in the intent file, or the fastest flight case when it has none.

| Quantity | Default criterion |
|---|---|
| Pitch stability | Static margin between 8% and 15% of the mean aerodynamic chord in every mass case, checked at the most aft CG with power on. From the library |
| Elevator | Trims at the slowest speed with the most forward CG and at the fastest speed with the most aft CG, using no more than one third of the throw. From the library |
| Rotation at takeoff | The elevator can raise the nose at liftoff speed with the most forward CG at takeoff mass |
| Cruise stall margin | Cruise lift coefficient is no more than 0.7 of the maximum. From the library |
| Tip stall | No section exceeds its local maximum lift at the stall speed limit, and the tip Reynolds number stays above the library's floor. From the library |
| Ailerons | Full deflection reaches a wingtip helix angle of 0.07 at the slowest speed |
| Rudder | Holds the sideslip set by the brief's crosswind at approach speed, and at least 10°, with travel left |
| Yaw and roll stability | Vertical-tail volume of at least 0.030, from the library. Yaw stability of at least 0.04 per radian. The dihedral effect has the stabilising sign. All in every mass case |
| V-tail | Each surface stays within its travel under combined full pitch and yaw commands from trim. The rudder case is checked together with pitch trim |

Nose ballast, up to a stated limit, is an allowed way to bring the CG into range unless the designer locks it out.

Static stability, control effectiveness and dynamic response answer different questions. A demanded roll time, gust recovery or hand-launch trajectory requires appropriate inertia, dynamics and control evidence. When those are unavailable, report the static result and leave the dynamic requirement open.

At acceptance, return the precise geometry and control revisions, the mass cases and CG range, flight cases and margins, source assumptions, and unresolved checks. CAD supplies the root sections, placement references and swept control envelopes the designer uses to build the fuselage. The Optimizer never needs to author the fuselage to do this work. What the handoff should say about the drag and pitching moment the fuselage is allowed is open, to be answered on the first run through the process.

**The agent** proposes placements, shapes and control layouts inside its delegation. The designer may lock a root, hinge, outer curve or sightline and delegate the surrounding region. The agent shows a geometry clearance, an authority estimate and a dynamic requirement as three separate findings.

## 7. Re-evaluate the designer's assembly and reconcile CFD

**The Optimizer provides** the re-evaluation and the bookkeeping for corrections. SW-Agent exports the returned assembly, WingBuilder checks clearance through control travel on the export, and the new positions go into the workbook. The existing STEP import and fuselage drag check are reused here. Hold the returned CAD revision fixed. Check geometry readback, refresh the masses and mass cases from the workbook, and repeat comparable flight cases from the pre-fuselage baseline. Explain deltas in CG, stability, pitch-control travel, power and reserve along with changes in model coverage.

CFD returns a studied quantity or delta, its reference convention, numerical uncertainty or sensitivity, and applicable conditions. Record any correction applied in the Optimizer, identify the estimate it replaces, and rerun the dependent cases. Do not apply a cruise drag correction across a takeoff or high-incidence case without supporting evidence.

**The agent** decides when CFD is warranted. CFD is appropriate when the affected physics is local or beyond the quick model's useful range. Examples include root-junction separation and pod wake at the tail. The agent supplies the actual geometry, pose, case, speed and appropriate trim or matched-lift comparison. If a change exceeds an accepted margin, it returns the finding to the designer before exploring a new candidate.

## 8. Agent operations and evidence on the board

**The Optimizer provides** operations that let the agent inspect available models and valid ranges, request partial-input sizing, check power sets, evaluate a fixed candidate, screen airfoils, run matched sweeps, refine a bounded family, compare saved revisions, and re-evaluate with updated data. These are proposed operations rather than commitments to particular command names. Placement edits and new WingBuilder definitions create candidate revisions through WingBuilder and the workbook.

Every study returns a concise verdict with supporting detail on demand. Include the target, predicted value and unit, margin, governing case, source revisions and assumptions. Every result also carries a list of what it checked and a list of what it did not, and the agent's stored prediction, with a mismatch of more than about a factor of two flagged unless the intent file sets another threshold. Return geometry overlays and spanwise data where useful. Preserve solver status, meshes, inputs and provenance for replay. A plot or comparison on the board must remain attached to the actual geometry, model settings, propulsion data and mass case used in that run. A change in any of them invalidates the relevant reuse decision. Keep the previous result as a historical comparison.

**The agent** combines these operations to answer the designer's question and never edits a result. It keeps to its delegation's limits itself. It never marks a revision accepted: that is recorded only from the designer's own act.

## 9. Checks on the Optimizer itself

There is no reference aircraft to compare with, so the Optimizer is checked in ways that need none.

- **Textbook cases with exact answers.** Examples are the induced drag of an elliptically loaded wing, a thrust moment against force times arm, and the CG of a set of point masses. These run whenever the model changes.
- **Airfoil data.** Section predictions are compared with UIUC tunnel data where a measurement covers the profile and Reynolds number.
- **A CFD spot check on the candidate.** Before the designer accepts the lifting surfaces, the agent runs one or two flight cases of that candidate in the CFD workbench and the results are compared. At the configuration decision this is one of the second routes and is not required. A disagreement larger than the stated model band is flagged, and it becomes evidence for the band.
- **A second route.** The Optimizer can rerun a deciding comparison at a finer setting so the agent can show whether the ranking holds.

## 10. Acceptance scenarios for the redesign

1. **Start from a blank aircraft.** With a payload, mission and uncertain construction mass, produce a sizing chart, a sizing range and the requirement that sets the size. One line of the chart matches a hand calculation. The agent can explain why it requests the first wing geometry.
2. **Preserve curvature.** Submit two wings with identical gross dimensions and different curvature between stations. Overlay each solver mesh on its corresponding loft and quantify the readback error. A small performance difference is an acceptable result if the model actually evaluated both shapes.
3. **Exercise hardware placement.** Move the same battery between two allowed stations. Total battery mass stays fixed; CG and the affected trim cases update. Relocate an otherwise identical motor above the CG and verify the expected change in thrust moment against an independent force-times-arm calculation.
4. **Expose a control limit.** Use a V-tail whose individual pitch and yaw commands fit but whose combined command exceeds travel on one side. Name the physical surface, limiting case and amount. Keep CAD interference as a separate check.
5. **Fail a complete mission honestly.** Exercise two candidates with valid cruise points, one of which misses landing reserve. Then exercise a required hand launch that meets its screening criteria. The results distinguish an energy failure, an open departure verdict and a demonstrated pass.
6. **Replay the accepted revision.** Save, rebuild and re-evaluate a candidate with non-default and curved geometry. Check that the geometry scored is the geometry displayed and exported. A new prop table or an iteration-limited smoke run cannot contaminate a later full result.
7. **Close the fuselage loop.** Import a human CAD revision that moves the battery and alters a fairing. Re-evaluate as drawn, compare to the accepted baseline, attach a CFD delta for one supported flight case and verify that only the applicable aircraft predictions change.
8. **Refuse an invented part.** Ask for a power set the catalog cannot supply. The result reports the gap and the nearest entries and evaluates nothing made up. A discontinued part is flagged as one.
9. **Report a range honestly.** Evaluate a candidate that passes at the expected printed mass and fails at the high end. The verdict names that input. Two candidates whose ranking reverses inside the range are reported as indistinguishable, with the reversal point.
10. **Reproduce known answers.** The textbook cases pass within a stated tolerance. A prediction entered far from the result is flagged as a mismatch.
11. **Respect the spar rule.** A wing whose WingBuilder spar check failed cannot be marked eligible for acceptance, whatever its performance.
12. **Take off from the ground.** Compute the ground roll for a candidate and compare it with a hand calculation. Move the CG forward until the elevator can no longer rotate the aircraft at liftoff speed, and check that the result names that limit.
13. **Exercise the control defaults.** Run left and right aileron and rudder commands on a symmetric aircraft with independently specified expected moment signs. The roll, sideslip and yaw-stability criteria each return a verdict with its governing case.
14. **Mark stale results.** Change one requirement in the intent file. Every verdict tied to that requirement's ID is marked stale, and the others are not.
15. **Flag a must that does not plainly pass.** Evaluate a candidate whose must-requirement passes only at the expected values, and one with a required hand launch. Both are flagged for the designer and neither is passed or failed. The designer's approval is stored in the design record.

The outcome is a reproducible answer to a designer's question, with enough evidence to adopt a candidate, keep an alternative, or identify the next discriminating study.

## Open questions

- The aerodynamic method, to be settled when that code is built.
- Where motor, battery and ESC data for the catalog comes from.
- What the fuselage handoff says about allowed drag and pitching moment.
- The model error bands, which are assumed until the checks in section 9 supply them.
