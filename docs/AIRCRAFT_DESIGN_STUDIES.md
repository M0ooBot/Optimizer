# Optimizer studies for agent-led aircraft design

**Working draft · 2026-11-26.** This is the proposed role of the Optimizer in the aircraft design harness. The designer and agent may begin with a rulebook and no geometry. Later the same tool must evaluate a particular wing, tail and human-designed fuselage as drawn. This document defines what questions the tool should answer, the inputs and evidence needed, and how to tell whether the answers are usable. [Harness design process](https://github.com/M0ooBot/Engineering-Harness/blob/main/AIRCRAFT_DESIGN_PROCESS.md)

## 1. Contract and current starting point

The agent proposes an aircraft, chooses the next experiment and explains trade-offs. The designer controls study freedoms and adopts revisions. The Optimizer supplies deterministic sizing, performance calculations, comparison and bounded numerical search. It must support a plain evaluation that holds the entire geometry fixed. A designer asking whether their wing trims has authorized analysis of that wing.

The existing library already has mass and propulsion models, propeller performance tables, fixed-design speed sweeps, longitudinal trim, stability checks, airfoil studies and a gradient-based optimizer. The existing search starts from a Python aircraft definition and can run a scripted sequence of discrete studies. Those are useful calculations and search options to expose independently. The stages below describe the target interface. [Existing concept](OPTIMIZATION_CONCEPT.md) · [Model details](MODEL_DETAILS.md) · [Reviewed solver revision](https://github.com/M0ooBot/Optimizer/blob/1c6301660c08fbee26e771fab6d461602099fd58/src/planeopt/solve.py)

| Reuse or adapt | Capability the harness still needs |
|---|---|
| Fixed-geometry pitch trim, mass estimates and propulsion operating points | Candidate inputs from WingBuilder and Systems, arbitrary installation poses, loading-case sweeps and an agent-callable evaluation |
| Airfoil comparisons and numerical optimization | Provisional profile screening from an estimated operating range, assessment of the as-built profile, and bounded studies on editable regions of an arbitrary supported wing |
| Existing mission objectives and speed sweeps | Concept sizing from partial inputs, separate launch methods, climb and complete mission energy accounting |
| Static stability and pitch-control travel checks | Physical control definitions, mixed commands, roll and yaw authority studies, with separate evidence for time-response requirements |
| Run artifacts and diagnostics | Results linked to design, geometry, data and models across tools, comparison views, and condition-specific CFD updates |

The prior code review recorded concrete regressions to carry forward, including incompatible checkpoint reuse, disagreement between solver and final sweep limits, changed propeller data retaining a fingerprint, and regenerated geometry ignoring saved design variables. These matter directly to the new interface. [Optimizer issue tracker](https://github.com/M0ooBot/Optimizer/issues)

## 2. Design data crossing the tool boundary

A study identifies the aircraft revision, mission or scoring rule, flight and loading cases, objective or comparison question, authorized changes, hard limits and expected output. Every physical quantity crossing tool boundaries has units and a coordinate convention. Persist the transformation to the fixed aircraft datum used by the board and CAD.

- **Geometry.** Read the actual section profiles, placements and planes, planform and skin guides, incidence, twist and control definitions from the candidate. Link the geometry revision used to create the analysis mesh. An early layout may consist of lifting surfaces, a provisional body estimate and dimensioned equipment envelopes. A later evaluation may reference the designer's CAD assembly.
- **Parts and poses.** Systems owns the shared parts list and mass ledger. The aircraft workbook references the same component IDs. A candidate places instances in 3D with an orientation, shape envelope and local mass centre where known. Two identical motors use two instances of one component definition. Payload and battery positions may vary by loading case. When the designer places the real component in CAD, use that pose for the assembled revision.
- **Propeller performance.** Package the propeller tables and their evaluation in the Optimizer. Keep the part's identity, mass and physical dimensions in the shared parts record with a reference to its performance data. Record which table, units, source, test or fit range and override revision were used. A flight-condition query combines that table with motor, ESC and battery data to return thrust or required RPM, current, voltage, electrical power and applicable margins.
- **Installation.** Accept a motor axis, thrust application point and propeller disc at any candidate pose. The board or geometry service checks swept-disc clearance. Use the force direction and moment arm about the case's CG in aircraft equilibrium. A free-air propeller table can support a free-air prediction while a wake-interaction claim remains open. Representing a pose geometrically does not itself establish an installation correction.
- **Early mass estimates.** Use process-specific evidence where available: skin mass per unit surface area, rib and joint estimates, spars, hardware and local allowances. Count each item once. Carry estimate ranges, print assumptions and provenance. Recalculate mass and its location when shape changes; a rigid translation moves the mass centre while preserving component mass. Sliced estimates and measured masses can subsequently replace assumptions for new revisions.

## 3. First call with no aircraft definition

Given a partial brief, estimate ranges for flying mass, wing area, installed power and usable energy. Accept known payload and hardware, mission segments, launch or runway conditions, span or packed-size limits, and bounded estimates for missing inputs. A designer may bring dimensions of a payload box before knowing its final station.

For each assumed takeoff mass and maximum lift coefficient, estimate the minimum wing area for the stall requirement. Combine this with applicable launch, climb and cruise power requirements. Estimate energy by segment; update battery and airframe mass assumptions, then repeat or expose the range over which the mass and energy budgets close. Report the controlling condition and sensitivity to uncertain inputs. If a segment cannot yet be evaluated, show its missing contribution and withhold any dependent mission verdict.

The result should support a concrete agent decision such as requesting two WingBuilder concepts at comparable wing area and with different mounting positions. It must also support the outcome that several concepts remain credible, or that the declared mass, energy and operating limits conflict. The choice of layout stays with the agent and designer.

## 4. Profiles and curved lifting surfaces

**Choosing a provisional profile.** Estimate required lift and speed and chord ranges before detailed geometry exists. Research a documented reference profile for the surface's role. Check data at relevant root and tip Reynolds numbers, drag across the working lift range, pitching moment and available thickness at the spar location. Tails and canards deserve their own operating assumptions. Record the exact coordinates, data source, coverage, and the condition that would cause the selection to be reopened. Where measured data is absent, label predictions and extrapolations. Match print-related changes such as a blunt trailing edge to the profile actually passed to the analysis.

**Geometry ownership.** Airfoil Converter prepares and places profiles. WingBuilder defines lifting surfaces and the SolidWorks bridge rebuilds shortlisted lofts. The analysis adapter must preserve dihedral, sweep, varying twist, orientation and offsets along curved surfaces, including what happens between specified stations. Return the derived analysis mesh and its discrepancy from the source at a stated sampling resolution. A finer mesh or alternative method is appropriate when a design decision turns on a bend the current analysis poorly resolves. [WingBuilder study specification](https://github.com/M0ooBot/Wing-Builder/blob/main/docs/SHAPE_OPTIMIZATION.md)

**Freedom to change shape.** The agent can create sections or change guides in an editable region. It can try two wings with identical projected span, area, and root and tip chords and still produce materially different sweep, section loading and curvature. A particular numerical study may free only onset of a bend or a few guide coordinates. That choice belongs to that study. No refit to a fixed wing family may silently replace the geometry that was submitted.

**Comparing profiles and configurations.** Re-trim shortlisted profiles on the aircraft, taking their effect on required trim, mass estimates, clearance and mission performance into account. Hold geometry fixed during a profile-only comparison and identify any change that forces a separate design revision. For a configuration trade, give each candidate its own reasonable dimensions and controls under common mission and hardware constraints. The existing greedy order of discrete studies is one optional search strategy. The agent can also compare coherent alternatives directly and reopen coupled choices.

## 5. Fixed-geometry studies and mission performance

For a flight condition and loading, compute the CG, find the equilibrium within the supported model, and report its residuals and applicable limits. Existing longitudinal trim solves for angle of attack and pitch-control deflection at specified airspeed, weight and longitudinal CG. Its propulsion match supplies the thrust required to offset drag. The proposed installed-propulsion model must also account for thrust direction and pitching moment. A raised pusher and a nose puller need separate installation inputs even when their propeller tables match.

Every numerical path claiming a legal operating point must apply consistent voltage, current, placard, flight-envelope, stability and control-travel limits as applicable. Record a converged point that exceeds a limit as a violation. Preserve a distinction between failure to converge, a substantiated violation, and a condition the model cannot assess. Show the governing case and margin, including very small margins that could reverse under the stated uncertainty.

| Mission segment | Inputs and result expected |
|---|---|
| Hand launch, when specified | Assess a realistic range of release speed, attitude, height and motor state. Determine whether the model supports a transition into controlled climb within the available height and control limits. A steady trimmed point at climb speed alone leaves the launch verdict open. |
| Ground takeoff, when specified | Assess runway length, surface, wind, slope and obstacle requirement where applicable. Report acceleration, rotation/liftoff, distance and margin. Identify unsupported ground or transition effects. |
| Climb | Start from a compatible departure state. Assess required thrust and power, current, climb capability, time and energy to the target altitude. |
| Transit, on station and return | Use the appropriate loading, wind, achievable ground progress, time or distance, aerodynamic condition and electrical loads for each segment. Distinguish time on station from total endurance. |
| Approach and landing | Check applicable stall margin, trim and control authority over the required loading conditions. Record any additional landing requirement that needs a separate model or test. |

Honor the brief's launch logic: hand launch, runway departure, both, or either. Sum segment energy and apply the required reserve once. For powered phases, preserve the flight condition and installation assumptions. An energy-per-distance objective needs positive progress over the intended ground track. A score based on completed missions must respect mission feasibility before ranking payload or other rewards.

## 6. Settling surfaces and controls ahead of the fuselage

Run placement and shape studies against a fixed datum and explicit body estimates. Compare wing position and incidence, tail or canard arm and placement, and surface shape with mass and CG updated for every revision. The designer may lock a root, hinge, outer curve or sightline and delegate the surrounding region. Return viable ranges as well as points so the human knows what freedom remains for packaging.

A control definition belongs to a specific parent surface and revision. It includes its outline, hinge axis, inboard and outboard boundaries, positive direction, travel limits and mixing. Downstream outputs identify the physical left or right surface, commanded and achieved deflection, and remaining travel after trim. Geometry tools check swept volumes on a CAD rebuild. The analysis model assesses forces and moments where it has coverage.

| Study | Required comparison |
|---|---|
| Elevator or movable canard | Trim deflection across speed and CG, and remaining nose-up and nose-down authority from each trimmed point |
| Ailerons | Rolling and yawing moments for both command directions, differential where relevant, and limiting outboard section conditions |
| Rudder | Yawing authority at the specified sideslip or command case and the condition that limits it |
| V-tail or elevons | Per-side travel under simultaneous commands, including baseline trim. Independent single-axis passes do not establish mixed-command feasibility. |
| Flaps, when present | Lift, drag, pitching-moment change, revised trim and remaining control authority for the applicable segment |

Static stability, control effectiveness and dynamic response answer different questions. A demanded roll time, gust recovery or hand-launch trajectory requires appropriate inertia, dynamics and control evidence. When those are unavailable, report the static result and leave the dynamic requirement open.

At acceptance, return the precise geometry and control revisions, the loading and CG range, flight cases and margins, source assumptions, and unresolved checks. CAD supplies the root sections, placement references and swept control envelopes the designer uses to build the fuselage. The Optimizer never needs to author the fuselage to do this work.

## 7. Re-evaluate the designer's assembly and reconcile CFD

Hold the returned CAD revision fixed. Check geometry readback, refresh mass and loading from Systems, and repeat comparable conditions from the pre-fuselage baseline. Explain deltas in CG, stability, pitch-control travel, power and reserve along with changes in model coverage. If a change exceeds an accepted margin, return the finding to the designer before exploring a new candidate.

CFD is appropriate when the affected physics is local or beyond the quick model's useful range. Examples include root-junction separation and pod wake at the tail. Supply the actual geometry, pose, case, speed and appropriate trim or matched-lift comparison. CFD returns a studied quantity or delta, its reference convention, numerical uncertainty or sensitivity, and applicable conditions. Record any correction applied in the Optimizer, identify the estimate it replaces, and rerun the dependent cases. Do not apply a cruise drag correction across a launch or high-incidence case without supporting evidence.

## 8. Agent operations and evidence on the board

The interface should let a coordinator inspect available models and valid ranges, request partial-input sizing, evaluate a fixed candidate, screen airfoils or propulsion choices, run matched sweeps, refine a bounded family, compare saved revisions, and re-evaluate with updated data. These are proposed operations rather than commitments to particular command names. Placement edits and new WingBuilder definitions create candidate revisions through their owning tools.

Every study returns a concise verdict with supporting detail on demand. Include the target, predicted value and unit, margin, governing case, source revisions, assumptions and checks outside the model. Return geometry overlays and spanwise data where useful. Preserve solver status, meshes, inputs and provenance for replay. A plot or comparison on the board must remain attached to the actual geometry, model settings, propeller data and mass case used in that run. A change in any of them invalidates the relevant reuse decision. Keep the previous result as a historical comparison.

## 9. Acceptance scenarios for the redesign

1. **Start from a blank aircraft.** With a payload, mission and uncertain construction mass, produce a sizing range and the controlling assumption. The agent can explain why it requests the first wing geometry.
2. **Preserve curvature.** Submit two wings with identical gross dimensions and different curvature between stations. Overlay each solver mesh on its corresponding loft and quantify the readback error. A small performance difference is an acceptable result if the model actually evaluated both shapes.
3. **Exercise hardware placement.** Move the same battery between two allowed stations. Total battery mass stays fixed; CG and the affected trim cases update. Relocate an otherwise identical motor above the CG and verify the expected change in thrust moment against an independent force-times-arm calculation.
4. **Expose a control limit.** Use a mixed surface whose individual pitch and yaw commands fit but whose combined command exceeds travel on one side. Name the physical surface, limiting case and amount. Keep CAD interference as a separate check.
5. **Fail a complete mission honestly.** Exercise two candidates with valid cruise points, one of which misses landing reserve. Then exercise a required launch case with missing evidence. The results distinguish energy failure, unresolved departure and a demonstrated pass.
6. **Replay the accepted revision.** Save, rebuild and re-evaluate a candidate with non-default and curved geometry. Check that the geometry scored is the geometry displayed and exported. A new prop table or an iteration-limited smoke run cannot contaminate a later full result.
7. **Close the fuselage loop.** Import a human CAD revision that moves the battery and alters a fairing. Re-evaluate as drawn, compare to the accepted baseline, attach a CFD delta for one supported condition and verify that only the applicable aircraft predictions change.

The outcome is a reproducible answer to a designer's question, with enough evidence to adopt a candidate, keep an alternative, or identify the next discriminating study.
