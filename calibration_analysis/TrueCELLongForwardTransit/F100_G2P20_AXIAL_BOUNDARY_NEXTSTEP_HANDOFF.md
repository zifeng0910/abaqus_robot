# F100 G2P20 Axial Boundary Next-Step Handoff

## Scope

This handoff closes the fixed-trajectory L12/L24/L36 axial-domain campaign and starts the CEL axial-boundary root-cause qualification. No new full robot dynamics solve is authorized until the input audit, force-partition cross-check, and Abaqus documentation review are complete.

## Prior campaign state

- L12, L24, and L36 all completed as one-step prescribed GLOBAL V plus spatial VR replays to 40 ms.
- The trajectories are effectively identical: L24 maximum position/orientation errors were `6.10998e-7 mm` and `5.36263e-5 deg`; L36 values were `5.90339e-7 mm` and `5.65169e-5 deg`. Both p99 axial-speed and rocking-speed errors passed their gates.
- L24 used 88 x 20 x 20 Eulerian cells over the doubled axial domain; L36 used 132 x 20 x 20. Cross-sectional spacing and physics were retained.
- The previous primary classification was `AXIAL_BOUNDARY_SENSITIVITY_UNRESOLVED`.

## Required impulse table

| Case | C1 (N s) | C2 (N s) | C3 (N s) | C4 (N s) |
|---|---:|---:|---:|---:|
| L12 | -6.0801e-08 | -6.5336e-07 | -1.6452e-05 | -1.3698e-06 |
| L24 | -3.6646e-08 | -7.5419e-06 | +6.4695e-07 | -1.4021e-06 |
| L36 | -2.3239e-08 | -7.0170e-07 | -1.0775e-05 | -3.1962e-06 |

C3 is non-monotonic, changing from negative to positive and back to negative as the domain is enlarged. The prescribed trajectories are identical, so this is not replay drift. L12's negative C3 term must not be described as validated physical recoil.

## Force interpretation gate

The reported robot-CEL quantity is `whole robot General Contact - direct robot-wall contribution`. It is an inferred residual, not a direct fluid traction. The cycle balance closes algebraically, but Abaqus did not provide a valid independent Eulerian robot-fluid traction, face mass-flux, or native pressure history. Overlapping head/tail contacts also prevent unique end-force attribution. Until those semantics are resolved, do not promote the inferred CEL residual to a physical recoil force.

## Next work

1. Complete `F100_G2P20_L12_L24_L36_INPUT_AUDIT.md` from the actual decks and case identities.
2. Complete the documented force-partition cross-check using the existing histories only.
3. Record the installed Abaqus 2025 documentation search and any supported Eulerian boundary semantics; do not invent keyword syntax.
4. Build at most one cheap straight-pipe reflection harness only after steps 1-3. Do not launch L48/L72/L100 or another full robot run.

