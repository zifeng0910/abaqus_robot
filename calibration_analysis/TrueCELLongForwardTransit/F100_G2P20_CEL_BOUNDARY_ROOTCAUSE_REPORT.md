# F100 G2P20 CEL Boundary Root-Cause Qualification

## Primary classification

`FLUID_FORCE_PARTITION_SEMANTICS_UNRESOLVED`

## Decision

The axial-domain campaign demonstrates real domain sensitivity in the inferred residual force, but it does not identify whether the cause is boundary reflection, initialization, or another CEL force-partition effect. The prescribed trajectories are matched, the initial fluid recipe is consistent, and the global robot balance closes. The required physical fluid-force quantity is not independently available in the current Abaqus output. Stop before any new full robot dynamics solve.

## Campaign summary

| Case | Axial domain | C1 (N s) | C2 (N s) | C3 (N s) | C4 (N s) | Kinematic gate |
|---|---|---:|---:|---:|---:|---|
| L12 | `[-6,6] mm` | -6.0801e-08 | -6.5336e-07 | -1.6452e-05 | -1.3698e-06 | archived pass |
| L24 | `[-12,12] mm` | -3.6646e-08 | -7.5419e-06 | +6.4695e-07 | -1.4021e-06 | pass |
| L36 | `[-18,18] mm` | -2.3239e-08 | -7.0170e-07 | -1.0775e-05 | -3.1962e-06 | pass |

L24 and L36 replay errors are below `0.001 mm`, `0.02 deg`, `0.05 mm/s` p99 axial speed, and `0.5 rad/s` p99 rocking speed. C3 is non-monotonic across the three lengths, so the L12 negative C3 term is not validated physical recoil. The prior classification `AXIAL_BOUNDARY_SENSITIVITY_UNRESOLVED` remains the correct campaign-level description, but the root-cause qualification now classifies the force semantics as unresolved.

## Input and initialization audit

The detailed audit is in `F100_G2P20_L12_L24_L36_INPUT_AUDIT.md`. The actual decks show:

- EC3D8R Eulerian meshes of 44/88/132 by 20 by 20 cells, with unchanged cross-sectional spacing `0.075 mm` and axial spacing `0.272727 mm`.
- Full-domain `FLUID_CEL_INIT_ALL` volume-fraction initialization, zero Eulerian nodal velocity, zero robot-fluid overlap, and consistent filled volumes and masses.
- The same density, viscosity, USUP EOS, `c0=100000 mm/s`, explicit scale factor, bulk viscosity, no mass scaling, and no restart read.
- Only the pipe reference point is encastre constrained. The axial Eulerian ends are unconstrained in the decks, with no explicit pressure, inflow/outflow, non-reflecting, acoustic, reservoir, or reference-pressure keyword.
- Robot-fluid and wall-fluid contact inclusions use the Eulerian material surface. Robot-only terminal stops are not fluid outlets.

No initialization inconsistency was found in the supplied input or case identities.

## Prescribed-motion force cross-check

For each cycle, the existing histories satisfy the signed balance

`reaction = inertia - magnetic - wall - inferred_fluid`

where `inferred_fluid = whole General Contact - direct robot-wall`. The cycle closure residuals are:

- L24: `3.51e-13`, `-7.93e-13`, `7.41e-13`, `-3.20e-13 N s`.
- L36: `-1.94e-13`, `5.12e-13`, `2.55e-12`, `3.80e-13 N s`.

The reported force-closure p99 residual/inertia ratios are `1.05193e-4` (L24) and `1.19734e-4` (L36). This verifies the internal prescribed-motion bookkeeping. It does not verify the physical partition because the independent Eulerian traction is missing. The L12/L24 history comparison has cycle correlations `-0.3786, 0.0227, 0.2263, 0.0417` and normalized RMS differences `1.5458, 1.0020, 0.9852, 1.1126`; these are domain-sensitivity diagnostics, not a substitute for direct traction validation.

## Official documentation review

The installed 2025 tree at `I:\SIMULIA\EstProducts\2025` was searched before proposing a boundary keyword. It contains CAADoc and runtime resources, but no local Abaqus Analysis User's Guide corpus or searchable Eulerian/CEL boundary reference with a documented pressure outlet, non-reflecting outlet, reservoir, or acoustic termination for this EC3D8R Explicit setup. Therefore no alternative keyword syntax or boundary semantics is asserted here. The deck's “unconstrained Eulerian material-flow boundary” comment is treated as model metadata, not as an independently documented non-reflecting condition.

## Harness gate

The straight-pipe reflection harness was not run. The force-partition cross-check reached the stop condition first: the current model has no independent fluid traction or face-flux observable, and the installed documentation set did not provide a verified alternative boundary formulation. Running a harness with the same ambiguous end treatment would only reproduce the ambiguity and would not validate a candidate boundary. No L48/L72/L100 replay and no wall-fluid-off replay were launched.

## Recommendation

Do not use the inferred robot-CEL residual as a validated physical recoil or make a larger-domain claim. Obtain a documented fluid-force and open-boundary formulation first. If the Abaqus 2025 documentation package can be installed or supplied, revisit one cheap straight-pipe harness with native pressure and face-flux outputs. If a suitable pressure/non-reflecting boundary is unavailable in Abaqus CEL, evaluate a CFD solver with explicit open-boundary control while retaining Abaqus for rigid-body and contact mechanics.

