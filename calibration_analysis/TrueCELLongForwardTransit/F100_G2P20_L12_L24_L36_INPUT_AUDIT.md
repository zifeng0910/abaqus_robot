# F100 G2P20 L12/L24/L36 Input and Initialization Audit

## Cases and source identity

All three `case_identity.json` files report `status=SOLVED`, `dynamics_run_count=1`, `fluid_mode=TRUE_CEL`, `FAST_SURROGATE_FLUID=OFF`, `ReducedHydro=OFF`, `restart_read=false`, and the same source commit `a6dbd1422e898e708833bd427e4faa8860ee700e`. The shared V/VR driver SHA256 is `4EF12426EE062B576DB9AF02E12ECA750EEBC3AC8A8AE1C948C093EDAF66C8D8`. Input hashes are L12 `C439C4A76ED84F18BF7AAD04350BC8572E2A3CACDC08B94DB1503093D71728A1`, L24 `EB2BCB1B246253E1C70540D81D86B4FE433C677522982E1935BF55D2901C23D3`, and L36 `E0DE815D34EAF97A08209828D70ABA232B4C813AFD1EDF2D7A18D321D8D3095B`.

## Eulerian mesh and extent

The actual fluid part in each INP is `EC3D8R` with 20 x 20 cells in the two cross-sectional directions and spacing `0.075 x 0.075 mm`. Axial spacing is `0.272727272727 mm`. The actual fluid counts are:

| Case | Axial cells | Fluid elements | Fluid nodes | Axial extent declared by deck |
|---|---:|---:|---:|---|
| L12 | 44 | 17,600 | 19,845 | `[-6, 6] mm` |
| L24 | 88 | 35,200 | 39,249 | `[-12, 12] mm` |
| L36 | 132 | 52,800 | 58,653 | `[-18, 18] mm` |

The decks retain the same cross-sectional mesh and axial spacing. L24 and L36 comments and element counts show only axial extension changes. The generated node coordinates are in the rotated Abaqus frame; their raw coordinate ranges therefore do not equal the canonical axial interval. The canonical intervals above are the deck construction parameters and case identities.

The wall helper is extended with each domain. The robot-only low/high terminal stops are separate rigid surfaces and are explicitly excluded from the CEL material contact; they cannot serve as fluid outlets.

## Initial conditions and material

Each deck contains `*Initial Conditions, type=VOLUME FRACTION` with the actual assembly elset `FLUID_CEL_INIT_ALL, Fluid_Eulerian-1_WATER, 1.0`. Case identities report full-domain convex-hull/cell-exclusion initialization, no local axial shortcut, zero robot-fluid initial overlap, and resolution bound `0.0327872 mm`. The initialized fluid volumes/masses are:

| Case | Filled volume (mm3) | Fluid mass (mg) | Initial positive EVF minimum | Initial EVF maximum |
|---|---:|---:|---:|---:|
| L12 | 12.1868182 | 12.1953490 | 0.008 | 1.0 |
| L24 | 26.2268182 | 26.2451769 | 0.008 | 1.0 |
| L36 | 40.2668182 | 40.2950050 | 0.008 | 1.0 |

The decks also initialize all Eulerian nodes with zero velocity. Material data are the same: density `1.0007e-09 tonne/mm3`, viscosity `7.1e-10 N s/mm2`, USUP EOS, and `c0=100000 mm/s` in the case identity. No initialization inconsistency was found in the supplied labels or reported volumes.

## Boundary conditions and contacts

The only explicit `*Boundary` in the step is `RP_PIPE, ENCASTRE`. The deck comments state `axial CEL ends: unconstrained Eulerian material-flow boundaries; no explicit pressure outlet`. No Eulerian end node set is constrained, and no pressure, inflow, outflow, non-reflecting, acoustic, reservoir, or reference-pressure keyword was found in these decks. The robot-only terminal stops are not fluid boundary conditions.

The actual `*Contact Inclusions` are:

- robot exterior surface with the Eulerian material surface;
- wall-helper surface with the Eulerian material surface.

The pipe-solid/robot direct wall pair is present in the general contact inventory, while the robot-fluid and wall-fluid assignments use `PROP_CEL_FLUID_HARD`. The selected softened contact law is retained for robot-wall only. No separate Eulerian face traction output was accepted by datacheck; the requested second-surface route was rejected.

## Numerical controls and output limits

All three decks use `*Dynamic, Explicit, SCALE FACTOR=0.4`, bulk viscosity, no mass scaling, and `*Restart, write, number interval=10, time marks=YES`; no restart read is present. History output is at `1e-5 s`; field output is sparse relative to the wave transit and contains velocity, EVF, and water stress. Native `PRESS` and face mass-flux histories are absent. The reported pressure is therefore only `p_proxy=-trace(S_water)/3`.

## Audit conclusion

The byte/semantic audit finds no evidence that L24 or L36 changed the initialization recipe, cross-sectional mesh, EOS, density, viscosity, contacts, timestep control, or prescribed V/VR driver. The unresolved issue is the open Eulerian axial end treatment and the lack of an independently defined fluid-force partition, not a discovered initial-condition inconsistency.

