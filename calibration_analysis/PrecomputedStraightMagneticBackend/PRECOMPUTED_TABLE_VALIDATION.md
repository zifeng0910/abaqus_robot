# Precomputed Straight Magnetic Backend Validation

## Architecture

The selected representation is `FIELD_GRADIENT_TABLE`, not a pose/load or trajectory table. The production branch is the analytic, point-dipole-equivalent `ROBOT_LOCAL_ELLIPTIC_ROCKING` implementation in `magpylib_socket_server_fast.py`.

The normalized lookup is a separable tensor product:

- `B_unit(phase)`, sampled directly from the live production implementation.
- `grad_unit(s_eff)`, where `s_eff = robot_s - driver_s(t)` and `driver_s(t) = driver_s0 + 6 mm/s * t`.
- Runtime reconstruction uses the current Abaqus `UR` axis-angle rotation, `F = m dot grad(B)`, and `M = m x B`.
- The existing 1 ms smooth ramp is applied in compiled Fortran.
- No robot position, rotation, contact state, or trajectory is precomputed.
- The table is loaded once per Abaqus MPI process at initialization. No Python, TCP, subprocess, or Magpylib call occurs during table-backed dynamics.
- Leaving the `s_eff` range raises a clear stop diagnostic; long extrapolation is not allowed.

`MAGNETIC_BACKEND` and `FLUID_MODE` are independent in `straight_tube_case_config.json`. The validated short A/B uses `FLUID_MODE=NONE`; the same table is reusable with `FAST_SURROGATE` and later `TRUE_CEL`.

## Table

| Quantity | Value |
|---|---:|
| Phase range/resolution | 0..360 deg / 1 deg |
| Phase stations | 361 including wrap |
| Effective axial range/resolution | -70..70 mm / 0.25 mm |
| Effective axial stations | 561 |
| Logical `(s_eff, phase)` shape | 561 x 361 |
| Stored size | 172,996 bytes |
| Precomputation wallclock | 1.325 s |

The active analytic production model has no radial dependence. Direct center and +/-0.25 mm RouteA `n/b` checks measured zero variation, so no `(u,v)` dimensions are needed.

## Offline Validation

The 64-state random validation exceeded the requested minimum of 50 states. Maximum relative errors against the live production implementation were:

| Quantity | Maximum error | Gate |
|---|---:|---:|
| B | 0.0009175% | 0.5% |
| Gradient | 0.0009483% | diagnostic |
| Force | 0.0009483% | 3% |
| Torque | 0.006954% | 1% |

All gates passed at `ds=0.25 mm`; the one permitted refinement to 0.10 mm was not needed. B0 and G normalization is exact by construction in the production analytic branch.

## Dynamic A/B

Both cases used the same non-CEL rigid robot, tube, General Contact, initial pose, 120 Hz, 10 mT, 0.15 mT gradient setting, 45 mm gradient length, and fixed `dt=2e-7 s`. Only the magnetic backend differed. Duration was exactly 0.5 ms.

| Comparison | Measured maximum | Gate |
|---|---:|---:|
| Axial position difference | 2.461e-11 mm | 0.01 mm |
| Rotation-vector difference | 3.295e-6 deg | 0.2 deg |
| Force error across the two dynamic trajectories | 0.0008977% | 3% |
| Torque error across the two dynamic trajectories | 0.007618% | 1% |

The dynamic A/B passed.

## CPU Benchmark

| CPUs | Explicit main loop | Complete job | U/UR consistency |
|---:|---:|---:|---|
| 1 | 5 s | 18.762 s | reference |
| 4 | 8 s | 19.785 s | exact at exported points |
| 6 | 22 s | 33.829 s | exact at exported points |

`cpus=1` is selected. This small model is MPI-overhead dominated. The table was loaded independently by every rank and all runs completed, confirming process-safe read-only use.

A two-cycle 120 Hz screen is 16.667 ms. Linear extrapolation from the selected run gives 166.7 s of solver time, or about 180 s including observed fixed job overhead. Plan on approximately **3 minutes** for `PRECOMPUTED_TABLE + FAST_SURROGATE` on one CPU. This is an estimate; the benchmark itself used `FLUID_MODE=NONE`.

No long CEL job and no forward/reverse motion screen was run in this task.
