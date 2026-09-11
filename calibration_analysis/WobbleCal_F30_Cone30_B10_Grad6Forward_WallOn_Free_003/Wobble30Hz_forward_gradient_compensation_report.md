# Wobble 30 Hz forward-gradient compensation

## Scope

This is the zero-cost replay and the single approved dynamic probe. The
Wall-ON free-translation baseline is `WobbleCal_F30_Cone30_B10_WallOn_Free_003`.
Only the analytic-gradient amplitude is changed for the new probe. All other
parameters remain unchanged (B=10 mT, spin=30 Hz, cone=30°, Bias=40°,
phase=248°, mass=10 mg, D055, same geometry/contact and trajectory).

## Gradient implementation audit

Production code (`J:\\magpy\\magpylib_socket_server.py`) uses
`analytic_gradient_b_t` as a localized axial-gradient field amplitude in tesla,
with Gaussian length scale `analytic_gradient_length_mm` (45 mm). It computes
`dB/ds = -(delta_s/L^2) exp[-delta_s^2/(2L^2)] G * 1000` in T/m and applies
the resulting force directly; the gradient is not a constant T/mm. Thus G=3 mT
means `0.003 T` amplitude, not `0.003 T/mm`.

## Replay result on the same trajectory

The replay uses the canonical centerline arclength increasing in +s as the
sole forward convention and the production robot moment/UR rotation mapping.
The sign gate is unambiguous:

| replay | mean Ft (µN) | median Ft (µN) | min..max (µN) | positive fraction | Ft impulse (N·s) |
|---|---:|---:|---:|---:|---:|
| G3, + sign | 3.533 | 4.159 | 0.000..4.500 | 100% | 1.040e-8 |
| G3, − sign | −3.533 | −4.159 | −4.500..0.000 | 0% | −1.040e-8 |
| G6, + sign | 7.066 | 8.319 | 0.000..9.000 | 100% | 2.080e-8 |
| G6, − sign | −7.066 | −8.319 | −9.000..0.000 | 0% | −2.080e-8 |

The measured production G3 telemetry on this run has mean Ft=2.873 µN,
median=2.743 µN, and is positive on every row. The replay therefore reproduces
the production sign and confirms that +G is the physical +s direction for this
trajectory. No `abs(Ft)` or CSV direction flip is used.

![gradient replay](gradient_sign_replay.png)

## Baseline dynamic reference

For the 3 ms Wall-ON free baseline, robot arclength changed from 13.814214 to
13.498125 mm (`Delta s=-0.316089 mm`) despite magnetic Ft remaining positive.
Peak COM speed was 511.5 mm/s, peak CPRESS 8.987 MPa, and exact minimum wall
gap −2.727 µm (node 55 / wall element 442 at 2.900 ms). This establishes that
the backward motion is produced by wall interaction/impulse, not by a negative
analytic gradient sign.

## Dynamic probe status

`WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003` was submitted with
`--analytic-gradient-b-t 0.006` (6 mT = 0.006 T), port 65496, and the same
3 ms deck. Datacheck/packaging completed without new node-face or edge-edge
overclosures. It completed successfully (31 frames, 3.000 ms).

| metric | G3 baseline | G6 (+6 mT) | change |
|---|---:|---:|---:|
| robot Δs (mm) | −0.316089 | −0.316767 | essentially unchanged |
| peak COM speed (mm/s) | 511.48 | 513.66 | +0.4% |
| final COM speed (mm/s) | 294.21 | 298.13 | +1.3% |
| UR1 peak-to-peak (rad) | 0.8501 | 0.8471 | unchanged |
| UR2 peak-to-peak (rad) | 0.4750 | 0.4960 | +4.4% |
| UR3 peak-to-peak (rad) | 0.7072 | 0.7033 | unchanged |
| CPRESS max (MPa) | 8.987 | 0.450 | −95.0% |
| exact minimum gap (µm) | −2.727 | −0.533 | 80.5% less penetration |
| measured Ft mean (µN) | 2.873 | 5.738 | ×2.00, positive all rows |

The fixed-camera GIF includes the corrected Abaqus-frame centerline/wall cloud,
robot trajectory, and per-frame red **B** and green **+s gradient** arrows:

`WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003_CEL_FSI_bend_validation.gif`

The stronger gradient substantially softens contact severity, but does **not**
recover forward transport: Δs remains −0.317 mm and peak speed remains about
514 mm/s. Magnetic Ft is positive throughout both runs, so the remaining
reverse displacement is a wall-contact impulse effect, not a sign error.

## Decision gate

The G6 gate is therefore **partially passed** (contact severity improved) but
the forward-motion gate is failed. Because reverse displacement did not improve
by ≥50%, a G9 probe is not justified under the stated rule. The next work should
target wall-impulse direction/geometry or a physically motivated forward-bias
redesign, not more scalar gradient magnitude. No long run or other parameter
change is authorized by this report.
