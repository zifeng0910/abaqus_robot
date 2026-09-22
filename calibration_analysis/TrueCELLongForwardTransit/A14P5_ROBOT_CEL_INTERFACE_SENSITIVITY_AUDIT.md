# A14P5 TAIL robot-interface discretization sensitivity (offline)

**Primary classification: ROBOT_CEL_NEGATIVE_LOAD_STRONGLY_INTERFACE_SENSITIVE**

**Independent-fluid-evidence classification: CONTROL_VOLUME_UNAVAILABLE**

These are classifications of the *existing numerical force ledger*, not proof
that the adverse load is physically real or a proven node-count artifact. No
Abaqus process, new dynamics case, mesh, physics change or ODB write was used.
The reproducible analysis is `scripts/audit_a14p5_interface_sensitivity.py`.

## Provenance and contact definitions

The original coarse `TRUECEL_A14P5_DIAG_2CYCLES` has no pair-isolated native
contact history. Its output-only `TRUECEL_A14P5_FORCE_FLUX_DIAG_2CYCLES`
successor adds contact history requests. Every value, including timestamps,
in the native RP U, UR, V and VR channels matches the original **exactly**
(maximum absolute difference zero); the input diff changes job identifiers and
output requests, not the physical cards. Therefore the successor is the coarse
*force-output proxy*, and is identified as such in every table.

The input's General Contact inclusions are robot-wall, robot-CEL and wall-CEL.
Only the first two act directly on the robot. Thus native robot whole-surface
`CFT` minus native robot/wall pair `CFT` is the **inferred robot-CEL resultant**
under these explicit inclusions; it is not a directly output pair-specific
fluid traction. Nodal `CNORMF+CSHEARF` and `CPRESS` are **whole-General-Contact**
fields and cannot be regionally split into wall and CEL. Frame-wise nodal axial
sum agrees with interpolated native whole-surface `CFT` to within
2.81e-7 of the peak native magnitude (maximum over the common window). This
supports the resultant extraction, not the physical origin of that load.

## Common-window ledger

All impulses below use 8.333333 ms as the start and source-native timestamps,
with identical endpoints for both runs. Units are 1e-8 N s. `S_J` uses the
coarse inferred fluid impulse as denominator.

| Endpoint ms | Coarse Jfluid | Refined Jfluid | S_J | Coarse Jmag / Jwall / delta p | Refined Jmag / Jwall / delta p |
|---:|---:|---:|---:|---|---|
| 9.000 | -3.387 | -4.023 | 0.188 | +1.769 / 0 / -1.642 | +1.767 / 0 / -2.315 |
| 9.500 | -9.289 | -1.486 | 0.840 | +3.084 / -1.112 / -7.143 | +3.086 / -0.782 / +0.937 |
| 10.000 | -8.821 | -3.090 | 0.650 | +4.383 / -1.560 / -5.666 | +4.387 / -1.112 / +0.502 |
| 10.300 | -5.809 | -10.196 | 0.755 | +5.163 / -1.560 / -1.888 | +5.165 / -1.121 / -5.845 |
| 10.682021 | -5.633 | -13.112 | 1.328 | +6.155 / -1.560 / -0.738 | +6.156 / -1.121 / -7.783 |

At the last common time, residual/normalized closure is +3.000e-9 N s / 2.247%
coarse and +2.933e-9 N s / 1.438% refined. The magnetic impulse differs by
about 0.017%, while the inferred fluid impulse differs by **133%** relative to
coarse. The refined run reaches zero axial speed at that time, whereas coarse
still has +13.161 mm/s. Displacement since cycle 2 began is +0.016748 vs
+0.023887 mm (refined minus coarse -0.007138 mm). The separate full-cycle-2
displacements from the prior audit are -0.17716 vs -0.26473 mm and should not
be substituted for these pre-reversal window displacements. The changed
inferred-fluid impulse (-7.479e-8 N s) has the sign and approximate size of
the changed momentum (-7.045e-8 N s); magnetic and wall differences and
closure account for the remainder.

At common magnetic phase, `A14P5_INTERFACE_COMMON_PHASE_STATES.csv` gives
position, speed, angle, angular speed and all three axial forces at six common
times. This is **not** a state match: already at 8.500 ms the angles are
-4.186/-5.062 degrees and angular speeds 222.68/208.61 rad/s, with speeds
13.327/8.427 mm/s (coarse/refined). Requiring phase within 25 us, angle within
0.15 degree and angular speed within 5 rad/s before 9.300 ms yields no pairs;
no state-controlled force contrast is asserted.

## Localization and time ordering

TAIL/MID/HEAD are the same body-coordinate partitions in both meshes:
`s < -0.6`, `-0.6 <= s <= +0.6`, `s > +0.6` mm. TAIL physical exterior area
is 2.2072717568 mm2 **in both**, while TAIL face count is **438 -> 1134**
and surface node count **466 -> 814**. Maximum actively loaded TAIL nodes in
the 8.333-9.000 ms window are 223 -> 460, using whole-GC nodal force norm
greater than 1e-8 N. TAIL whole-GC axial nodal impulse over the final common
window changes from -9.261e-7 to -4.524e-6 N s; MID and HEAD contain large
compensating positive contributions. These large *regional* values cancel in
the native whole resultant. At fixed area, the TAIL regional magnitude grows
4.89x versus 2.59x as many faces and 1.75x as many surface nodes. It does
not simply scale with count, and includes any wall contribution; the spatial
change is suspicious but is **not** a proof of a counting bug. The force map
and force/area plot preserve identical physical bins and scales.

The largest negative instantaneous inferred forces in the window are
-0.002043 N at 9.160995 ms coarse and -0.002306 N at 9.723990 ms refined:
13% larger in magnitude and 0.563 ms later, showing changed waveform as well
as changed impulse. An operational 50-us trailing average first below
-5e-5 N occurs at 8.452869 vs 8.397617 ms (refined 55.3 us earlier). This
onset has an intrinsic 50-us trailing delay and approximately 0.00009-ms
native-sample half-width; it is *not* a unique physical trigger. Raw native
velocity differentiation first remains below -1000 mm/s2 for 30 us at
8.437977 vs 9.360102 ms, but this criterion is especially sensitive to
oscillation/threshold and cannot establish a causal onset. The EVF-weighted
absolute pressure metric in a translating +/-1.8-mm axial slab first differs
from its cycle-start value by >max(50%, 1e-4 N/mm2) at 9.025056 vs
8.950178 ms. Field-frame half-width is about 0.0125 ms, plus threshold
uncertainty. The 9.5-ms metric is 5.49e-4 vs 2.31e-3 N/mm2. Pressure metric
is not pressure traction and these onset definitions are not commensurate.

## Independent fluid evidence and CPRESS limits

The fixed-coordinate, fixed-scale CEL figure shows pressure, axial velocity
and EVF near both robot silhouettes at nearest frames to 8.5, 9.0, 9.5,
10.0, 10.3 and 10.682 ms. Shared absolute pressure and speed color limits
are 0.02116 N/mm2 and 597.55 mm/s. The pressure field differs strongly, but
the extracted field is cell-centered water-phase pressure, not pressure mapped
onto the robot boundary. With unknown interface pressure and the sloped robot
normal `n_s`, neither a scalar TAIL-minus-HEAD pressure nor a fabricated
`integral(-p*n_s*dA)` establishes the axial reaction. Classification:
**DIRECT_PRESSURE_TRACTION_UNAVAILABLE** and
**DIRECT_VISCOUS_TRACTION_UNAVAILABLE**.

The existing CPRESS summary contains only whole-GC frame peak, positive-value
count and value count; no peak coordinates, pair identity, actual active-patch
area or element-face mapping. Its full-run peaks 47.686/75.920 N/mm2 occur
at **1.850/1.875 ms**, outside the pre-reversal interval. In this interval
the maxima are 9.545/19.069 N/mm2. It is impossible from these saved values
to locate the peak, tell wall from CEL, calculate high-pressure physical area
or integrate pressure to match the native resultant. These tests remain
**unavailable**, not failed or passed.

A local Eulerian control volume would need the fluid momentum inventory,
advective and stress fluxes on **all** enclosing faces, fluid/solid cut-cell
bookkeeping and moving-boundary terms if the CV follows the robot. The saved
NPZ has EVF, cell-centered pressure and nodal velocity but no native face
flux, local face viscous stress or validated moving-interface pressure map;
the saved `boundary_stress_private.npz` is restricted to distant axial
*domain ends* in the coarse case. A pressure-only cell-centered rectangle
would omit unknown shear and cut-cell transfer and cannot independently infer
the robot reaction. No defensible two-size local CV impulse or percentage
difference is available: **CONTROL_VOLUME_UNAVAILABLE**. Therefore the
negative contact-inferred load is **not independently physically validated**.

The 11.675-ms event remains a **GENERAL CONTACT REGION FLAG TRANSITION**;
zero direct wall-pair force and increasing positive TAIL gap exclude calling
it confirmed wall recontact. No conclusion here invokes an external magnet
or transverse magnetic trap. Next numerically controlled work, if authorized,
would isolate and resolve the robot-CEL interface before tuning magnetic
control, contact or physical geometry.

## Deliverables

The common-window CSV is `A14P5_ROBOT_CEL_INTERFACE_SENSITIVITY.csv`;
common-phase states are in `A14P5_INTERFACE_COMMON_PHASE_STATES.csv` and
machine-readable limitations/metrics in `A14P5_ROBOT_CEL_INTERFACE_SENSITIVITY.json`.
The five figures are `A14P5_COARSE_VS_REFINE_FLUID_FORCE_HISTORY.png`,
`A14P5_COARSE_VS_REFINE_FLUID_IMPULSE.png`,
`A14P5_ROBOT_INTERFACE_AXIAL_FORCE_MAP.png`,
`A14P5_INTERFACE_FORCE_PER_AREA.png` and
`A14P5_COARSE_VS_REFINE_CEL_FIELDS.png`. A control-volume figure was not
created because no defensible estimate exists.
