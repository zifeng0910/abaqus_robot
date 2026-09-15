# Straight-pipe prescribed-wobble contact torque audit

## Decision

**Primary classification: `PRESCRIBED_MOTION_MASKED_GEOMETRIC_BRIDGE` (CASE C).**

The freely driven straight baseline is not primarily stalled by normal-contact damping capture. It reaches a persistent opposing-wall support state in which the inferred contact torque almost exactly balances the available magnetic torque. The controlled prescribed replay does not produce clean hit-separate wobble: it remains in opposing HEAD/TAIL support for the complete 2.3-16.667 ms analysis window and forces the robot through geometrically incompatible wall positions.

Do not run the proposed `zeta=0.20` magnetic case. Do not start a geometry sweep. The next controlled change should first reduce the commanded tilt trajectory into the current geometry's feasible envelope.

## Historical prescribed-motion search

No pre-Magpylib prescribed/programmed-motion implementation was found in the current repository or its reachable Git history. Searches covered Abaqus `*BOUNDARY` and `*AMPLITUDE` definitions, prescribed rotational DOFs, angular velocity, VDISP/VUAMP, programmed/kinematic wobble terms, and GIF-generation paths. Existing VUAMP sources are load interfaces, not a recovered prescribed-rotation case. Consequently, the remembered historical animation cannot be assigned a commit, geometry, clearance, friction, damping, or motion law and is not quantitatively comparable.

## Zero-new-solve baseline audit

The completed `PROD_LOCAL30_G0_STRAIGHT_CTRL` ODB was reused without rerunning it. Its history output contains per-increment RP pose and velocity plus whole-surface `CFN/CFS/CFT` resultants, but no node-level contact force and contact-point position. Therefore `sum((r_i-COM) x F_i)` is not directly recoverable, especially during opposing contacts whose net force can cancel while retaining a contact couple.

Contact torque was instead inferred from rigid-body angular momentum balance using measured RP angular velocity, the mass-scaled L2400 inertia, recorded magnetic torque, and the unchanged Reduced-Hydro law. A 10.1 us Savitzky-Golay derivative suppresses increment-scale differentiation chatter; sub-window impact peaks must be interpreted by impulse rather than instantaneous amplitude.

From 2.3 to 33.333 ms:

- both-end 20 um bridge fraction: **100%** in the dense reconstructed interval;
- magnetic torque RMS during the bridge: **0.00289493 N mm**;
- inferred contact torque RMS: **0.00289492 N mm**;
- median transverse inferred-contact/magnetic ratio: **0.999999**;
- 95th percentile ratio: **1.000041**;
- contact events over the complete cycle: **353**.

This is a near-static torque balance, not a loss of magnetic-load transmission. Node-level normal/tangential relative velocity and event restitution are not identifiable during simultaneous opposing contact from the archived whole-surface resultants.

## Contact damping audit

Existing same-contact first-impact evidence at `zeta=0.50`, tangent fraction 0, gives effective normal restitution **0.621897**, a reopened gap up to **68.016 um**, and no contact interval approaching 0.10 ms. That evidence does not meet `NORMAL_CONTACT_DAMPING_CAPTURE`: rebound is not quenched and the gap does reopen. Friction and contact damping can still influence later motion, but there is no evidence basis for authorizing a `zeta=0.20` long magnetic run.

## Single prescribed replay

`STRAIGHT_PRESCRIBED_WOBBLE_CONTACT_AUDIT` retained the current L=2.40 mm, D=0.815 mm robot, straight tube, mesh, wall, General Contact, `mu=0.03`, `zeta=0.50`, tangent fraction 0, Reduced-Hydro coefficients, and direct `dt=1e-7 s`. Translation remained free. Magpylib loads were removed and RP rotations were prescribed using a 30 Hz, nominal 30 deg conical trajectory with a 1 ms quintic startup ramp. The run completed 166,670 increments successfully in 326 s solver wall time.

After startup (2.3-16.667 ms):

- actual finite-rotation tilt range: **24.441-29.858 deg**, mean **26.552 deg**;
- both-end 20 um support fraction: **100%**;
- contact-active fraction: **100%**;
- minimum reconstructed geometric gap: **-246.511 um**;
- transverse reaction moment RMS: **11,506.84 N mm**;
- baseline available transverse magnetic torque RMS: **0.00295388 N mm**.

The apparent RMS ratio is about **3.90 million**, and the reaction exceeds the magnetic torque throughout the bridge window. However, this number is not a physical torque requirement: Abaqus issued a deep-penetration warning, and the prescribed path activated penalty-scale reaction. It is only diagnostic evidence that the imposed trajectory is geometrically incompatible. CASE A is not assigned because its prerequisite, clean hit-separate prescribed wobble, was not met.

## Interpretation

The old visual impression of wobble is insufficient evidence of valid impact-separation. A rotational boundary condition can rotate the RP while the robot remains supported by opposing walls and while constraint reactions become arbitrarily large. In the current geometry, the free magnetic case settles into a contact couple that consumes the available torque; forcing a much larger cone does not release it, but drives both ends through the wall constraint.

The evidence therefore ranks causes as:

1. **Geometric opposing support under the requested trajectory: decisive.**
2. **Finite magnetic torque versus the established contact couple: limiting in the free case.**
3. **Normal damping capture: not supported by current restitution/reopening evidence.**
4. **Historical prescribed-case comparability: unresolved because the original case was not found.**

## Output limitations

- Baseline contact torque is an angular-momentum-balance inversion, not direct node-level `r x F` summation.
- Normal and frictional contact torque cannot be separated from available ODB output.
- GIF normal velocity is the time derivative of reconstructed minimum gap, explicitly labeled as a proxy.
- Whole-surface contact resultants cannot reconstruct simultaneous opposing HEAD/TAIL forces.
- The replay reaction ratio is penalty-dominated and must not be used to size a magnetic actuator.

## Exactly one next step

Perform a zero-solve geometric feasibility calculation for the current exact robot surface and straight tube to define the maximum nonpenetrating local tilt as a function of azimuth and COM radial offset. Use that envelope to formulate a lower-amplitude prescribed diagnostic before changing L/D. This recommendation is analysis only; no additional dynamic case is authorized here.
