# Magnetic drive frame and wobble root-cause audit

## Executive conclusion

The production magnetic field is not a 30 deg cone about the robot-local tube axis. It is a 30 deg cone about an axis built at the prescribed driver arc, first tilted 40 deg from the production DXF tangent toward projected global +Z. That fitted cone axis is 119.275 deg from the polarity-correct robot-local tube axis. Over a complete 30 Hz command cycle, the production field has only 0.000005 local turns about the tube.

The controlled full-cycle Abaqus comparison isolates this error. `FRAME_CURRENT` produces only -0.0538 robot-axis turns and a short, wall-constrained arc. With the same L2.4 carrier, contact, hydrodynamics, B0, frequency and G=0, `FRAME_LOCAL30` and `FRAME_LOCAL40` produce 1.4360 and 1.4199 robot-axis turns, respectively, with repeated wall-sector changes and true separations. The primary classification is therefore:

**`MAGNETIC_DRIVE_FRAME_IS_PRIMARY_ROOT_CAUSE`**

The static implementation classification is **`MAGNETIC_DRIVE_FRAME_MISALIGNED`**. Abaqus UR mapping is valid and the body-fixed magnetic moment is correctly aligned with HEAD-to-TAIL, so neither is the common-mode cause. A large first impact in `FRAME_LOCAL30` also activates the separately gated 5 ms ramp diagnostic; its final result is reported below.

This is a one-cycle mechanism diagnostic, not evidence of a settled periodic limit cycle. The previous 6.667 ms screen covered only 0.20 command cycles (72 deg) and was intrinsically incapable of validating steady wobble or cycle closure.

## Scope and controlled cases

All dynamic cases use L=2.40 mm, D=0.815 mm, a rigid 1,588-node/7,302-element coarse robot, nominal 0.10 mm mesh, B0=10 mT, f=30 Hz, G=0, mu=0.03, zeta=0.50, baseline Reduced-Hydrodynamics, direct dt=0.1 us, 100 us field output, dense RP/contact history, and duration 33.333333 ms. Only the magnetic-frame definition differs between the three mandatory cases.

| Case | Magnetic frame | Cone | Startup ramp |
|---|---|---:|---:|
| `FRAME_CURRENT` | Unmodified production architecture | 30 deg, Bias=40 deg | 1 ms |
| `FRAME_LOCAL30` | Screening-only COM-local transported tube frame | 30 deg | 1 ms |
| `FRAME_LOCAL40` | Screening-only COM-local transported tube frame | 40 deg | 1 ms |
| `FRAME_LOCAL30_RAMP` | Gated repeat of `FRAME_LOCAL30` | 30 deg | 5 ms |

No broad parameter sweep was run. The contact-assisted optional case was not activated because both local-frame cases already show repeated true separations.

## Exact production forcing path

The traced path is:

`VUAMP/VUFORC -> SOCKET_POSE -> socket server -> analytic magnetic model -> mapped global force/torque -> Abaqus RP`

For the constant-frequency production configuration,

```text
phi(t) = phase0 + sense * 2*pi*f*t
```

At the selected field arc, the server constructs the production DXF tangent `t_dxf` and the projection of global +Z normal to it, `n_z`. With nonzero Bias beta,

```text
k  = normalize(cos(beta)*t_dxf + sin(beta)*n_z)
e1 = normalize(n_z - (n_z dot k)*k)
e2 = normalize(k cross e1)
B  = B0*[cos(alpha)*k + sin(alpha)*(cos(phi)*e1 + sin(phi)*e2)]
```

Here `alpha=cone=30 deg`, `beta=Bias=40 deg`, `phase0=248 deg`, and `sense=+1`. Thus `cone` is the half-angle about the already biased axis `k`; it is not the angle to the tube axis.

The global magnetic moment and loads are

```text
m_global = R(UR) * [polarity * |m| * a_body]
T_Nmm    = 1000 * (m_global cross B)
F_N      = m_global dot grad(B)
```

The smooth startup scale multiplies both force and torque. At t=0 it is zero, rises with `x^2*(3-2*x)`, and reaches one at the requested ramp time.

### Parameter semantics

| Parameter | Implemented meaning and frame |
|---|---|
| `f` | Command phase rate in Hz; global scalar, integrated into `phi(t)` |
| `B0` | Magnitude of the analytic rotating field; 10 mT in these diagnostics |
| `cone` | Half-angle `alpha` about production axis `k`, not about the robot-local tube axis |
| `Bias` | Angle `beta` from production DXF tangent toward projected global +Z |
| `phase` | Initial azimuth in the `e1/e2` basis; 248 deg |
| `sense` | Sign of phase accumulation; +1 here |
| `G` | Independent analytic spatial derivative amplitude used for force, not added to B; zero in frame diagnostics |
| `Lgradient` | Gaussian derivative length scale in arc length; production default 45 mm; inactive when G=0 |
| robot pose | Rotates body moment through Rodrigues `R(UR)`; should select the local field arc when follow-robot is active |

There is an additional production logic defect: `_driver_pose()` populates `robot_arc_mm` only when `adaptive_lead_mm>0`. The present setting is zero, so `--analytic-follow-robot` cannot select robot COM arc and falls back to driver arc. At t=0, robot COM arc is 13.7999606 mm while driver arc is 18.899960626 mm.

The production DXF tangent polarity is also opposite the requested continuously oriented tube direction `c_hat`, whose sign is chosen by `c_hat dot a0 > 0`.

## Zero-dynamics field audit

At the authoritative initial COM, the polarity-correct transported frame is

```text
c_hat = [ 0.99012918,  0.07930163,  0.11556579]
e1    = [-0.11868692,  0.03582724,  0.99228515]
e2    = [ 0.07454942, -0.99620664,  0.04488566]
c_hat dot a0 = 0.97292785
```

The one-cycle production B fit gives:

| Quantity | Result |
|---|---:|
| Fitted cone axis, global Abaqus | `[-0.573220, -0.175682, 0.800346]` |
| Fitted axis in `(c,e1,e2)` | `[-0.489001, 0.855911, 0.168207]` |
| Fitted half-angle | 30.00004 deg |
| Cone-fit RMS angular residual | 0.01947 deg |
| Fitted axis vs robot-local `c_hat` | 119.27497 deg |
| `theta_B_t` min/mean/max | 89.2927 / 116.9937 / 149.2927 deg |
| Local B azimuth winding | 0.000005 turns |
| Effective local azimuth rate | 0.00015 Hz |
| `|B|` | 10 mT to numerical precision |
| `B_parallel` range | -8.598 to +0.123 mT |
| `B_perp` range | 5.107 to 10.000 mT |

The field-only animation is [`CurrentProductionField_LocalTubeFrame_1cycle.gif`](CurrentProductionField_LocalTubeFrame_1cycle.gif). The fitted orbit and component diagnostics are in [`figures/field_tip_local_orbits.png`](figures/field_tip_local_orbits.png) and [`figures/field_local_frame_diagnostics.png`](figures/field_local_frame_diagnostics.png).

Along the already-completed L2.3 trajectory, with no Abaqus rerun, production B remains 80.567-98.273 deg from the local tube axis. Across 0.20 command cycles its local winding is -0.09413 turns and mean local phase rate is -14.12 Hz, opposite in sign and different in magnitude from the +30 Hz command.

## Pose and magnetic-moment audit

Abaqus Explicit RP `UR1..UR3` are total axis-angle rotation-vector components in radians for this rigid body. The server applies the Rodrigues exponential map, not XYZ Euler angles.

Against 1,588 rigid nodes over 68 actual ODB frames and rotations up to 121.56 deg, the maximum server-to-Kabsch orientation error is 4.18e-6 deg and Kabsch node RMS is 4.77e-7 mm. Synthetic 0, 15, 30, 45, 60 and 90 deg quaternion comparisons have zero numerical error. Classification: **`POSE_MAPPING_VALID`**.

The normalized body moment is `[0.965033, -0.110154, 0.237862]` in the Abaqus initial body representation and has dot product 0.999957 with HEAD-to-TAIL `a0`. The implemented moment direction therefore matches the intended long-axis magnetization.

At t=0, the full-amplitude field would be 86.416 deg from the moment. Its raw torque is `[0.000389, -0.009155, -0.005819] N mm`, or `[-0.001013, -0.006149, +0.008888] N mm` in `(tube-axis,e1,e2)` components. The actual t=0 torque is zero because of the 1 ms smooth startup ramp. The large radial components still produce a strong acquisition transient as the ramp reaches full amplitude.

## Full-cycle dynamic results

| Metric | `FRAME_CURRENT` | `FRAME_LOCAL30` | `FRAME_LOCAL40` | `FRAME_LOCAL30_RAMP` |
|---|---:|---:|---:|---:|
| Command winding | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| Local B winding | 0.00061 | 1.00852 | 1.01177 | 1.00738 |
| Robot-axis winding | -0.05377 | 1.43601 | 1.41991 | 1.42276 |
| Circular mean robot-minus-B phase | +5.10 deg | -2.20 deg | -1.99 deg | -3.57 deg |
| Maximum directed tilt | 38.62 deg | 37.70 deg | 40.17 deg | 35.57 deg |
| Final directed tilt | 36.50 deg | 25.80 deg | 34.14 deg | 24.49 deg |
| Contact events | 154 | 32 | 139 | 33 |
| True separations >=20 us | 52 | 32 | 110 | 30 |
| Contact-side transitions >30 deg | 99 | 19 | 65 | 17 |
| Both-end support fraction | 0.9400 | 0.0665 | 0.4970 | 0.0480 |
| Longest continuous contact | 4081.3 us | 2.2 us | 834.7 us | 725.5 us |
| Longest opposing bridge | 29.920 ms | 0.500 ms | 7.600 ms | 0.520 ms |
| First contact time | 0.916 ms | 1.672 ms | 1.489 ms | 3.188 ms |
| Peak contact force | 7.212 N | 4.436 N | 6.076 N | 2.850 N |
| Peak contact time | 1.291 ms | 1.672 ms | 1.490 ms | 3.188 ms |
| Arc displacement `Delta s` | -3.629 mm | +1.744 mm | +4.791 mm | +1.234 mm |

`FRAME_CURRENT` travels only a short arc/cluster in `(a dot e1, a dot e2)`, then remains wall constrained. It reaches the lower centerline projection bound; this is real reverse motion, not endpoint taper. Telemetry confirms B=10 mT, drive scale=1 and nonzero magnetic torque after the 1 ms ramp.

Both local-frame cases trace a complete large orbit, show repeated contact-sector rotation and true separations, and visually exhibit the requested rotate-hit-separate-continue/change-sector sequence. Their 1.42-1.44 turns per one field cycle include acquisition overshoot; the trajectories do not close exactly and must not be called steady periodic solutions.

The synchronized primary artifact is [`Current_vs_Local30_vs_Local40_OneCycle.gif`](Current_vs_Local30_vs_Local40_OneCycle.gif). Individual dual views and orbit plots are stored in each case directory. The shared winding comparison is [`figures/dynamic_winding_comparison.png`](figures/dynamic_winding_comparison.png).

## Gated startup-ramp diagnostic

`FRAME_LOCAL30` first contacts the wall at 1.671897 ms with a 4.43579 N peak at 1.672297 ms, after rapidly moving from the initial 13.36 deg tilt into the 30-40 deg orbit band. Because the subsequent trajectory has a valid full orbit, this satisfies the user-defined gate for exactly one `FRAME_LOCAL30_RAMP` case with a 5 ms smooth ramp.

The gated case completed successfully through the full 33.333333 ms cycle. Across 335 status samples, the reported total-energy balance magnitude remained at or below `1.573e-9`, only `2.53e-6` of the `6.228e-4` peak kinetic energy, supporting numerical stability under the prescribed direct time increment. Relative to the otherwise identical 1 ms-ramp case, the 5 ms ramp delayed first contact from 1.671897 to 3.187694 ms and reduced the peak contact force from 4.43579 to 2.84982 N, a 35.8% reduction. It retained the intended mechanism: robot-axis winding is 1.42276 turns, local-B winding is 1.00738 turns, the circular mean phase lag is -3.57 deg, and there are 17 contact-sector transitions with 30 true separations.

The ramp therefore softens and delays the initial wall slam without changing the primary frame conclusion. It does not eliminate acquisition contact: the longest continuous contact increases from 2.2 us to 725.5 us, while both-end support remains low (0.0480) and the longest opposing bridge remains short (0.520 ms). Its one-cycle arc displacement is +1.234 mm versus +1.744 mm for the 1 ms ramp. The result supports treating startup shaping as a secondary transient-control parameter after the production magnetic frame is corrected, not as an alternative root-cause explanation.

The ramp-case dual view is [`cases/FRAME_LOCAL30_RAMP/FRAME_LOCAL30_RAMP_DualView_33p333.gif`](cases/FRAME_LOCAL30_RAMP/FRAME_LOCAL30_RAMP_DualView_33p333.gif), and its axis orbit is reported in the same case directory.

## Answers to the 14 required questions

1. **What axis does production B rotate around?** About the biased driver-arc axis `k=cos(40 deg)t_dxf+sin(40 deg)n_z`, fitted globally as `[-0.573220,-0.175682,0.800346]`, not about robot-local `c_hat`.
2. **What does `cone=30 deg` mean?** It is the half-angle around that already biased axis.
3. **Why was tube-relative field tilt about 82-97 deg in the prior trajectory?** The field basis is built at driver arc, `analytic-follow-robot` is disabled by the `adaptive_lead=0` logic path, tangent polarity is opposite `c_hat`, and Bias tilts the axis another 40 deg.
4. **Does local B wind 360 deg per 30 Hz cycle?** Production: no, approximately 0.000005 turns in the fixed-COM audit. Local diagnostic: yes, about 1.01 turns along the moving trajectories.
5. **Is Abaqus UR converted correctly?** Yes; maximum actual ODB orientation error is 4.18e-6 deg.
6. **Is the body-fixed magnetic moment direction correct?** Yes; dot product with HEAD-to-TAIL is 0.999957.
7. **How large is the t=0 torque mismatch?** Raw full-amplitude angle is 86.416 deg and torque components are `[-0.001013,-0.006149,+0.008888] N mm` locally. Applied torque at exactly t=0 is zero due to the ramp.
8. **Does `FRAME_CURRENT` make a closed axis orbit?** No; it is a short wall-constrained arc with -0.0538 net turns.
9. **Does `FRAME_LOCAL30`?** It completes a clear orbit with 1.4360 turns, but does not close exactly after acquisition.
10. **Does `FRAME_LOCAL40`?** It completes a clear orbit with 1.4199 turns, but also does not close exactly.
11. **Which cases repeatedly rotate contact sector?** `FRAME_LOCAL30`, `FRAME_LOCAL40`, and the gated ramp repeat do. Current has many impact pulses while translating along/against the wall, but its robot axis does not wind; event count alone is not wobble evidence.
12. **Does any case visually show rotational wall-impact wobble?** Yes, all three local-tangent cases show it over the diagnostic cycle. This is not yet a multi-cycle steady-state claim.
13. **Was the old 0.20-cycle screen intrinsically too short?** Yes. It was useful for early transient boundaries but cannot establish a full revolution, wall-sector exchange, cycle closure or steady wobble.
14. **Primary cause?** Magnetic drive-frame construction is primary. Pose mapping and magnetization are valid. The gated ramp confirms startup is a secondary acquisition issue: it reduces peak force by 35.8% while retaining the local-tangent orbit, but does not remove wall contact.

## Decision boundary and next experiment

The evidence justifies replacing the production field-frame construction before any renewed length/friction/zeta/hydrodynamic sweep. It does not justify silently replacing production with the screening wrapper, selecting 30 vs 40 deg as a winner, or claiming a stable limit cycle from one period.

After the ramp result is incorporated, the next physics validation should be a limited multi-cycle run of an explicitly reviewed local-tangent field implementation, with propulsion gradient still separated from the rotational drive. The acceptance gate should require repeatable cycle-to-cycle axis winding, phase lag, wall-sector sequence, and bounded energy, not displacement alone.
