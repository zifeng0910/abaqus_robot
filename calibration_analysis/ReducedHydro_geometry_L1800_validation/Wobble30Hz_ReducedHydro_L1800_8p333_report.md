# Reduced-Hydrodynamics L=1.800 mm geometry validation

## 1. Decision

`L1800_PARTIAL_WOBBLE_SURVIVAL_NO_PERSISTENT_BRIDGE`

Shortening the already validated zeta=0.50 robot removes the persistent opposing-wall bridge, but it does not yet produce sustained post-impact wobble. The robot advances in local true-axis phase, then slows and reverses late while reaching a near-transverse posture.

## 2. Corrected baseline geometry identity

The validated zeta=0.50 INP already contains a 0.90 transverse PCA shrink. Direct mesh audit gives a baseline PCA envelope of **2.921463527 x 0.814538836 mm**, not the legacy pre-shrink 2.931414686 x 0.903074896 mm description. This candidate changes only the PCA axial span to **1.800000000 mm** and leaves the transverse coordinates unchanged. Job labels and Socket metadata now use D0815 to avoid perpetuating the old identity error.

## 3. Frozen physics and scaling

Contact remains SmoothWall114, HARD penalty, mu=0.03, zeta=0.50 and zero tangent damping. The 30 Hz, 10 mT, cone 30 deg, Bias 40 deg, phase 248 deg, sense +1, G=6 mT and L=45 mm drive is unchanged. Reduced-Hydro remains Cparallel=4e-9, Cperp=1.2e-8 N s/mm, Kspin=1e-9 and Kwobble=3e-9 N mm s. Constant density and magnetization give mass **6.161278579 mg** and explicit moment **0.000719639311 A m2**, both scaled by 0.616129547. COM and RP are unchanged.

## 4. Eligibility and run completion

Old jam-pose replay found no opposing bridge at 20, 5 or 0 um and a minimum replayed gap of +135.885 um. Datacheck had zero errors and no initial overclosure. The one permitted dynamic job completed to **8.333 ms** in **83,330** direct increments at fixed dt=1.0e-07 s. There was no retry or timestep change.

## 5. Contact and exact gap

First physical contact began at **0.9487 ms**, lasted **1.4 us**, peaked at **5.443692 N**, and had effective separable restitution 0.567105. There were **18** events; the longest lasted **2.6 us**. The minimum exact node-to-triangle gap was **-1.357477 um**. This is 0.357 um beyond the previous preferred -1 um penetration target, so contact remains provisional even though the run is numerically valid.

## 6. Opposing-wall bridge

The first opposing bridge occurred at **4.9986 ms** and the longest lasted **0.2588 ms**. No episode exceeded 0.5 ms and none persisted to the end. The baseline longest bridge was 5.5604 ms. Axial shortening therefore removed the previously dominant geometric jam mechanism.

## 7. True-axis phase and tilt

After first contact, local true-axis phase advanced **38.9525 deg**; total advance peaked at **69.5904 deg** and ended at 64.0791 deg. Late 0.75 ms window rates ranged from -20.460 to 7.819 Hz and became negative. Maximum tilt was **89.9999 deg**, final tilt **81.6053 deg**, and the 7-8.333 ms mean was 78.9710 deg. This is dynamic tumbling/return, not the old 33 deg bridged lock, but it is not sustained wobble.

## 8. Plateau boundary

The longest qualifying support interval is **0.5001 ms**. Only 4 center samples qualify, in two isolated pairs near 4.0649 and 8.0177 ms. The result is only 0.1 us above the 0.5 ms definition and is reported as a marginal threshold crossing, not a long persistent plateau.

## 9. Torque and energy

Pre/post-impact magnetic torque RMS is **0.00387303/0.00224082 N mm**; it does not collapse. Full-window magnetic/hydrodynamic torque RMS is 0.00248234/2.10408249e-06 N mm. `Thydro dot omega` is never positive above tolerance. Maximum |ETOTAL| is 2.284610e-10 J (0.00399% of peak |ALLWK|), final rigid/ODB kinetic-energy mismatch is 3.286e-14 J, and ALLAE remains zero. No numerical energy growth is detected.

## 10. Translation

Canonical displacement is **-0.484862 mm** and final Vt is **-60.868744 mm/s** under increasing-s FORWARD. Translation remains diagnostic, not the wobble gate.

## 11. Baseline comparison

| Metric | Frozen baseline | L=1.800 mm |
|---|---:|---:|
| Contact events | 204 | 18 |
| Longest contact (us) | 415.2 | 2.6 |
| Longest bridge (ms) | 5.5604 | 0.2588 |
| Persistent bridge | True | False |
| Post-impact phase advance (deg) | 47.5699 | 38.9525 |
| Maximum/final tilt (deg) | 34.33/33.39 | 90.00/81.61 |
| Canonical delta s (mm) | -2.7635 | -0.4849 |

The geometry change reduces contact events by 91.2% and removes the long bridge, but does not increase net post-impact phase advance over this quarter-cycle.

## 12. Exactly one next step

**Audit the existing L=1.800 mm trajectory's rotational response and Reduced-Hydro rotational coefficients without another Abaqus run.** The next question is why the free, non-bridged robot approaches 90 deg tilt and reverses phase late. Do not change magnetic forcing, contact damping, or geometry again until that audit separates expected tumble dynamics from a rotational-model limitation. As before, 8.333 ms is only one quarter of a 30 Hz cycle and cannot establish frequency lock.
