# Reduced-Hydrodynamics zeta=0.50 wall-on validation at 8.333 ms

## 1. Frozen model identity

The single dynamic job was `Wobble_F30_G6L45_ReducedHydro_Zeta050_WallOn_Free_0083`. It retained the 2.931414686 x 0.903074896 mm rigid robot, approximately 9.99997 mg mass, original inertia, RP, initial pose, and HEAD/TAIL definition. The magnetic drive remained 30 Hz, 10 mT, cone 30 deg, Bias 40 deg, phase 248 deg, sense +1, G=6 mT and L=45 mm through the production Magpylib Socket; Abaqus performed no electromagnetic solve. Reduced-Hydro remained body-following and anisotropic with Cparallel=4e-9, Cperp=1.2e-8 N s/mm, Kspin=1e-9 and Kwobble=3e-9 N mm s. Contact remained SmoothWall114, HARD penalty, mu=0.03, zeta=0.50 and tangent fraction 0.

Input identity was audited against the approved 1.3 ms source. The source SHA-256 is `0a4338fb29b9f85a2c991f576edb92e0b6243a3de1374f894720825cff302d82` and the 8.333 ms candidate SHA-256 is `ac84bb2f03ce0f0fe609305998a3ff2647231586ade43e58d86d2ddec13b0fa0`. Only duration, field cadence and an output comment changed.

## 2. Why zeta=0.50 is provisional

The normal damping fraction is a provisional engineering contact model, not an experimentally calibrated restitution law. Its short-window first-impact gate gave e_n=0.621897, 1.7 us duration, 5.490071 N peak force, -0.513129 um minimum gap, reopening and closed momentum/energy budgets. This validation therefore freezes zeta=0.50; it does not claim experimental restitution calibration.

## 3. Run completion and timestep

The job completed successfully to **8.333 ms** with **83,330** increments at fixed direct `dt=1.0e-07 s`. Abaqus reported 13 warnings, no initial overclosure, no Explicit failure and no automatic retry. The ODB contains 335 field frames at 25 us cadence and 83,331 RP/history samples. This is one quarter of a 30 Hz period, so it cannot establish periodic lock.

## 4. True-axis reconstruction

The robot axis was reconstructed at every increment by applying rigid U/UR to the reference exterior-node HEAD and TAIL groups. The reported phase and tilt do not use global UR component ranges. Maximum true tilt was **34.3251 deg** and the final tilt was **33.3909 deg**.

## 5. Canonical local frame

COM was projected continuously onto the authoritative `TRUE_centerline`; increasing canonical s is FORWARD. The local tangent and continuous e1/e2 frame define both robot and magnetic phases. Previous-s continuity and endpoint extrapolation were retained. The comparison cases were reprocessed with the same definitions.

## 6. First impact

First physical contact began at **1.0036 ms** and lasted **1.7000 us**. Peak wall resultant was **5.490071 N**, impulse magnitude was **5.337964e-06 N s**, and exact minimum gap was **-0.513129 um**. The separable normal velocities were -887.578 and 551.982 mm/s, giving effective e_n=0.621897, consistent with the accepted short probe.

## 7. Repeated contact events

There were **204** solver-active events. The longest was event 204, from 7.8991 to 8.3142 ms, lasting **415.2000 us** with a low 0.246285 N peak; it represents sustained multipoint contact rather than a clean separable impact. Event impulse, node, wall triangle, gap and pre/post normal velocity are catalogued in `zeta050_8p333_contact_events.csv`.

## 8. Exact gap

Rigid-body geometry was reconstructed for all near-wall intervals. Exact node-to-triangle distance used 1 us coarse coverage and 0.1 us backfill wherever gap <=50 um or solver contact was active. The full-run minimum exact gap was **-1.145295 um**. This gate uses impulse, gap and duration, not CPRESS peak.

## 9. Opposing-wall bridge

The first opposing-sector bridge occurred at **1.7712 ms**. A transient bridge ran from 1.7712 to 1.8637 ms; the persistent episode ran from **2.7727 ms through 8.3330 ms**, lasting **5.5604 ms**. Maximum near-wall node counts at 20/5/0 um were 9/6/3, with two separated wall-normal clusters at every threshold. Persistent bridge is therefore **present**.

## 10. True-axis phase progression

After first impact, the unwrapped true-axis local phase advanced only **47.5699 deg** by 8.333 ms, below the approximately 90 deg strong-survival reference. It advanced intermittently, reversed over some windows, and approached zero/negative short-window rate near the end rather than sustaining wobble progression.

## 11. Phase plateau

The longest qualifying 0.5 ms phase-support window was **0.5001 ms**, just above the formal threshold, while the local field phase continued. Thus a >0.5 ms phase plateau is present. Late 0.75 ms robot-rate windows ranged from -5.541 to 30.997 Hz, including reversal; these are short-window indicators only.

## 12. Tilt evolution

Tilt increased from approximately 12.7 deg to a high 30-34 deg posture. After the persistent bridge formed it remained near that slanted state, ending at 33.3909 deg. Early variations do not offset the combined late bridge, phase plateau and high-tilt gate.

## 13. Magnetic torque

Magnetic torque RMS was **0.00692476 N mm** before first impact and **0.00950673 N mm** after it. It did not collapse; post-impact RMS increased. Full-window |Tmag| RMS was 0.00923438 N mm. A persisting drive together with stalled robot phase supports a contact/geometric stall rather than loss of magnetic forcing.

## 14. Hydro torque

Full-window |Thydro| RMS was 8.60486408e-07 N mm, far below |Tmag|, and `Thydro dot omega` never became positive above tolerance: maximum **0.000e+00 W**, positive fraction **0.000**. The Reduced-Hydro rotational term is therefore dissipative throughout. No coefficient was refit.

## 15. Energy consistency

ODB assembly energies were converted from the active N-mm-s unit system (N mm) to joules before reporting. Maximum absolute Abaqus ETOTAL residual was **4.558492e-11 J**, artificial energy ALLAE remained zero, and the final independently reconstructed rigid KE differed from ODB ALLKE by only **4.137741e-13 J**. Magnetic work ended at 5.165303e-06 J, hydro work at -1.233517e-08 J, and inferred contact work at -4.324413e-06 J. There is no numerical energy-growth failure.

## 16. Canonical translation

Translation is diagnostic, not the success gate. Canonical displacement was **-2.763453 mm** and final tangent velocity was **-392.562550 mm/s**. Both are negative, so motion was backward under the authoritative increasing-s FORWARD definition; no legacy robot_arc sign was used.

## 17. CEL versus Reduced-Hydro

| Metric | CEL | RH old contact | RH zeta0.50 |
|---|---:|---:|---:|
| First contact (ms) | 1.5001 | 6.1500 | 1.0036 |
| Contact events | 24 | 1 | 204 |
| Longest contact (ms) | 1.2329 | 0.0500 | 0.4152 |
| Minimum gap (um) | -2.0594 | -0.6313 | -1.1453 |
| Maximum tilt (deg) | 34.2822 | 32.8343 | 34.3251 |
| Post-impact phase advance (deg) | 39.8441 | 9.9884 | 47.5699 |
| Longest plateau (ms) | 0.0000 | 0.0000 | 0.5001 |
| Bridge onset (ms) | 2.2500 | not detected | 1.7712 |
| Longest bridge (ms) | 2.8830 | 0.0000 | 5.5604 |
| Canonical delta s (mm) | -2.0255 | -2.8327 | -2.7635 |
| Final Vt (mm/s) | -304.3447 | -411.7654 | -392.5626 |

The older cases have only 50 us geometry snapshots, so their micro-impact event counts are lower bounds. Under the unified local-frame geometry definitions, zeta=0.50 does not show a clear wobble-survival improvement over CEL: both reach a high tilted state and form long opposing bridges.

## 18. Visual GIF interpretation

`Wobble_F30_ReducedHydro_Zeta050_WallOn_8p333.gif` uses a fixed camera and shows the full robot, translucent wall, HEAD/TAIL, canonical centerline, local frame, magnetic vector, true axis and live audit quantities. It visibly follows the sequence **turn, become slanted, remain geometrically constrained**. The synchronized CEL comparison shows no clear qualitative rescue of post-impact wobble by deleting CEL and using zeta=0.50 contact.

## 19. Decision

`GEOMETRY_CONTACT_STILL_JAMS_WITH_REDUCED_HYDRO`

The run is implementation-valid and energetically stable. Contact rebound is now provisionally reasonable and magnetic torque survives, yet the 2.931 x 0.903 mm robot develops a 5.5604 ms opposing-wall bridge, a >0.5 ms phase plateau and a persistent approximately 33 deg tilted posture. The current limiting mechanism is therefore robot-to-pipe geometry, not CEL, missing magnetic torque or the accepted normal damping. This result does **not** permit extension to 16.667 ms and does not establish 30 Hz lock.

## 20. Exactly one next step

**Begin modifying the robot geometry dimensions, prioritizing a shorter robot length.** Do not continue tuning contact damping, Reduced-Hydro coefficients or magnetic-field parameters before that geometry change.
