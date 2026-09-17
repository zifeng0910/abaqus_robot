# FAST precomputed F120 dual-end no-rebound screen

Final classification: `FAST_MOTION_TOPOLOGY_FAILED`

This is one non-CEL `FAST_SURROGATE_FLUID` dynamics run. Magnetic loads came from the validated local `PRECOMPUTED_TABLE`; dynamics made zero Python, TCP, socket, or Magpylib calls. General Contact remained the robot-wall architecture. Its single normal-law change was HARD to progressive SCALE FACTOR (`r=5%`, geometric scale 10, initial scale 0.01), with mu=0.03 and zeta=1.0 unchanged.

- delta_s: -0.08340030 mm
- rocking: -12.8678 to +13.3278 deg
- HEAD episodes: 0; max R_v: None
- TAIL episodes: 13; max R_v: 59.04272435528249
- longest both-end bridge: 0.000000 ms
- magnetic lookups: 83334; socket calls: 0
- failed gates: TAIL_single_episode_per_half_cycle, HEAD_R_v_below_0p2, TAIL_R_v_below_0p2, overall_plus_s_progression
- primary failure reason: CONTACT_STIFFNESS_TOO_HIGH
- exactly one next FAST_SCREEN modification: Change only the General Contact geometric stiffness scale from 10.0 to 2.0; keep initial scale 0.01 and all other settings fixed.

Performance: 1455.775 s Abaqus wallclock for 16.666668 ms physical time on one CPU (87.3465 s/ms), with 83,334 magnetic table lookups and zero socket calls. Dense 10 us contact-field ODB output dominated I/O. See `FAST_PRECOMP_F120_DUALEND_NOREBOUND_Performance.md` for comparison context and the increment-CSV diagnostic limitation.
