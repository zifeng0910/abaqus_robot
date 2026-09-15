# Local Head-Tail Rocking 5 Hz Validation

## Decision

**LOCAL_ROCKING_MODE_RECOVERED**

1. RouteA is `theta_rock(t)=10 deg*r(t)*sin(2*pi*f*t)`, initially at 0 deg, about the transported binormal. Its legacy 0.03 s display uses 666.667 Hz and a 2 ms smooth ramp; the physical reference is 5 Hz.
2. The source establishes `(c_hat,n_rock,b_rock)=(T,N,T x N)`. In the current production PT gauge, `n_rock=cos(chi)e1+sin(chi)e2` and `b_rock=c_hat x n_rock`, with `chi=-61.372847575963 deg`; it is not the unrotated `(e1,e2)` pair.
3. The unchanged straight geometry is exactly feasible: solver ID `1.334690480 mm`, analytic 10 deg margin `0.115316535 mm`, and static exact minimum gap `69.198 um`.
4. The primary target changed because RouteA is alternating planar HEAD-TAIL rocking, while the former 30 deg cone has a 360 deg cross-section orbit and is a different motion family.
5. `ROBOT_LOCAL_ROCKING` uses `B=B0[cos(alpha_B)c_hat+sin(alpha_B)n_rock]`, `alpha_B=10 deg*sin(2*pi*5t)`, with continuous COM projection and the existing PT frame.
6. Dynamic `|B|` remains 10 mT; maximum sampled magnitude error is `1.516e-11 T`.
7. The field remains in one local plane and has zero commanded 360 deg winding, as established by the offline one-cycle regression.
8. No robot UR is prescribed. The response comes from the validated finite moment, `m x B`, rigid-body dynamics, Reduced-Hydro loads, and one-wall General Contact.
9. `theta_rock` spans `-10.250..10.251 deg` (5 Hz fitted amplitude `10.012 deg`).
10. `theta_cross` RMS/max is `3.614e-05/8.644e-05 deg`; PCA line-likeness is `1.000000000` (descriptive, not a hard gate).
11. HEAD-only/TAIL-only support fractions are `0.0000/0.0000`; support-indicator sign changes: `0`.
12. Both-end support fraction is `0.0000`, longest bridge `0.000 ms`; it is not the majority state.
13. Minimum exact surface gap is `64.877 um`; deep penetration is `False`. Peak kinetic energy is `3.459e-07 N mm`, and maximum absolute `ETOTAL` residual is `6.825e-16 N mm`.
14. `RouteA_vs_MagneticRocking.gif` provides the normalized-phase visual gate. The scalar classification does not override direct visual review.
15. Exactly one next step: **Transfer this identical ROBOT_LOCAL_ROCKING input, contact model, and Reduced-Hydro model to the curved tube for one turn-section validation.**

The legacy RouteA GIF is used only for topology/amplitude/phase comparison because it contains CEL deep-penetration diagnostics and prescribed RP motion. The HighEndStop reference remains secondary and only frames 0-27 are retained as its trusted early window.
