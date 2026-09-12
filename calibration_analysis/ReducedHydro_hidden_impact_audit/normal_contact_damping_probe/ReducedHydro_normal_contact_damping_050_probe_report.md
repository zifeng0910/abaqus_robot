# Reduced-Hydro normal-contact damping 0.50 probe

**Status:** `PROVISIONAL_NORMAL_DAMPING_PROBE`

**Decision:** `PROVISIONAL_NORMAL_CONTACT_DAMPING_ACCEPTABLE_FOR_WOBBLE_TEST`

**Long-run status:** The short-window gate permits a later 8.333 ms validation, but no long run was launched in this work.

## 1. Why 0.20 failed

The ζ=0.20, tangent=0 case was numerically stable and dissipative, but its contact-point normal restitution was 0.825871, above the 0.70 upper warning limit. It therefore remained too elastic for the provisional screening objective.

## 2. Why 0.50 was selected

ζ=0.50 was chosen as one higher-information probe. A local two-point empirical extrapolation suggested restitution near 0.58, inside the provisional 0.30-0.65 design window. This was candidate selection only: it is not an experimental calibration, theoretical restitution law, or Abaqus-defined mapping.

## 3. Why this is now a one-variable comparison

The formal reference is ζ=0.20, tangent fraction=0 and the candidate is ζ=0.50, tangent fraction=0. Only the normal critical damping fraction changed. Magnetic forcing, Reduced-Hydro coefficients, robot geometry, mass, inertia and pose, wall, friction coefficient 0.03, General Contact, hard pressure-overclosure, penalty enforcement, surfaces, output cadence, duration, and direct timestep were frozen.

## 4. Datacheck and dt stability

The independent Abaqus 2025 datacheck passed with 0 errors, the same 13 warnings as the reference, and no initial overclosure. The sole 1.3 ms dynamic run completed successfully with fixed `*Dynamic, Explicit, DIRECT` `dt=1e-7 s`. No retry, timestep change, damping sweep, or 8.333 ms run was performed. Energy remained finite and showed no blow-up.

## 5. Contact-point normal kinematics

Kinematics use the actual contacting node and active wall-triangle normal with `Vcontact = VCOM + omega x rcontact`. Incoming normal speed was -887.577908 mm/s and outgoing normal speed was +551.982120 mm/s. Tangential speed changed from 617.501 to 26.998 mm/s. The post-impact gap slope was +551.979979 mm/s, so the body opened away from the wall rather than remaining attached.

## 6. Restitution

The contact-point result is `e_n=0.621897`. The independent gap-slope estimate is 0.621902. Relative to ζ=0.20, restitution fell by 0.203974 absolute or 24.698%. The candidate lies inside the provisional 0.30-0.65 numerical window; this window is not a measured material range.

## 7. Contact force and impulse

First contact ran from 1.0036 to 1.0052 ms, lasting 1.7 µs. Peak wall-force magnitude was 5.490071 N. Wall impulse magnitude was 5.337964354e-06 N s with vector [1.519258424e-06, 3.433546004e-06, -3.794269248e-06] N s.

## 8. Momentum closure

For the 0.95-1.10 ms audit window, `m DeltaV = Jmag + Jhydro + Jwall` closed with relative residual 3.758e-08; the clean pre-impact residual was 2.272e-08. Robot-side plus wall-side CFT had maximum action/reaction error 0.000e+00 N. These are below the preferred `1e-4` residual criterion.

## 9. Energy transfer

Across first-impact shoulders, Delta Ktrans=1.423204 µJ, Delta Krot=-2.352690 µJ, and Delta Ktotal=-0.929486 µJ. Wmag=0.007262 µJ and Whydro=-4.536010e-06 µJ. Direct contact work decomposes into Wnormal=-0.898948 µJ and combined Wtangential=-0.037795 µJ, giving total Wcontact=-0.936743 µJ. The independent rigid-body energy balance gives -0.936743 µJ. Total kinetic energy retained was 63.497%. Contact is strictly dissipative.

## 10. Penetration

The exact minimum first-event node-to-triangle gap was -0.513129 µm. The minimum over every 0.1 µs increment of the entire 1.3 ms trajectory was -0.928826 µm, in the second event. Both satisfy the preferred `>-1 µm` criterion.

## 11. Sticking / reopening

The longest solver-active interval was 1.8 µs, far below 0.10 ms. Persistent sticking did not occur. After the first event, the exact gap reopened to as much as 68.016449 µm before the next event.

## 12. Second impact

A second impact occurred at 1.2551-1.2568 ms for 1.8 µs. Its peak wall force was 5.827337 N, impulse was 6.405741943e-06 N s, and contact-point restitution was 0.664311.

## 13. 0.20 versus 0.50

| Metric | ζ=0.20, tangent=0 | ζ=0.50, tangent=0 |
|---|---:|---:|
| Normal restitution | 0.825871 | 0.621897 |
| Gap-slope restitution | 0.825883 | 0.621902 |
| COM recoil (mm/s) | 601.242 | 533.797 |
| First duration (µs) | 1.7 | 1.7 |
| First peak force (N) | 5.686401 | 5.490071 |
| First wall impulse (N s) | 6.012417017e-06 | 5.337964354e-06 |
| Wnormal (µJ) | -0.459107 | -0.898948 |
| Total KE retained | 80.539% | 63.497% |

The candidate COM recoil vector was [151.925554, 343.354672, -379.427314] mm/s, with magnitude 533.797 mm/s. This is 11.218% below ζ=0.20. Recoil is reported only as a secondary diagnostic because rotation, moment arm, tangential motion, friction, and contact geometry also contribute.

## 14. Decision

The result is **`PROVISIONAL_NORMAL_CONTACT_DAMPING_ACCEPTABLE_FOR_WOBBLE_TEST`**. The primary restitution is inside 0.30-0.65; momentum closes; total contact work is non-positive; the entire-window penetration is preferred-level; the gap reopens; no contact interval approaches 0.10 ms; and fixed `dt=1e-7 s` remained stable. Thus the same ζ=0.50, tangent=0 setting is eligible for a separately authorized 8.333 ms wobble validation. This does not establish experimental realism.

## 15. Exactly one next step

Run one 8.333 ms Reduced-Hydro sustained-wobble validation with exactly the same ζ=0.50, tangent=0 contact, magnetic drive, hydrodynamic coefficients, robot, wall, friction, and timestep. Evaluate true-axis local phase progression, sustained wobble, phase plateau, actual tilt, contact timeline, and absence of a persistent opposing-wall bridge. Do not perform an automatic damping sweep.

## Output traceability

All reported and plotted quantities trace to public CSV files in this directory. ODB, solver binaries, private extraction arrays, and raw Socket/Hydro telemetry are excluded from Git.
