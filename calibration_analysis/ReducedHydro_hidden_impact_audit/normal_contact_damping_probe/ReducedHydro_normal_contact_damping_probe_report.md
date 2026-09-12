# Reduced-Hydro normal-contact damping probe

**Status:** `PROVISIONAL_CONTACT_DISSIPATION_PROBE`

**Decision:** `NORMAL_CONTACT_STILL_TOO_ELASTIC`
**Long-run gate:** **Do not run the 8.333 ms wobble test with this candidate.**

## 1. Previous force-transfer issue closed

Commit `962fdff34bc2859f3b50406396d7c883c505d6d7` closed the Socket/CLOAD, unit, mass, RP-COM, hidden-impact, momentum, and first-impact energy audits. The accepted conclusion remains `FORCE_TRANSFER_VALID_HIDDEN_WALL_IMPACT_CONFIRMED`; this probe does not reopen those questions.

## 2. Why the current rebound requires contact calibration

The baseline contact-point normal restitution is 0.944097, and its first-impact COM recoil is 640.315 mm/s. The collision is energetically consistent but nearly elastic, so it can distort later wobble dynamics even though the force chain is correct.

## 3. Current Abaqus contact property

The baseline uses hard pressure-overclosure, General Contact penalty enforcement, friction coefficient 0.03, and critical contact damping fraction **0.055**. No penalty stiffness scaling or soft pressure-overclosure is used.

## 4. Tangent-fraction audit

The baseline input does not specify `TANGENT FRACTION`; Abaqus/Explicit 2025 therefore applies its default **1.0**. Consequently, baseline CFS can include both Coulomb friction and tangential viscous contact damping. Available CFS output cannot separate those two contributions.

## 5. Current effective normal restitution

Using the active wall triangle, true robot contact node, rigid-body COM velocity, angular velocity, and `V_contact = V_COM + omega x r_c`:

- incoming normal velocity: -887.577908 mm/s
- outgoing normal velocity: +837.959840 mm/s
- contact-point restitution: 0.944097
- independent gap-slope proxy: 0.944104
- incoming/outgoing tangential speed: 617.501 / 147.652 mm/s
- normal work: -0.145818 µJ
- combined tangential work: -0.049298 µJ

## 6. Experimental target search

The repository contains no credible experimental restitution coefficient, high-speed rebound measurement, or contact-compliance dataset for this robot-wall pair. This run is therefore not an experimental calibration.

## 7. Why 0.20 is provisional if no experiment exists

The candidate damping fraction 0.20 is a single engineering probe intended to test whether moderate subcritical normal damping suppresses the near-elastic rebound. The design window 0.3-0.65 is a numerical screening window, not a measured material property.

## 8. Datacheck / timestep stability

Abaqus 2025 datacheck passed with 0 errors, the same 13 warnings as baseline, and no initial overclosure. The single 1.3 ms dynamic run completed successfully at fixed direct `dt = 1e-7 s`. Energy remained finite through both impacts; no timestep adjustment or retry was performed.

## 9. First-impact normal kinematics

The candidate first contact starts at 1.0036 ms and ends at 1.0052 ms, for 1.7 µs. Peak wall force is 5.686401 N and wall impulse is 6.012417017e-06 N s. The incoming normal velocity is unchanged at -887.577908 mm/s, while outgoing velocity falls to +733.025105 mm/s, giving `e_n = 0.825871` and a gap-slope proxy of 0.825883.

## 10. Momentum closure

For 0.95-1.10 ms, `m DeltaV = Jmag + Jhydro + Jwall` closes with relative residual **7.374e-09**. The clean pre-impact residual is 2.272e-08. Robot-side plus wall-side CFT has a maximum action/reaction error of 0.000e+00 N. These checks exceed the preferred numerical requirement by a wide margin.

## 11. Energy transfer

Across the first-impact shoulders:

- Delta Ktrans = 1.805782 µJ
- Delta Krot = -2.301331 µJ
- Delta Ktotal = -0.495549 µJ
- Wmag = +0.007159 µJ
- Whydro = -4.846040e-06 µJ
- inferred net Wcontact = -0.502703 µJ
- total kinetic-energy retained = 80.539%

Contact does not create net energy. The contact-point decomposition gives Wnormal = -0.459107 µJ and combined Wtangential = -0.043595 µJ; their sum agrees with the independent rigid-body energy balance to rounding.

## 12. Gap / penetration

The first-impact minimum exact node-to-triangle gap is -0.555829 µm. The minimum over both 1.3 ms contact events is -0.917839 µm, occurring in the second event. Both remain inside the preferred `>-1 µm` penetration gate.

## 13. Contact duration and second impact

No solver-active contact interval exceeds 0.10 ms, so there is no persistent sticking. A second event remains: 1.2648-1.2665 ms, duration 1.8 µs, peak force 5.803292 N, and impulse 6.765676478e-06 N s.

## 14. Baseline versus candidate

| Metric | Baseline 0.055, tangent 1 | Probe 0.20, tangent 0 |
|---|---:|---:|
| Normal restitution | 0.944097 | 0.825871 |
| Gap-slope restitution | 0.944104 | 0.825883 |
| COM recoil (mm/s) | 640.315 | 601.242 |
| Recoil reduction | - | 6.10% |
| Peak first force (N) | 5.891618 | 5.686401 |
| First duration (µs) | 1.7 | 1.7 |
| Minimum first-event gap (µm) | -0.582271 | -0.555829 |
| Normal work (µJ) | -0.145818 | -0.459107 |

Only two contact parameters changed: critical damping fraction `0.055 -> 0.20` and tangent fraction `implicit 1.0 -> explicit 0.0`. Magnetic forcing, Reduced-Hydro coefficients, robot geometry/mass/inertia/pose, wall geometry, friction, General Contact, hard pressure-overclosure, surfaces, penalty formulation, duration, and timestep are identical.

## 15. Decision

The candidate is **`NORMAL_CONTACT_STILL_TOO_ELASTIC`**. Although dissipation increased and all numerical validity gates pass, `e_n = 0.825871` is above 0.70 and recoil remains 601.242 mm/s, well above 400 mm/s. The recoil reduction is only 6.10%, far below the preferred 40-50%. It is therefore ineligible for the long wobble test.

## 16. Exactly one next step

Do not run 8.333 ms and do not launch an automatic damping sweep. The next step is to choose one stronger normal-only damping candidate using these measured restitution/work results, document that choice as another provisional probe, and obtain explicit user approval before any new dynamic run.

## Output traceability

All plotted quantities come from the CSV files in this directory. Private NPZ extraction arrays, ODB, solver binaries, and raw Socket/Hydro telemetry are intentionally excluded from Git.
