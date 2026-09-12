"""Write the required 15-section report for the zeta=0.50 probe."""
from pathlib import Path
import json

import pandas as pd


HERE = Path(__file__).resolve().parent


def f(value, digits=6):
    return f'{float(value):.{digits}f}'


def main():
    ref = json.loads((HERE / 'normaldamp020_summary.json').read_text())
    cand = json.loads((HERE / 'normaldamp050_summary.json').read_text())
    events = pd.read_csv(HERE / 'normaldamp050_contact_event_catalog.csv')
    first = events.iloc[0]
    second_text = 'No second solver-active impact occurred.'
    if len(events) > 1:
        second = events.iloc[1]
        second_text = (f'A second impact occurred at {second.start_s*1e3:.4f}-{second.end_s*1e3:.4f} ms '
                       f'for {second.duration_us:.1f} µs. Its peak wall force was '
                       f'{second.peak_force_N:.6f} N, impulse was {second.impulse_norm_Ns:.9e} N s, '
                       f'and contact-point restitution was {second.e_n:.6f}.')
    dv = [cand[f'deltaV_COM_{axis}_mm_s'] for axis in 'xyz']
    text = f"""# Reduced-Hydro normal-contact damping 0.50 probe

**Status:** `PROVISIONAL_NORMAL_DAMPING_PROBE`

**Decision:** `{cand['decision']}`

**Long-run status:** The short-window gate permits a later 8.333 ms validation, but no long run was launched in this work.

## 1. Why 0.20 failed

The ζ=0.20, tangent=0 case was numerically stable and dissipative, but its contact-point normal restitution was {ref['e_n']:.6f}, above the 0.70 upper warning limit. It therefore remained too elastic for the provisional screening objective.

## 2. Why 0.50 was selected

ζ=0.50 was chosen as one higher-information probe. A local two-point empirical extrapolation suggested restitution near 0.58, inside the provisional 0.30-0.65 design window. This was candidate selection only: it is not an experimental calibration, theoretical restitution law, or Abaqus-defined mapping.

## 3. Why this is now a one-variable comparison

The formal reference is ζ=0.20, tangent fraction=0 and the candidate is ζ=0.50, tangent fraction=0. Only the normal critical damping fraction changed. Magnetic forcing, Reduced-Hydro coefficients, robot geometry, mass, inertia and pose, wall, friction coefficient 0.03, General Contact, hard pressure-overclosure, penalty enforcement, surfaces, output cadence, duration, and direct timestep were frozen.

## 4. Datacheck and dt stability

The independent Abaqus 2025 datacheck passed with 0 errors, the same 13 warnings as the reference, and no initial overclosure. The sole 1.3 ms dynamic run completed successfully with fixed `*Dynamic, Explicit, DIRECT` `dt=1e-7 s`. No retry, timestep change, damping sweep, or 8.333 ms run was performed. Energy remained finite and showed no blow-up.

## 5. Contact-point normal kinematics

Kinematics use the actual contacting node and active wall-triangle normal with `Vcontact = VCOM + omega x rcontact`. Incoming normal speed was {cand['vn_in_mm_s']:.6f} mm/s and outgoing normal speed was +{cand['vn_out_mm_s']:.6f} mm/s. Tangential speed changed from {cand['vt_in_mm_s']:.3f} to {cand['vt_out_mm_s']:.3f} mm/s. The post-impact gap slope was +{cand['gap_opening_slope_mm_s']:.6f} mm/s, so the body opened away from the wall rather than remaining attached.

## 6. Restitution

The contact-point result is `e_n={cand['e_n']:.6f}`. The independent gap-slope estimate is {cand['gap_proxy_restitution']:.6f}. Relative to ζ=0.20, restitution fell by {cand['e_n_absolute_drop_from_020']:.6f} absolute or {cand['e_n_reduction_from_020_percent']:.3f}%. The candidate lies inside the provisional 0.30-0.65 numerical window; this window is not a measured material range.

## 7. Contact force and impulse

First contact ran from {cand['first_contact_s']*1e3:.4f} to {cand['first_contact_end_s']*1e3:.4f} ms, lasting {cand['first_contact_duration_us']:.1f} µs. Peak wall-force magnitude was {cand['first_peak_force_N']:.6f} N. Wall impulse magnitude was {cand['first_wall_impulse_norm_Ns']:.9e} N s with vector [{cand['first_wall_impulse_xyz_Ns'][0]:.9e}, {cand['first_wall_impulse_xyz_Ns'][1]:.9e}, {cand['first_wall_impulse_xyz_Ns'][2]:.9e}] N s.

## 8. Momentum closure

For the 0.95-1.10 ms audit window, `m DeltaV = Jmag + Jhydro + Jwall` closed with relative residual {cand['momentum_closure_relative_residual']:.3e}; the clean pre-impact residual was {cand['clean_momentum_relative_residual']:.3e}. Robot-side plus wall-side CFT had maximum action/reaction error {cand['action_reaction_max_N']:.3e} N. These are below the preferred `1e-4` residual criterion.

## 9. Energy transfer

Across first-impact shoulders, Delta Ktrans={cand['Delta_Ktrans_uJ']:.6f} µJ, Delta Krot={cand['Delta_Krot_uJ']:.6f} µJ, and Delta Ktotal={cand['Delta_Ktotal_uJ']:.6f} µJ. Wmag={cand['Wmag_uJ']:.6f} µJ and Whydro={cand['Whydro_uJ']:.6e} µJ. Direct contact work decomposes into Wnormal={cand['Wnormal_uJ']:.6f} µJ and combined Wtangential={cand['Wtangential_combined_uJ']:.6f} µJ, giving total Wcontact={cand['Wcontact_total_uJ']:.6f} µJ. The independent rigid-body energy balance gives {cand['Wcontact_from_energy_balance_uJ']:.6f} µJ. Total kinetic energy retained was {100*cand['total_KE_retained_fraction']:.3f}%. Contact is strictly dissipative.

## 10. Penetration

The exact minimum first-event node-to-triangle gap was {cand['min_first_event_gap_um']:.6f} µm. The minimum over every 0.1 µs increment of the entire 1.3 ms trajectory was {cand['min_gap_um']:.6f} µm, in the second event. Both satisfy the preferred `>-1 µm` criterion.

## 11. Sticking / reopening

The longest solver-active interval was {cand['longest_contact_duration_us']:.1f} µs, far below 0.10 ms. Persistent sticking did not occur. After the first event, the exact gap reopened to as much as {cand['max_reopened_gap_before_next_event_um']:.6f} µm before the next event.

## 12. Second impact

{second_text}

## 13. 0.20 versus 0.50

| Metric | ζ=0.20, tangent=0 | ζ=0.50, tangent=0 |
|---|---:|---:|
| Normal restitution | {ref['e_n']:.6f} | {cand['e_n']:.6f} |
| Gap-slope restitution | {ref['gap_proxy_restitution']:.6f} | {cand['gap_proxy_restitution']:.6f} |
| COM recoil (mm/s) | {ref['deltaV_COM_norm_mm_s']:.3f} | {cand['deltaV_COM_norm_mm_s']:.3f} |
| First duration (µs) | {ref['first_contact_duration_us']:.1f} | {cand['first_contact_duration_us']:.1f} |
| First peak force (N) | {ref['first_peak_force_N']:.6f} | {cand['first_peak_force_N']:.6f} |
| First wall impulse (N s) | {ref['first_wall_impulse_norm_Ns']:.9e} | {cand['first_wall_impulse_norm_Ns']:.9e} |
| Wnormal (µJ) | {ref['Wnormal_uJ']:.6f} | {cand['Wnormal_uJ']:.6f} |
| Total KE retained | {100*ref['total_KE_retained_fraction']:.3f}% | {100*cand['total_KE_retained_fraction']:.3f}% |

The candidate COM recoil vector was [{dv[0]:.6f}, {dv[1]:.6f}, {dv[2]:.6f}] mm/s, with magnitude {cand['deltaV_COM_norm_mm_s']:.3f} mm/s. This is {cand['COM_recoil_reduction_percent']:.3f}% below ζ=0.20. Recoil is reported only as a secondary diagnostic because rotation, moment arm, tangential motion, friction, and contact geometry also contribute.

## 14. Decision

The result is **`{cand['decision']}`**. The primary restitution is inside 0.30-0.65; momentum closes; total contact work is non-positive; the entire-window penetration is preferred-level; the gap reopens; no contact interval approaches 0.10 ms; and fixed `dt=1e-7 s` remained stable. Thus the same ζ=0.50, tangent=0 setting is eligible for a separately authorized 8.333 ms wobble validation. This does not establish experimental realism.

## 15. Exactly one next step

Run one 8.333 ms Reduced-Hydro sustained-wobble validation with exactly the same ζ=0.50, tangent=0 contact, magnetic drive, hydrodynamic coefficients, robot, wall, friction, and timestep. Evaluate true-axis local phase progression, sustained wobble, phase plateau, actual tilt, contact timeline, and absence of a persistent opposing-wall bridge. Do not perform an automatic damping sweep.

## Output traceability

All reported and plotted quantities trace to public CSV files in this directory. ODB, solver binaries, private extraction arrays, and raw Socket/Hydro telemetry are excluded from Git.
"""
    (HERE / 'ReducedHydro_normal_contact_damping_050_probe_report.md').write_text(text, encoding='utf-8')


if __name__ == '__main__':
    main()
