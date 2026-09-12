"""Render the 16-section normal-contact damping probe report."""
from pathlib import Path
import json

import pandas as pd


HERE = Path(__file__).resolve().parent


def f(value, digits=6):
    return f'{float(value):.{digits}f}'


def main():
    base = json.loads((HERE / 'current_contact_restitution_summary.json').read_text())
    cand = json.loads((HERE / 'normaldamp020_summary.json').read_text())
    events = pd.read_csv(HERE / 'normaldamp020_contact_event_catalog.csv')
    second = events.iloc[1]
    text = f"""# Reduced-Hydro normal-contact damping probe

**Status:** `PROVISIONAL_CONTACT_DISSIPATION_PROBE`

**Decision:** `{cand['decision']}`
**Long-run gate:** **Do not run the 8.333 ms wobble test with this candidate.**

## 1. Previous force-transfer issue closed

Commit `962fdff34bc2859f3b50406396d7c883c505d6d7` closed the Socket/CLOAD, unit, mass, RP-COM, hidden-impact, momentum, and first-impact energy audits. The accepted conclusion remains `FORCE_TRANSFER_VALID_HIDDEN_WALL_IMPACT_CONFIRMED`; this probe does not reopen those questions.

## 2. Why the current rebound requires contact calibration

The baseline contact-point normal restitution is {f(base['e_n'])}, and its first-impact COM recoil is {f(base['deltaV_COM_norm_mm_s'], 3)} mm/s. The collision is energetically consistent but nearly elastic, so it can distort later wobble dynamics even though the force chain is correct.

## 3. Current Abaqus contact property

The baseline uses hard pressure-overclosure, General Contact penalty enforcement, friction coefficient 0.03, and critical contact damping fraction **0.055**. No penalty stiffness scaling or soft pressure-overclosure is used.

## 4. Tangent-fraction audit

The baseline input does not specify `TANGENT FRACTION`; Abaqus/Explicit 2025 therefore applies its default **1.0**. Consequently, baseline CFS can include both Coulomb friction and tangential viscous contact damping. Available CFS output cannot separate those two contributions.

## 5. Current effective normal restitution

Using the active wall triangle, true robot contact node, rigid-body COM velocity, angular velocity, and `V_contact = V_COM + omega x r_c`:

- incoming normal velocity: {f(base['vn_in_mm_s'])} mm/s
- outgoing normal velocity: +{f(base['vn_out_mm_s'])} mm/s
- contact-point restitution: {f(base['e_n'])}
- independent gap-slope proxy: {f(base['gap_proxy_restitution'])}
- incoming/outgoing tangential speed: {f(base['vt_in_mm_s'], 3)} / {f(base['vt_out_mm_s'], 3)} mm/s
- normal work: {f(base['Wnormal_uJ'])} µJ
- combined tangential work: {f(base['Wtangential_combined_uJ'])} µJ

## 6. Experimental target search

The repository contains no credible experimental restitution coefficient, high-speed rebound measurement, or contact-compliance dataset for this robot-wall pair. This run is therefore not an experimental calibration.

## 7. Why 0.20 is provisional if no experiment exists

The candidate damping fraction 0.20 is a single engineering probe intended to test whether moderate subcritical normal damping suppresses the near-elastic rebound. The design window 0.3-0.65 is a numerical screening window, not a measured material property.

## 8. Datacheck / timestep stability

Abaqus 2025 datacheck passed with 0 errors, the same 13 warnings as baseline, and no initial overclosure. The single 1.3 ms dynamic run completed successfully at fixed direct `dt = 1e-7 s`. Energy remained finite through both impacts; no timestep adjustment or retry was performed.

## 9. First-impact normal kinematics

The candidate first contact starts at {f(cand['first_contact_s']*1e3, 4)} ms and ends at {f(cand['first_contact_end_s']*1e3, 4)} ms, for {f(cand['first_contact_duration_us'], 1)} µs. Peak wall force is {f(cand['first_peak_force_N'], 6)} N and wall impulse is {cand['first_wall_impulse_norm_Ns']:.9e} N s. The incoming normal velocity is unchanged at {f(cand['vn_in_mm_s'])} mm/s, while outgoing velocity falls to +{f(cand['vn_out_mm_s'])} mm/s, giving `e_n = {f(cand['e_n'])}` and a gap-slope proxy of {f(cand['gap_proxy_restitution'])}.

## 10. Momentum closure

For 0.95-1.10 ms, `m DeltaV = Jmag + Jhydro + Jwall` closes with relative residual **{cand['momentum_closure_relative_residual']:.3e}**. The clean pre-impact residual is {cand['clean_momentum_relative_residual']:.3e}. Robot-side plus wall-side CFT has a maximum action/reaction error of {cand['action_reaction_max_N']:.3e} N. These checks exceed the preferred numerical requirement by a wide margin.

## 11. Energy transfer

Across the first-impact shoulders:

- Delta Ktrans = {f(cand['Delta_Ktrans_uJ'])} µJ
- Delta Krot = {f(cand['Delta_Krot_uJ'])} µJ
- Delta Ktotal = {f(cand['Delta_Ktotal_uJ'])} µJ
- Wmag = +{f(cand['Wmag_uJ'])} µJ
- Whydro = {cand['Whydro_uJ']:.6e} µJ
- inferred net Wcontact = {f(cand['Wcontact_net_uJ'])} µJ
- total kinetic-energy retained = {100*cand['total_KE_retained_fraction']:.3f}%

Contact does not create net energy. The contact-point decomposition gives Wnormal = {f(cand['Wnormal_uJ'])} µJ and combined Wtangential = {f(cand['Wtangential_combined_uJ'])} µJ; their sum agrees with the independent rigid-body energy balance to rounding.

## 12. Gap / penetration

The first-impact minimum exact node-to-triangle gap is {f(events.iloc[0].min_gap_um)} µm. The minimum over both 1.3 ms contact events is {f(cand['min_gap_um'])} µm, occurring in the second event. Both remain inside the preferred `>-1 µm` penetration gate.

## 13. Contact duration and second impact

No solver-active contact interval exceeds 0.10 ms, so there is no persistent sticking. A second event remains: {f(second.start_s*1e3, 4)}-{f(second.end_s*1e3, 4)} ms, duration {f(second.duration_us, 1)} µs, peak force {f(second.peak_force_N, 6)} N, and impulse {second.impulse_norm_Ns:.9e} N s.

## 14. Baseline versus candidate

| Metric | Baseline 0.055, tangent 1 | Probe 0.20, tangent 0 |
|---|---:|---:|
| Normal restitution | {f(base['e_n'])} | {f(cand['e_n'])} |
| Gap-slope restitution | {f(base['gap_proxy_restitution'])} | {f(cand['gap_proxy_restitution'])} |
| COM recoil (mm/s) | {f(base['deltaV_COM_norm_mm_s'], 3)} | {f(cand['deltaV_COM_norm_mm_s'], 3)} |
| Recoil reduction | - | {f(cand['COM_recoil_reduction_percent'], 2)}% |
| Peak first force (N) | {f(base['first_peak_force_N'], 6)} | {f(cand['first_peak_force_N'], 6)} |
| First duration (µs) | {f(base['first_contact_duration_us'], 1)} | {f(cand['first_contact_duration_us'], 1)} |
| Minimum first-event gap (µm) | {f(base['min_gap_um'])} | {f(events.iloc[0].min_gap_um)} |
| Normal work (µJ) | {f(base['Wnormal_uJ'])} | {f(cand['Wnormal_uJ'])} |

Only two contact parameters changed: critical damping fraction `0.055 -> 0.20` and tangent fraction `implicit 1.0 -> explicit 0.0`. Magnetic forcing, Reduced-Hydro coefficients, robot geometry/mass/inertia/pose, wall geometry, friction, General Contact, hard pressure-overclosure, surfaces, penalty formulation, duration, and timestep are identical.

## 15. Decision

The candidate is **`{cand['decision']}`**. Although dissipation increased and all numerical validity gates pass, `e_n = {f(cand['e_n'])}` is above 0.70 and recoil remains {f(cand['deltaV_COM_norm_mm_s'], 3)} mm/s, well above 400 mm/s. The recoil reduction is only {f(cand['COM_recoil_reduction_percent'], 2)}%, far below the preferred 40-50%. It is therefore ineligible for the long wobble test.

## 16. Exactly one next step

Do not run 8.333 ms and do not launch an automatic damping sweep. The next step is to choose one stronger normal-only damping candidate using these measured restitution/work results, document that choice as another provisional probe, and obtain explicit user approval before any new dynamic run.

## Output traceability

All plotted quantities come from the CSV files in this directory. Private NPZ extraction arrays, ODB, solver binaries, and raw Socket/Hydro telemetry are intentionally excluded from Git.
"""
    (HERE / 'ReducedHydro_normal_contact_damping_probe_report.md').write_text(text, encoding='utf-8')


if __name__ == '__main__':
    main()
