"""Write the geometry-candidate comparison and decision report."""
from pathlib import Path
import json

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "ReducedHydro_zeta050_8p333_validation"
DECISION = "L1800_PARTIAL_WOBBLE_SURVIVAL_NO_PERSISTENT_BRIDGE"


def main():
    new = json.loads((HERE / "L1800_8p333_summary.json").read_text())
    old = json.loads((BASE / "zeta050_8p333_summary.json").read_text())
    identity = json.loads((HERE / "L1800_input_identity.json").read_text())
    events = pd.read_csv(HERE / "L1800_8p333_contact_events.csv")
    phase = pd.read_csv(HERE / "L1800_8p333_local_phase.csv")
    bridge = pd.read_csv(HERE / "L1800_8p333_bridge_timeline.csv")
    energy = pd.read_csv(HERE / "L1800_8p333_energy.csv")
    torque = pd.read_csv(HERE / "L1800_8p333_magnetic_hydro_torque.csv")
    windows = pd.read_csv(HERE / "L1800_8p333_phase_windows.csv")

    comparison = pd.DataFrame([
        {"case": "Frozen baseline", "pca_length_mm": identity["base_pca_axial_span_mm"], **old},
        {"case": "L=1.800 mm", "pca_length_mm": identity["candidate_pca_axial_span_mm"], **new},
    ])
    columns = ["case", "job", "pca_length_mm", "contact_event_count", "longest_contact_us",
               "min_exact_gap_um", "first_bridge_s", "longest_bridge_ms", "persistent_bridge",
               "max_tilt_deg", "late_tilt_deg", "postimpact_phase_advance_deg",
               "longest_phase_plateau_ms", "canonical_delta_s_mm", "final_Vt_mm_s"]
    comparison[columns].to_csv(HERE / "baseline_vs_L1800_summary.csv", index=False)

    plateau_times = bridge.loc[bridge.phase_plateau.astype(bool), "time_s"].to_numpy()
    phase_peak = float(phase.phi_robot_advance_deg.max())
    late = windows[windows.center_s >= .007]
    etotal_ratio = float(energy.Abaqus_ETOTAL_J.abs().max() / energy.Abaqus_ALLWK_J.abs().max())
    ke_error = float(abs(energy.Ktotal_rigid_J.iloc[-1] - energy.Abaqus_ALLKE_J.iloc[-1]))
    first = events.iloc[0]
    report = f"""# Reduced-Hydrodynamics L=1.800 mm geometry validation

## 1. Decision

`{DECISION}`

Shortening the already validated zeta=0.50 robot removes the persistent opposing-wall bridge, but it does not yet produce sustained post-impact wobble. The robot advances in local true-axis phase, then slows and reverses late while reaching a near-transverse posture.

## 2. Corrected baseline geometry identity

The validated zeta=0.50 INP already contains a 0.90 transverse PCA shrink. Direct mesh audit gives a baseline PCA envelope of **{identity['base_pca_axial_span_mm']:.9f} x {identity['base_pca_transverse_diameter_mm']:.9f} mm**, not the legacy pre-shrink 2.931414686 x 0.903074896 mm description. This candidate changes only the PCA axial span to **{identity['candidate_pca_axial_span_mm']:.9f} mm** and leaves the transverse coordinates unchanged. Job labels and Socket metadata now use D0815 to avoid perpetuating the old identity error.

## 3. Frozen physics and scaling

Contact remains SmoothWall114, HARD penalty, mu=0.03, zeta=0.50 and zero tangent damping. The 30 Hz, 10 mT, cone 30 deg, Bias 40 deg, phase 248 deg, sense +1, G=6 mT and L=45 mm drive is unchanged. Reduced-Hydro remains Cparallel=4e-9, Cperp=1.2e-8 N s/mm, Kspin=1e-9 and Kwobble=3e-9 N mm s. Constant density and magnetization give mass **{identity['candidate_mass_mg']:.9f} mg** and explicit moment **{identity['candidate_magnetic_moment_Am2']:.12f} A m2**, both scaled by {identity['mass_volume_ratio']:.9f}. COM and RP are unchanged.

## 4. Eligibility and run completion

Old jam-pose replay found no opposing bridge at 20, 5 or 0 um and a minimum replayed gap of +135.885 um. Datacheck had zero errors and no initial overclosure. The one permitted dynamic job completed to **{new['completed_s']*1e3:.3f} ms** in **{new['increments']:,}** direct increments at fixed dt={new['dt_s']:.1e} s. There was no retry or timestep change.

## 5. Contact and exact gap

First physical contact began at **{first.start_s*1e3:.4f} ms**, lasted **{first.duration_us:.1f} us**, peaked at **{first.peak_force_N:.6f} N**, and had effective separable restitution {first.effective_restitution:.6f}. There were **{new['contact_event_count']}** events; the longest lasted **{new['longest_contact_us']:.1f} us**. The minimum exact node-to-triangle gap was **{new['min_exact_gap_um']:.6f} um**. This is 0.357 um beyond the previous preferred -1 um penetration target, so contact remains provisional even though the run is numerically valid.

## 6. Opposing-wall bridge

The first opposing bridge occurred at **{new['first_bridge_s']*1e3:.4f} ms** and the longest lasted **{new['longest_bridge_ms']:.4f} ms**. No episode exceeded 0.5 ms and none persisted to the end. The baseline longest bridge was {old['longest_bridge_ms']:.4f} ms. Axial shortening therefore removed the previously dominant geometric jam mechanism.

## 7. True-axis phase and tilt

After first contact, local true-axis phase advanced **{new['postimpact_phase_advance_deg']:.4f} deg**; total advance peaked at **{phase_peak:.4f} deg** and ended at {phase.phi_robot_advance_deg.iloc[-1]:.4f} deg. Late 0.75 ms window rates ranged from {late.robot_phase_rate_Hz.min():.3f} to {late.robot_phase_rate_Hz.max():.3f} Hz and became negative. Maximum tilt was **{new['max_tilt_deg']:.4f} deg**, final tilt **{new['late_tilt_deg']:.4f} deg**, and the 7-8.333 ms mean was {pd.read_csv(HERE / 'L1800_8p333_true_axis.csv').query('time_s >= 0.007').tilt_deg.mean():.4f} deg. This is dynamic tumbling/return, not the old 33 deg bridged lock, but it is not sustained wobble.

## 8. Plateau boundary

The longest qualifying support interval is **{new['longest_phase_plateau_ms']:.4f} ms**. Only {len(plateau_times)} center samples qualify, in two isolated pairs near {plateau_times.min()*1e3:.4f} and {plateau_times.max()*1e3:.4f} ms. The result is only 0.1 us above the 0.5 ms definition and is reported as a marginal threshold crossing, not a long persistent plateau.

## 9. Torque and energy

Pre/post-impact magnetic torque RMS is **{new['preimpact_Tmag_RMS_Nmm']:.8f}/{new['postimpact_Tmag_RMS_Nmm']:.8f} N mm**; it does not collapse. Full-window magnetic/hydrodynamic torque RMS is {np.sqrt(np.mean(torque.Tmag_norm_Nmm**2)):.8f}/{np.sqrt(np.mean(torque.Thydro_norm_Nmm**2)):.8e} N mm. `Thydro dot omega` is never positive above tolerance. Maximum |ETOTAL| is {energy.Abaqus_ETOTAL_J.abs().max():.6e} J ({etotal_ratio*100:.5f}% of peak |ALLWK|), final rigid/ODB kinetic-energy mismatch is {ke_error:.3e} J, and ALLAE remains zero. No numerical energy growth is detected.

## 10. Translation

Canonical displacement is **{new['canonical_delta_s_mm']:.6f} mm** and final Vt is **{new['final_Vt_mm_s']:.6f} mm/s** under increasing-s FORWARD. Translation remains diagnostic, not the wobble gate.

## 11. Baseline comparison

| Metric | Frozen baseline | L=1.800 mm |
|---|---:|---:|
| Contact events | {old['contact_event_count']} | {new['contact_event_count']} |
| Longest contact (us) | {old['longest_contact_us']:.1f} | {new['longest_contact_us']:.1f} |
| Longest bridge (ms) | {old['longest_bridge_ms']:.4f} | {new['longest_bridge_ms']:.4f} |
| Persistent bridge | {old['persistent_bridge']} | {new['persistent_bridge']} |
| Post-impact phase advance (deg) | {old['postimpact_phase_advance_deg']:.4f} | {new['postimpact_phase_advance_deg']:.4f} |
| Maximum/final tilt (deg) | {old['max_tilt_deg']:.2f}/{old['late_tilt_deg']:.2f} | {new['max_tilt_deg']:.2f}/{new['late_tilt_deg']:.2f} |
| Canonical delta s (mm) | {old['canonical_delta_s_mm']:.4f} | {new['canonical_delta_s_mm']:.4f} |

The geometry change reduces contact events by {(1-new['contact_event_count']/old['contact_event_count'])*100:.1f}% and removes the long bridge, but does not increase net post-impact phase advance over this quarter-cycle.

## 12. Exactly one next step

**Audit the existing L=1.800 mm trajectory's rotational response and Reduced-Hydro rotational coefficients without another Abaqus run.** The next question is why the free, non-bridged robot approaches 90 deg tilt and reverses phase late. Do not change magnetic forcing, contact damping, or geometry again until that audit separates expected tumble dynamics from a rotational-model limitation. As before, 8.333 ms is only one quarter of a 30 Hz cycle and cannot establish frequency lock.
"""
    (HERE / "Wobble30Hz_ReducedHydro_L1800_8p333_report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
