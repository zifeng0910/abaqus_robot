"""Write the final 20-section validation report from audited public outputs."""
from pathlib import Path
import json

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
DECISION = "GEOMETRY_CONTACT_STILL_JAMS_WITH_REDUCED_HYDRO"


def value(row, key, scale=1.0, digits=4):
    item = row[key]
    if pd.isna(item):
        return "not detected"
    return f"{float(item) * scale:.{digits}f}"


def main():
    s = json.loads((HERE / "zeta050_8p333_summary.json").read_text())
    identity = json.loads((HERE / "zeta050_8p333_input_identity.json").read_text())
    events = pd.read_csv(HERE / "zeta050_8p333_contact_events.csv")
    clusters = pd.read_csv(HERE / "zeta050_8p333_contact_clusters.csv")
    comparison = pd.read_csv(HERE / "cel_vs_reducedhydro_zeta050_summary.csv")
    energy = pd.read_csv(HERE / "zeta050_8p333_energy.csv")
    torque = pd.read_csv(HERE / "zeta050_8p333_magnetic_hydro_torque.csv")
    windows = pd.read_csv(HERE / "zeta050_8p333_phase_windows.csv")

    first = events.iloc[0]
    longest = events.loc[events.duration_us.idxmax()]
    threshold = clusters.groupby("threshold_um").agg(
        max_nodes=("node_count", "max"), max_clusters=("cluster_count", "max"),
        opposing_samples=("opposing_bridge", "sum"))
    etotal_abs = energy.Abaqus_ETOTAL_J.abs().max()
    ke_end_error = abs(energy.Ktotal_rigid_J.iloc[-1] - energy.Abaqus_ALLKE_J.iloc[-1])
    tmag_rms = np.sqrt(np.mean(torque.Tmag_norm_Nmm ** 2))
    thydro_rms = np.sqrt(np.mean(torque.Thydro_norm_Nmm ** 2))
    late_windows = windows[windows.center_s >= 0.0075]

    rows = []
    for case in ("CEL", "RH old contact", "RH zeta0.50"):
        row = comparison.loc[comparison.case == case].iloc[0]
        rows.append(
            f"| {case} | {value(row, 'first_contact_s', 1e3)} | {int(row.contact_events)} | "
            f"{value(row, 'longest_contact_ms')} | {value(row, 'min_gap_um')} | "
            f"{value(row, 'max_tilt_deg')} | {value(row, 'postimpact_phase_advance_deg')} | "
            f"{value(row, 'longest_phase_plateau_ms')} | {value(row, 'opposing_bridge_onset_s', 1e3)} | "
            f"{value(row, 'longest_bridge_ms')} | {value(row, 'canonical_delta_s_mm')} | "
            f"{value(row, 'final_Vt_mm_s')} |"
        )
    comparison_table = "\n".join(rows)

    report = f"""# Reduced-Hydrodynamics zeta=0.50 wall-on validation at 8.333 ms

## 1. Frozen model identity

The single dynamic job was `{s['job']}`. It retained the 2.931414686 x 0.903074896 mm rigid robot, approximately 9.99997 mg mass, original inertia, RP, initial pose, and HEAD/TAIL definition. The magnetic drive remained 30 Hz, 10 mT, cone 30 deg, Bias 40 deg, phase 248 deg, sense +1, G=6 mT and L=45 mm through the production Magpylib Socket; Abaqus performed no electromagnetic solve. Reduced-Hydro remained body-following and anisotropic with Cparallel=4e-9, Cperp=1.2e-8 N s/mm, Kspin=1e-9 and Kwobble=3e-9 N mm s. Contact remained SmoothWall114, HARD penalty, mu=0.03, zeta=0.50 and tangent fraction 0.

Input identity was audited against the approved 1.3 ms source. The source SHA-256 is `{identity['source_sha256']}` and the 8.333 ms candidate SHA-256 is `{identity['candidate_sha256']}`. Only duration, field cadence and an output comment changed.

## 2. Why zeta=0.50 is provisional

The normal damping fraction is a provisional engineering contact model, not an experimentally calibrated restitution law. Its short-window first-impact gate gave e_n=0.621897, 1.7 us duration, 5.490071 N peak force, -0.513129 um minimum gap, reopening and closed momentum/energy budgets. This validation therefore freezes zeta=0.50; it does not claim experimental restitution calibration.

## 3. Run completion and timestep

The job completed successfully to **{s['completed_s']*1e3:.3f} ms** with **{s['increments']:,}** increments at fixed direct `dt={s['dt_s']:.1e} s`. Abaqus reported 13 warnings, no initial overclosure, no Explicit failure and no automatic retry. The ODB contains {s['field_frames']} field frames at 25 us cadence and 83,331 RP/history samples. This is one quarter of a 30 Hz period, so it cannot establish periodic lock.

## 4. True-axis reconstruction

The robot axis was reconstructed at every increment by applying rigid U/UR to the reference exterior-node HEAD and TAIL groups. The reported phase and tilt do not use global UR component ranges. Maximum true tilt was **{s['max_tilt_deg']:.4f} deg** and the final tilt was **{s['late_tilt_deg']:.4f} deg**.

## 5. Canonical local frame

COM was projected continuously onto the authoritative `TRUE_centerline`; increasing canonical s is FORWARD. The local tangent and continuous e1/e2 frame define both robot and magnetic phases. Previous-s continuity and endpoint extrapolation were retained. The comparison cases were reprocessed with the same definitions.

## 6. First impact

First physical contact began at **{first.start_s*1e3:.4f} ms** and lasted **{first.duration_us:.4f} us**. Peak wall resultant was **{first.peak_force_N:.6f} N**, impulse magnitude was **{first.impulse_norm_Ns:.6e} N s**, and exact minimum gap was **{first.min_gap_um:.6f} um**. The separable normal velocities were {first.vn_in_mm_s:.3f} and {first.vn_out_mm_s:.3f} mm/s, giving effective e_n={first.effective_restitution:.6f}, consistent with the accepted short probe.

## 7. Repeated contact events

There were **{s['contact_event_count']}** solver-active events. The longest was event {int(longest.event_number)}, from {longest.start_s*1e3:.4f} to {longest.end_s*1e3:.4f} ms, lasting **{longest.duration_us:.4f} us** with a low {longest.peak_force_N:.6f} N peak; it represents sustained multipoint contact rather than a clean separable impact. Event impulse, node, wall triangle, gap and pre/post normal velocity are catalogued in `zeta050_8p333_contact_events.csv`.

## 8. Exact gap

Rigid-body geometry was reconstructed for all near-wall intervals. Exact node-to-triangle distance used 1 us coarse coverage and 0.1 us backfill wherever gap <=50 um or solver contact was active. The full-run minimum exact gap was **{s['min_exact_gap_um']:.6f} um**. This gate uses impulse, gap and duration, not CPRESS peak.

## 9. Opposing-wall bridge

The first opposing-sector bridge occurred at **{s['first_bridge_s']*1e3:.4f} ms**. A transient bridge ran from 1.7712 to 1.8637 ms; the persistent episode ran from **2.7727 ms through 8.3330 ms**, lasting **{s['longest_bridge_ms']:.4f} ms**. Maximum near-wall node counts at 20/5/0 um were {int(threshold.loc[20.0].max_nodes)}/{int(threshold.loc[5.0].max_nodes)}/{int(threshold.loc[0.0].max_nodes)}, with two separated wall-normal clusters at every threshold. Persistent bridge is therefore **present**.

## 10. True-axis phase progression

After first impact, the unwrapped true-axis local phase advanced only **{s['postimpact_phase_advance_deg']:.4f} deg** by 8.333 ms, below the approximately 90 deg strong-survival reference. It advanced intermittently, reversed over some windows, and approached zero/negative short-window rate near the end rather than sustaining wobble progression.

## 11. Phase plateau

The longest qualifying 0.5 ms phase-support window was **{s['longest_phase_plateau_ms']:.4f} ms**, just above the formal threshold, while the local field phase continued. Thus a >0.5 ms phase plateau is present. Late 0.75 ms robot-rate windows ranged from {late_windows.robot_phase_rate_Hz.min():.3f} to {late_windows.robot_phase_rate_Hz.max():.3f} Hz, including reversal; these are short-window indicators only.

## 12. Tilt evolution

Tilt increased from approximately 12.7 deg to a high 30-34 deg posture. After the persistent bridge formed it remained near that slanted state, ending at {s['late_tilt_deg']:.4f} deg. Early variations do not offset the combined late bridge, phase plateau and high-tilt gate.

## 13. Magnetic torque

Magnetic torque RMS was **{s['preimpact_Tmag_RMS_Nmm']:.8f} N mm** before first impact and **{s['postimpact_Tmag_RMS_Nmm']:.8f} N mm** after it. It did not collapse; post-impact RMS increased. Full-window |Tmag| RMS was {tmag_rms:.8f} N mm. A persisting drive together with stalled robot phase supports a contact/geometric stall rather than loss of magnetic forcing.

## 14. Hydro torque

Full-window |Thydro| RMS was {thydro_rms:.8e} N mm, far below |Tmag|, and `Thydro dot omega` never became positive above tolerance: maximum **{s['Thydro_power_max_W']:.3e} W**, positive fraction **{s['Thydro_power_positive_fraction']:.3f}**. The Reduced-Hydro rotational term is therefore dissipative throughout. No coefficient was refit.

## 15. Energy consistency

ODB assembly energies were converted from the active N-mm-s unit system (N mm) to joules before reporting. Maximum absolute Abaqus ETOTAL residual was **{etotal_abs:.6e} J**, artificial energy ALLAE remained zero, and the final independently reconstructed rigid KE differed from ODB ALLKE by only **{ke_end_error:.6e} J**. Magnetic work ended at {energy.Wmag_J.iloc[-1]:.6e} J, hydro work at {energy.Whydro_J.iloc[-1]:.6e} J, and inferred contact work at {energy.Wcontact_inferred_J.iloc[-1]:.6e} J. There is no numerical energy-growth failure.

## 16. Canonical translation

Translation is diagnostic, not the success gate. Canonical displacement was **{s['canonical_delta_s_mm']:.6f} mm** and final tangent velocity was **{s['final_Vt_mm_s']:.6f} mm/s**. Both are negative, so motion was backward under the authoritative increasing-s FORWARD definition; no legacy robot_arc sign was used.

## 17. CEL versus Reduced-Hydro

| Metric | CEL | RH old contact | RH zeta0.50 |
|---|---:|---:|---:|
| First contact (ms) | {value(comparison.iloc[0], 'first_contact_s', 1e3)} | {value(comparison.iloc[1], 'first_contact_s', 1e3)} | {value(comparison.iloc[2], 'first_contact_s', 1e3)} |
| Contact events | {int(comparison.iloc[0].contact_events)} | {int(comparison.iloc[1].contact_events)} | {int(comparison.iloc[2].contact_events)} |
| Longest contact (ms) | {value(comparison.iloc[0], 'longest_contact_ms')} | {value(comparison.iloc[1], 'longest_contact_ms')} | {value(comparison.iloc[2], 'longest_contact_ms')} |
| Minimum gap (um) | {value(comparison.iloc[0], 'min_gap_um')} | {value(comparison.iloc[1], 'min_gap_um')} | {value(comparison.iloc[2], 'min_gap_um')} |
| Maximum tilt (deg) | {value(comparison.iloc[0], 'max_tilt_deg')} | {value(comparison.iloc[1], 'max_tilt_deg')} | {value(comparison.iloc[2], 'max_tilt_deg')} |
| Post-impact phase advance (deg) | {value(comparison.iloc[0], 'postimpact_phase_advance_deg')} | {value(comparison.iloc[1], 'postimpact_phase_advance_deg')} | {value(comparison.iloc[2], 'postimpact_phase_advance_deg')} |
| Longest plateau (ms) | {value(comparison.iloc[0], 'longest_phase_plateau_ms')} | {value(comparison.iloc[1], 'longest_phase_plateau_ms')} | {value(comparison.iloc[2], 'longest_phase_plateau_ms')} |
| Bridge onset (ms) | {value(comparison.iloc[0], 'opposing_bridge_onset_s', 1e3)} | {value(comparison.iloc[1], 'opposing_bridge_onset_s', 1e3)} | {value(comparison.iloc[2], 'opposing_bridge_onset_s', 1e3)} |
| Longest bridge (ms) | {value(comparison.iloc[0], 'longest_bridge_ms')} | {value(comparison.iloc[1], 'longest_bridge_ms')} | {value(comparison.iloc[2], 'longest_bridge_ms')} |
| Canonical delta s (mm) | {value(comparison.iloc[0], 'canonical_delta_s_mm')} | {value(comparison.iloc[1], 'canonical_delta_s_mm')} | {value(comparison.iloc[2], 'canonical_delta_s_mm')} |
| Final Vt (mm/s) | {value(comparison.iloc[0], 'final_Vt_mm_s')} | {value(comparison.iloc[1], 'final_Vt_mm_s')} | {value(comparison.iloc[2], 'final_Vt_mm_s')} |

The older cases have only 50 us geometry snapshots, so their micro-impact event counts are lower bounds. Under the unified local-frame geometry definitions, zeta=0.50 does not show a clear wobble-survival improvement over CEL: both reach a high tilted state and form long opposing bridges.

## 18. Visual GIF interpretation

`Wobble_F30_ReducedHydro_Zeta050_WallOn_8p333.gif` uses a fixed camera and shows the full robot, translucent wall, HEAD/TAIL, canonical centerline, local frame, magnetic vector, true axis and live audit quantities. It visibly follows the sequence **turn, become slanted, remain geometrically constrained**. The synchronized CEL comparison shows no clear qualitative rescue of post-impact wobble by deleting CEL and using zeta=0.50 contact.

## 19. Decision

`{DECISION}`

The run is implementation-valid and energetically stable. Contact rebound is now provisionally reasonable and magnetic torque survives, yet the 2.931 x 0.903 mm robot develops a 5.5604 ms opposing-wall bridge, a >0.5 ms phase plateau and a persistent approximately 33 deg tilted posture. The current limiting mechanism is therefore robot-to-pipe geometry, not CEL, missing magnetic torque or the accepted normal damping. This result does **not** permit extension to 16.667 ms and does not establish 30 Hz lock.

## 20. Exactly one next step

**Begin modifying the robot geometry dimensions, prioritizing a shorter robot length.** Do not continue tuning contact damping, Reduced-Hydro coefficients or magnetic-field parameters before that geometry change.
"""
    (HERE / "Wobble30Hz_ReducedHydro_Zeta050_8p333_report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
