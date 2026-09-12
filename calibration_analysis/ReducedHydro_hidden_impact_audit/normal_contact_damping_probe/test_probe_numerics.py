"""Regression checks for the one-run normal-contact damping probe."""
from pathlib import Path
import json

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent


def main():
    base = json.loads((HERE / 'current_contact_restitution_summary.json').read_text())
    cand = json.loads((HERE / 'normaldamp020_summary.json').read_text())
    momentum = pd.read_csv(HERE / 'normaldamp020_momentum_closure.csv').set_index('window')
    energy = pd.read_csv(HERE / 'normaldamp020_energy_summary.csv').iloc[0]
    events = pd.read_csv(HERE / 'normaldamp020_contact_event_catalog.csv')
    gap = pd.read_csv(HERE / 'normaldamp020_gap_history.csv')
    history = pd.read_csv(HERE / 'normaldamp020_first_impact_history.csv')
    action_reaction = pd.read_csv(HERE / 'normaldamp020_contact_action_reaction.csv').iloc[0]

    assert cand['decision'] == 'NORMAL_CONTACT_STILL_TOO_ELASTIC'
    assert cand['eligible_for_8p333ms'] is False
    assert cand['run_completed_successfully'] and cand['dt_stable']
    assert cand['dt_s'] == 1e-7 and cand['warning_count'] == 13
    assert cand['e_n'] > .70 and cand['deltaV_COM_norm_mm_s'] > 400
    assert 0 < cand['COM_recoil_reduction_percent'] < 40
    assert cand['e_n'] < base['e_n']
    assert np.isclose(cand['e_n'], cand['gap_proxy_restitution'], rtol=2e-5)
    assert momentum.loc['first_impact_audit_0p95_to_1p10ms', 'relative_residual'] < 1e-4
    assert momentum.loc['clean_0_to_0p95ms', 'relative_residual'] < 1e-4
    assert action_reaction.max_robot_plus_wall_CFT_error_N < 1e-6
    assert action_reaction.max_CFT_minus_CFN_plus_CFS_error_N < 1e-6
    assert energy.Wcontact_net_J <= 0
    assert abs(energy.energy_identity_error_J) < 1e-18
    assert cand['min_gap_um'] > -1
    assert not cand['persistent_contact_over_0p10ms']
    assert len(events) == 2 and events.duration_us.max() < 100
    assert gap.gap_um.min() == cand['min_gap_um']
    assert set(['vn_mm_s', 'gap_um', 'Fwall_norm_N', 'Wnormal_cumulative_J']).issubset(history.columns)
    required = [
        'contact_damping_baseline_vs_020.csv', 'normaldamp020_contact_force_history.csv',
        'normaldamp020_momentum_closure.csv', 'normaldamp020_energy_summary.csv',
        'ReducedHydro_normal_contact_damping_probe_report.md',
        'ReducedHydro_FirstImpact_Baseline_vs_NormalDamp020.gif',
    ]
    assert all((HERE / name).is_file() and (HERE / name).stat().st_size > 0 for name in required)
    print('PASS: candidate numerics, decision gates, and required artifacts')


if __name__ == '__main__':
    main()
