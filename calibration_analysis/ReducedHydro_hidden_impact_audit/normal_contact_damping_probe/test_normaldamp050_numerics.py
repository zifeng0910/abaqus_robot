"""Regression checks for the single zeta=0.50 normal-only contact probe."""
from pathlib import Path
import json

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent


def main():
    ref = json.loads((HERE / 'normaldamp020_summary.json').read_text())
    cand = json.loads((HERE / 'normaldamp050_summary.json').read_text())
    momentum = pd.read_csv(HERE / 'normaldamp050_momentum_closure.csv').set_index('window')
    energy = pd.read_csv(HERE / 'normaldamp050_energy_summary.csv').iloc[0]
    events = pd.read_csv(HERE / 'normaldamp050_contact_event_catalog.csv')
    gap = pd.read_csv(HERE / 'normaldamp050_gap_history.csv')
    history = pd.read_csv(HERE / 'normaldamp050_first_impact_history.csv')
    comparison = pd.read_csv(HERE / 'normaldamp020_vs_050_summary.csv')
    action = pd.read_csv(HERE / 'normaldamp050_contact_action_reaction.csv').iloc[0]

    assert cand['decision'] == 'PROVISIONAL_NORMAL_CONTACT_DAMPING_ACCEPTABLE_FOR_WOBBLE_TEST'
    assert cand['eligible_for_8p333ms'] is True
    assert cand['run_completed_successfully'] and cand['dt_stable']
    assert cand['dt_s'] == 1e-7 and cand['warning_count'] == 13
    assert .30 <= cand['e_n'] <= .65 and cand['e_n'] < ref['e_n']
    assert np.isclose(cand['e_n'], cand['gap_proxy_restitution'], rtol=2e-5)
    assert cand['e_n_reduction_from_020_percent'] > 0
    assert momentum.loc['first_impact_audit_0p95_to_1p10ms', 'relative_residual'] < 1e-4
    assert momentum.loc['clean_0_to_0p95ms', 'relative_residual'] < 1e-4
    assert action.max_robot_plus_wall_CFT_error_N < 1e-6
    assert action.max_CFT_minus_CFN_plus_CFS_error_N < 1e-6
    assert energy.Wnormal_J <= 0 and energy.Wcontact_total_J <= 0
    assert abs(energy.contact_work_balance_difference_J) < 1e-12
    assert cand['min_first_event_gap_um'] > -1 and cand['min_gap_um'] > -1
    assert not cand['persistent_contact_over_0p10ms']
    assert cand['gap_reopened_after_first_event']
    assert cand['longest_contact_duration_us'] < 100
    assert len(events) == 2 and comparison.shape[0] == 2
    assert len(gap) == 13001 and np.isclose(gap.gap_um.min(), cand['min_gap_um'])
    assert set(['vn_mm_s', 'gap_um', 'Fwall_norm_N', 'Wnormal_cumulative_J']).issubset(history.columns)
    required = [
        'normaldamp020_vs_050_summary.csv', 'normaldamp050_first_impact_history.csv',
        'normaldamp050_contact_force_history.csv', 'normaldamp050_normal_velocity.csv',
        'normaldamp050_gap_history.csv', 'normaldamp050_energy_summary.csv',
        'normaldamp050_momentum_closure.csv', 'normaldamp050_contact_event_catalog.csv',
        'normaldamp020_vs_050_restitution.png', 'normaldamp020_vs_050_gap.png',
        'normaldamp020_vs_050_wall_force.png', 'normaldamp020_vs_050_energy.png',
        'normaldamp020_vs_050_COM_recoil.png',
        'ReducedHydro_FirstImpact_NormalDamp020_vs_050.gif',
        'ReducedHydro_normal_contact_damping_050_probe_report.md',
    ]
    assert all((HERE / name).is_file() and (HERE / name).stat().st_size > 0 for name in required)
    report = (HERE / 'ReducedHydro_normal_contact_damping_050_probe_report.md').read_text(encoding='utf-8')
    assert all(f'## {number}.' in report for number in range(1, 16))
    print('PASS: zeta=0.50 numerics, gates, traceability, and required artifacts')


if __name__ == '__main__':
    main()
