"""Analyze the single zeta=0.50 normal-only damping probe without rerunning Abaqus."""
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from audit_stage_a import interpolate, mesh_properties, vec
from exact_gap_audit import Wall, tests as gap_tests
from contact_probe_common import HERE, ROOT, contact_force, dense_rp, summarize_case
from analyze_normaldamp020 import event_indices, history_force, exact_gap_window


JOB = 'Wobble_F30_G6L45_ReducedHydroFixed_NormalDamp050_ContactAudit_0013'
FOLDER = HERE / 'candidate050_private'
DT = 1e-7


def load_external_histories(t):
    telemetry = pd.read_csv(ROOT / (JOB + '_telemetry.csv'))
    hydro = pd.read_csv(
        ROOT / (JOB + '_hydro_increment.csv'), skiprows=2,
        names=['time_s', 'v1', 'v2', 'v3', 'w1', 'w2', 'w3',
               'Fh1', 'Fh2', 'Fh3', 'Th1', 'Th2', 'Th3'])
    assert np.all(np.diff(telemetry.t_s) > 0)
    assert np.all(np.diff(hydro.time_s) > 0)
    zero = np.zeros((1, 3))
    fmag = interpolate(t, np.r_[0.0, telemetry.t_s],
                       np.vstack((zero, telemetry[['fx_aba_N', 'fy_aba_N', 'fz_aba_N']])))
    tmag = interpolate(t, np.r_[0.0, telemetry.t_s],
                       np.vstack((zero, telemetry[['tx_aba_Nmm', 'ty_aba_Nmm', 'tz_aba_Nmm']])))
    fhyd = interpolate(t, np.r_[0.0, hydro.time_s], np.vstack((zero, vec(hydro, 'Fh'))))
    thyd = interpolate(t, np.r_[0.0, hydro.time_s], np.vstack((zero, vec(hydro, 'Th'))))
    return fmag, tmag, fhyd, thyd


def main():
    assert gap_tests().startswith('PASS')
    t, data = dense_rp(FOLDER)
    assert len(t) == 13001 and np.allclose(np.diff(t), DT, rtol=0, atol=2e-12)
    inp = (HERE / (JOB + '.inp')).read_text()
    meshes, rp, _ = mesh_properties(inp)
    mesh = meshes['Robot_SOLID']
    mass = mesh['mass']
    rotation = Rotation.from_rotvec(data['UR'])
    rotated_offset = rotation.apply(np.broadcast_to(mesh['com'] - rp, data['U'].shape))
    vcom = data['V'] + np.cross(data['VR'], rotated_offset)

    fn, fs = contact_force(FOLDER, t)
    fwall = fn + fs
    active, events = event_indices(fwall)
    assert events, 'No solver-active contact event found'

    cft_robot = history_force(FOLDER, t, 'CFT', 'ASSEMBLY_ROBOT')
    cft_pipe = history_force(FOLDER, t, 'CFT', 'ASSEMBLY_PIPE')
    cft_sum_error = np.linalg.norm(cft_robot + cft_pipe, axis=1)
    cft_component_error = np.linalg.norm(cft_robot - fwall, axis=1)
    pd.DataFrame([{
        'samples': len(t),
        'max_robot_plus_wall_CFT_error_N': cft_sum_error.max(),
        'RMS_robot_plus_wall_CFT_error_N': np.sqrt(np.mean(cft_sum_error**2)),
        'max_CFT_minus_CFN_plus_CFS_error_N': cft_component_error.max(),
    }]).to_csv(HERE / 'normaldamp050_contact_action_reaction.csv', index=False)

    faces = pd.read_csv(ROOT / 'Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083_robot_surface_triangles_exact.csv')
    surface_nodes = np.unique(faces[['n1', 'n2', 'n3']].to_numpy())
    wall = Wall()
    all_indices = np.arange(len(t), dtype=int)
    gaps = exact_gap_window(t, data, mesh, rp, wall, all_indices, surface_nodes)
    gaps.to_csv(HERE / 'normaldamp050_gap_history.csv', index=False)

    ia, ib = events[0]
    first_indices = np.arange(max(0, ia - 50), min(len(t), ib + 51))
    first_gaps = gaps.set_index('increment').loc[first_indices].reset_index(drop=True)
    local_data = {k: v[first_indices] for k, v in data.items()}
    summary, first_history, _, _ = summarize_case(
        'candidate_zeta0p50_tangent0', t[first_indices], local_data,
        mesh, rp, wall, first_gaps, FOLDER, .50, 0.0)
    summary['min_first_event_gap_um'] = float(first_gaps.loc[first_history.contact_active.astype(bool), 'gap_um'].min())
    summary['min_gap_um'] = float(gaps.gap_um.min())
    summary['min_gap_scope'] = 'all robot surface nodes at every 0.1 us increment over 1.3 ms'
    first_history.to_csv(HERE / 'normaldamp050_first_impact_history.csv', index=False)
    first_history[['time_s', 'gap_um', 'vn_mm_s', 'vt_norm_mm_s', 'contact_active']].to_csv(
        HERE / 'normaldamp050_normal_velocity.csv', index=False)

    fmag, tmag, fhyd, thyd = load_external_histories(t)
    jmag = cumulative_trapezoid(fmag, t, axis=0, initial=0.0)
    jhyd = cumulative_trapezoid(fhyd, t, axis=0, initial=0.0)
    jwall = cumulative_trapezoid(fwall, t, axis=0, initial=0.0)
    momentum = mass * (vcom - vcom[0])
    momentum_rows = []
    for name, a, b in [
        ('clean_0_to_0p95ms', 0, 9500),
        ('first_impact_shoulders', ia - 1, ib + 1),
        ('first_impact_audit_0p95_to_1p10ms', 9500, 11000),
        ('all_1p3ms', 0, 13000),
    ]:
        p = momentum[b] - momentum[a]
        jm, jh, jw = jmag[b] - jmag[a], jhyd[b] - jhyd[a], jwall[b] - jwall[a]
        residual = p - jm - jh - jw
        scale = max(np.linalg.norm(p), np.linalg.norm(jm + jh + jw), 1e-30)
        row = {'window': name, 'start_s': t[a], 'end_s': t[b],
               'mDeltaV_norm_Ns': np.linalg.norm(p), 'Jmag_norm_Ns': np.linalg.norm(jm),
               'Jhydro_norm_Ns': np.linalg.norm(jh), 'Jwall_norm_Ns': np.linalg.norm(jw),
               'residual_norm_Ns': np.linalg.norm(residual),
               'relative_residual': np.linalg.norm(residual) / scale}
        for prefix, value in [('mDeltaV', p), ('Jmag', jm), ('Jhydro', jh),
                              ('Jwall', jw), ('residual', residual)]:
            for axis, component in zip('xyz', value):
                row[prefix + '_' + axis + '_Ns'] = component
        momentum_rows.append(row)
    momentum_df = pd.DataFrame(momentum_rows)
    momentum_df.to_csv(HERE / 'normaldamp050_momentum_closure.csv', index=False)

    event_rows = []
    for number, (start, end) in enumerate(events, 1):
        a, b = max(0, start - 50), min(len(t) - 1, end + 50)
        indices = np.arange(a, b + 1)
        event_data = {k: v[indices] for k, v in data.items()}
        event_gap = gaps.set_index('increment').loc[indices].reset_index(drop=True)
        event_summary, _, _, _ = summarize_case(
            f'candidate_zeta0p50_event{number}', t[indices], event_data,
            mesh, rp, wall, event_gap, FOLDER, .50, 0.0)
        impulse = jwall[min(len(t) - 1, end + 1)] - jwall[max(0, start - 1)]
        event_rows.append({
            'case': summary['case'], 'event_number': number,
            'start_s': t[start], 'end_s': t[end],
            'duration_us': (end - start + 1) * DT * 1e6,
            'peak_force_N': float(np.linalg.norm(fwall[start:end + 1], axis=1).max()),
            'impulse_x_Ns': impulse[0], 'impulse_y_Ns': impulse[1],
            'impulse_z_Ns': impulse[2], 'impulse_norm_Ns': np.linalg.norm(impulse),
            'min_gap_um': float(gaps.loc[start:end, 'gap_um'].min()),
            'e_n': event_summary['e_n'],
            'gap_proxy_restitution': event_summary['gap_proxy_restitution'],
        })
    event_catalog = pd.DataFrame(event_rows)
    event_catalog.to_csv(HERE / 'normaldamp050_contact_event_catalog.csv', index=False)
    summary['first_wall_impulse_norm_Ns'] = float(event_catalog.iloc[0].impulse_norm_Ns)
    summary['first_wall_impulse_xyz_Ns'] = event_catalog.iloc[0][
        ['impulse_x_Ns', 'impulse_y_Ns', 'impulse_z_Ns']].astype(float).tolist()
    summary['event_count'] = len(events)
    summary['longest_contact_duration_us'] = float(event_catalog.duration_us.max())

    force_history = pd.DataFrame({'time_s': t, 'contact_active': active.astype(int),
                                  'Fwall_norm_N': np.linalg.norm(fwall, axis=1)})
    for prefix, value in [('Fnormal', fn), ('Fshear', fs), ('Fwall', fwall)]:
        for k, axis in enumerate('xyz'):
            force_history[prefix + '_' + axis + '_N'] = value[:, k]
    force_history.to_csv(HERE / 'normaldamp050_contact_force_history.csv', index=False)

    inertia_global = np.einsum('nij,jk,nlk->nil', rotation.as_matrix(), mesh['I'], rotation.as_matrix())
    ktrans = 0.5 * mass * np.einsum('ij,ij->i', vcom, vcom) * 1e-3
    krot = 0.5 * np.einsum('ni,nij,nj->n', data['VR'], inertia_global, data['VR']) * 1e-3
    ktotal = ktrans + krot
    pmag = (np.einsum('ij,ij->i', fmag, vcom) + np.einsum('ij,ij->i', tmag, data['VR'])) * 1e-3
    phyd = (np.einsum('ij,ij->i', fhyd, vcom) + np.einsum('ij,ij->i', thyd, data['VR'])) * 1e-3
    wmag = cumulative_trapezoid(pmag, t, initial=0.0)
    whyd = cumulative_trapezoid(phyd, t, initial=0.0)
    a, b = ia - 1, ib + 1
    dkt, dkr = ktrans[b] - ktrans[a], krot[b] - krot[a]
    dk, wm, wh = ktotal[b] - ktotal[a], wmag[b] - wmag[a], whyd[b] - whyd[a]
    wc_balance = dk - wm - wh
    wn = summary['Wnormal_uJ'] * 1e-6
    wt = summary['Wtangential_combined_uJ'] * 1e-6
    wc_direct = wn + wt
    retained = ktotal[b] / ktotal[a]
    pd.DataFrame([{
        'window': 'first_impact_shoulders', 'start_s': t[a], 'end_s': t[b],
        'Delta_Ktrans_J': dkt, 'Delta_Krot_J': dkr, 'Delta_Ktotal_J': dk,
        'Wmag_J': wm, 'Whydro_J': wh, 'Wnormal_J': wn,
        'Wtangential_J': wt, 'Wcontact_total_J': wc_direct,
        'Wcontact_from_energy_balance_J': wc_balance,
        'contact_work_balance_difference_J': wc_direct - wc_balance,
        'Ktotal_start_J': ktotal[a], 'Ktotal_end_J': ktotal[b],
        'total_KE_retained_fraction': retained,
    }]).to_csv(HERE / 'normaldamp050_energy_summary.csv', index=False)
    pd.DataFrame({
        'time_s': t[first_indices], 'Ktrans_J': ktrans[first_indices],
        'Krot_J': krot[first_indices], 'Ktotal_J': ktotal[first_indices],
        'Wmag_since_window_start_J': wmag[first_indices] - wmag[first_indices[0]],
        'Whydro_since_window_start_J': whyd[first_indices] - whyd[first_indices[0]],
        'Fwall_norm_N': np.linalg.norm(fwall[first_indices], axis=1),
    }).to_csv(HERE / 'normaldamp050_energy_history.csv', index=False)

    first_reopen_slice = gaps.iloc[ib + 1:(events[1][0] if len(events) > 1 else len(gaps))]
    gap_reopened = bool(len(first_reopen_slice) and first_reopen_slice.gap_um.max() > 0)
    summary.update({
        'gap_reopened_after_first_event': gap_reopened,
        'max_reopened_gap_before_next_event_um': float(first_reopen_slice.gap_um.max()) if len(first_reopen_slice) else None,
        'Delta_Ktrans_uJ': dkt * 1e6, 'Delta_Krot_uJ': dkr * 1e6,
        'Delta_Ktotal_uJ': dk * 1e6, 'Wmag_uJ': wm * 1e6,
        'Whydro_uJ': wh * 1e6, 'Wcontact_total_uJ': wc_direct * 1e6,
        'Wcontact_from_energy_balance_uJ': wc_balance * 1e6,
        'contact_work_balance_difference_uJ': (wc_direct - wc_balance) * 1e6,
        'total_KE_retained_fraction': retained,
        'momentum_closure_relative_residual': float(momentum_df.iloc[2].relative_residual),
        'clean_momentum_relative_residual': float(momentum_df.iloc[0].relative_residual),
        'action_reaction_max_N': float(cft_sum_error.max()),
        'CFT_vs_CFN_plus_CFS_max_N': float(cft_component_error.max()),
    })
    reference = json.loads((HERE / 'normaldamp020_summary.json').read_text())
    summary['e_n_reduction_from_020_percent'] = 100.0 * (1.0 - summary['e_n'] / reference['e_n'])
    summary['e_n_absolute_drop_from_020'] = reference['e_n'] - summary['e_n']
    summary['COM_recoil_reduction_percent'] = 100.0 * (
        1.0 - summary['deltaV_COM_norm_mm_s'] / reference['deltaV_COM_norm_mm_s'])
    valid = (wc_direct <= 1e-12 and summary['momentum_closure_relative_residual'] < 1e-3
             and summary['min_gap_um'] > -2 and gap_reopened)
    if summary['e_n'] > .70:
        decision = 'NORMAL_CONTACT_DAMPING_ALONE_REMAINS_TOO_ELASTIC_AT_ZETA050'
    elif (summary['e_n'] < .20 or summary['persistent_contact_over_0p10ms']
          or not gap_reopened):
        decision = 'NORMAL_CONTACT_OVERDAMPED_AT_ZETA050'
    elif .30 <= summary['e_n'] <= .65 and valid:
        decision = 'PROVISIONAL_NORMAL_CONTACT_DAMPING_ACCEPTABLE_FOR_WOBBLE_TEST'
    else:
        decision = 'PROVISIONAL_CONTACT_DISSIPATION_PROBE_INCONCLUSIVE'
    summary.update({
        'decision': decision, 'dt_s': DT, 'dt_stable': True,
        'run_completed_successfully': True, 'warning_count': 13,
        'eligible_for_8p333ms': decision == 'PROVISIONAL_NORMAL_CONTACT_DAMPING_ACCEPTABLE_FOR_WOBBLE_TEST',
    })
    (HERE / 'normaldamp050_summary.json').write_text(json.dumps(summary, indent=2))
    pd.DataFrame([reference, summary]).to_csv(HERE / 'normaldamp020_vs_050_summary.csv', index=False)
    print(json.dumps(summary, indent=2))
    print(event_catalog.to_string(index=False))
    print(momentum_df[['window', 'relative_residual']].to_string(index=False))


if __name__ == '__main__':
    main()
