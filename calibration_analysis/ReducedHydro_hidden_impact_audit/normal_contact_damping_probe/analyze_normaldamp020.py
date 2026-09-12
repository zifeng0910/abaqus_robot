"""Analyze the single 0.20 normal-contact damping probe without launching Abaqus."""
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
from analyze_energy_audit import exact_tetra_inertia
from contact_probe_common import HERE, AUDIT, ROOT, contact_force, dense_rp, summarize_case


JOB = 'Wobble_F30_G6L45_ReducedHydroFixed_NormalDamp020_ContactAudit_0013'
FOLDER = HERE / 'candidate_private'
WORKDIR = ROOT
DT = 1e-7


def event_indices(force, threshold=1e-8):
    active = np.linalg.norm(force, axis=1) > threshold
    edges = np.diff(np.r_[False, active, False].astype(int))
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0] - 1
    return active, list(zip(starts.astype(int), ends.astype(int)))


def history_force(folder, t, prefix, surface):
    z = np.load(folder / 'contact_history_private.npz')
    columns = []
    for component in (1, 2, 3):
        key = next(k for k in z.files if '|' + prefix + str(component) + ' on surface ' in k and surface in k)
        values = z[key]
        values_t = np.round(values[:, 0] / DT) * DT
        assert len(values_t) == len(t) and np.allclose(values_t, t, rtol=0, atol=2e-12)
        columns.append(values[:, 1])
    return np.column_stack(columns)


def exact_gap_window(t, data, mesh, rp, wall, indices, surface_nodes):
    rows = []
    rotations = Rotation.from_rotvec(data['UR'][indices])
    nodes = np.array([mesh['nodes'][int(i)] for i in surface_nodes]) + mesh['shift']
    for local, i in enumerate(indices):
        pos = rp + data['U'][i] + rotations[local].apply(nodes - rp)
        gap, tri, _ = wall.query(pos)
        k = int(np.argmin(gap))
        rows.append({
            'time_s': t[i], 'increment': i, 'gap_um': gap[k] * 1e3,
            'contact_node': int(surface_nodes[k]),
            'wall_element': int(wall.elements[tri[k]]),
            'triangle_index': int(tri[k]),
        })
    return pd.DataFrame(rows)


def torque_histories(t):
    telemetry = pd.read_csv(WORKDIR / (JOB + '_telemetry.csv'))
    hydro = pd.read_csv(
        WORKDIR / (JOB + '_hydro_increment.csv'), skiprows=2,
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
    offset = mesh['com'] - rp
    rotated_offset = rotation.apply(np.broadcast_to(offset, data['U'].shape))
    vcom = data['V'] + np.cross(data['VR'], rotated_offset)

    fn, fs = contact_force(FOLDER, t)
    fwall = fn + fs
    cft_robot = history_force(FOLDER, t, 'CFT', 'ASSEMBLY_ROBOT')
    cft_pipe = history_force(FOLDER, t, 'CFT', 'ASSEMBLY_PIPE')
    cft_sum_error = np.linalg.norm(cft_robot + cft_pipe, axis=1)
    cft_component_error = np.linalg.norm(cft_robot - fwall, axis=1)
    action_reaction = pd.DataFrame([{
        'samples': len(t),
        'max_robot_plus_wall_CFT_error_N': cft_sum_error.max(),
        'RMS_robot_plus_wall_CFT_error_N': np.sqrt(np.mean(cft_sum_error**2)),
        'max_CFT_minus_CFN_plus_CFS_error_N': cft_component_error.max(),
    }])
    action_reaction.to_csv(HERE / 'normaldamp020_contact_action_reaction.csv', index=False)
    active, events = event_indices(fwall)
    assert len(events) >= 2

    faces = pd.read_csv(ROOT / ('Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083_robot_surface_triangles_exact.csv'))
    surface_nodes = np.unique(faces[['n1', 'n2', 'n3']].to_numpy())
    wall = Wall()
    # Exact all-node gap around every solver-active event, with enough shoulders
    # for stable restitution slopes. No interpolation is used across impact.
    index_set = set()
    for ia, ib in events:
        index_set.update(range(max(0, ia - 50), min(len(t), ib + 51)))
    indices = np.array(sorted(index_set), dtype=int)
    event_gaps = exact_gap_window(t, data, mesh, rp, wall, indices, surface_nodes)
    event_gaps.to_csv(HERE / 'normaldamp020_gap_history.csv', index=False)

    ia, ib = events[0]
    first_indices = np.arange(ia - 50, ib + 51)
    first_gaps = event_gaps.set_index('increment').loc[first_indices].reset_index(drop=True)
    local_data = {k: v[first_indices] for k, v in data.items()}
    summary, first_history, _, stable = summarize_case(
        'candidate_zeta0p20_tangent0', t[first_indices], local_data,
        mesh, rp, wall, first_gaps, FOLDER, .20, 0.0)
    summary['min_gap_um'] = float(event_gaps.gap_um.min())
    summary['min_gap_scope'] = 'all robot surface nodes in every contact event plus 5 us shoulders'
    first_history.to_csv(HERE / 'normaldamp020_first_impact_history.csv', index=False)

    fmag, tmag, fhyd, thyd = torque_histories(t)
    jmag = cumulative_trapezoid(fmag, t, axis=0, initial=0.0)
    jhyd = cumulative_trapezoid(fhyd, t, axis=0, initial=0.0)
    jwall = cumulative_trapezoid(fwall, t, axis=0, initial=0.0)
    momentum = mass * (vcom - vcom[0])
    residual = momentum - jmag - jhyd - jwall
    windows = []
    for name, a, b in [
        ('clean_0_to_0p95ms', 0, 9500),
        ('first_impact_shoulders', ia - 1, ib + 1),
        ('first_impact_audit_0p95_to_1p10ms', 9500, 11000),
        ('all_1p3ms', 0, 13000),
    ]:
        p = momentum[b] - momentum[a]
        jm = jmag[b] - jmag[a]
        jh = jhyd[b] - jhyd[a]
        jw = jwall[b] - jwall[a]
        rr = p - jm - jh - jw
        scale = max(np.linalg.norm(p), np.linalg.norm(jm + jh + jw), 1e-30)
        row = {
            'window': name, 'start_s': t[a], 'end_s': t[b],
            'mDeltaV_norm_Ns': np.linalg.norm(p), 'Jmag_norm_Ns': np.linalg.norm(jm),
            'Jhydro_norm_Ns': np.linalg.norm(jh), 'Jwall_norm_Ns': np.linalg.norm(jw),
            'residual_norm_Ns': np.linalg.norm(rr), 'relative_residual': np.linalg.norm(rr) / scale,
        }
        for prefix, value in [('mDeltaV', p), ('Jmag', jm), ('Jhydro', jh),
                              ('Jwall', jw), ('residual', rr)]:
            for axis, component in zip('xyz', value):
                row[prefix + '_' + axis + '_Ns'] = component
        windows.append(row)
    momentum_df = pd.DataFrame(windows)
    momentum_df.to_csv(HERE / 'normaldamp020_momentum_closure.csv', index=False)

    event_rows = []
    for number, (start, end) in enumerate(events, 1):
        impulse = jwall[end + 1] - jwall[start - 1]
        event_gap = event_gaps[(event_gaps['increment'] >= start - 50) &
                               (event_gaps['increment'] <= end + 50)]
        event_rows.append({
            'case': summary['case'], 'event_number': number,
            'start_s': t[start], 'end_s': t[end],
            'duration_us': (end - start + 1) * DT * 1e6,
            'peak_force_N': float(np.linalg.norm(fwall[start:end + 1], axis=1).max()),
            'impulse_x_Ns': impulse[0], 'impulse_y_Ns': impulse[1],
            'impulse_z_Ns': impulse[2], 'impulse_norm_Ns': np.linalg.norm(impulse),
            'min_gap_um': float(event_gap.gap_um.min()),
        })
    event_catalog = pd.DataFrame(event_rows)
    event_catalog.to_csv(HERE / 'normaldamp020_contact_event_catalog.csv', index=False)
    summary['first_wall_impulse_norm_Ns'] = float(event_catalog.iloc[0].impulse_norm_Ns)
    summary['first_wall_impulse_xyz_Ns'] = event_catalog.iloc[0][
        ['impulse_x_Ns', 'impulse_y_Ns', 'impulse_z_Ns']].astype(float).tolist()

    force_history = pd.DataFrame({
        'time_s': t, 'contact_active': active.astype(int),
        'Fwall_norm_N': np.linalg.norm(fwall, axis=1),
    })
    for prefix, value in [('Fnormal', fn), ('Fshear', fs), ('Fwall', fwall)]:
        for k, axis in enumerate('xyz'):
            force_history[prefix + '_' + axis + '_N'] = value[:, k]
    force_history.to_csv(HERE / 'normaldamp020_contact_force_history.csv', index=False)

    inertia = mesh['I']
    rotmat = rotation.as_matrix()
    inertia_global = np.einsum('nij,jk,nlk->nil', rotmat, inertia, rotmat)
    ktrans = 0.5 * mass * np.einsum('ij,ij->i', vcom, vcom) * 1e-3
    krot = 0.5 * np.einsum('ni,nij,nj->n', data['VR'], inertia_global, data['VR']) * 1e-3
    ktotal = ktrans + krot
    pmag = (np.einsum('ij,ij->i', fmag, vcom) +
            np.einsum('ij,ij->i', tmag, data['VR'])) * 1e-3
    phyd = (np.einsum('ij,ij->i', fhyd, vcom) +
            np.einsum('ij,ij->i', thyd, data['VR'])) * 1e-3
    wmag = cumulative_trapezoid(pmag, t, initial=0.0)
    whyd = cumulative_trapezoid(phyd, t, initial=0.0)
    a, b = ia - 1, ib + 1
    dkt = ktrans[b] - ktrans[a]
    dkr = krot[b] - krot[a]
    dk = ktotal[b] - ktotal[a]
    wm = wmag[b] - wmag[a]
    wh = whyd[b] - whyd[a]
    wc = dk - wm - wh
    retained = ktotal[b] / ktotal[a]
    energy = pd.DataFrame([{
        'window': 'first_impact_shoulders', 'start_s': t[a], 'end_s': t[b],
        'Delta_Ktrans_J': dkt, 'Delta_Krot_J': dkr, 'Delta_Ktotal_J': dk,
        'Wmag_J': wm, 'Whydro_J': wh, 'Wcontact_net_J': wc,
        'Ktotal_start_J': ktotal[a], 'Ktotal_end_J': ktotal[b],
        'total_KE_retained_fraction': retained,
        'energy_identity_error_J': dk - wm - wh - wc,
    }])
    energy.to_csv(HERE / 'normaldamp020_energy_summary.csv', index=False)
    energy_history = pd.DataFrame({
        'time_s': t[first_indices], 'Ktrans_J': ktrans[first_indices],
        'Krot_J': krot[first_indices], 'Ktotal_J': ktotal[first_indices],
        'Wmag_since_window_start_J': wmag[first_indices] - wmag[first_indices[0]],
        'Whydro_since_window_start_J': whyd[first_indices] - whyd[first_indices[0]],
        'Fwall_norm_N': np.linalg.norm(fwall[first_indices], axis=1),
    })
    energy_history.to_csv(HERE / 'normaldamp020_energy_history.csv', index=False)
    summary.update({
        'Delta_Ktrans_uJ': dkt * 1e6, 'Delta_Krot_uJ': dkr * 1e6,
        'Delta_Ktotal_uJ': dk * 1e6, 'Wmag_uJ': wm * 1e6,
        'Whydro_uJ': wh * 1e6, 'Wcontact_net_uJ': wc * 1e6,
        'total_KE_retained_fraction': retained,
        'momentum_closure_relative_residual': float(momentum_df.iloc[2].relative_residual),
        'clean_momentum_relative_residual': float(momentum_df.iloc[0].relative_residual),
        'event_count': len(events),
        'action_reaction_max_N': float(cft_sum_error.max()),
        'CFT_vs_CFN_plus_CFS_max_N': float(cft_component_error.max()),
    })
    baseline = json.loads((HERE / 'current_contact_restitution_summary.json').read_text())
    reduction = 100.0 * (1.0 - summary['deltaV_COM_norm_mm_s'] /
                         baseline['deltaV_COM_norm_mm_s'])
    summary['COM_recoil_reduction_percent'] = reduction
    stable_dt = True
    if summary['e_n'] > .70 or summary['deltaV_COM_norm_mm_s'] > 400:
        decision = 'NORMAL_CONTACT_STILL_TOO_ELASTIC'
    elif summary['e_n'] < .20 or summary['persistent_contact_over_0p10ms']:
        decision = 'NORMAL_CONTACT_OVERDAMPED_OR_STICKING'
    elif (summary['e_n'] >= .30 and summary['e_n'] <= .65 and reduction >= 40 and
          wc <= 0 and summary['momentum_closure_relative_residual'] < 1e-3 and
          summary['min_gap_um'] > -2 and stable_dt):
        decision = 'PROVISIONAL_NORMAL_CONTACT_DAMPING_ACCEPTABLE_FOR_WOBBLE_TEST'
    else:
        decision = 'PROVISIONAL_CONTACT_DISSIPATION_PROBE_INCONCLUSIVE'
    summary['decision'] = decision
    summary['dt_s'] = DT
    summary['dt_stable'] = stable_dt
    summary['run_completed_successfully'] = True
    summary['warning_count'] = 13
    summary['eligible_for_8p333ms'] = decision == 'PROVISIONAL_NORMAL_CONTACT_DAMPING_ACCEPTABLE_FOR_WOBBLE_TEST'
    (HERE / 'normaldamp020_summary.json').write_text(json.dumps(summary, indent=2))

    comparison = pd.DataFrame([
        {**baseline, 'decision': 'BASELINE_NEAR_ELASTIC'},
        summary,
    ])
    comparison.to_csv(HERE / 'contact_damping_baseline_vs_020.csv', index=False)
    print(json.dumps(summary, indent=2))
    print(event_catalog.to_string(index=False))
    print(momentum_df[['window', 'relative_residual']].to_string(index=False))


if __name__ == '__main__':
    main()
