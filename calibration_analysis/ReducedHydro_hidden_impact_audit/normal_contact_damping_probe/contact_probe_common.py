"""Shared contact-point kinematics for baseline and one damping candidate."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent
ROOT = AUDIT.parents[2]


def _get(z, variable):
    return max((z[k] for k in z.files if k == variable or k.startswith(variable + ' (Repeated:')), key=len)


def dense_rp(folder):
    z = np.load(folder / 'rp_history_private.npz')
    t0 = _get(z, 'V1')[:, 0]
    t = np.round(t0 / 1e-7) * 1e-7
    data = {}
    for prefix in ('U', 'UR', 'V', 'VR'):
        data[prefix] = np.column_stack([
            np.interp(t0, _get(z, prefix + str(i))[:, 0], _get(z, prefix + str(i))[:, 1])
            for i in (1, 2, 3)])
    return t, data


def contact_force(folder, t):
    z = np.load(folder / 'contact_history_private.npz')
    def force(prefix, surface):
        cols = []
        for i in (1, 2, 3):
            key = next(k for k in z.files if '|' + prefix + str(i) + ' on surface ' in k and surface in k)
            a = z[key]
            at = np.round(a[:, 0] / 1e-7) * 1e-7
            assert np.all(np.diff(at) > 0)
            cols.append(np.interp(t, at, a[:, 1]))
        return np.column_stack(cols)
    return force('CFN', 'ASSEMBLY_ROBOT'), force('CFS', 'ASSEMBLY_ROBOT')


def first_event(force, threshold=1e-8):
    active = np.linalg.norm(force, axis=1) > threshold
    edge = np.diff(np.r_[False, active, False].astype(int))
    starts = np.where(edge == 1)[0]
    ends = np.where(edge == -1)[0] - 1
    assert len(starts), 'No solver-active contact event found'
    return active, int(starts[0]), int(ends[0]), list(zip(starts.astype(int), ends.astype(int)))


def stable_windows(ia, ib, count=5):
    pre = np.arange(ia - count, ia)
    post = np.arange(ib + 1, ib + count + 1)
    assert len(pre) == count and len(post) == count
    return pre, post


def kinematics(t, data, mesh, rp, node_label, wall, gaps):
    node0 = np.asarray(mesh['nodes'][int(node_label)]) + mesh['shift']
    rotations = Rotation.from_rotvec(data['UR'])
    xcom = rp + data['U'] + rotations.apply(np.broadcast_to(mesh['com'] - rp, data['U'].shape))
    xcontact = rp + data['U'] + rotations.apply(np.broadcast_to(node0 - rp, data['U'].shape))
    rcontact = xcontact - xcom
    vcom = data['V'] + np.cross(data['VR'], rotations.apply(np.broadcast_to(mesh['com'] - rp, data['U'].shape)))
    vcontact = vcom + np.cross(data['VR'], rcontact)
    triangle = gaps.triangle_index.to_numpy(int)
    n = wall.normals[triangle]
    vn = np.einsum('ij,ij->i', vcontact, n)
    vt = vcontact - vn[:, None] * n
    return xcom, xcontact, vcom, vcontact, n, vn, vt


def summarize_case(case_name, t, data, mesh, rp, wall, gaps, folder, damping, tangent_fraction):
    fn, fs = contact_force(folder, t)
    force = fn + fs
    active, ia, ib, events = first_event(force)
    pre, post = stable_windows(ia, ib)
    node = int(gaps.iloc[ia].contact_node)
    xcom, xc, vcom, vc, normal, vn, vt = kinematics(t, data, mesh, rp, node, wall, gaps)
    gap_mm = gaps.gap_um.to_numpy(float) * 1e-3
    gap_pre = np.polyfit(t[pre] - t[pre].mean(), gap_mm[pre], 1)[0]
    gap_post = np.polyfit(t[post] - t[post].mean(), gap_mm[post], 1)[0]
    vn_in = float(vn[pre].mean())
    vn_out = float(vn[post].mean())
    en = vn_out / (-vn_in)
    egap = gap_post / (-gap_pre)
    power_n = np.einsum('ij,ij->i', fn, vn[:, None] * normal) * 1e-3
    power_t = np.einsum('ij,ij->i', fs, vt) * 1e-3
    wn = cumulative_trapezoid(power_n, t, initial=0.0)
    wt = cumulative_trapezoid(power_t, t, initial=0.0)
    a, b = ia - 1, ib + 1
    dv = vcom[b] - vcom[a]
    event_rows = []
    for start, end in events:
        impulse = np.trapz(force[start - 1:end + 2], t[start - 1:end + 2], axis=0)
        event_rows.append({
            'case': case_name, 'event_number': len(event_rows) + 1,
            'start_s': t[start], 'end_s': t[end], 'duration_us': (end - start + 1) * .1,
            'peak_force_N': np.linalg.norm(force[start:end + 1], axis=1).max(),
            'impulse_x_Ns': impulse[0], 'impulse_y_Ns': impulse[1], 'impulse_z_Ns': impulse[2],
            'impulse_norm_Ns': np.linalg.norm(impulse),
        })
    summary = {
        'case': case_name, 'critical_damping_fraction': damping,
        'tangent_fraction': tangent_fraction, 'contact_node': node,
        'first_contact_s': t[ia], 'first_contact_end_s': t[ib],
        'first_contact_duration_us': (ib - ia + 1) * .1,
        'first_peak_force_N': np.linalg.norm(force[ia:ib + 1], axis=1).max(),
        'vn_in_mm_s': vn_in, 'vn_out_mm_s': vn_out, 'e_n': en,
        'vt_in_mm_s': np.linalg.norm(vt[pre], axis=1).mean(),
        'vt_out_mm_s': np.linalg.norm(vt[post], axis=1).mean(),
        'gap_closing_slope_mm_s': gap_pre, 'gap_opening_slope_mm_s': gap_post,
        'gap_proxy_restitution': egap,
        'Wnormal_uJ': (wn[b] - wn[a]) * 1e6,
        'Wtangential_combined_uJ': (wt[b] - wt[a]) * 1e6,
        'work_decomposition_limit': ('CFS combines Coulomb friction and tangential contact damping; '
                                     'the two cannot be separated from available output'),
        'deltaV_COM_x_mm_s': dv[0], 'deltaV_COM_y_mm_s': dv[1], 'deltaV_COM_z_mm_s': dv[2],
        'deltaV_COM_norm_mm_s': np.linalg.norm(dv),
        'min_gap_um': gaps.gap_um.min(),
        'persistent_contact_over_0p10ms': any((end - start + 1) * 1e-7 > 1e-4 for start, end in events),
    }
    history = pd.DataFrame({
        'time_s': t, 'gap_um': gaps.gap_um, 'contact_node': gaps.contact_node,
        'wall_element': gaps.wall_element, 'triangle_index': gaps.triangle_index,
        'n_away_x': normal[:, 0], 'n_away_y': normal[:, 1], 'n_away_z': normal[:, 2],
        'Vcontact_x_mm_s': vc[:, 0], 'Vcontact_y_mm_s': vc[:, 1], 'Vcontact_z_mm_s': vc[:, 2],
        'vn_mm_s': vn, 'vt_norm_mm_s': np.linalg.norm(vt, axis=1),
        'Vcom_x_mm_s': vcom[:, 0], 'Vcom_y_mm_s': vcom[:, 1], 'Vcom_z_mm_s': vcom[:, 2],
        'Fnormal_x_N': fn[:, 0], 'Fnormal_y_N': fn[:, 1], 'Fnormal_z_N': fn[:, 2],
        'Fshear_x_N': fs[:, 0], 'Fshear_y_N': fs[:, 1], 'Fshear_z_N': fs[:, 2],
        'Fwall_norm_N': np.linalg.norm(force, axis=1), 'contact_active': active.astype(int),
        'Wnormal_cumulative_J': wn, 'Wtangential_combined_cumulative_J': wt,
    })
    return summary, history, pd.DataFrame(event_rows), (pre, post)
