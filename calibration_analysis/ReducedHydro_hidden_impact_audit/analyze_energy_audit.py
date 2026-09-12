"""Read-only first-impact rigid-body energy and external-work audit.

Units: mm, tonne, s, N and N mm. Kinetic energy is converted from N mm to J.
No solver is launched and no force/work term is fitted to close the balance.
"""
from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.spatial.transform import Rotation

from audit_stage_a import HERE, ROOT, blocks, dense_history, interpolate, mesh_properties, vec


JOB = 'Wobble_F30_G6L45_ReducedHydroFixed_ContactAudit_0012'
WINDOW = (0.00098, 0.00108)


def exact_tetra_inertia(mesh):
    """Continuum inertia of constant-density linear tetrahedra about mesh COM."""
    rho = mesh['mass'] / mesh['vol'].sum()
    second = np.zeros((3, 3))
    for verts, volume in zip(mesh['verts'], mesh['vol']):
        # Integral(x x^T)dV = V/20 * (sum(v)sum(v)^T + sum(v v^T)).
        s = verts.sum(axis=0)
        second += volume / 20.0 * (np.outer(s, s) + verts.T @ verts)
    com = mesh['com']
    second_com = rho * second - mesh['mass'] * np.outer(com, com)
    return np.eye(3) * np.trace(second_com) - second_com


def load_torque(t):
    tel = pd.read_csv(ROOT / (JOB + '_telemetry.csv'))
    hyd = pd.read_csv(
        ROOT / (JOB + '_hydro_increment.csv'), skiprows=2,
        names=['time_s', 'v1', 'v2', 'v3', 'w1', 'w2', 'w3',
               'Fh1', 'Fh2', 'Fh3', 'Th1', 'Th2', 'Th3'])
    assert np.all(np.diff(tel.t_s) > 0) and np.all(np.diff(hyd.time_s) > 0)
    zero3 = np.zeros((1, 3))
    tmag = interpolate(
        t, np.r_[0.0, tel.t_s],
        np.vstack((zero3, tel[['tx_aba_Nmm', 'ty_aba_Nmm', 'tz_aba_Nmm']])))
    thyd = interpolate(t, np.r_[0.0, hyd.time_s], np.vstack((zero3, vec(hyd, 'Th'))))
    return tmag, thyd


def delta(a, i0, i1):
    return float(a[i1] - a[i0])


def main():
    text = (HERE / (JOB + '.inp')).read_text()
    meshes, rp, _ = mesh_properties(text)
    mesh = meshes['Robot_SOLID']
    mass = mesh['mass']
    inertia_lumped = mesh['I']
    inertia_exact = exact_tetra_inertia(mesh)

    t0, data, _ = dense_history(HERE / 'diagnostic')
    t = np.round(t0 / 1e-7) * 1e-7
    assert np.all(np.diff(t) > 0)
    d = pd.DataFrame(data)
    ur = vec(d, 'UR')
    omega = vec(d, 'VR')
    velocity = vec(d, 'V')
    offset = mesh['com'] - rp
    rot = Rotation.from_rotvec(ur).as_matrix()
    rotated_offset = np.einsum('nij,j->ni', rot, offset)
    vcom = velocity + np.cross(omega, rotated_offset)

    itrans = 0.5 * mass * np.einsum('ni,ni->n', vcom, vcom) * 1e-3

    def rotational_ke(inertia):
        ig = np.einsum('nij,jk,nlk->nil', rot, inertia, rot)
        return 0.5 * np.einsum('ni,nij,nj->n', omega, ig, omega) * 1e-3

    krot_lumped = rotational_ke(inertia_lumped)
    krot_exact = rotational_ke(inertia_exact)
    ktot_lumped = itrans + krot_lumped
    ktot_exact = itrans + krot_exact

    dense = np.load(HERE / 'diagnostic_dense_private.npz')
    assert np.allclose(t, dense['t'], rtol=0, atol=1e-12)
    fmag, fhyd = dense['F'], dense['H']
    tmag, thyd = load_torque(t)
    pmag_trans = np.einsum('ni,ni->n', fmag, vcom) * 1e-3
    phyd_trans = np.einsum('ni,ni->n', fhyd, vcom) * 1e-3
    pmag_rot = np.einsum('ni,ni->n', tmag, omega) * 1e-3
    phyd_rot = np.einsum('ni,ni->n', thyd, omega) * 1e-3
    work = {}
    for name, power in [('mag_trans', pmag_trans), ('mag_rot', pmag_rot),
                        ('hydro_trans', phyd_trans), ('hydro_rot', phyd_rot)]:
        work[name] = cumulative_trapezoid(power, t, initial=0.0)
    applied_work = sum(work.values())
    contact_net_work = (ktot_lumped - ktot_lumped[0]) - applied_work

    ah = np.load(HERE / 'diagnostic' / 'assembly_history_private.npz')
    allke = ah['ALLKE']
    allke_t = np.round(allke[:, 0] / 1e-7) * 1e-7
    allke_j = allke[:, 1] * 1e-3
    sampled_lumped = np.interp(allke_t, t, ktot_lumped)
    sampled_exact = np.interp(allke_t, t, ktot_exact)
    allke_abs = np.maximum(np.abs(allke_j), 1e-30)
    lumped_err = np.abs(sampled_lumped - allke_j)
    exact_err = np.abs(sampled_exact - allke_j)
    crosscheck = pd.DataFrame({
        'time_s': allke_t, 'ALLKE_J': allke_j,
        'Krigid_lumped_J': sampled_lumped, 'Krigid_continuum_J': sampled_exact,
        'lumped_abs_error_J': lumped_err, 'lumped_relative_error': lumped_err / allke_abs,
        'continuum_abs_error_J': exact_err, 'continuum_relative_error': exact_err / allke_abs,
    })
    crosscheck.to_csv(HERE / 'first_impact_ALLKE_crosscheck.csv', index=False)

    ilo = int(round(WINDOW[0] / 1e-7))
    ihi = int(round(WINDOW[1] / 1e-7))
    assert t[ilo] == WINDOW[0] and t[ihi] == WINDOW[1]
    force = np.linalg.norm(dense['Fwall'], axis=1)
    active = force > 1e-8
    active_idx = np.where(active & (t >= WINDOW[0]) & (t <= WINDOW[1]))[0]
    ia, ib = int(active_idx[0]), int(active_idx[-1])
    # Half-step shoulders capture the trapezoidal impulse/work around first/last nonzero samples.
    ipre, ipost = ia - 1, ib + 1

    timeline = pd.DataFrame({
        'time_s': t[ilo:ihi + 1],
        'Ktrans_J': itrans[ilo:ihi + 1],
        'Krot_J': krot_lumped[ilo:ihi + 1],
        'Ktotal_J': ktot_lumped[ilo:ihi + 1],
        'Wmag_trans_J': work['mag_trans'][ilo:ihi + 1] - work['mag_trans'][ilo],
        'Wmag_rot_J': work['mag_rot'][ilo:ihi + 1] - work['mag_rot'][ilo],
        'Whydro_trans_J': work['hydro_trans'][ilo:ihi + 1] - work['hydro_trans'][ilo],
        'Whydro_rot_J': work['hydro_rot'][ilo:ihi + 1] - work['hydro_rot'][ilo],
        'Wcontact_net_J': contact_net_work[ilo:ihi + 1] - contact_net_work[ilo],
        'Fwall_norm_N': force[ilo:ihi + 1],
        'omega_norm_rad_s': np.linalg.norm(omega[ilo:ihi + 1], axis=1),
        'Vcom_norm_mm_s': np.linalg.norm(vcom[ilo:ihi + 1], axis=1),
    })
    timeline.to_csv(HERE / 'first_impact_energy_timeline.csv', index=False)

    windows = []
    for name, a, b in [
        ('audit_0p98_to_1p08ms', ilo, ihi),
        ('first_impact_shoulders', ipre, ipost),
        ('field_endpoints_1p0000_to_1p0501ms', 10000, 10501),
    ]:
        vals = {
            'window': name, 'start_s': t[a], 'end_s': t[b],
            'Delta_Ktrans_J': delta(itrans, a, b),
            'Delta_Krot_J': delta(krot_lumped, a, b),
            'Delta_Ktotal_J': delta(ktot_lumped, a, b),
            'Wmag_trans_J': delta(work['mag_trans'], a, b),
            'Wmag_rot_J': delta(work['mag_rot'], a, b),
            'Whydro_trans_J': delta(work['hydro_trans'], a, b),
            'Whydro_rot_J': delta(work['hydro_rot'], a, b),
            'Wcontact_net_J': delta(contact_net_work, a, b),
            'Ktrans_start_J': itrans[a], 'Ktrans_end_J': itrans[b],
            'Krot_start_J': krot_lumped[a], 'Krot_end_J': krot_lumped[b],
        }
        vals['Wapplied_J'] = sum(vals[x] for x in (
            'Wmag_trans_J', 'Wmag_rot_J', 'Whydro_trans_J', 'Whydro_rot_J'))
        vals['energy_identity_error_J'] = vals['Delta_Ktotal_J'] - vals['Wapplied_J'] - vals['Wcontact_net_J']
        windows.append(vals)
    pd.DataFrame(windows).to_csv(HERE / 'first_impact_energy_windows.csv', index=False)

    event = windows[1]
    positive_external_work = max(event['Wapplied_J'], 0.0)
    transferred_lower_bound = max(event['Delta_Ktrans_J'] - positive_external_work, 0.0)
    summary = {
        'decision': 'ROTATION_TO_TRANSLATION_WITH_NET_CONTACT_DISSIPATION',
        'first_force_nonzero_ms': t[ia] * 1e3,
        'last_force_nonzero_ms': t[ib] * 1e3,
        'energy_window_start_ms': t[ipre] * 1e3,
        'energy_window_end_ms': t[ipost] * 1e3,
        'Delta_Ktrans_uJ': event['Delta_Ktrans_J'] * 1e6,
        'Delta_Krot_uJ': event['Delta_Krot_J'] * 1e6,
        'Delta_Ktotal_uJ': event['Delta_Ktotal_J'] * 1e6,
        'Wmag_uJ': (event['Wmag_trans_J'] + event['Wmag_rot_J']) * 1e6,
        'Whydro_uJ': (event['Whydro_trans_J'] + event['Whydro_rot_J']) * 1e6,
        'Wcontact_net_uJ': event['Wcontact_net_J'] * 1e6,
        'rotation_to_translation_lower_bound_uJ': transferred_lower_bound * 1e6,
        'translation_gain_from_preimpact_rotation_lower_bound_fraction': (
            transferred_lower_bound / max(event['Delta_Ktrans_J'], 1e-30)),
        'total_KE_retained_fraction': event['Ktrans_end_J'] + event['Krot_end_J'],
        'ALLKE_lumped_max_abs_error_uJ': float(lumped_err.max() * 1e6),
        'ALLKE_lumped_max_relative_error': float((lumped_err / allke_abs).max()),
        'ALLKE_continuum_max_relative_error': float((exact_err / allke_abs).max()),
        'inertia_lumped_tonne_mm2': inertia_lumped.tolist(),
        'inertia_continuum_tonne_mm2': inertia_exact.tolist(),
    }
    summary['total_KE_retained_fraction'] = (
        (event['Ktrans_end_J'] + event['Krot_end_J']) /
        max(event['Ktrans_start_J'] + event['Krot_start_J'], 1e-30))

    def assembly_delta(name, ta, tb):
        arr = ah[name]
        at = np.round(arr[:, 0] / 1e-7) * 1e-7
        a = int(np.where(at == ta)[0][0])
        b = int(np.where(at == tb)[0][0])
        return float((arr[b, 1] - arr[a, 1]) * 1e-3)

    abaqus = {
        'window': 'Abaqus_sparse_1p0000_to_1p0501ms',
        'start_s': 0.001, 'end_s': 0.0010501,
        'Delta_ALLKE_J': assembly_delta('ALLKE', 0.001, 0.0010501),
        'Delta_ALLWK_J': assembly_delta('ALLWK', 0.001, 0.0010501),
        'Delta_ALLVD_J': assembly_delta('ALLVD', 0.001, 0.0010501),
        'Delta_ALLFD_J': assembly_delta('ALLFD', 0.001, 0.0010501),
        'Delta_ALLPW_J': assembly_delta('ALLPW', 0.001, 0.0010501),
        'Delta_ETOTAL_J': assembly_delta('ETOTAL', 0.001, 0.0010501),
    }
    # Abaqus Explicit total-energy identity for this rigid-contact model.
    abaqus['contact_dissipation_effective_J'] = (
        abaqus['Delta_ALLVD_J'] + abaqus['Delta_ALLFD_J'] - abaqus['Delta_ALLPW_J'])
    abaqus['identity_error_J'] = (
        abaqus['Delta_ALLKE_J'] - abaqus['Delta_ALLWK_J'] +
        abaqus['contact_dissipation_effective_J'] - abaqus['Delta_ETOTAL_J'])
    dense_endpoint = windows[2]
    abaqus['dense_Delta_Ktotal_J'] = dense_endpoint['Delta_Ktotal_J']
    abaqus['dense_Wapplied_J'] = dense_endpoint['Wapplied_J']
    abaqus['dense_contact_dissipation_J'] = -dense_endpoint['Wcontact_net_J']
    pd.DataFrame([abaqus]).to_csv(HERE / 'first_impact_abaqus_energy_crosscheck.csv', index=False)
    summary.update({
        'Abaqus_1p0000_to_1p0501ms_contact_dissipation_uJ': abaqus['contact_dissipation_effective_J'] * 1e6,
        'dense_1p0000_to_1p0501ms_contact_dissipation_uJ': -dense_endpoint['Wcontact_net_J'] * 1e6,
        'Abaqus_energy_identity_error_uJ': abaqus['identity_error_J'] * 1e6,
    })
    if summary['Wcontact_net_uJ'] > 1e-3:
        summary['decision'] = 'ENERGY_AUDIT_AMBIGUOUS_POSITIVE_CONTACT_WORK'
    (HERE / 'first_impact_energy_summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(pd.DataFrame(windows).to_string(index=False))


if __name__ == '__main__':
    main()
