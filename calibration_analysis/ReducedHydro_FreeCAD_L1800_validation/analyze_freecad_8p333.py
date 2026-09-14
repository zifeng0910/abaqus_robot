"""Unified true-axis, contact, bridge, torque and energy audit for the 8.333 ms run."""
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import maximum_filter1d, minimum_filter1d
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT = REPO.parent
AUDIT = REPO / 'calibration_analysis' / 'ReducedHydro_hidden_impact_audit'
PROBE = AUDIT / 'normal_contact_damping_probe'
sys.path.insert(0, str(AUDIT))
sys.path.insert(0, str(PROBE))
from audit_stage_a import mesh_properties
from exact_gap_audit import Wall, tests as gap_tests
from contact_probe_common import contact_force, dense_rp
from analyze_normaldamp020 import event_indices

JOB = 'Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_WallOn_Free_0083'
DT = 1e-7
AXIS_SEED = np.array([0.9647382600216, -0.1188742372140, 0.2348382536499])


def unit(v):
    v = np.asarray(v, float)
    return v / np.maximum(np.linalg.norm(v, axis=-1, keepdims=True), 1e-30)


def interp_vector(frame, t, columns, time_col):
    assert np.all(np.diff(frame[time_col]) > 0)
    return np.column_stack([np.interp(t, frame[time_col], frame[c]) for c in columns])


def canonical_projection(points, curve):
    cv = curve[['x_mm', 'y_mm', 'z_mm']].to_numpy(float)
    arc = curve.arclength_mm.to_numpy(float)
    edge = np.diff(cv, axis=0)
    edge2 = np.einsum('ij,ij->i', edge, edge)
    tangents = unit(edge)
    rows = []
    previous_segment = None
    previous_s = None
    for point in points:
        if previous_segment is None:
            candidates = np.arange(len(edge))
        else:
            candidates = np.arange(max(0, previous_segment - 8), min(len(edge), previous_segment + 9))
        a = cv[candidates]
        u = np.einsum('ij,ij->i', point - a, edge[candidates]) / edge2[candidates]
        for q, seg in enumerate(candidates):
            if seg == 0 and u[q] < 0:
                continue
            if seg == len(edge) - 1 and u[q] > 1:
                continue
            u[q] = np.clip(u[q], 0, 1)
        projected = a + u[:, None] * edge[candidates]
        s = arc[candidates] + u * (arc[candidates + 1] - arc[candidates])
        cost = np.einsum('ij,ij->i', point - projected, point - projected)
        if previous_s is not None:
            cost += (.02 * (s - previous_s)) ** 2
        k = int(np.argmin(cost))
        seg = int(candidates[k])
        # A large residual triggers a global search instead of trusting continuity.
        if cost[k] > 1.0:
            candidates = np.arange(len(edge))
            a = cv[:-1]
            u = np.einsum('ij,ij->i', point - a, edge) / edge2
            u[1:-1] = np.clip(u[1:-1], 0, 1)
            u[0] = min(u[0], 1)
            u[-1] = max(u[-1], 0)
            projected = a + u[:, None] * edge
            s = arc[:-1] + u * np.diff(arc)
            cost = np.einsum('ij,ij->i', point - projected, point - projected)
            if previous_s is not None:
                cost += (.02 * (s - previous_s)) ** 2
            k = int(np.argmin(cost)); seg = k
        rows.append((s[k], *projected[k], *tangents[seg], seg))
        previous_segment, previous_s = seg, s[k]
    return pd.DataFrame(rows, columns=['s_mm', 'proj_x', 'proj_y', 'proj_z',
                                       't_x', 't_y', 't_z', 'segment'])


def local_frame(tangent):
    e1 = np.empty_like(tangent)
    for i, th in enumerate(tangent):
        candidate = np.array([0., 0., 1.]) - th[2] * th
        if np.linalg.norm(candidate) < 1e-8:
            candidate = np.array([0., 1., 0.]) - th[1] * th
        candidate /= np.linalg.norm(candidate)
        if i and np.dot(candidate, e1[i - 1]) < 0:
            candidate *= -1
        e1[i] = candidate
    return e1, unit(np.cross(tangent, e1))


def robust_slope(time, phase):
    x = time - time.mean()
    slope, intercept = np.polyfit(x, phase, 1)
    residual = phase - (slope * x + intercept)
    mad = np.median(np.abs(residual - np.median(residual)))
    keep = np.ones(len(x), bool) if mad < 1e-14 else np.abs(residual - np.median(residual)) <= 4 * 1.4826 * mad
    if keep.sum() >= max(3, int(.7 * len(x))):
        slope = np.polyfit(x[keep], phase[keep], 1)[0]
    return float(slope)


def phase_windows(t, robot_phase, field_phase, width_s=.00075, stride_s=.00005):
    rows = []
    half = width_s / 2
    for center in np.arange(t[0] + half, t[-1] - half + DT, stride_s):
        mask = (t >= center - half) & (t <= center + half)
        rows.append({
            'center_s': center, 'start_s': center - half, 'end_s': center + half,
            'samples': int(mask.sum()),
            'robot_phase_rate_Hz': robust_slope(t[mask], robot_phase[mask]) / (2 * np.pi),
            'field_local_phase_rate_Hz': robust_slope(t[mask], field_phase[mask]) / (2 * np.pi),
            'robot_phase_advance_deg': np.degrees(robot_phase[mask][-1] - robot_phase[mask][0]),
            'field_phase_advance_deg': np.degrees(field_phase[mask][-1] - field_phase[mask][0]),
        })
    return pd.DataFrame(rows)


def plateau_timeline(t, robot_phase, field_phase, window_s=.0005, tolerance_deg=10.):
    count = int(round(window_s / np.median(np.diff(t)))) + 1
    phase_range = maximum_filter1d(robot_phase, size=count, mode='nearest') - minimum_filter1d(
        robot_phase, size=count, mode='nearest')
    field_rate = np.gradient(field_phase, t)
    candidate = (phase_range < np.radians(tolerance_deg)) & (np.abs(field_rate) > np.radians(10) / window_s)
    edges = np.diff(np.r_[False, candidate, False].astype(int))
    starts, ends = np.where(edges == 1)[0], np.where(edges == -1)[0] - 1
    events = []
    for a, b in zip(starts, ends):
        start = max(t[0], t[a] - window_s / 2)
        end = min(t[-1], t[b] + window_s / 2)
        events.append((start, end))
    return candidate, events, np.degrees(phase_range)


def exact_state(index, data, mesh, rp, rotation, nodes, node_ids, wall):
    positions = rp + data['U'][index] + rotation[index].apply(nodes - rp)
    gaps, triangles, closest = wall.query(positions)
    k = int(np.argmin(gaps))
    return positions, gaps, triangles, closest, k


def angular_clusters(angles, max_gap_deg=45.):
    if not len(angles):
        return []
    angles = np.mod(np.asarray(angles), 2 * np.pi)
    order = np.argsort(angles); values = angles[order]
    gaps = np.diff(np.r_[values, values[0] + 2 * np.pi])
    cut = int(np.argmax(gaps))
    rotated = np.r_[values[cut + 1:], values[:cut + 1] + 2 * np.pi]
    breaks = np.where(np.diff(rotated) > np.radians(max_gap_deg))[0] + 1
    groups = np.split(rotated, breaks)
    return [float(np.mod(np.angle(np.mean(np.exp(1j * group))), 2 * np.pi)) for group in groups if len(group)]


def longest_true_interval(t, flag):
    edges = np.diff(np.r_[False, flag, False].astype(int))
    starts, ends = np.where(edges == 1)[0], np.where(edges == -1)[0] - 1
    events = [(float(t[a]), float(t[b]), float(t[b] - t[a] + np.median(np.diff(t)))) for a, b in zip(starts, ends)]
    return events, max((event[2] for event in events), default=0.)


def analyze_new_case():
    assert gap_tests().startswith('PASS')
    t, data = dense_rp(HERE / 'candidate_private')
    assert len(t) == 83331 and np.allclose(np.diff(t), DT, rtol=0, atol=2e-12)
    meshes, rp, _ = mesh_properties((HERE / f'{JOB}.inp').read_text())
    mesh = meshes['Robot_SOLID']; rotation = Rotation.from_rotvec(data['UR'])
    all_nodes = np.array([mesh['nodes'][i] for i in sorted(mesh['nodes'])]) + mesh['shift']
    all_ids = np.array(sorted(mesh['nodes']))
    projection = all_nodes @ (AXIS_SEED / np.linalg.norm(AXIS_SEED))
    low, high = np.quantile(projection, [.10, .90])
    head_ids, tail_ids = all_ids[projection <= low], all_ids[projection >= high]
    head0 = all_nodes[projection <= low].mean(axis=0); tail0 = all_nodes[projection >= high].mean(axis=0)
    tail = rp + data['U'] + rotation.apply(np.broadcast_to(tail0 - rp, data['U'].shape))
    head = rp + data['U'] + rotation.apply(np.broadcast_to(head0 - rp, data['U'].shape))
    axis = unit(tail - head)
    com_offset = rotation.apply(np.broadcast_to(mesh['com'] - rp, data['U'].shape))
    com = rp + data['U'] + com_offset
    vcom = data['V'] + np.cross(data['VR'], com_offset)

    curve = pd.read_csv(ROOT / 'CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_TRUE_centerline_odb.csv')
    projected = canonical_projection(com, curve)
    tangent = projected[['t_x', 't_y', 't_z']].to_numpy()
    e1, e2 = local_frame(tangent)
    aperp = axis - np.einsum('ij,ij->i', axis, tangent)[:, None] * tangent
    phi_robot = np.unwrap(np.arctan2(np.einsum('ij,ij->i', aperp, e2),
                                     np.einsum('ij,ij->i', aperp, e1)))
    axis_dot_tangent = np.einsum('ij,ij->i', axis, tangent)
    directed_tilt = np.degrees(np.arccos(np.clip(axis_dot_tangent, -1, 1)))
    tilt = np.degrees(np.arccos(np.clip(np.abs(axis_dot_tangent), 0, 1)))

    telemetry = pd.read_csv(ROOT / f'{JOB}_telemetry.csv')
    B = interp_vector(telemetry, t, ['Bx_aba_T', 'By_aba_T', 'Bz_aba_T'], 't_s')
    fmag = interp_vector(telemetry, t, ['fx_aba_N', 'fy_aba_N', 'fz_aba_N'], 't_s')
    tmag = interp_vector(telemetry, t, ['tx_aba_Nmm', 'ty_aba_Nmm', 'tz_aba_Nmm'], 't_s')
    phi_B = np.unwrap(np.arctan2(np.einsum('ij,ij->i', B, e2), np.einsum('ij,ij->i', B, e1)))
    phi_cmd = np.interp(t, telemetry.t_s, telemetry.instantaneous_phase_rad)
    phase_diff = np.unwrap(phi_robot - phi_B)
    Bhat = unit(B)
    theta_B = np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i', Bhat, tangent), -1, 1)))
    robot_B_misalignment = np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i', axis, Bhat), -1, 1)))
    windows = phase_windows(t, phi_robot, phi_B)
    windows.to_csv(HERE / 'freecad_8p333_phase_windows.csv', index=False)
    plateau, plateau_events, phase_range_deg = plateau_timeline(t, phi_robot, phi_B)

    hydro = pd.read_csv(ROOT / f'{JOB}_hydro_increment.csv', skiprows=2,
                        names=['time_s', 'v1', 'v2', 'v3', 'w1', 'w2', 'w3',
                               'Fh1', 'Fh2', 'Fh3', 'Th1', 'Th2', 'Th3'])
    fhyd = interp_vector(hydro, t, ['Fh1', 'Fh2', 'Fh3'], 'time_s')
    thyd = interp_vector(hydro, t, ['Th1', 'Th2', 'Th3'], 'time_s')
    omega = data['VR']
    tmag_axis = np.einsum('ij,ij->i', tmag, axis)
    tmag_tilt = np.linalg.norm(tmag - tmag_axis[:, None] * axis, axis=1)
    thydro_power = np.einsum('ij,ij->i', thyd, omega) * 1e-3
    pmag = (np.einsum('ij,ij->i', fmag, vcom) + np.einsum('ij,ij->i', tmag, omega)) * 1e-3
    phydro = (np.einsum('ij,ij->i', fhyd, vcom) + np.einsum('ij,ij->i', thyd, omega)) * 1e-3
    wmag = cumulative_trapezoid(pmag, t, initial=0.)
    whydro = cumulative_trapezoid(phydro, t, initial=0.)

    fn, fs = contact_force(HERE / 'candidate_private', t); fwall = fn + fs
    active, events = event_indices(fwall)
    faces = pd.read_csv(HERE / 'freecad_robot_surface_triangles_exact.csv')
    node_ids = np.unique(faces[['n1', 'n2', 'n3']].to_numpy()).astype(int)
    nodes = np.array([mesh['nodes'][int(i)] for i in node_ids]) + mesh['shift']
    wall = Wall()

    def query_indices(indices, collect_clusters=False):
        cluster_rows = []
        indices = np.asarray(indices, dtype=int)
        base = nodes - rp
        for offset in range(0, len(indices), 24):
            block = indices[offset:offset + 24]
            matrices = rotation[block].as_matrix()
            positions = (rp + data['U'][block, None, :]
                         + np.einsum('bij,nj->bni', matrices, base))
            gaps, triangles, closest = wall.query(positions.reshape(-1, 3))
            gaps = gaps.reshape(len(block), len(nodes))
            triangles = triangles.reshape(len(block), len(nodes))
            closest = closest.reshape(len(block), len(nodes), 3)
            for row, i in enumerate(block):
                k = int(np.argmin(gaps[row]))
                exact[int(i)] = (float(gaps[row, k]), int(node_ids[k]), int(triangles[row, k]))
                if not collect_clusters or gaps[row, k] > .02:
                    continue
                for threshold_um in (20., 5., 0.):
                    selected = np.where(gaps[row] * 1e3 <= threshold_um)[0]
                    normals = wall.normals[triangles[row, selected]] if len(selected) else np.empty((0, 3))
                    angles = np.arctan2(normals @ e2[i], normals @ e1[i]) if len(selected) else np.array([])
                    centers = angular_clusters(angles)
                    separation = max((np.degrees(np.arccos(np.clip(np.cos(a - b), -1, 1)))
                                      for q, a in enumerate(centers) for b in centers[q + 1:]), default=0.)
                    opposing = bool(len(centers) >= 2 and separation >= 120.)
                    ids = node_ids[selected]
                    cluster_rows.append({
                        'time_s': t[i], 'increment': int(i), 'threshold_um': threshold_um,
                        'node_count': len(ids), 'cluster_count': len(centers),
                        'sector_centers_deg': ';'.join(f'{np.degrees(a):.3f}' for a in centers),
                        'max_sector_separation_deg': separation,
                        'head_involved': bool(np.intersect1d(ids, head_ids).size),
                        'tail_involved': bool(np.intersect1d(ids, tail_ids).size),
                        'participating_nodes': ';'.join(str(int(x)) for x in ids),
                        'opposing_bridge': opposing,
                    })
        return cluster_rows

    coarse = np.arange(0, len(t), 10, dtype=int)
    if coarse[-1] != len(t) - 1:
        coarse = np.r_[coarse, len(t) - 1]
    exact = {}
    cluster_rows = query_indices(coarse, collect_clusters=True)
    dense = set()
    for start, end in events:
        dense.update(range(max(0, start - 50), min(len(t), end + 51)))
    cluster_rows.extend(query_indices(sorted(dense - set(coarse)), collect_clusters=True))
    gap_rows = [{'time_s': t[i], 'increment': i, 'gap_um': value[0] * 1e3,
                 'contact_node': value[1], 'wall_triangle': value[2],
                 'sampling': 'dense_0p1us' if i in dense else 'coarse_1us'}
                for i, value in sorted(exact.items())]
    gap_df = pd.DataFrame(gap_rows)
    gap_df.to_csv(HERE / 'freecad_8p333_exact_gap.csv', index=False)

    clusters = pd.DataFrame(cluster_rows)
    clusters.to_csv(HERE / 'freecad_8p333_contact_clusters.csv', index=False)
    bridge_flag = np.zeros(len(t), bool)
    at20 = clusters[clusters.threshold_um == 20.]
    bridge_samples = dict(zip(at20.increment.to_numpy(int), at20.opposing_bridge.to_numpy(bool)))
    for a, z in zip(coarse[:-1], coarse[1:]):
        bridge_flag[a:z] = bridge_samples.get(int(a), False)
    bridge_flag[coarse[-1]] = bridge_samples.get(int(coarse[-1]), False)
    for i in dense:
        if i in bridge_samples: bridge_flag[i] = bridge_samples[i]
    bridge_events, longest_bridge = longest_true_interval(t, bridge_flag)
    pd.DataFrame({'time_s': t, 'opposing_bridge': bridge_flag.astype(int),
                  'phase_plateau': plateau.astype(int), 'tilt_deg': tilt,
                  'field_phase_continues': (np.abs(np.gradient(phi_B, t)) > np.radians(10) / .0005).astype(int)}
                 ).to_csv(HERE / 'freecad_8p333_bridge_timeline.csv', index=False)

    event_rows = []
    for number, (start, end) in enumerate(events, 1):
        indices = [i for i in range(start, end + 1) if i in exact]
        min_i = min(indices, key=lambda i: exact[i][0])
        node = exact[min_i][1]; tri = exact[min_i][2]
        node0 = np.asarray(mesh['nodes'][node]) + mesh['shift']
        xnode = rp + data['U'] + rotation.apply(np.broadcast_to(node0 - rp, data['U'].shape))
        rcontact = xnode - com
        vcontact = vcom + np.cross(omega, rcontact)
        normal = wall.normals[tri]
        vn = vcontact @ normal
        pre = np.arange(max(0, start - 5), start)
        post = np.arange(end + 1, min(len(t), end + 6))
        vn_in = float(vn[pre].mean()) if len(pre) else np.nan
        vn_out = float(vn[post].mean()) if len(post) else np.nan
        impulse = np.trapz(fwall[max(0, start - 1):min(len(t), end + 2)],
                           t[max(0, start - 1):min(len(t), end + 2)], axis=0)
        event_rows.append({
            'event_number': number, 'start_s': t[start], 'end_s': t[end],
            'duration_us': (end - start + 1) * DT * 1e6,
            'peak_force_N': np.linalg.norm(fwall[start:end + 1], axis=1).max(),
            'impulse_x_Ns': impulse[0], 'impulse_y_Ns': impulse[1], 'impulse_z_Ns': impulse[2],
            'impulse_norm_Ns': np.linalg.norm(impulse), 'contact_node': node,
            'wall_triangle': tri, 'min_gap_um': exact[min_i][0] * 1e3,
            'contact_x_mm':xnode[min_i,0],'contact_y_mm':xnode[min_i,1],'contact_z_mm':xnode[min_i,2],
            'vn_in_mm_s': vn_in, 'vn_out_mm_s': vn_out,
            'effective_restitution': vn_out / -vn_in if vn_in < 0 and vn_out > 0 else np.nan,
            'separable_rebound': bool(vn_in < 0 and vn_out > 0),
        })
    event_catalog = pd.DataFrame(event_rows)
    event_catalog.to_csv(HERE / 'freecad_8p333_contact_events.csv', index=False)

    phase_data = pd.DataFrame({
        'time_s': t, 'phi_robot_rad': phi_robot, 'phi_robot_advance_deg': np.degrees(phi_robot - phi_robot[0]),
        'phi_B_local_rad': phi_B, 'phi_cmd_rad': phi_cmd,
        'phi_robot_minus_B_rad': phase_diff, 'phase_window_range_deg': phase_range_deg,
        'phase_plateau': plateau.astype(int),
    })
    phase_data.to_csv(HERE / 'freecad_8p333_local_phase.csv', index=False)
    pd.DataFrame({'time_s':t,'theta_B_deg':theta_B,'B_local_phase_rad':phi_B,
                  'robot_B_misalignment_deg':robot_B_misalignment,
                  'B_x_T':B[:,0],'B_y_T':B[:,1],'B_z_T':B[:,2]}).to_csv(
                      HERE/'freecad_8p333_field_orientation.csv',index=False)
    pd.DataFrame({'time_s':t,'directed_tilt_deg':directed_tilt,'folded_tilt_deg':tilt,
                  'axis_dot_tangent':axis_dot_tangent}).to_csv(
                      HERE/'freecad_8p333_directed_tilt.csv',index=False)
    axis_data = pd.DataFrame({
        'time_s': t, 'head_x_mm': head[:, 0], 'head_y_mm': head[:, 1], 'head_z_mm': head[:, 2],
        'tail_x_mm': tail[:, 0], 'tail_y_mm': tail[:, 1], 'tail_z_mm': tail[:, 2],
        'axis_x': axis[:, 0], 'axis_y': axis[:, 1], 'axis_z': axis[:, 2],
        'tangent_x': tangent[:, 0], 'tangent_y': tangent[:, 1], 'tangent_z': tangent[:, 2],
        'e1_x': e1[:, 0], 'e1_y': e1[:, 1], 'e1_z': e1[:, 2],
        'e2_x': e2[:, 0], 'e2_y': e2[:, 1], 'e2_z': e2[:, 2], 'tilt_deg': tilt,
        'directed_tilt_deg':directed_tilt,'axis_dot_tangent':axis_dot_tangent,
    })
    axis_data.to_csv(HERE / 'freecad_8p333_true_axis.csv', index=False)
    torque_data = pd.DataFrame({
        'time_s': t, 'Tmag_x_Nmm': tmag[:, 0], 'Tmag_y_Nmm': tmag[:, 1], 'Tmag_z_Nmm': tmag[:, 2],
        'Tmag_norm_Nmm': np.linalg.norm(tmag, axis=1), 'Tmag_axis_Nmm': tmag_axis,
        'Tmag_tilt_Nmm': tmag_tilt, 'Thydro_x_Nmm': thyd[:, 0],
        'Thydro_y_Nmm': thyd[:, 1], 'Thydro_z_Nmm': thyd[:, 2],
        'Thydro_norm_Nmm': np.linalg.norm(thyd, axis=1), 'Thydro_dot_omega_W': thydro_power,
        'Tmag_to_Thydro_ratio': np.linalg.norm(tmag, axis=1) / np.maximum(np.linalg.norm(thyd, axis=1), 1e-30),
    })
    torque_data.to_csv(HERE / 'freecad_8p333_magnetic_hydro_torque.csv', index=False)
    translation = pd.DataFrame({
        'time_s': t, 'canonical_s_mm': projected.s_mm,
        'delta_s_mm': projected.s_mm - projected.s_mm.iloc[0],
        'projection_residual_mm': np.linalg.norm(com - projected[['proj_x', 'proj_y', 'proj_z']].to_numpy(), axis=1),
        'Vt_mm_s': np.einsum('ij,ij->i', vcom, tangent),
    })
    translation.to_csv(HERE / 'freecad_8p333_canonical_translation.csv', index=False)

    rotmat = rotation.as_matrix()
    inertia_global = np.einsum('nij,jk,nlk->nil', rotmat, mesh['I'], rotmat)
    ktrans = .5 * mesh['mass'] * np.einsum('ij,ij->i', vcom, vcom) * 1e-3
    krot = .5 * np.einsum('ni,nij,nj->n', omega, inertia_global, omega) * 1e-3
    assembly = np.load(HERE / 'candidate_private' / 'assembly_history_private.npz')
    energy = pd.DataFrame({'time_s': t, 'Ktrans_J': ktrans, 'Krot_J': krot,
                           'Ktotal_rigid_J': ktrans + krot, 'Wmag_J': wmag, 'Whydro_J': whydro,
                           'Wcontact_inferred_J': (ktrans + krot) - (ktrans[0] + krot[0]) - wmag - whydro})
    for key in assembly.files:
        values = assembly[key]
        # ODB energies use the active N-mm-s system, so N mm must be multiplied by 1e-3 for joules.
        energy['Abaqus_' + key + '_J'] = np.interp(t, values[:, 0], values[:, 1]) * 1e-3
    energy.to_csv(HERE / 'freecad_8p333_energy.csv', index=False)

    first_contact = float(event_catalog.start_s.iloc[0])
    first_index = int(round(first_contact / DT))
    post_phase_advance = float(np.degrees(phi_robot[-1] - phi_robot[first_index]))
    pre = t < first_contact
    post = t > event_catalog.end_s.iloc[0]
    plateau_longest = max((end - start for start, end in plateau_events), default=0.)
    runtime = {
        'job': JOB, 'completed_s': float(t[-1]), 'increments': len(t) - 1,
        'dt_s': DT, 'dt_stable': True, 'warning_count': 13,
        'field_frames': 335, 'field_interval_s': 2.5e-5,
        'first_contact_s': first_contact, 'contact_event_count': len(events),
        'longest_contact_us': float(event_catalog.duration_us.max()),
        'min_exact_gap_um': float(gap_df.gap_um.min()),
        'first_bridge_s': bridge_events[0][0] if bridge_events else None,
        'longest_bridge_ms': longest_bridge * 1e3,
        'persistent_bridge': bool(longest_bridge > .0005),
        'max_tilt_deg': float(tilt.max()), 'late_tilt_deg': float(tilt[-1]),
        'max_directed_tilt_deg':float(directed_tilt.max()),'final_directed_tilt_deg':float(directed_tilt[-1]),
        'axis_crossed_90deg':bool(np.any(axis_dot_tangent<0)),
        # Canonical increasing-s points opposite production a0 initially.  A
        # reversal is therefore a sign change relative to the initial polarity,
        # not the first negative sample.
        'first_axis_reversal_s': (float(t[np.where(
            np.signbit(axis_dot_tangent) != np.signbit(axis_dot_tangent[0]))[0][0]])
            if np.any(np.signbit(axis_dot_tangent) != np.signbit(axis_dot_tangent[0])) else None),
        'theta_B_min_deg':float(theta_B.min()),'theta_B_max_deg':float(theta_B.max()),
        'robot_B_misalignment_max_deg':float(robot_B_misalignment.max()),
        'postimpact_phase_advance_deg': post_phase_advance,
        'longest_phase_plateau_ms': plateau_longest * 1e3,
        'phase_plateau_over_0p5ms': bool(plateau_longest > .0005),
        'preimpact_Tmag_RMS_Nmm': float(np.sqrt(np.mean(np.linalg.norm(tmag[pre], axis=1) ** 2))),
        'postimpact_Tmag_RMS_Nmm': float(np.sqrt(np.mean(np.linalg.norm(tmag[post], axis=1) ** 2))),
        'Tmag_collapsed': bool(np.sqrt(np.mean(np.linalg.norm(tmag[post], axis=1) ** 2)) <
                               .5 * np.sqrt(np.mean(np.linalg.norm(tmag[pre], axis=1) ** 2))),
        'Thydro_power_max_W': float(thydro_power.max()),
        'Thydro_power_positive_fraction': float(np.mean(thydro_power > 1e-15)),
        'Abaqus_ETOTAL_min_J': float(energy.Abaqus_ETOTAL_J.min()),
        'Abaqus_ETOTAL_max_J': float(energy.Abaqus_ETOTAL_J.max()),
        'Abaqus_ALLCD_end_J': float(energy.Abaqus_ALLCD_J.iloc[-1]),
        'canonical_delta_s_mm': float(translation.delta_s_mm.iloc[-1]),
        'final_Vt_mm_s': float(translation.Vt_mm_s.iloc[-1]),
        'head_node_ids': [int(x) for x in head_ids], 'tail_node_ids': [int(x) for x in tail_ids],
        'centerline': 'CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_TRUE_centerline_odb.csv',
    }
    (HERE / 'freecad_8p333_summary.json').write_text(json.dumps(runtime, indent=2))
    pd.DataFrame([runtime]).drop(columns=['head_node_ids', 'tail_node_ids']).to_csv(
        HERE / 'freecad_8p333_runtime_identity.csv', index=False)
    return runtime, dict(t=t, data=data, mesh=mesh, rp=rp, rotation=rotation,
                         axis=axis, head=head, tail=tail, tangent=tangent, e1=e1, e2=e2,
                         com=com, phase=phi_robot, phi_B=phi_B, tilt=tilt,
                         bridge=bridge_flag, active=active, fwall=fwall,
                         translation=translation, gap=gap_df, wall=wall, node_ids=node_ids, nodes=nodes)


def old_case_summary(label, job, rp_frame, gap_path, mesh, rp):
    t = rp_frame.time_s.to_numpy(float)
    U = rp_frame[['U1', 'U2', 'U3']].to_numpy(float)
    UR = rp_frame[['UR1', 'UR2', 'UR3']].to_numpy(float)
    rotation = Rotation.from_rotvec(UR)
    all_nodes = np.array([mesh['nodes'][i] for i in sorted(mesh['nodes'])]) + mesh['shift']
    all_ids = np.array(sorted(mesh['nodes']))
    q = all_nodes @ (AXIS_SEED / np.linalg.norm(AXIS_SEED)); lo, hi = np.quantile(q, [.10, .90])
    axis0 = unit(all_nodes[q >= hi].mean(axis=0) - all_nodes[q <= lo].mean(axis=0))
    axis = rotation.apply(np.broadcast_to(axis0, U.shape))
    com = rp + U + rotation.apply(np.broadcast_to(mesh['com'] - rp, U.shape))
    curve = pd.read_csv(ROOT / 'CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_TRUE_centerline_odb.csv')
    projected = canonical_projection(com, curve)
    tangent = projected[['t_x', 't_y', 't_z']].to_numpy(); e1, e2 = local_frame(tangent)
    aperp = axis - np.einsum('ij,ij->i', axis, tangent)[:, None] * tangent
    phase = np.unwrap(np.arctan2(np.einsum('ij,ij->i', aperp, e2), np.einsum('ij,ij->i', aperp, e1)))
    tilt = np.degrees(np.arccos(np.clip(np.abs(np.einsum('ij,ij->i', axis, tangent)), 0, 1)))
    faces = pd.read_csv(HERE / 'freecad_robot_surface_triangles_exact.csv')
    node_ids = np.unique(faces[['n1', 'n2', 'n3']].to_numpy()).astype(int)
    nodes = np.array([mesh['nodes'][int(i)] for i in node_ids]) + mesh['shift']
    wall = Wall(); exact_gap = []; bridge = []
    for i in range(len(t)):
        positions = rp + U[i] + rotation[i].apply(nodes - rp)
        gaps, triangles, _ = wall.query(positions)
        exact_gap.append(float(gaps.min()))
        selected = np.where(gaps <= .02)[0]
        normals = wall.normals[triangles[selected]] if len(selected) else np.empty((0, 3))
        angles = np.arctan2(normals @ e2[i], normals @ e1[i]) if len(selected) else np.array([])
        centers = angular_clusters(angles)
        separation = max((np.degrees(np.arccos(np.clip(np.cos(a - b), -1, 1)))
                          for q, a in enumerate(centers) for b in centers[q + 1:]), default=0.)
        bridge.append(len(centers) >= 2 and separation >= 120.)
    exact_gap = np.asarray(exact_gap)
    contact = exact_gap <= 0
    contact_events, longest_contact = longest_true_interval(t, contact)
    bridge_events, longest_bridge = longest_true_interval(t, np.asarray(bridge, bool))
    first = contact_events[0][0] if contact_events else np.nan
    start = int(np.argmin(np.abs(t - first))) if np.isfinite(first) else 0
    _, plateau_events, _ = plateau_timeline(t, phase, 2 * np.pi * 30. * t)
    longest_plateau = max((b - a for a, b in plateau_events), default=0.)
    V = (rp_frame[['V1', 'V2', 'V3']].to_numpy(float)
         if set(['V1', 'V2', 'V3']).issubset(rp_frame.columns) else np.gradient(com, t, axis=0))
    return {'case': label, 'job': job, 'first_contact_s': first,
            'contact_events': len(contact_events), 'longest_contact_ms': longest_contact * 1e3,
            'min_gap_um': float(exact_gap.min() * 1e3), 'max_tilt_deg': float(tilt.max()),
            'postimpact_phase_advance_deg': float(np.degrees(phase[-1] - phase[start])),
            'longest_phase_plateau_ms': longest_plateau * 1e3,
            'opposing_bridge_onset_s': bridge_events[0][0] if bridge_events else np.nan,
            'longest_bridge_ms': longest_bridge * 1e3,
            'canonical_delta_s_mm': float(projected.s_mm.iloc[-1] - projected.s_mm.iloc[0]),
            'final_Vt_mm_s': float(np.dot(V[-1], tangent[-1])),
            'contact_source': '50 us exact-gap snapshots; micro-impacts may be missed',
            'time_s': t, 'phase_rad': phase, 'tilt_deg': tilt}


def comparison(runtime, new_data):
    meshes, rp, _ = mesh_properties((HERE / f'{JOB}.inp').read_text()); mesh = meshes['Robot_SOLID']
    rh_job = 'Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083'
    rh_rp = pd.read_csv(ROOT / f'{rh_job}_rp_fields.csv')
    rh = old_case_summary('RH old contact', rh_job, rh_rp,
                          ROOT / f'{rh_job}_exact_wall_penetration.csv', mesh, rp)
    cel_job = 'Wobble_F30_G6L45_WallOn_Free_0083_WobbleSurvival'
    raw = json.loads((ROOT / f'{cel_job}_rp_fields.json').read_text())['rows']
    cel_rp = pd.DataFrame({'time_s': [row['time'] for row in raw]})
    for prefix in ('U', 'UR'):
        values = np.asarray([row[prefix] for row in raw], float)
        for i in range(3): cel_rp[prefix + str(i + 1)] = values[:, i]
    cel = old_case_summary('CEL', cel_job, cel_rp,
                           ROOT / f'{cel_job}_exact_wall_penetration.csv', mesh, rp)
    new = {'case': 'RH zeta0.50', 'job': JOB, 'first_contact_s': runtime['first_contact_s'],
           'contact_events': runtime['contact_event_count'],
           'longest_contact_ms': runtime['longest_contact_us'] / 1e3,
           'min_gap_um': runtime['min_exact_gap_um'], 'max_tilt_deg': runtime['max_tilt_deg'],
           'postimpact_phase_advance_deg': runtime['postimpact_phase_advance_deg'],
           'longest_phase_plateau_ms': runtime['longest_phase_plateau_ms'],
           'opposing_bridge_onset_s': runtime['first_bridge_s'],
           'longest_bridge_ms': runtime['longest_bridge_ms'],
           'canonical_delta_s_mm': runtime['canonical_delta_s_mm'],
           'final_Vt_mm_s': runtime['final_Vt_mm_s'],
           'contact_source': '0.1 us resultant history plus exact dense near-wall geometry',
           'time_s': new_data['t'], 'phase_rad': new_data['phase'], 'tilt_deg': new_data['tilt']}
    public_columns = [key for key in new if key not in ('time_s', 'phase_rad', 'tilt_deg')]
    pd.DataFrame([{key: case[key] for key in public_columns} for case in (cel, rh, new)]).to_csv(
        HERE / 'cel_vs_reducedhydro_zeta050_summary.csv', index=False)
    np.savez_compressed(HERE / 'comparison_series_private.npz',
                        cel_t=cel['time_s'], cel_phase=cel['phase_rad'], cel_tilt=cel['tilt_deg'],
                        rh_t=rh['time_s'], rh_phase=rh['phase_rad'], rh_tilt=rh['tilt_deg'])


if __name__ == '__main__':
    result, series = analyze_new_case()
    comparison(result, series)
    print(json.dumps(result, indent=2))

