"""Fixed-camera slow-motion comparison of baseline and normal-damping probe."""
from pathlib import Path
import json
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from audit_stage_a import mesh_properties
from exact_gap_audit import Wall
from contact_probe_common import HERE, AUDIT, ROOT, dense_rp


JOB = 'Wobble_F30_G6L45_ReducedHydroFixed_NormalDamp020_ContactAudit_0013'
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'],
    'font.size': 7,
    'axes.titlesize': 8,
    'axes.labelsize': 6,
    'xtick.labelsize': 5,
    'ytick.labelsize': 5,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
})


def case(folder, inp, force_csv, summary):
    t, data = dense_rp(folder)
    meshes, rp, _ = mesh_properties(inp.read_text())
    mesh = meshes['Robot_SOLID']
    rot = Rotation.from_rotvec(data['UR'])
    off = mesh['com'] - rp
    vcom = data['V'] + np.cross(data['VR'], rot.apply(np.broadcast_to(off, data['U'].shape)))
    force = pd.read_csv(force_csv)
    force_t = force.time_s.to_numpy()
    assert np.all(np.diff(force_t) > 0)
    if 'Fwall_norm_N' in force:
        fnorm = force.Fwall_norm_N.to_numpy()
    else:
        fnorm = np.linalg.norm(force[['Ftotal_robot_x_N', 'Ftotal_robot_y_N', 'Ftotal_robot_z_N']], axis=1)
    assert len(force_t) == len(t) and np.allclose(force_t, t, rtol=0, atol=2e-12)
    return dict(t=t, data=data, mesh=mesh, rp=rp, rot=rot, vcom=vcom,
                fnorm=fnorm, summary=summary)


def main():
    wall = Wall()
    sb = json.loads((HERE / 'current_contact_restitution_summary.json').read_text())
    sc = json.loads((HERE / 'normaldamp020_summary.json').read_text())
    baseline = case(
        AUDIT / 'diagnostic', AUDIT / 'Wobble_F30_G6L45_ReducedHydroFixed_ContactAudit_0012.inp',
        AUDIT / 'first_impact_contact_force_history.csv', sb)
    candidate = case(
        HERE / 'candidate_private', HERE / (JOB + '.inp'),
        HERE / 'normaldamp020_contact_force_history.csv', sc)
    face_table = pd.read_csv(ROOT / 'Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083_robot_surface_triangles_exact.csv')
    face_ids = face_table[['n1', 'n2', 'n3']].to_numpy(int)
    node_ids = np.unique(face_ids)

    # Sparse full-window context plus dense frames around the microsecond impact.
    frame_times = np.unique(np.r_[
        np.linspace(0, .00095, 12),
        np.linspace(.00098, .000999, 8),
        np.arange(.001000, .0010101, .0000002),
        np.linspace(.001015, .00130, 14),
    ])
    frame_indices = np.rint(frame_times / 1e-7).astype(int)

    first_i = int(round(sb['first_contact_s'] / 1e-7))
    mesh = baseline['mesh']
    nodes0 = np.array([mesh['nodes'][int(i)] for i in node_ids]) + mesh['shift']
    pos0 = baseline['rp'] + baseline['data']['U'][first_i] + baseline['rot'][first_i].apply(nodes0 - baseline['rp'])
    g0, tri0, _ = wall.query(pos0)
    contact_center = wall.centers[tri0[np.argmin(g0)]]
    wall_mask = np.linalg.norm(wall.centers - contact_center, axis=1) < 3.2
    wall_polys = wall.tri[wall_mask]
    all_points = np.vstack((wall_polys.reshape(-1, 3), pos0))
    center = all_points.mean(axis=0)
    span = np.ptp(all_points, axis=0)
    half = max(span.max() * .58, 2.0)

    fig = plt.figure(figsize=(7.2047, 3.1496), constrained_layout=True)
    axes = [fig.add_subplot(1, 2, j + 1, projection='3d') for j in range(2)]

    def draw(frame_number):
        i = frame_indices[frame_number]
        for ax, c, title, color in zip(
                axes, [baseline, candidate], ['Baseline: ζ=0.055, tangent=1', 'Probe: ζ=0.20, tangent=0'],
                ['#5B6573', '#D55E00']):
            ax.clear()
            mesh = c['mesh']
            all_nodes = {label: np.asarray(coord) + mesh['shift'] for label, coord in mesh['nodes'].items()}
            transformed = {
                label: c['rp'] + c['data']['U'][i] + c['rot'][i].apply(point - c['rp'])
                for label, point in all_nodes.items()}
            robot_polys = np.array([[transformed[int(label)] for label in row] for row in face_ids])
            surface_pos = np.array([transformed[int(label)] for label in node_ids])
            gap, triangle, closest = wall.query(surface_pos)
            k = int(np.argmin(gap))
            node = surface_pos[k]
            tri = int(triangle[k])
            normal = wall.normals[tri]
            xcom = c['rp'] + c['data']['U'][i] + c['rot'][i].apply(mesh['com'] - c['rp'])
            rcontact = node - xcom
            vcontact = c['vcom'][i] + np.cross(c['data']['VR'][i], rcontact)
            vn = float(np.dot(vcontact, normal))
            ax.add_collection3d(Poly3DCollection(wall_polys, facecolor='#56B4E9', edgecolor='none', alpha=.10))
            ax.add_collection3d(Poly3DCollection(robot_polys, facecolor=color, edgecolor='#30343B', linewidth=.15, alpha=.82))
            if c['fnorm'][i] > 1e-8:
                ax.add_collection3d(Poly3DCollection([wall.tri[tri]], facecolor='#F0E442', edgecolor='#9A7D00', linewidth=.8, alpha=.85))
            ax.scatter(*node, s=18, color='#CC3311', depthshade=False)
            ax.quiver(*node, *(normal * .55), color='#0072B2', linewidth=1.2, arrow_length_ratio=.18)
            scale = 4e-4
            ax.quiver(*node, *(vcontact * scale), color='#009E73', linewidth=1.2, arrow_length_ratio=.18)
            ax.set(xlim=(center[0] - half, center[0] + half),
                   ylim=(center[1] - half, center[1] + half),
                   zlim=(center[2] - half, center[2] + half))
            ax.set_box_aspect((1, 1, 1))
            ax.view_init(elev=19, azim=-57)
            ax.set_title(title, pad=2)
            ax.set_xlabel('x (mm)', labelpad=-3)
            ax.set_ylabel('y (mm)', labelpad=-3)
            ax.set_zlabel('z (mm)', labelpad=-4)
            separated = c['t'][i] > c['summary']['first_contact_end_s']
            restitution = f'  e_n={c["summary"]["e_n"]:.3f}' if separated else ''
            ax.text2D(.02, .97,
                      f't={c["t"][i]*1e3:.4f} ms\ngap={gap[k]*1e3:+.3f} µm\n'
                      f'|Fwall|={c["fnorm"][i]:.3f} N\nv_n={vn:+.1f} mm s-1{restitution}',
                      transform=ax.transAxes, va='top')
        fig.suptitle('Reduced-Hydro first-impact comparison | fixed camera', fontsize=8)

    draw(int(np.argmin(np.abs(frame_indices - first_i))))
    poster = HERE / 'ReducedHydro_FirstImpact_Baseline_vs_NormalDamp020_poster'
    fig.savefig(str(poster) + '.png', dpi=600, bbox_inches='tight')
    fig.savefig(str(poster) + '.tiff', dpi=600, bbox_inches='tight')
    fig.savefig(str(poster) + '.pdf', bbox_inches='tight')
    fig.savefig(str(poster) + '.svg', bbox_inches='tight')
    animation = FuncAnimation(fig, draw, frames=len(frame_indices), interval=90, repeat=True)
    animation.save(HERE / 'ReducedHydro_FirstImpact_Baseline_vs_NormalDamp020.gif',
                   writer=PillowWriter(fps=10))
    plt.close(fig)


if __name__ == '__main__':
    main()
