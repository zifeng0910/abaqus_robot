"""Fixed-camera first-impact animation for the strict 0.20 versus 0.50 comparison."""
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
from contact_probe_common import HERE, ROOT, dense_rp


mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'],
    'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 6,
    'xtick.labelsize': 5, 'ytick.labelsize': 5,
    'svg.fonttype': 'none', 'pdf.fonttype': 42,
})


def load_case(zeta):
    tag = f'normaldamp0{zeta}'
    job = f'Wobble_F30_G6L45_ReducedHydroFixed_NormalDamp0{zeta}_ContactAudit_0013'
    folder = HERE / ('candidate_private' if zeta == 20 else 'candidate050_private')
    t, data = dense_rp(folder)
    meshes, rp, _ = mesh_properties((HERE / f'{job}.inp').read_text())
    mesh = meshes['Robot_SOLID']
    rot = Rotation.from_rotvec(data['UR'])
    off = mesh['com'] - rp
    vcom = data['V'] + np.cross(data['VR'], rot.apply(np.broadcast_to(off, data['U'].shape)))
    force = pd.read_csv(HERE / f'{tag}_contact_force_history.csv')
    assert len(force) == len(t) and np.allclose(force.time_s, t, rtol=0, atol=2e-12)
    summary = json.loads((HERE / f'{tag}_summary.json').read_text())
    return dict(t=t, data=data, mesh=mesh, rp=rp, rot=rot, vcom=vcom,
                fnorm=force.Fwall_norm_N.to_numpy(), summary=summary)


def main():
    wall = Wall()
    ref, cand = load_case(20), load_case(50)
    faces = pd.read_csv(ROOT / 'Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083_robot_surface_triangles_exact.csv')
    face_ids = faces[['n1', 'n2', 'n3']].to_numpy(int)
    node_ids = np.unique(face_ids)
    frame_times = np.unique(np.r_[np.linspace(.000995, .001000, 5),
                                  np.arange(.001000, .0010101, .0000002),
                                  np.linspace(.001012, .001040, 8)])
    frame_indices = np.rint(frame_times / 1e-7).astype(int)

    first_i = int(round(ref['summary']['first_contact_s'] / 1e-7))
    nodes0 = np.array([ref['mesh']['nodes'][int(i)] for i in node_ids]) + ref['mesh']['shift']
    pos0 = ref['rp'] + ref['data']['U'][first_i] + ref['rot'][first_i].apply(nodes0 - ref['rp'])
    gap0, tri0, _ = wall.query(pos0)
    contact_center = wall.centers[tri0[np.argmin(gap0)]]
    wall_polys = wall.tri[np.linalg.norm(wall.centers - contact_center, axis=1) < 3.2]
    cloud = np.vstack((wall_polys.reshape(-1, 3), pos0))
    center = cloud.mean(axis=0)
    half = max(np.ptp(cloud, axis=0).max() * .58, 2.0)

    fig = plt.figure(figsize=(7.2047, 3.1496), constrained_layout=True)
    axes = [fig.add_subplot(1, 2, i + 1, projection='3d') for i in range(2)]

    def draw(frame_number):
        i = frame_indices[frame_number]
        for ax, case, title, color in zip(
                axes, [ref, cand], ['Reference: ζ=0.20, tangent=0', 'Candidate: ζ=0.50, tangent=0'],
                ['#5B6573', '#D55E00']):
            ax.clear()
            mesh = case['mesh']
            transformed = {
                label: case['rp'] + case['data']['U'][i]
                + case['rot'][i].apply(np.asarray(coord) + mesh['shift'] - case['rp'])
                for label, coord in mesh['nodes'].items()}
            robot_polys = np.array([[transformed[int(label)] for label in row] for row in face_ids])
            surface_pos = np.array([transformed[int(label)] for label in node_ids])
            gap, triangle, _ = wall.query(surface_pos)
            k = int(np.argmin(gap))
            node, tri = surface_pos[k], int(triangle[k])
            normal = wall.normals[tri]
            xcom = (case['rp'] + case['data']['U'][i]
                    + case['rot'][i].apply(mesh['com'] - case['rp']))
            vcontact = case['vcom'][i] + np.cross(case['data']['VR'][i], node - xcom)
            vn = float(np.dot(vcontact, normal))
            ax.add_collection3d(Poly3DCollection(wall_polys, facecolor='#56B4E9', edgecolor='none', alpha=.10))
            ax.add_collection3d(Poly3DCollection(robot_polys, facecolor=color,
                                                 edgecolor='#30343B', linewidth=.15, alpha=.82))
            if case['fnorm'][i] > 1e-8:
                ax.add_collection3d(Poly3DCollection([wall.tri[tri]], facecolor='#F0E442',
                                                     edgecolor='#9A7D00', linewidth=.8, alpha=.85))
            ax.scatter(*node, s=18, color='#CC3311', depthshade=False)
            ax.quiver(*node, *(normal * .55), color='#0072B2', linewidth=1.2, arrow_length_ratio=.18)
            ax.quiver(*node, *(vcontact * 4e-4), color='#009E73', linewidth=1.2, arrow_length_ratio=.18)
            ax.set(xlim=(center[0] - half, center[0] + half),
                   ylim=(center[1] - half, center[1] + half),
                   zlim=(center[2] - half, center[2] + half))
            ax.set_box_aspect((1, 1, 1))
            ax.view_init(elev=19, azim=-57)
            ax.set_title(title, pad=2)
            ax.set_xlabel('x (mm)', labelpad=-3)
            ax.set_ylabel('y (mm)', labelpad=-3)
            ax.set_zlabel('z (mm)', labelpad=-4)
            separated = case['t'][i] > case['summary']['first_contact_end_s']
            restitution = f'\ne_n={case["summary"]["e_n"]:.3f}' if separated else ''
            ax.text2D(.02, .97,
                      f't={case["t"][i]*1e3:.4f} ms\ngap={gap[k]*1e3:+.3f} µm\n'
                      f'|Fwall|={case["fnorm"][i]:.3f} N\nv_n={vn:+.1f} mm s^-1{restitution}',
                      transform=ax.transAxes, va='top')
        fig.suptitle('Reduced-Hydro first impact | fixed camera and one-variable comparison', fontsize=8)

    draw(int(np.argmin(np.abs(frame_indices - first_i))))
    poster = HERE / 'ReducedHydro_FirstImpact_NormalDamp020_vs_050_poster'
    fig.savefig(str(poster) + '.png', dpi=600, bbox_inches='tight')
    fig.savefig(str(poster) + '.tiff', dpi=600, bbox_inches='tight')
    fig.savefig(str(poster) + '.pdf', bbox_inches='tight')
    fig.savefig(str(poster) + '.svg', bbox_inches='tight')
    animation = FuncAnimation(fig, draw, frames=len(frame_indices), interval=90, repeat=True)
    animation.save(HERE / 'ReducedHydro_FirstImpact_NormalDamp020_vs_050.gif',
                   writer=PillowWriter(fps=10))
    plt.close(fig)


if __name__ == '__main__':
    main()
