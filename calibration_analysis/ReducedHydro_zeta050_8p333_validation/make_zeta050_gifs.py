"""Fixed-camera 3-D animations for the zeta=0.50 long validation."""
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


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT = REPO.parent
AUDIT = REPO / "calibration_analysis" / "ReducedHydro_hidden_impact_audit"
sys.path.insert(0, str(AUDIT))
from audit_stage_a import mesh_properties
from exact_gap_audit import Wall

JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_WallOn_Free_0083"
OLD_JOB = "Wobble_F30_G6L45_WallOn_Free_0083_WobbleSurvival"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 6,
    "xtick.labelsize": 5,
    "ytick.labelsize": 5,
})


def transform_polys(mesh, rp, displacement, rotation, face_ids):
    labels = np.unique(face_ids)
    base = {int(i): np.asarray(mesh["nodes"][int(i)]) + mesh["shift"] for i in labels}
    moved = {label: rp + displacement + rotation.apply(point - rp) for label, point in base.items()}
    return np.asarray([[moved[int(label)] for label in face] for face in face_ids])


def configure(ax, center, half, title):
    ax.set(xlim=(center[0] - half, center[0] + half),
           ylim=(center[1] - half, center[1] + half),
           zlim=(center[2] - half, center[2] + half))
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=19, azim=-57)
    ax.set_title(title, pad=2)
    ax.set_xlabel("x (mm)", labelpad=-3)
    ax.set_ylabel("y (mm)", labelpad=-3)
    ax.set_zlabel("z (mm)", labelpad=-4)


def main():
    field = pd.read_csv(HERE / "candidate_private" / "rp_fields.csv")
    axis = pd.read_csv(HERE / "zeta050_8p333_true_axis.csv")
    phase = pd.read_csv(HERE / "zeta050_8p333_local_phase.csv")
    windows = pd.read_csv(HERE / "zeta050_8p333_phase_windows.csv")
    gap = pd.read_csv(HERE / "zeta050_8p333_exact_gap.csv")
    bridge = pd.read_csv(HERE / "zeta050_8p333_bridge_timeline.csv")
    translation = pd.read_csv(HERE / "zeta050_8p333_canonical_translation.csv")
    events = pd.read_csv(HERE / "zeta050_8p333_contact_events.csv")
    centerline = pd.read_csv(ROOT / "CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_TRUE_centerline_odb.csv")
    face_table = pd.read_csv(ROOT / "Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083_robot_surface_triangles_exact.csv")
    face_ids = face_table[["n1", "n2", "n3"]].to_numpy(int)
    meshes, rp, _ = mesh_properties((HERE / f"{JOB}.inp").read_text())
    mesh = meshes["Robot_SOLID"]
    wall = Wall()

    frame_rows = np.unique(np.rint(np.linspace(0, len(field) - 1, 125)).astype(int))
    t_field = field.time_s.to_numpy(float)
    dense_index = np.clip(np.rint(t_field / 1e-7).astype(int), 0, len(axis) - 1)
    rotations = Rotation.from_rotvec(field[["UR1", "UR2", "UR3"]].to_numpy(float))
    initial_polys = transform_polys(mesh, rp, field.loc[0, ["U1", "U2", "U3"]].to_numpy(float), rotations[0], face_ids)
    center = initial_polys.reshape(-1, 3).mean(axis=0)
    wall_mask = np.linalg.norm(wall.centers - center, axis=1) < 4.0
    wall_polys = wall.tri[wall_mask]
    half = 2.35
    curve = centerline[["x_mm", "y_mm", "z_mm"]].to_numpy(float)
    curve_mask = np.linalg.norm(curve - center, axis=1) < 4.5
    curve = curve[curve_mask]

    event_intervals = events[["start_s", "end_s"]].to_numpy(float)

    fig = plt.figure(figsize=(5.2, 4.2), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")

    def draw_main(frame_number):
        row = frame_rows[frame_number]
        i = dense_index[row]
        t = t_field[row]
        ax.clear()
        robot = transform_polys(mesh, rp, field.loc[row, ["U1", "U2", "U3"]].to_numpy(float), rotations[row], face_ids)
        ax.add_collection3d(Poly3DCollection(wall_polys, facecolor="#56B4E9", edgecolor="none", alpha=0.10))
        ax.add_collection3d(Poly3DCollection(robot, facecolor="#0072B2", edgecolor="#30343B", linewidth=0.12, alpha=0.82))
        ax.plot(curve[:, 0], curve[:, 1], curve[:, 2], color="#5B6573", lw=0.9, ls="--")
        head = axis.loc[i, ["head_x_mm", "head_y_mm", "head_z_mm"]].to_numpy(float)
        tail = axis.loc[i, ["tail_x_mm", "tail_y_mm", "tail_z_mm"]].to_numpy(float)
        midpoint = 0.5 * (head + tail)
        tangent = axis.loc[i, ["tangent_x", "tangent_y", "tangent_z"]].to_numpy(float)
        e1 = axis.loc[i, ["e1_x", "e1_y", "e1_z"]].to_numpy(float)
        e2 = axis.loc[i, ["e2_x", "e2_y", "e2_z"]].to_numpy(float)
        robot_axis = axis.loc[i, ["axis_x", "axis_y", "axis_z"]].to_numpy(float)
        phi_b = phase.phi_B_local_rad.iloc[i]
        bvec = np.cos(phi_b) * e1 + np.sin(phi_b) * e2
        ax.scatter(*head, s=18, color="#D55E00", label="HEAD", depthshade=False)
        ax.scatter(*tail, s=18, color="#009E73", label="TAIL", depthshade=False)
        for vector, color, label, scale in ((tangent, "#5B6573", "Local tangent", 0.8),
                                             (e1, "#CC79A7", "Local e1", 0.6),
                                             (e2, "#F0E442", "Local e2", 0.6),
                                             (robot_axis, "#0072B2", "True axis", 1.0),
                                             (bvec, "#D55E00", "B", 0.8)):
            ax.quiver(*midpoint, *(vector * scale), color=color, linewidth=1.1, arrow_length_ratio=0.16, label=label)
        active = np.any((t >= event_intervals[:, 0]) & (t <= event_intervals[:, 1])) if len(event_intervals) else False
        g = np.interp(t, gap.time_s, gap.gap_um)
        rate = np.interp(t, windows.center_s, windows.robot_phase_rate_Hz)
        delta_s = np.interp(t, translation.time_s, translation.delta_s_mm)
        ax.text2D(0.02, 0.98,
                  f"t={t*1e3:.3f} ms\ntrue tilt={axis.tilt_deg.iloc[i]:.1f} deg\n"
                  f"local phase={np.degrees(phase.phi_robot_rad.iloc[i]):.1f} deg\n"
                  f"window rate={rate:.1f} Hz\nexact gap={g:+.2f} um\n"
                  f"contact={'ACTIVE' if active else 'off'} | bridge={'YES' if bridge.opposing_bridge.iloc[i] else 'no'}\n"
                  f"canonical delta s={delta_s:+.3f} mm",
                  transform=ax.transAxes, va="top")
        configure(ax, center, half, "Reduced-Hydro zeta=0.50 | fixed camera")

    draw_main(0)
    animation = FuncAnimation(fig, draw_main, frames=len(frame_rows), interval=80, repeat=True)
    animation.save(HERE / "Wobble_F30_ReducedHydro_Zeta050_WallOn_8p333.gif", writer=PillowWriter(fps=12))
    plt.close(fig)

    old_raw = json.loads((ROOT / f"{OLD_JOB}_rp_fields.json").read_text())["rows"]
    old_t = np.asarray([row["time"] for row in old_raw], float)
    old_u = np.asarray([row["U"] for row in old_raw], float)
    old_ur = Rotation.from_rotvec(np.asarray([row["UR"] for row in old_raw], float))
    sync_times = np.linspace(0, min(old_t[-1], t_field[-1]), 100)
    old_indices = np.searchsorted(old_t, sync_times).clip(0, len(old_t) - 1)
    new_indices = np.searchsorted(t_field, sync_times).clip(0, len(t_field) - 1)

    fig = plt.figure(figsize=(7.2, 3.4), constrained_layout=True)
    axes = [fig.add_subplot(1, 2, j + 1, projection="3d") for j in range(2)]

    def draw_comparison(frame_number):
        cases = ((axes[0], old_u[old_indices[frame_number]], old_ur[old_indices[frame_number]], "CEL", "#5B6573"),
                 (axes[1], field.loc[new_indices[frame_number], ["U1", "U2", "U3"]].to_numpy(float),
                  rotations[new_indices[frame_number]], "Reduced-Hydro zeta=0.50", "#0072B2"))
        for panel, displacement, rotation, title, color in cases:
            panel.clear()
            robot = transform_polys(mesh, rp, displacement, rotation, face_ids)
            panel.add_collection3d(Poly3DCollection(wall_polys, facecolor="#56B4E9", edgecolor="none", alpha=0.10))
            panel.add_collection3d(Poly3DCollection(robot, facecolor=color, edgecolor="#30343B", linewidth=0.12, alpha=0.82))
            panel.plot(curve[:, 0], curve[:, 1], curve[:, 2], color="#5B6573", lw=0.8, ls="--")
            panel.text2D(0.02, 0.97, f"t={sync_times[frame_number]*1e3:.3f} ms", transform=panel.transAxes, va="top")
            configure(panel, center, half, title)
        fig.suptitle("CEL versus Reduced-Hydro post-impact motion | synchronized fixed camera", fontsize=8)

    draw_comparison(0)
    animation = FuncAnimation(fig, draw_comparison, frames=len(sync_times), interval=85, repeat=True)
    animation.save(HERE / "CEL_vs_ReducedHydro_Zeta050_8p333.gif", writer=PillowWriter(fps=12))
    plt.close(fig)


if __name__ == "__main__":
    main()
