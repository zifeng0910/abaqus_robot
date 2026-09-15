"""Publication figures and synchronized GIFs for the three-cycle result."""
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
sys.path[:0] = [str(REPO / "calibration_analysis" / "ReducedHydro_hidden_impact_audit")]
from exact_gap_audit import Wall

BLUE = "#0072B2"; ORANGE = "#D55E00"; GREEN = "#009E73"
GRAY = "#65717E"; RED = "#CC3311"; PURPLE = "#882255"
PERIOD = 1.0 / 30.0
COM_X = 1.25606694407031

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.linewidth": 0.8,
    "axes.spines.right": False, "axes.spines.top": False,
    "legend.frameon": False, "svg.fonttype": "none", "pdf.fonttype": 42,
})


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def orbit_figure(data):
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.65), sharex=True, sharey=True, constrained_layout=True)
    colors = (GRAY, BLUE, ORANGE)
    for cycle, axis_plot, color in zip((1, 2, 3), axes, colors):
        subset = data[data.cycle == cycle]
        axis_plot.plot(subset.q1, subset.q2, color=color, lw=1.15)
        axis_plot.scatter(subset.q1.iloc[0], subset.q2.iloc[0], s=14, color="black", zorder=3, label="Start")
        axis_plot.scatter(subset.q1.iloc[-1], subset.q2.iloc[-1], s=18, facecolor="white", edgecolor=color, zorder=3, label="End")
        axis_plot.axhline(0, color="#D6DADF", lw=0.6); axis_plot.axvline(0, color="#D6DADF", lw=0.6)
        axis_plot.set(xlim=(-1.02, 1.02), ylim=(-1.02, 1.02), aspect="equal",
                      xlabel="a dot e1", title=f"Cycle {cycle}")
    axes[0].set_ylabel("a dot e2")
    axes[0].legend(loc="lower left", fontsize=6)
    fig.suptitle("Robot-axis local orbit | identical axes", fontsize=9)
    save(fig, "Cycle1_Cycle2_Cycle3_axis_orbit")


def cycle_gif(data):
    phase_grid = np.linspace(0.0, 2.0 * np.pi, 121)
    cycles = []
    for cycle in (1, 2, 3):
        subset = data[data.cycle == cycle].copy()
        phase = np.unwrap(subset.B_phase_rad.to_numpy())
        phase -= phase[0]
        order = np.argsort(phase)
        phase = phase[order]
        q1 = subset.q1.to_numpy()[order]; q2 = subset.q2.to_numpy()[order]
        unique = np.r_[True, np.diff(phase) > 1.0e-12]
        cycles.append((np.interp(phase_grid, phase[unique], q1[unique]),
                       np.interp(phase_grid, phase[unique], q2[unique])))
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.65), sharex=True, sharey=True, constrained_layout=True)
    colors = (GRAY, BLUE, ORANGE)

    def draw(frame):
        for cycle, axis_plot, color, (q1, q2) in zip((1, 2, 3), axes, colors, cycles):
            axis_plot.clear()
            axis_plot.plot(q1[:frame + 1], q2[:frame + 1], color=color, lw=1.2)
            axis_plot.scatter(q1[frame], q2[frame], s=24, color=color)
            axis_plot.axhline(0, color="#D6DADF", lw=0.6); axis_plot.axvline(0, color="#D6DADF", lw=0.6)
            axis_plot.set(xlim=(-1.02, 1.02), ylim=(-1.02, 1.02), aspect="equal",
                          xlabel="a dot e1", title=f"Cycle {cycle}")
        axes[0].set_ylabel("a dot e2")
        fig.suptitle(f"Phase-normalized local orbit | B phase {np.degrees(phase_grid[frame]):.0f} deg", fontsize=9)
        return []

    draw(0)
    animation = FuncAnimation(fig, draw, frames=len(phase_grid), interval=55, repeat=True)
    animation.save(OUT / "Cycle1_vs_Cycle2_vs_Cycle3_LocalOrbit.gif", writer=PillowWriter(fps=18))
    plt.close(fig)


def dual_view_gif(data):
    gap = pd.read_csv(OUT / "exact_gap_timeseries.csv")
    wall = Wall()
    frame_index = np.unique(np.rint(np.linspace(0, len(data) - 1, 181)).astype(int))
    centers = data[["center_x", "center_y", "center_z"]].to_numpy(float)
    com_all = data[["com_x", "com_y", "com_z"]].to_numpy(float)
    selected_center = centers[frame_index]
    lumen_radius, _, _ = wall.query(selected_center)
    midpoint = np.mean(com_all, axis=0)
    wall_mask = np.linalg.norm(wall.centers - midpoint, axis=1) < 7.5
    wall_polys = wall.tri[wall_mask]
    global_min = np.min(np.vstack((com_all, centers)), axis=0) - 1.8
    global_max = np.max(np.vstack((com_all, centers)), axis=0) + 1.8
    half = max(float(np.max(global_max - global_min)) / 2.0, 3.0)
    global_center = 0.5 * (global_min + global_max)
    fig = plt.figure(figsize=(7.2, 3.5), constrained_layout=True)
    left = fig.add_subplot(1, 2, 1, projection="3d")
    right = fig.add_subplot(1, 2, 2)

    def draw(frame):
        index = frame_index[frame]
        row = data.iloc[index]
        com = row[["com_x", "com_y", "com_z"]].to_numpy(float)
        center = row[["center_x", "center_y", "center_z"]].to_numpy(float)
        axis = row[["axis_x", "axis_y", "axis_z"]].to_numpy(float)
        tangent = row[["t_x", "t_y", "t_z"]].to_numpy(float)
        e1 = row[["e1_x", "e1_y", "e1_z"]].to_numpy(float)
        e2 = row[["e2_x", "e2_y", "e2_z"]].to_numpy(float)
        b = row[["Bx_T", "By_T", "Bz_T"]].to_numpy(float)
        head = com - COM_X * axis; tail = com + (2.4 - COM_X) * axis
        left.clear(); right.clear()
        left.add_collection3d(Poly3DCollection(wall_polys, facecolor="#76B7D5", edgecolor="none", alpha=0.08))
        left.plot([head[0], tail[0]], [head[1], tail[1]], [head[2], tail[2]], color=BLUE, lw=5, alpha=0.8)
        left.scatter(*head, color=RED, s=20, depthshade=False, label="HEAD")
        left.scatter(*tail, color=GREEN, s=20, depthshade=False, label="TAIL")
        left.scatter(*center, color="black", s=10, depthshade=False)
        for vector, color, label, scale in ((tangent, GRAY, "t", 0.9), (e1, PURPLE, "e1", 0.6),
                                             (e2, GREEN, "e2", 0.6), (b / np.linalg.norm(b), ORANGE, "B", 0.9)):
            left.quiver(*center, *(scale * vector), color=color, linewidth=1.1, arrow_length_ratio=0.18, label=label)
        left.set(xlim=(global_center[0] - half, global_center[0] + half),
                 ylim=(global_center[1] - half, global_center[1] + half),
                 zlim=(global_center[2] - half, global_center[2] + half))
        left.set_box_aspect((1, 1, 1)); left.view_init(elev=20, azim=-56)
        left.set_title("Fixed global tube view"); left.set_xlabel("x (mm)"); left.set_ylabel("y (mm)"); left.set_zlabel("z (mm)")
        left.legend(loc="upper left", fontsize=5, ncol=2)
        theta = np.linspace(0, 2 * np.pi, 361); radius = abs(float(lumen_radius[frame]))
        right.plot(radius * np.cos(theta), radius * np.sin(theta), color=GRAY, lw=1.1)
        projected = []
        for point in (head, tail):
            delta = point - center
            projected.append((np.dot(delta, e1), np.dot(delta, e2)))
        right.plot([projected[0][0], projected[1][0]], [projected[0][1], projected[1][1]], color=BLUE, lw=4, alpha=0.8)
        right.scatter(*projected[0], color=RED, s=22); right.scatter(*projected[1], color=GREEN, s=22)
        right.arrow(0, 0, 0.55 * np.dot(b / np.linalg.norm(b), e1), 0.55 * np.dot(b / np.linalg.norm(b), e2),
                    width=0.008, head_width=0.06, color=ORANGE, length_includes_head=True)
        right.arrow(0, 0, 0.55 * row.q1, 0.55 * row.q2, width=0.006, head_width=0.05,
                    color=BLUE, length_includes_head=True)
        right.axhline(0, color="#D6DADF", lw=0.6); right.axvline(0, color="#D6DADF", lw=0.6)
        nearest = int(np.argmin(np.abs(gap.time_s.to_numpy() - row.time_s)))
        sector = gap.nearest_wall_sector_deg.iloc[nearest]
        right.plot([0, radius * np.cos(np.radians(sector))], [0, radius * np.sin(np.radians(sector))],
                   color=RED if row.contact_active else GRAY, lw=0.9, ls=":" )
        right.set(xlim=(-0.85, 0.85), ylim=(-0.85, 0.85), aspect="equal",
                  xlabel="e1 (mm)", ylabel="e2 (mm)", title="COM-following local cross-section")
        right.text(0.02, 0.98, f"Cycle {int(row.cycle)} | t={row.time_s * 1e3:.2f} ms\n"
                   f"robot phase={np.degrees(row.robot_phase_rad) % 360:.1f} deg\n"
                   f"B phase={np.degrees(row.B_phase_rad) % 360:.1f} deg\n"
                   f"tilt={row.tilt_deg:.1f} deg | contact={'ON' if row.contact_active else 'off'}\n"
                   f"sector={sector:.0f} deg | delta s={row.delta_s_mm:+.3f} mm",
                   transform=right.transAxes, va="top", fontsize=6)
        fig.suptitle("PROD_LOCAL30_G0_3CYCLE | production local-tangent field", fontsize=9)
        return []

    draw(0)
    fig.savefig(OUT / "PROD_LOCAL30_G0_3CYCLE_DualView_preview.png", dpi=300, bbox_inches="tight")
    animation = FuncAnimation(fig, draw, frames=len(frame_index), interval=65, repeat=True)
    animation.save(OUT / "PROD_LOCAL30_G0_3CYCLE_DualView.gif", writer=PillowWriter(fps=15))
    plt.close(fig)


def main():
    data = pd.read_csv(OUT / "three_cycle_pose.csv")
    orbit_figure(data)
    cycle_gif(data)
    dual_view_gif(data)


if __name__ == "__main__":
    main()
