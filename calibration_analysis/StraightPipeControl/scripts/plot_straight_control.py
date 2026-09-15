"""Render the straight-control diagnostics and curved-versus-straight comparison."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
BASE = OUT.parent / "ProductionLocalFrameValidation"
BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
GRAY = "#65717E"
LIGHT = "#D6DADF"
RED = "#CC3311"
PERIOD = 1.0 / 30.0

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.linewidth": 0.8,
    "axes.spines.right": False, "axes.spines.top": False,
    "legend.frameon": False, "svg.fonttype": "none", "pdf.fonttype": 42,
})


def frame_indices(length, count=181):
    return np.unique(np.rint(np.linspace(0, length - 1, count)).astype(int))


def local_projection(row, prefix):
    center = row[["center_x", "center_y", "center_z"]].to_numpy(float)
    e1 = row[["e1_x", "e1_y", "e1_z"]].to_numpy(float)
    e2 = row[["e2_x", "e2_y", "e2_z"]].to_numpy(float)
    point = row[[f"{prefix}_x", f"{prefix}_y", f"{prefix}_z"]].to_numpy(float) - center
    return np.dot(point, e1), np.dot(point, e2)


def draw_local(ax, row, sector, radius, title, trail=None):
    theta = np.linspace(0.0, 2.0 * np.pi, 361)
    ax.plot(radius * np.cos(theta), radius * np.sin(theta), color=GRAY, lw=1.1)
    head = local_projection(row, "head")
    tail = local_projection(row, "tail")
    ax.plot([head[0], tail[0]], [head[1], tail[1]], color=BLUE, lw=4.5, alpha=0.82)
    ax.scatter(*head, color=RED, s=22, zorder=3, label="Head")
    ax.scatter(*tail, color=GREEN, s=22, zorder=3, label="Tail")
    b = row[["Bx_T", "By_T", "Bz_T"]].to_numpy(float)
    b /= max(np.linalg.norm(b), 1.0e-15)
    e1 = row[["e1_x", "e1_y", "e1_z"]].to_numpy(float)
    e2 = row[["e2_x", "e2_y", "e2_z"]].to_numpy(float)
    ax.arrow(0, 0, 0.55 * np.dot(b, e1), 0.55 * np.dot(b, e2), width=0.008,
             head_width=0.06, color=ORANGE, length_includes_head=True, label="B")
    ax.arrow(0, 0, 0.55 * row.q1, 0.55 * row.q2, width=0.006,
             head_width=0.05, color=BLUE, length_includes_head=True, label="Axis")
    if trail is not None and len(trail):
        ax.plot(trail[:, 0], trail[:, 1], color=BLUE, lw=0.8, alpha=0.45)
    ax.plot([0, radius * np.cos(np.radians(sector))], [0, radius * np.sin(np.radians(sector))],
            color=RED if row.contact_active else GRAY, lw=0.9, ls=":")
    ax.axhline(0, color=LIGHT, lw=0.6); ax.axvline(0, color=LIGHT, lw=0.6)
    ax.set(xlim=(-0.85, 0.85), ylim=(-0.85, 0.85), aspect="equal",
           xlabel="e1 (mm)", ylabel="e2 (mm)", title=title)


def dual_view_gif(data, gap, identity):
    selected = frame_indices(len(data))
    radius = float(identity["wall_radius_mm"])
    tangent = np.asarray(identity["straight_tangent_aba"], float)
    e1 = np.asarray(identity["initial_e1_aba"], float)
    e2 = np.asarray(identity["initial_e2_aba"], float)
    center = np.asarray(identity["initial_center_aba_mm"], float)
    axial = np.linspace(-4.0, 4.0, 42)
    theta = np.linspace(0.0, 2.0 * np.pi, 32)
    surface = (center[None, None, :] + axial[:, None, None] * tangent +
               radius * (np.cos(theta)[None, :, None] * e1 + np.sin(theta)[None, :, None] * e2))
    fig = plt.figure(figsize=(7.2, 3.45), constrained_layout=True)
    left = fig.add_subplot(1, 2, 1, projection="3d")
    right = fig.add_subplot(1, 2, 2)

    def draw(frame):
        index = selected[frame]
        row = data.iloc[index]
        nearest = int(np.argmin(np.abs(gap.time_s.to_numpy() - row.time_s)))
        sector = float(gap.nearest_wall_sector_deg.iloc[nearest])
        left.clear(); right.clear()
        left.plot_surface(surface[:, :, 0], surface[:, :, 1], surface[:, :, 2],
                          color="#76B7D5", alpha=0.10, linewidth=0, shade=False)
        head = row[["head_x", "head_y", "head_z"]].to_numpy(float)
        tail = row[["tail_x", "tail_y", "tail_z"]].to_numpy(float)
        local_center = row[["center_x", "center_y", "center_z"]].to_numpy(float)
        b = row[["Bx_T", "By_T", "Bz_T"]].to_numpy(float); b /= np.linalg.norm(b)
        left.plot([head[0], tail[0]], [head[1], tail[1]], [head[2], tail[2]], color=BLUE, lw=5)
        left.scatter(*head, color=RED, s=20); left.scatter(*tail, color=GREEN, s=20)
        for vector, color, label, scale in ((tangent, GRAY, "t", 0.9), (e1, GREEN, "e1", 0.6),
                                             (e2, RED, "e2", 0.6), (b, ORANGE, "B", 0.9)):
            left.quiver(*local_center, *(scale * vector), color=color, linewidth=1.1,
                        arrow_length_ratio=0.18, label=label)
        cloud = np.vstack((surface.reshape(-1, 3), data[["head_x", "head_y", "head_z"]].to_numpy(),
                           data[["tail_x", "tail_y", "tail_z"]].to_numpy()))
        mid = 0.5 * (cloud.min(axis=0) + cloud.max(axis=0)); half = 4.4
        left.set(xlim=(mid[0]-half, mid[0]+half), ylim=(mid[1]-half, mid[1]+half), zlim=(mid[2]-half, mid[2]+half))
        left.set_box_aspect((1, 1, 1)); left.view_init(elev=20, azim=-56)
        left.set_title("Fixed global straight-tube view")
        left.set_xlabel("x (mm)"); left.set_ylabel("y (mm)"); left.set_zlabel("z (mm)")
        left.legend(loc="upper left", fontsize=5, ncol=2)
        trail = data.loc[:index, ["q1", "q2"]].to_numpy(float)
        draw_local(right, row, sector, radius, "COM-following local cross-section", trail=trail)
        right.text(0.02, 0.98,
                   f"t={row.time_s*1e3:.2f} ms | phase={row.cycle_fraction*360:.0f} deg\n"
                   f"tilt={row.tilt_deg:.1f} deg | contact={'ON' if row.contact_active else 'off'}\n"
                   f"sector={sector:.0f} deg | delta s={row.delta_s_mm:+.3f} mm",
                   transform=right.transAxes, va="top", fontsize=6)
        fig.suptitle("PROD_LOCAL30_G0_STRAIGHT_CTRL | one 30 Hz cycle", fontsize=9)
        return []

    draw(0)
    fig.savefig(OUT / "PROD_LOCAL30_G0_STRAIGHT_CTRL_DualView_preview.png", dpi=300, bbox_inches="tight")
    FuncAnimation(fig, draw, frames=len(selected), interval=65, repeat=True).save(
        OUT / "PROD_LOCAL30_G0_STRAIGHT_CTRL_DualView.gif", writer=PillowWriter(fps=15))
    plt.close(fig)


def local_orbit_gif(data, gap, identity):
    selected = frame_indices(len(data), 151)
    radius = float(identity["wall_radius_mm"])
    fig, ax = plt.subplots(figsize=(4.2, 3.8), constrained_layout=True)

    def draw(frame):
        index = selected[frame]
        row = data.iloc[index]
        nearest = int(np.argmin(np.abs(gap.time_s.to_numpy() - row.time_s)))
        sector = float(gap.nearest_wall_sector_deg.iloc[nearest])
        ax.clear()
        trail = data.loc[:index, ["q1", "q2"]].to_numpy(float)
        draw_local(ax, row, sector, radius, "Straight control | local orbit", trail=trail)
        ax.text(0.02, 0.98,
                f"t={row.time_s*1e3:.2f} ms | tilt={row.tilt_deg:.1f} deg\n"
                f"contact={'ON' if row.contact_active else 'off'} | delta s={row.delta_s_mm:+.3f} mm",
                transform=ax.transAxes, va="top", fontsize=6)
        return []

    draw(0)
    FuncAnimation(fig, draw, frames=len(selected), interval=65, repeat=True).save(
        OUT / "PROD_LOCAL30_G0_STRAIGHT_CTRL_LocalOrbit.gif", writer=PillowWriter(fps=15))
    plt.close(fig)


def phase_resample(data):
    command = data.time_s.to_numpy(float) / PERIOD * 2.0 * np.pi
    if np.any(np.diff(command) <= 0.0):
        raise ValueError("Animation interpolation requires strictly increasing time")
    phase = np.linspace(0.0, 2.0 * np.pi, 151)
    columns = {}
    for name in ("q1", "q2", "tilt_deg", "delta_s_mm", "contact_active"):
        columns[name] = np.interp(phase, command, data[name].to_numpy(float))
    return phase, columns


def comparison_gif(straight, curved):
    curved = curved[curved.cycle == 1].copy()
    phase, s = phase_resample(straight)
    _, c = phase_resample(curved)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35), sharex=True, sharey=True, constrained_layout=True)

    def draw(frame):
        for ax, values, title, color in ((axes[0], c, "Curved pipe | cycle 1", GRAY),
                                          (axes[1], s, "Straight pipe | control", BLUE)):
            ax.clear()
            ax.plot(values["q1"][:frame+1], values["q2"][:frame+1], color=color, lw=1.2)
            ax.scatter(values["q1"][frame], values["q2"][frame], color=RED if values["contact_active"][frame] > 0.5 else color, s=25)
            ax.axhline(0, color=LIGHT, lw=0.6); ax.axvline(0, color=LIGHT, lw=0.6)
            ax.set(xlim=(-0.75, 0.75), ylim=(-0.75, 0.75), aspect="equal",
                   xlabel="axis dot e1", title=title)
            ax.text(0.02, 0.98,
                    f"tilt={values['tilt_deg'][frame]:.1f} deg\n"
                    f"delta s={values['delta_s_mm'][frame]:+.3f} mm\n"
                    f"contact={'ON' if values['contact_active'][frame] > 0.5 else 'off'}",
                    transform=ax.transAxes, va="top", fontsize=6)
        axes[0].set_ylabel("axis dot e2")
        fig.suptitle(f"Curved versus straight | command phase {np.degrees(phase[frame]):.0f} deg", fontsize=9)
        return []

    draw(0)
    FuncAnimation(fig, draw, frames=len(phase), interval=65, repeat=True).save(
        OUT / "Curved_vs_Straight_Comparison.gif", writer=PillowWriter(fps=15))
    plt.close(fig)


def static_summary(straight, curved, metrics):
    curved = curved[curved.cycle == 1]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.55), constrained_layout=True)
    axes[0].plot(curved.time_s * 1e3, curved.tilt_deg, color=GRAY, lw=1.0, label="Curved")
    axes[0].plot(straight.time_s * 1e3, straight.tilt_deg, color=BLUE, lw=1.0, label="Straight")
    axes[0].set(xlabel="Time (ms)", ylabel="Directed tilt (deg)", title="a  Tilt response")
    axes[0].legend(fontsize=6)
    axes[1].plot(curved.q1, curved.q2, color=GRAY, lw=1.0, label="Curved")
    axes[1].plot(straight.q1, straight.q2, color=BLUE, lw=1.0, label="Straight")
    axes[1].axhline(0, color=LIGHT, lw=0.6); axes[1].axvline(0, color=LIGHT, lw=0.6)
    axes[1].set(xlabel="axis dot e1", ylabel="axis dot e2", title="b  Local-axis orbit", aspect="equal")
    curved_metrics = json.loads((BASE / "three_cycle_summary.json").read_text())["cycles"][0]
    labels = ["Axis\nwinding", "Mean\ntilt", "Contact\nevents", "20 um\nbridge"]
    curved_values = [1.0, 1.0, 1.0, 1.0]
    straight_values = [metrics["robot_axis_winding"] / curved_metrics["robot_axis_winding"],
                       metrics["tilt_mean_deg"] / curved_metrics["tilt_mean_deg"],
                       metrics["contact_event_count"] / curved_metrics["contact_event_count"],
                       metrics["longest_20um_bridge_ms"] / curved_metrics["longest_bridge_ms"]]
    if min(curved_values + straight_values) <= 0.0:
        raise ValueError("Log-ratio panel requires strictly positive metrics")
    x = np.arange(len(labels)); width = 0.36
    log_floor = 0.3
    axes[2].bar(x-width/2, np.asarray(curved_values)-log_floor, width, bottom=log_floor,
                color=GRAY, label="Curved")
    axes[2].bar(x+width/2, np.asarray(straight_values)-log_floor, width, bottom=log_floor,
                color=BLUE, label="Straight")
    axes[2].set_xticks(x, labels); axes[2].set_ylabel("Ratio to curved cycle 1")
    axes[2].set_yscale("log"); axes[2].set_ylim(0.3, 60.0)
    axes[2].set_yticks([0.5, 1.0, 10.0, 50.0], ["0.5", "1", "10", "50"])
    axes[2].axhline(1.0, color=LIGHT, lw=0.7)
    axes[2].set_title("c  Cycle-level comparison"); axes[2].legend(fontsize=6)
    fig.savefig(OUT / "Curved_vs_Straight_Summary.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / "Curved_vs_Straight_Summary.pdf", bbox_inches="tight")
    fig.savefig(OUT / "Curved_vs_Straight_Summary.svg", bbox_inches="tight")
    plt.close(fig)


def main():
    straight = pd.read_csv(OUT / "straight_pose.csv")
    gap = pd.read_csv(OUT / "straight_exact_gap_timeseries.csv")
    curved = pd.read_csv(BASE / "three_cycle_pose.csv")
    metrics = json.loads((OUT / "straight_metrics.json").read_text())
    identity = json.loads((OUT / "case" / "PROD_LOCAL30_G0_STRAIGHT_CTRL" / "case_identity.json").read_text())
    dual_view_gif(straight, gap, identity)
    local_orbit_gif(straight, gap, identity)
    comparison_gif(straight, curved)
    static_summary(straight, curved, metrics)


if __name__ == "__main__":
    main()
