"""Python-only figures for continuous local-field validation."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import pandas as pd

from continuous_segment_reference import ContinuousSegmentReference
from regress_production_local_frame import DXF, TRANSFORM, VENDORED, build, load_module, RP0


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
BLUE = "#0072B2"
ORANGE = "#D55E00"
GRAY = "#65717E"
GREEN = "#009E73"

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


def continuity_figure():
    data = pd.read_csv(OUT / "temporal_continuity.csv")
    time_ms = data.time_s.to_numpy() * 1e3
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 5.1), sharex=True, constrained_layout=True)
    axes[0].plot(time_ms, data.old_B_step_deg, color=GRAY, lw=0.7, label="Old nearest vertex")
    axes[0].plot(time_ms, data.new_B_step_deg, color=BLUE, lw=0.8, label="Continuous segment")
    axes[0].set(ylabel="B direction step (deg)", title="Continuous projection suppresses vertex-switch field jumps")
    axes[0].legend(ncol=2, loc="upper right")
    axes[1].plot(time_ms, data.old_t_step_deg, color=GRAY, lw=0.7)
    axes[1].plot(time_ms, data.new_t_step_deg, color=BLUE, lw=0.8)
    axes[1].set(ylabel="Tangent step (deg)")
    axes[2].plot(time_ms, data.old_s_mm, color=GRAY, lw=0.8)
    axes[2].plot(time_ms, data.new_s_mm, color=BLUE, lw=0.9)
    axes[2].set(xlabel="Time (ms)", ylabel="Projected s (mm)")
    for letter, axis in zip("abc", axes):
        axis.text(-0.065, 1.03, letter, transform=axis.transAxes, fontsize=8, fontweight="bold")
    save(fig, "Continuous_vs_nearest_vertex_temporal_continuity")


def field_gif():
    data = pd.read_csv(OUT / "field_one_cycle.csv")
    production = load_module("field_gif_production", VENDORED)
    model, _, _ = build(production, field_mode="ROBOT_LOCAL_TANGENT")
    reference = ContinuousSegmentReference(DXF, TRANSFORM, model.robot_axis_global,
                                            polarity=model.robot_polarity)
    curve = reference.origin_aba + reference.curve.dot(reference.R)
    center = data.loc[0, ["center_vx", "center_vy", "center_vz"]].to_numpy(float)
    near = curve[np.linalg.norm(curve - center, axis=1) < 6.0]
    frame_rows = np.unique(np.rint(np.linspace(0, len(data) - 1, 121)).astype(int))
    fig = plt.figure(figsize=(7.2, 3.35), constrained_layout=True)
    global_axis = fig.add_subplot(1, 2, 1, projection="3d")
    local_axis = fig.add_subplot(1, 2, 2)

    def draw(frame_number):
        index = frame_rows[frame_number]
        row = data.iloc[index]
        tangent = row[["t_x", "t_y", "t_z"]].to_numpy(float)
        e1 = row[["e1_x", "e1_y", "e1_z"]].to_numpy(float)
        e2 = row[["e2_x", "e2_y", "e2_z"]].to_numpy(float)
        b = row[["Bx_T", "By_T", "Bz_T"]].to_numpy(float)
        global_axis.clear(); local_axis.clear()
        global_axis.plot(near[:, 0], near[:, 1], near[:, 2], color=GRAY, lw=1.2)
        global_axis.scatter(*RP0, color="black", s=20, label="Robot COM", depthshade=False)
        global_axis.scatter(*center, color=GREEN, s=18, label="Projected center", depthshade=False)
        global_axis.plot([RP0[0], center[0]], [RP0[1], center[1]], [RP0[2], center[2]],
                         color=GRAY, ls="--", lw=0.8)
        global_axis.quiver(*center, *(1.6 * tangent), color=BLUE, linewidth=1.4,
                           arrow_length_ratio=0.16, label="Local tangent")
        global_axis.quiver(*center, *(1.6 * b / np.linalg.norm(b)), color=ORANGE,
                           linewidth=1.5, arrow_length_ratio=0.16, label="B")
        global_axis.set(xlim=(center[0] - 2.7, center[0] + 2.7),
                        ylim=(center[1] - 2.7, center[1] + 2.7),
                        zlim=(center[2] - 2.7, center[2] + 2.7))
        global_axis.set_box_aspect((1, 1, 1)); global_axis.view_init(elev=22, azim=-56)
        global_axis.set_title("Global curved-tube view")
        global_axis.set_xlabel("x (mm)"); global_axis.set_ylabel("y (mm)"); global_axis.set_zlabel("z (mm)")
        global_axis.legend(loc="upper left", fontsize=6)
        upto = data.iloc[:index + 1]
        trace1 = 1e3 * np.einsum("ij,ij->i", upto[["Bx_T", "By_T", "Bz_T"]].to_numpy(),
                                 upto[["e1_x", "e1_y", "e1_z"]].to_numpy())
        trace2 = 1e3 * np.einsum("ij,ij->i", upto[["Bx_T", "By_T", "Bz_T"]].to_numpy(),
                                 upto[["e2_x", "e2_y", "e2_z"]].to_numpy())
        b1, b2 = 1e3 * np.dot(b, e1), 1e3 * np.dot(b, e2)
        circle = np.linspace(0, 2 * np.pi, 361)
        local_axis.plot(5 * np.cos(circle), 5 * np.sin(circle), color="#B8C0C8", lw=0.8, ls="--")
        local_axis.plot(trace1, trace2, color=ORANGE, lw=1.2)
        local_axis.arrow(0, 0, b1, b2, width=0.045, head_width=0.32, color=ORANGE,
                         length_includes_head=True)
        local_axis.axhline(0, color="#D6DADF", lw=0.6); local_axis.axvline(0, color="#D6DADF", lw=0.6)
        local_axis.set(xlim=(-5.7, 5.7), ylim=(-5.7, 5.7), aspect="equal",
                       xlabel="B dot e1 (mT)", ylabel="B dot e2 (mT)",
                       title="View along local tangent")
        local_axis.text(0.03, 0.97, f"t = {row.time_s * 1e3:.3f} ms\n"
                        f"command phase = {(248.0 + row.time_s * 30 * 360) % 360:.1f} deg\n"
                        "|B| = 10.000 mT\nangle(B,t) = 30.000 deg",
                        transform=local_axis.transAxes, va="top")
        fig.suptitle("Production ROBOT_LOCAL_TANGENT | one 30 Hz cycle", fontsize=9)
        return []

    draw(0)
    fig.savefig(OUT / "Production_LocalTangent_Field_OneCycle_preview.png", dpi=300, bbox_inches="tight")
    animation = FuncAnimation(fig, draw, frames=len(frame_rows), interval=55, repeat=True)
    animation.save(OUT / "Production_LocalTangent_Field_OneCycle.gif", writer=PillowWriter(fps=18))
    plt.close(fig)


def main():
    continuity_figure()
    field_gif()


if __name__ == "__main__":
    main()
