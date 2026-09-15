"""Python-only field-frame figures and one-cycle production-field GIF."""
from pathlib import Path
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import pandas as pd

from audit_common import AUDIT, CURVE_CSV, RP0


DATA = AUDIT / "field_only"
FIGURES = AUDIT / "figures"
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 7, "axes.linewidth": 0.7, "axes.spines.right": False,
                     "axes.spines.top": False, "legend.frameon": False,
                     "svg.fonttype": "none", "pdf.fonttype": 42})
BLUE, ORANGE, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#6B7278"


def save(fig, name):
    fig.savefig(FIGURES / (name + ".png"), dpi=600, bbox_inches="tight")
    fig.savefig(FIGURES / (name + ".svg"), bbox_inches="tight")
    fig.savefig(FIGURES / (name + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def static_figures():
    current = pd.read_csv(DATA / "current_production_field_one_cycle.csv")
    local30 = pd.read_csv(DATA / "local_tangent_field_30deg_one_cycle.csv")
    local40 = pd.read_csv(DATA / "local_tangent_field_40deg_one_cycle.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45), constrained_layout=True)
    for label, table, color, ax in zip(("Current production", "Local tangent 30 deg", "Local tangent 40 deg"),
                                        (current, local30, local40), (ORANGE, BLUE, GREEN), axes):
        ax.plot(table.B_e1_T * 1000, table.B_e2_T * 1000, color=color, lw=1.2)
        ax.scatter(table.B_e1_T.iloc[0] * 1000, table.B_e2_T.iloc[0] * 1000, s=22, color="black", label="Start")
        ax.scatter(table.B_e1_T.iloc[-1] * 1000, table.B_e2_T.iloc[-1] * 1000, s=24, facecolors="none", edgecolors=color, label="End")
        ax.axhline(0, color="#D3D7DA", lw=.6); ax.axvline(0, color="#D3D7DA", lw=.6)
        ax.set_aspect("equal"); ax.set(xlabel="B e1 (mT)", ylabel="B e2 (mT)", title=label)
    axes[0].legend(fontsize=6)
    for letter, ax in zip("abc", axes):
        ax.text(-.18, 1.08, letter, transform=ax.transAxes, fontweight="bold", fontsize=8)
    save(fig, "field_tip_local_orbits")

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 3.7), sharex=True, constrained_layout=True)
    time_ms = current.time_s * 1000
    axes[0].plot(time_ms, current.theta_B_t_deg, color=ORANGE, label="Current production")
    axes[0].plot(time_ms, local30.theta_B_t_deg, color=BLUE, label="Local tangent 30 deg")
    axes[0].plot(time_ms, local40.theta_B_t_deg, color=GREEN, label="Local tangent 40 deg")
    axes[0].set(ylabel="Field-to-tube angle (deg)"); axes[0].legend(ncol=3, fontsize=6)
    for label, table, color in (("Current production", current, ORANGE), ("Local tangent 30 deg", local30, BLUE), ("Local tangent 40 deg", local40, GREEN)):
        winding = (table.psi_B_unwrapped_rad - table.psi_B_unwrapped_rad.iloc[0]) / (2 * math.pi)
        axes[1].plot(time_ms, winding, color=color, label=label)
    axes[1].plot(time_ms, time_ms / time_ms.iloc[-1], color=GREY, lw=.8, ls="--", label="30 Hz command")
    axes[1].set(xlabel="Time (ms)", ylabel="Local azimuth winding")
    for letter, ax in zip("ab", axes): ax.text(-.08, 1.05, letter, transform=ax.transAxes, fontweight="bold", fontsize=8)
    save(fig, "field_local_frame_diagnostics")


def field_gif():
    table = pd.read_csv(DATA / "current_production_field_one_cycle.csv")
    curve = pd.read_csv(CURVE_CSV)[["x_mm", "y_mm", "z_mm"]].to_numpy(float)
    near = curve[np.linalg.norm(curve - RP0, axis=1) < 6.0]
    indices = np.linspace(0, len(table) - 1, 121).astype(int)
    fig = plt.figure(figsize=(7.2, 3.3), dpi=int("133"))
    ax = fig.add_subplot(1, 2, 1, projection="3d"); bx = fig.add_subplot(1, 2, 2)
    scale = 250.0
    def update(frame_index):
        i = indices[frame_index]; row = table.iloc[i]
        ax.cla(); bx.cla()
        ax.plot(near[:, 0], near[:, 1], near[:, 2], color="#AEB5BA", lw=5, alpha=.4)
        b = row[["Bx_T", "By_T", "Bz_T"]].to_numpy(float)
        ax.quiver(*RP0, *(scale * b), color=ORANGE, linewidth=2)
        ax.scatter(*RP0, color="black", s=16)
        bounds = np.vstack([near, RP0 - 3, RP0 + 3]); lo=bounds.min(axis=0); hi=bounds.max(axis=0)
        ax.set(xlim=(lo[0],hi[0]), ylim=(lo[1],hi[1]), zlim=(lo[2],hi[2])); ax.set_box_aspect(hi-lo)
        ax.view_init(elev=23, azim=-58); ax.set_axis_off(); ax.set_title("Global tube and production B")
        trace = table.iloc[:i+1]
        bx.plot(trace.B_e1_T * 1000, trace.B_e2_T * 1000, color=ORANGE, lw=1)
        bx.scatter(row.B_e1_T * 1000, row.B_e2_T * 1000, color=ORANGE, s=25)
        bx.axhline(0,color="#D3D7DA",lw=.6); bx.axvline(0,color="#D3D7DA",lw=.6)
        bx.set_aspect("equal"); bx.set(xlim=(-11,11),ylim=(-11,11),xlabel="B e1 (mT)",ylabel="B e2 (mT)",title="View down local tube axis")
        fig.suptitle("Current production field | t=%.3f ms | command phase=%.1f deg | local winding=%+.3f" %
                     (row.time_s*1000, math.degrees(row.command_phase_rad)%360,
                      (row.psi_B_unwrapped_rad-table.psi_B_unwrapped_rad.iloc[0])/(2*math.pi)), fontsize=8)
        return []
    animation = FuncAnimation(fig, update, frames=len(indices), interval=55, blit=False)
    animation.save(AUDIT / "CurrentProductionField_LocalTubeFrame_1cycle.gif", writer=PillowWriter(fps=18))
    plt.close(fig)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    static_figures(); field_gif()


if __name__ == "__main__":
    main()
