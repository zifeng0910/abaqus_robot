"""Generate the eight public validation figures from audited CSV outputs."""
from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
COLORS = {
    "robot": "#0072B2",
    "field": "#D55E00",
    "contact": "#CC3311",
    "bridge": "#882255",
    "hydro": "#009E73",
    "neutral": "#5B6573",
}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def save(fig, name):
    fig.savefig(HERE / f"{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(HERE / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(HERE / f"{name}.svg", bbox_inches="tight")
    fig.savefig(HERE / f"{name}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def event_spans(ax, events, alpha=0.12):
    for row in events.itertuples():
        ax.axvspan(row.start_s * 1e3, row.end_s * 1e3, color=COLORS["contact"], alpha=alpha, lw=0)


def main():
    phase = pd.read_csv(HERE / "zeta050_8p333_local_phase.csv")
    windows = pd.read_csv(HERE / "zeta050_8p333_phase_windows.csv")
    axis = pd.read_csv(HERE / "zeta050_8p333_true_axis.csv")
    gap = pd.read_csv(HERE / "zeta050_8p333_exact_gap.csv")
    events = pd.read_csv(HERE / "zeta050_8p333_contact_events.csv")
    bridge = pd.read_csv(HERE / "zeta050_8p333_bridge_timeline.csv")
    torque = pd.read_csv(HERE / "zeta050_8p333_magnetic_hydro_torque.csv")
    translation = pd.read_csv(HERE / "zeta050_8p333_canonical_translation.csv")
    summary = json.loads((HERE / "zeta050_8p333_summary.json").read_text())
    comp = np.load(HERE / "comparison_series_private.npz")

    time_ms = phase.time_s.to_numpy() * 1e3
    robot_advance = np.degrees(phase.phi_robot_rad - phase.phi_robot_rad.iloc[0])
    field_advance = np.degrees(phase.phi_B_local_rad - phase.phi_B_local_rad.iloc[0])

    fig, ax = plt.subplots(figsize=(3.54, 2.25), constrained_layout=True)
    ax.plot(time_ms, robot_advance, color=COLORS["robot"], lw=1.25, label="Robot true axis")
    ax.plot(time_ms, field_advance, color=COLORS["field"], lw=1.0, label="Local magnetic field")
    event_spans(ax, events)
    ax.set(xlabel="Time (ms)", ylabel="Unwrapped phase advance (deg)", title="Intermittent post-impact true-axis phase")
    ax.legend(loc="best")
    save(fig, "true_axis_phase_8p333")

    fig, axes = plt.subplots(2, 1, figsize=(3.54, 3.25), sharex=True, constrained_layout=True)
    axes[0].plot(time_ms, np.degrees(phase.phi_robot_minus_B_rad - phase.phi_robot_minus_B_rad.iloc[0]),
                 color=COLORS["robot"], lw=1.1)
    axes[0].set(ylabel="Robot - field (deg)", title="Local phase difference and short-window rate")
    axes[1].plot(windows.center_s * 1e3, windows.robot_phase_rate_Hz, color=COLORS["robot"], lw=1.1,
                 label="Robot")
    axes[1].plot(windows.center_s * 1e3, windows.field_local_phase_rate_Hz, color=COLORS["field"], lw=1.0,
                 label="Field")
    axes[1].axhline(30, color=COLORS["neutral"], ls="--", lw=0.8, label="30 Hz command")
    axes[1].set(xlabel="Time (ms)", ylabel="Window rate (Hz)")
    axes[1].legend(ncol=3, loc="best")
    save(fig, "phase_difference_8p333")

    fig, ax = plt.subplots(figsize=(3.54, 2.25), constrained_layout=True)
    ax.plot(axis.time_s * 1e3, axis.tilt_deg, color=COLORS["robot"], lw=1.2)
    event_spans(ax, events)
    ax.axhline(summary["late_tilt_deg"], color=COLORS["neutral"], ls="--", lw=0.8,
               label=f'Late tilt {summary["late_tilt_deg"]:.1f} deg')
    ax.set(xlabel="Time (ms)", ylabel="True tilt (deg)", title="Tilt settles near 33 deg after bridge formation")
    ax.legend(loc="best")
    save(fig, "tilt_angle_8p333")

    fig, axes = plt.subplots(2, 1, figsize=(3.54, 3.25), sharex=True, constrained_layout=True,
                             gridspec_kw={"height_ratios": [2.2, 1]})
    axes[0].plot(gap.time_s * 1e3, gap.gap_um, color=COLORS["neutral"], lw=0.75)
    axes[0].axhline(0, color="black", lw=0.7)
    axes[0].set_yscale("symlog", linthresh=2)
    axes[0].yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, _: f"{x:g}"))
    axes[0].set(ylabel="Exact gap (um; symlog)", title="Exact wall gap and solver-active contact")
    axes[0].text(0.99, 0.96, f'Minimum {summary["min_exact_gap_um"]:.3f} um',
                 transform=axes[0].transAxes, ha="right", va="top")
    for row in events.itertuples():
        axes[1].vlines((row.start_s + row.end_s) * 0.5e3, 0, row.peak_force_N,
                       color=COLORS["contact"], lw=1.2)
        axes[0].axvspan(row.start_s * 1e3, row.end_s * 1e3, color=COLORS["contact"], alpha=0.14, lw=0)
    axes[1].set(xlabel="Time (ms)", ylabel="Peak |Fwall| (N)")
    save(fig, "gap_and_contact_8p333")

    fig, axes = plt.subplots(2, 1, figsize=(3.54, 3.1), sharex=True, constrained_layout=True,
                             gridspec_kw={"height_ratios": [1, 2]})
    axes[0].fill_between(bridge.time_s * 1e3, 0, bridge.opposing_bridge, step="post",
                         color=COLORS["bridge"], alpha=0.7, label="Opposing bridge")
    axes[0].fill_between(bridge.time_s * 1e3, 0, bridge.phase_plateau, step="post",
                         color=COLORS["field"], alpha=0.35, label="Phase plateau")
    axes[0].set(ylim=(0, 1.05), ylabel="State", title="Persistent opposing bridge precedes tilted plateau")
    axes[0].legend(ncol=2, loc="upper right")
    axes[1].plot(bridge.time_s * 1e3, bridge.tilt_deg, color=COLORS["robot"], lw=1.0)
    axes[1].set(xlabel="Time (ms)", ylabel="True tilt (deg)")
    save(fig, "bridge_state_8p333")

    fig, axes = plt.subplots(2, 1, figsize=(3.54, 3.25), sharex=True, constrained_layout=True)
    axes[0].plot(torque.time_s * 1e3, torque.Tmag_norm_Nmm, color=COLORS["field"], lw=1.0, label="|Tmag|")
    axes[0].plot(torque.time_s * 1e3, torque.Thydro_norm_Nmm, color=COLORS["hydro"], lw=1.0,
                 label="|Thydro|")
    axes[0].set_yscale("symlog", linthresh=1e-6)
    axes[0].yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, _: f"{x:g}"))
    axes[0].set(ylabel="Torque (N mm; symlog)", title="Magnetic torque persists during geometric stall")
    axes[0].legend(loc="best")
    axes[1].plot(torque.time_s * 1e3, torque.Thydro_dot_omega_W * 1e6,
                 color=COLORS["hydro"], lw=1.0)
    axes[1].axhline(0, color="black", lw=0.7)
    axes[1].set(xlabel="Time (ms)", ylabel="Thydro dot omega (uW)")
    save(fig, "magnetic_vs_hydro_torque_8p333")

    fig, axes = plt.subplots(2, 1, figsize=(3.54, 3.25), sharex=True, constrained_layout=True)
    axes[0].plot(translation.time_s * 1e3, translation.delta_s_mm, color=COLORS["robot"], lw=1.1)
    axes[0].set(ylabel="Canonical delta s (mm)", title="Canonical centerline translation")
    axes[1].plot(translation.time_s * 1e3, translation.Vt_mm_s, color=COLORS["neutral"], lw=1.0)
    axes[1].axhline(0, color="black", lw=0.7)
    axes[1].set(xlabel="Time (ms)", ylabel="Vt (mm/s)")
    save(fig, "canonical_translation_8p333")

    fig, ax = plt.subplots(figsize=(3.54, 2.4), constrained_layout=True)
    for t_key, phase_key, label, color in (
            ("cel_t", "cel_phase", "CEL", COLORS["neutral"]),
            ("rh_t", "rh_phase", "RH old contact", COLORS["field"])):
        values = comp[phase_key]
        ax.plot(comp[t_key] * 1e3, np.degrees(values - values[0]), lw=1.0, label=label, color=color)
    ax.plot(time_ms, robot_advance, lw=1.3, label="RH zeta=0.50", color=COLORS["robot"])
    ax.set(xlabel="Time (ms)", ylabel="True-axis phase advance (deg)",
           title="Unified local-frame phase: no clear survival gain")
    ax.legend(loc="best")
    save(fig, "cel_vs_rh_zeta050_phase")


if __name__ == "__main__":
    main()
