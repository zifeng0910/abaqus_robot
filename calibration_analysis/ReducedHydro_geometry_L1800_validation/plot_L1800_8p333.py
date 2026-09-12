"""Publication-style figures for the L=1.800 mm geometry decision."""
from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "ReducedHydro_zeta050_8p333_validation"
BLUE, ORANGE, RED, GREEN, PURPLE, GREY = "#0072B2", "#D55E00", "#CC3311", "#009E73", "#882255", "#5B6573"
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
                     "xtick.labelsize": 7.2, "ytick.labelsize": 7.2, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": .8, "legend.frameon": False,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})


def save(fig, name):
    fig.savefig(HERE / f"{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(HERE / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(HERE / f"{name}.svg", bbox_inches="tight")
    fig.savefig(HERE / f"{name}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main():
    p = pd.read_csv(HERE / "L1800_8p333_local_phase.csv")
    a = pd.read_csv(HERE / "L1800_8p333_true_axis.csv")
    g = pd.read_csv(HERE / "L1800_8p333_exact_gap.csv")
    e = pd.read_csv(HERE / "L1800_8p333_contact_events.csv")
    b = pd.read_csv(HERE / "L1800_8p333_bridge_timeline.csv")
    w = pd.read_csv(HERE / "L1800_8p333_phase_windows.csv")
    q = pd.read_csv(HERE / "L1800_8p333_magnetic_hydro_torque.csv")
    en = pd.read_csv(HERE / "L1800_8p333_energy.csv")
    op = pd.read_csv(BASE / "zeta050_8p333_local_phase.csv")
    oa = pd.read_csv(BASE / "zeta050_8p333_true_axis.csv")
    old = json.loads((BASE / "zeta050_8p333_summary.json").read_text())
    new = json.loads((HERE / "L1800_8p333_summary.json").read_text())

    fig, axes = plt.subplots(2, 1, figsize=(3.54, 3.35), sharex=True, constrained_layout=True)
    axes[0].plot(p.time_s*1e3, p.phi_robot_advance_deg, color=BLUE, lw=1.15, label="Robot true axis")
    axes[0].plot(p.time_s*1e3, np.degrees(p.phi_B_local_rad-p.phi_B_local_rad.iloc[0]), color=ORANGE, lw=.9, label="Local field")
    axes[0].set(ylabel="Phase advance (deg)", title="Short robot advances, then reverses late")
    axes[0].legend(loc="upper left")
    axes[1].plot(a.time_s*1e3, a.tilt_deg, color=BLUE, lw=1.0)
    axes[1].set(xlabel="Time (ms)", ylabel="True tilt (deg)")
    save(fig, "L1800_phase_and_tilt")

    fig, axes = plt.subplots(3, 1, figsize=(3.54, 4.15), sharex=True, constrained_layout=True,
                             gridspec_kw={"height_ratios": [2, 1, 1]})
    axes[0].plot(g.time_s*1e3, g.gap_um, color=GREY, lw=.65); axes[0].axhline(0, color="black", lw=.6)
    axes[0].set_yscale("symlog", linthresh=2); axes[0].set(ylabel="Exact gap (um)", title="Contacts remain brief; opposing bridges are transient")
    axes[1].vlines((e.start_s+e.end_s)*.5e3, 0, e.peak_force_N, color=RED, lw=1.0); axes[1].set(ylabel="Peak |Fwall| (N)")
    axes[2].fill_between(b.time_s*1e3, 0, b.opposing_bridge, step="post", color=PURPLE, alpha=.75)
    axes[2].set(xlabel="Time (ms)", ylabel="Bridge", ylim=(0, 1.05))
    save(fig, "L1800_contact_and_bridge")

    fig, axes = plt.subplots(2, 1, figsize=(3.54, 3.3), sharex=True, constrained_layout=True)
    axes[0].plot(q.time_s*1e3, q.Tmag_norm_Nmm, color=ORANGE, lw=.9, label="|Tmag|")
    axes[0].plot(q.time_s*1e3, q.Thydro_norm_Nmm, color=GREEN, lw=.9, label="|Thydro|")
    axes[0].set_yscale("symlog", linthresh=1e-6); axes[0].set(ylabel="Torque (N mm)", title="Drive persists and hydro torque remains dissipative"); axes[0].legend()
    axes[1].plot(en.time_s*1e3, en.Abaqus_ETOTAL_J*1e9, color=GREY, lw=.8, label="ETOTAL")
    axes[1].plot(en.time_s*1e3, en.Abaqus_ALLAE_J*1e9, color=RED, lw=.8, label="ALLAE")
    axes[1].set(xlabel="Time (ms)", ylabel="Energy (nJ)"); axes[1].legend()
    save(fig, "L1800_torque_and_energy")

    fig, axes = plt.subplots(2, 1, figsize=(3.54, 3.4), sharex=True, constrained_layout=True)
    axes[0].plot(op.time_s*1e3, op.phi_robot_advance_deg, color=GREY, lw=.9, label="Baseline L=2.921 mm")
    axes[0].plot(p.time_s*1e3, p.phi_robot_advance_deg, color=BLUE, lw=1.1, label="L=1.800 mm")
    axes[0].set(ylabel="Phase advance (deg)", title="Axial shortening removes the long bridge, not late reversal"); axes[0].legend()
    axes[1].plot(oa.time_s*1e3, oa.tilt_deg, color=GREY, lw=.9)
    axes[1].plot(a.time_s*1e3, a.tilt_deg, color=BLUE, lw=1.1)
    axes[1].text(.98, .08, f"Longest bridge: {old['longest_bridge_ms']:.2f} to {new['longest_bridge_ms']:.2f} ms",
                 transform=axes[1].transAxes, ha="right", color=PURPLE)
    axes[1].set(xlabel="Time (ms)", ylabel="True tilt (deg)")
    save(fig, "baseline_vs_L1800_dynamics")


if __name__ == "__main__":
    main()
