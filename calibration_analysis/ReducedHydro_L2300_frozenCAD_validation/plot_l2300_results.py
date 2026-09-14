"""Create the publication-style L2300 validation figure set from public CSVs."""
from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
LONG = REPO/"calibration_analysis"/"ReducedHydro_zeta050_8p333_validation"
SHORT = REPO/"calibration_analysis"/"ReducedHydro_FreeCAD_L1800_validation"
BLUE = "#0072B2"; ORANGE = "#D55E00"; GREEN = "#009E73"; RED = "#CC3311"
GRAY = "#667085"; LIGHT = "#D9E2EC"
mpl.rcParams.update({"font.family":"sans-serif", "font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "font.size":7, "axes.titlesize":8, "axes.labelsize":7, "xtick.labelsize":6,
    "ytick.labelsize":6, "legend.fontsize":6, "axes.spines.top":False,
    "axes.spines.right":False, "axes.linewidth":.8, "legend.frameon":False,
    "svg.fonttype":"none", "pdf.fonttype":42})


def save(fig, name):
    fig.savefig(HERE/(name+".png"), dpi=300, bbox_inches="tight")
    fig.savefig(HERE/(name+".pdf"), bbox_inches="tight")
    fig.savefig(HERE/(name+".svg"), bbox_inches="tight")
    plt.close(fig)


def add_panel(ax, label):
    ax.text(-.11, 1.04, label, transform=ax.transAxes, fontweight="bold", fontsize=8)


def build_comparison(runtime):
    long = json.loads((LONG/"zeta050_8p333_summary.json").read_text())
    short = json.loads((SHORT/"freecad_8p333_summary.json").read_text())
    rows = [
        {"case":"Long L~2.9", "length_mm":2.92146, "outcome":"TOO LONG / JAM",
         "crossed_90deg":False, "longest_bridge_ms":long["longest_bridge_ms"],
         "contact_events":long["contact_event_count"], "delta_s_mm":long["canonical_delta_s_mm"]},
        {"case":"L1.8", "length_mm":1.8, "outcome":"TOO SHORT / TUMBLE",
         "crossed_90deg":True, "longest_bridge_ms":short["longest_bridge_ms"],
         "contact_events":short["contact_event_count"], "delta_s_mm":short["canonical_delta_s_mm"]},
        {"case":"L2.3", "length_mm":2.3, "outcome":runtime["classification"],
         "crossed_90deg":runtime["axis_crossed_90deg"],
         "longest_bridge_ms":runtime["longest_opposing_bridge_ms"],
         "contact_events":runtime["contact_event_count"], "delta_s_mm":runtime["canonical_delta_s_mm"]},
    ]
    frame = pd.DataFrame(rows)
    frame.to_csv(HERE/"Long_L1800_L2300_summary.csv", index=False)
    return frame


def main():
    runtime = json.loads((HERE/"L2300_dynamic_summary.json").read_text())
    tilt = pd.read_csv(HERE/"L2300_directed_tilt.csv")
    support = pd.read_csv(HERE/"L2300_wall_support_15_20_25um.csv")
    support_t = pd.read_csv(HERE/"L2300_wall_support_timeline.csv")
    sectors = pd.read_csv(HERE/"L2300_wall_sector_timeline.csv")
    gap = pd.read_csv(HERE/"L2300_exact_gap.csv")
    events = pd.read_csv(HERE/"L2300_contact_events.csv")
    bridge = pd.read_csv(HERE/"L2300_bridge_timeline.csv")
    phase = pd.read_csv(HERE/"L2300_phase.csv")
    torque = pd.read_csv(HERE/"L2300_torque.csv")
    energy = pd.read_csv(HERE/"L2300_energy.csv")
    translation = pd.read_csv(HERE/"L2300_translation.csv")
    regression = pd.read_csv(HERE/"L2300_threshold_aware_mesh_acceptance.csv")
    comparison = build_comparison(runtime)

    fig, ax = plt.subplots(figsize=(3.5, 2.75))
    status = regression.contact_threshold_status.eq("CONTACT_THRESHOLD_AMBIGUOUS") | regression.wall_support_threshold_status.eq("WALL_SUPPORT_THRESHOLD_AMBIGUOUS")
    ax.scatter(regression.min_gap_um[~status], regression.mesh_min_gap_um[~status], s=10, color=BLUE, alpha=.65, label="Resolved poses")
    ax.scatter(regression.min_gap_um[status], regression.mesh_min_gap_um[status], s=24, color=ORANGE, edgecolor="white", lw=.4, label="Threshold ambiguous")
    lo=min(regression.min_gap_um.min(),regression.mesh_min_gap_um.min());hi=max(regression.min_gap_um.max(),regression.mesh_min_gap_um.max())
    ax.plot([lo,hi],[lo,hi],color=GRAY,lw=.8,ls="--"); ax.fill_between([lo,hi],[lo-5,hi-5],[lo+5,hi+5],color=LIGHT,alpha=.45)
    ax.set(xlabel="Frozen CAD minimum gap (um)",ylabel="Solver-surface minimum gap (um)",title="All four Boolean mismatches lie inside threshold uncertainty")
    ax.legend(loc="upper left"); save(fig,"L2300_threshold_aware_surface_regression")

    fig, ax = plt.subplots(figsize=(3.5,2.5)); ms=tilt.time_s*1e3
    ax.plot(ms,tilt.directed_tilt_deg,color=BLUE,lw=1);ax.axhline(90,color=RED,ls="--",lw=.8,label="Tumble boundary")
    ax.set(xlabel="Time (ms)",ylabel="Directed tilt (deg)",title="Directed tilt remains on the initial-polarity side")
    ax.legend(); save(fig,"L2300_directed_tilt")

    fig, ax = plt.subplots(figsize=(3.5,2.5)); x=np.arange(len(support)); width=.19
    for j,(col,label,color) in enumerate((("HEAD_supported_fraction","HEAD",ORANGE),("TAIL_supported_fraction","TAIL",GREEN),
                                         ("either_end_supported_fraction","Either",BLUE),("both_end_supported_fraction","Both",GRAY))):
        ax.bar(x+(j-1.5)*width,support[col],width,color=color,label=label)
    ax.set_xticks(x,support.threshold_um.astype(int).astype(str));ax.set(xlabel="Gap threshold (um)",ylabel="Time fraction",ylim=(0,1),title="Wall support is robust to the 15-25 um threshold band")
    ax.legend(ncol=2);save(fig,"L2300_wall_support_fraction")

    fig, ax = plt.subplots(figsize=(3.5,2.7)); mask=sectors.HEAD_supported_20um.astype(bool)
    ax.scatter(np.degrees(np.unwrap(np.radians(sectors.robot_phase_deg)))[mask],sectors.HEAD_sector_deg[mask],s=5,color=ORANGE,label="HEAD",alpha=.65)
    mask=sectors.TAIL_supported_20um.astype(bool)
    ax.scatter(np.degrees(np.unwrap(np.radians(sectors.robot_phase_deg)))[mask],sectors.TAIL_sector_deg[mask],s=5,color=GREEN,label="TAIL",alpha=.65)
    ax.set(xlabel="Unwrapped robot phase (deg)",ylabel="Nearest wall sector (deg)",ylim=(0,360),title="Supported wall sector moves with local phase")
    ax.legend();save(fig,"L2300_wall_sector_vs_phase")

    fig, axes=plt.subplots(2,1,figsize=(5.2,3.7),sharex=True,gridspec_kw={"height_ratios":[2,1]})
    complete_gap=gap[gap.sampling.str.startswith("coarse")]
    for name,color in (("HEAD",ORANGE),("BODY",GRAY),("TAIL",GREEN)):
        axes[0].plot(complete_gap.time_s*1e3,complete_gap[name+"_gap_um"],color=color,lw=.65,label=name)
    axes[0].axhspan(15,25,color=LIGHT,alpha=.45);axes[0].axhline(0,color=RED,ls="--",lw=.7);axes[0].set(ylabel="Exact CAD gap (um)",ylim=(-5,100));axes[0].legend(ncol=3)
    for row in events.itertuples(): axes[1].axvspan(row.start_s*1e3,row.end_s*1e3,color=RED,alpha=.18)
    axes[1].step(bridge.time_s*1e3,bridge.opposing_bridge_20um,where="post",color=BLUE,lw=.8,label="Bridge 20 um")
    axes[1].set(xlabel="Time (ms)",ylabel="State",ylim=(-.05,1.1));axes[1].legend();axes[0].set_title("Exact CAD proximity and solver-contact windows remain distinct evidence")
    add_panel(axes[0],"a");add_panel(axes[1],"b");save(fig,"L2300_gap_contact_bridge")

    fig, axes=plt.subplots(2,1,figsize=(5.2,3.6),sharex=True)
    axes[0].plot(phase.time_s*1e3,phase.robot_phase_advance_deg,color=BLUE,lw=.9,label="Robot")
    axes[0].plot(phase.time_s*1e3,np.degrees(phase.B_local_phase_rad-phase.B_local_phase_rad.iloc[0]),color=ORANGE,lw=.8,label="Field")
    axes[0].set_ylabel("Phase advance (deg)");axes[0].legend()
    axes[1].plot(phase.time_s*1e3,phase.phase_window_range_deg,color=GREEN,lw=.8);axes[1].axhline(10,color=RED,ls="--",lw=.7)
    axes[1].set(xlabel="Time (ms)",ylabel="0.5 ms range (deg)");axes[0].set_title("True-axis phase progression over the quarter-cycle record")
    add_panel(axes[0],"a");add_panel(axes[1],"b");save(fig,"L2300_phase_progression")

    fig, ax=plt.subplots(figsize=(3.5,2.5));ax.semilogy(torque.time_s*1e3,np.maximum(torque.Tmag_norm_Nmm,1e-12),color=BLUE,lw=.8,label="Magnetic")
    ax.semilogy(torque.time_s*1e3,np.maximum(torque.Thydro_norm_Nmm,1e-12),color=ORANGE,lw=.8,label="Hydrodynamic")
    torque_ticks=[1e-11,1e-9,1e-7,1e-5,1e-3]
    ax.set_yticks(torque_ticks,["1e-11","1e-9","1e-7","1e-5","1e-3"])
    ax.set(xlabel="Time (ms)",ylabel="Torque norm (N mm)",title="Magnetic torque persists while hydrodynamics dissipate");ax.legend();save(fig,"L2300_torque")

    fig, axes=plt.subplots(2,1,figsize=(5.2,3.6),sharex=True)
    axes[0].plot(translation.time_s*1e3,translation.delta_s_mm,color=BLUE,lw=.9);axes[0].axvline(2,color=GRAY,ls="--",lw=.7)
    axes[0].set_ylabel("Delta s (mm)");axes[1].plot(translation.time_s*1e3,translation.Vt_mm_s,color=GREEN,lw=.7);axes[1].axhline(0,color=GRAY,lw=.6)
    axes[1].set(xlabel="Time (ms)",ylabel="Vt (mm/s)");axes[0].set_title("Canonical translation is a secondary success gate")
    add_panel(axes[0],"a");add_panel(axes[1],"b");save(fig,"L2300_translation")

    comparison["longest_phase_plateau_ms"]=[json.loads((LONG/"zeta050_8p333_summary.json").read_text())["longest_phase_plateau_ms"],
        json.loads((SHORT/"freecad_8p333_summary.json").read_text())["longest_phase_plateau_ms"],runtime["longest_phase_plateau_ms"]]
    comparison.to_csv(HERE/"Long_L1800_L2300_summary.csv",index=False)
    fig, axes=plt.subplots(1,4,figsize=(7.2,2.55));colors=[GRAY,ORANGE,BLUE]
    axes[0].bar(comparison.case,comparison.longest_bridge_ms,color=colors);axes[0].axhline(.5,color=RED,ls="--",lw=.7);axes[0].set_ylabel("Longest bridge (ms)")
    axes[1].bar(comparison.case,comparison.crossed_90deg.astype(int),color=colors);axes[1].set(ylabel="Crossed 90 deg",yticks=[0,1])
    axes[2].bar(comparison.case,comparison.longest_phase_plateau_ms,color=colors);axes[2].axhline(.5,color=RED,ls="--",lw=.7);axes[2].set_ylabel("Phase plateau (ms)")
    axes[3].bar(comparison.case,comparison.delta_s_mm,color=colors);axes[3].axhline(0,color=GRAY,lw=.6);axes[3].set_ylabel("Total Delta s (mm)")
    for i,ax in enumerate(axes): ax.tick_params(axis="x",rotation=25);add_panel(ax,chr(97+i))
    fig.suptitle("Bridge alone is insufficient for jam: posture and phase must also lock",fontsize=8);save(fig,"Long_L1800_L2300_summary")

    fig, axes=plt.subplots(2,1,figsize=(5.2,3.6),sharex=True)
    axes[0].plot(energy.time_s*1e3,energy.Ktotal_rigid_J*1e6,color=BLUE,lw=.8,label="Rigid KE")
    axes[0].plot(energy.time_s*1e3,energy.Abaqus_ALLKE_J*1e6,color=ORANGE,lw=.7,ls="--",label="Abaqus ALLKE")
    axes[0].set_ylabel("Energy (uJ)");axes[0].legend()
    axes[1].plot(energy.time_s*1e3,energy.Whydro_J*1e6,color=GREEN,lw=.8,label="Hydro work")
    axes[1].plot(energy.time_s*1e3,energy.Abaqus_ETOTAL_J*1e6,color=GRAY,lw=.8,label="ETOTAL")
    axes[1].set(xlabel="Time (ms)",ylabel="Energy/work (uJ)");axes[1].legend();axes[0].set_title("Energy remains bounded through repeated wall exchanges")
    add_panel(axes[0],"a");add_panel(axes[1],"b");save(fig,"L2300_energy")


if __name__ == "__main__":
    main()
