"""Publication figures and fixed-camera diagnostic animations (Python only)."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Circle, Polygon
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BASE_CASE = REPO / "calibration_analysis/StraightPipeControl/case/PROD_LOCAL30_G0_STRAIGHT_CTRL"
REPLAY_CASE = HERE / "case/STRAIGHT_PRESCRIBED_WOBBLE_CONTACT_AUDIT"
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
A0 = np.array([0.9647382600216, -0.1188742372140, 0.2348382536499]); A0 /= np.linalg.norm(A0)
RADIUS = 0.66734524
ROBOT_RADIUS = 0.4075
LENGTH = 2.4
COM_X = 1.25606694407031

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": .8, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "pdf.fonttype": 42,
    "svg.fonttype": "none", "savefig.facecolor": "white",
})
COL = {"mag": "#2474A6", "contact": "#C44E52", "hydro": "#4F786D",
       "head": "#D97706", "tail": "#2563A6", "bridge": "#E7A8A8"}
ANIMATION_RENDER_SCALE = 110


def get_array(z, name):
    return max((z[k] for k in z.files if k == name or k.startswith(name + " (Repeated:")), key=len)


def vec(z, prefix, t=None, surface=None):
    cols = []
    for i in (1, 2, 3):
        if surface:
            a = max((z[k] for k in z.files if "|" + prefix + str(i) + " on surface " in k and surface in k), key=len)
        else:
            a = get_array(z, prefix + str(i))
        if t is None: t = a[:, 0]
        cols.append(np.interp(t, a[:, 0], a[:, 1]))
    return t, np.column_stack(cols)


def save_figure(fig, stem):
    fig.savefig(HERE / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(HERE / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(HERE / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(HERE / f"{stem}.tiff", dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def torque_figures():
    base = pd.read_csv(HERE / "straight_bridge_torque_balance.csv")
    m = base.time_s >= .0023
    t = base.loc[m, "time_s"].to_numpy() * 1e3
    fig, ax = plt.subplots(4, 1, figsize=(7.1, 6.4), sharex=True,
                           gridspec_kw={"height_ratios": [2.0, 1.1, 1.1, .45]})
    ax[0].plot(t, base.loc[m, "Tmag_norm_Nmm"], color=COL["mag"], label=r"$|T_{mag}|$")
    ax[0].plot(t, base.loc[m, "Tcontact_inferred_norm_Nmm"], color=COL["contact"], lw=.8,
               label=r"$|T_{contact}|$ (balance inversion)")
    ax[0].plot(t, base.loc[m, "Thydro_norm_Nmm"], color=COL["hydro"], label=r"$|T_{hydro}|$")
    ax[0].set_ylabel("Torque (N mm)"); ax[0].legend(ncol=3, loc="upper right")
    ax[0].set_title("Persistent bridge balances the available magnetic torque")
    ax[1].plot(t, base.loc[m, "omega_wobble_rad_s"], color="#6B4C9A")
    ax[1].set_ylabel(r"$|\omega_{wobble}|$ (rad s$^{-1}$)")
    ax[2].plot(t, base.loc[m, "HEAD_gap_um"], color=COL["head"], label="HEAD")
    ax[2].plot(t, base.loc[m, "TAIL_gap_um"], color=COL["tail"], label="TAIL")
    ax[2].axhline(20, color="#777777", ls="--", lw=.7); ax[2].set_ylabel("Gap (um)"); ax[2].legend(ncol=2)
    ax[3].fill_between(t, 0, base.loc[m, "opposing_bridge"], color=COL["bridge"])
    ax[3].set_yticks([0, 1], ["open", "bridge"]); ax[3].set_xlabel("Time (ms)")
    fig.tight_layout(); save_figure(fig, "Straight_Bridge_Torque_Balance")

    rep = pd.read_csv(HERE / "prescribed_reaction_vs_magnetic.csv")
    m = rep.time_s >= .0023; t = rep.loc[m, "time_s"].to_numpy() * 1e3
    fig, ax = plt.subplots(3, 1, figsize=(7.1, 5.6), sharex=True,
                           gridspec_kw={"height_ratios": [2, 1.1, .5]})
    ax[0].semilogy(t, np.maximum(rep.loc[m, "required_transverse_Nmm"], 1e-8),
                   color=COL["contact"], label="Prescribed reaction moment")
    ax[0].semilogy(t, np.maximum(rep.loc[m, "available_transverse_Nmm"], 1e-8),
                   color=COL["mag"], label="Baseline magnetic moment")
    ax[0].set_ylabel("Transverse moment (N mm)"); ax[0].legend()
    ax[0].set_title("The imposed 30 deg path activates penalty-scale wall reaction")
    ax[1].plot(t, rep.loc[m, "HEAD_gap_um"], color=COL["head"], label="HEAD")
    ax[1].plot(t, rep.loc[m, "TAIL_gap_um"], color=COL["tail"], label="TAIL")
    ax[1].axhline(20, color="#777777", ls="--", lw=.7); ax[1].axhline(0, color="#222222", lw=.6)
    ax[1].set_ylabel("Geometric gap (um)"); ax[1].legend(ncol=2)
    ax[2].fill_between(t, 0, rep.loc[m, "opposing_bridge"], color=COL["bridge"])
    ax[2].set_yticks([0, 1], ["open", "bridge"]); ax[2].set_xlabel("Time (ms)")
    fig.tight_layout(); save_figure(fig, "Prescribed_Reaction_vs_Magnetic_Torque")


def pose(case, gap_path, torque_path, replay=False):
    rp = np.load(case / "private/rp_history_private.npz")
    t, u = vec(rp, "U"); _, ur = vec(rp, "UR", t)
    rot = Rotation.from_rotvec(ur); axis = rot.apply(np.broadcast_to(A0, u.shape)); com = RP0 + u
    gap = pd.read_csv(gap_path)
    torque = pd.read_csv(torque_path)
    times = np.linspace(0, min(.016667, t[-1]), 120)
    idx = np.searchsorted(t, times).clip(0, len(t)-1)
    ident = json.loads((BASE_CASE / "case_identity.json").read_text())
    tangent = np.asarray(ident["straight_tangent_aba"]); e1 = np.asarray(ident["initial_e1_aba"]); e2 = np.asarray(ident["initial_e2_aba"])
    rel = com[idx] - np.asarray(ident["initial_center_aba_mm"])
    out = {"time": t[idx], "com_s": rel @ tangent, "com_y": rel @ e1, "com_z": rel @ e2,
           "axis_s": axis[idx] @ tangent, "axis_y": axis[idx] @ e1, "axis_z": axis[idx] @ e2,
           "head_gap": np.interp(t[idx], gap.time_s, gap.HEAD_gap_um),
           "tail_gap": np.interp(t[idx], gap.time_s, gap.TAIL_gap_um)}
    out["gap_rate"] = np.gradient(np.minimum(out["head_gap"], out["tail_gap"])*1e-3, out["time"])
    for name in torque.columns:
        if name != "time_s": out[name] = np.interp(t[idx], torque.time_s, torque[name])
    tel = pd.read_csv(BASE_CASE / "PROD_LOCAL30_G0_STRAIGHT_CTRL_telemetry.csv")
    b = np.column_stack([np.interp(t[idx], tel.t_s, tel[c]) for c in ("Bx_aba_T", "By_aba_T", "Bz_aba_T")])
    out["B_y"] = b @ e1; out["B_z"] = b @ e2
    return out


def robot_patch(ax, s, r, ds, dr, color="#4A4A4A"):
    d = np.array([ds, dr]); d /= max(np.linalg.norm(d), 1e-12); p = np.array([-d[1], d[0]])
    head = np.array([s, r]) - COM_X*d; tail = np.array([s, r]) + (LENGTH-COM_X)*d
    poly = np.vstack((head+ROBOT_RADIUS*p, tail+ROBOT_RADIUS*p,
                      tail-ROBOT_RADIUS*p, head-ROBOT_RADIUS*p))
    ax.add_patch(Polygon(poly, closed=True, fc=color, ec="black", lw=.6, alpha=.92))
    ax.add_patch(Circle(head, ROBOT_RADIUS, fc="#D97706", ec="black", lw=.6))
    ax.text(*head, "H", ha="center", va="center", color="white", fontsize=6, weight="bold")
    ax.text(*tail, "T", ha="center", va="center", color="white", fontsize=6, weight="bold")


def draw_pose(ax_side, ax_cross, p, i, title, torque_prefix):
    for ax in (ax_side, ax_cross): ax.clear()
    ax_side.axhspan(-RADIUS, RADIUS, color="#B8D8E8", alpha=.18)
    ax_side.axhline(RADIUS, color="#4F91B3"); ax_side.axhline(-RADIUS, color="#4F91B3")
    radial = np.hypot(p["axis_y"][i], p["axis_z"][i]); sign = 1 if p["axis_y"][i] >= 0 else -1
    robot_patch(ax_side, p["com_s"][i], sign*np.hypot(p["com_y"][i], p["com_z"][i]), p["axis_s"][i], radial)
    ax_side.set(xlim=(-1.8, 1.8), ylim=(-.9, .9), aspect="equal", xlabel="Tube axis (mm)", ylabel="Radial (mm)")
    ax_side.set_title(title)
    ax_cross.add_patch(Circle((0,0), RADIUS, fc="#B8D8E8", alpha=.18, ec="#4F91B3"))
    cy, cz = p["com_y"][i], p["com_z"][i]
    hy, hz = cy-COM_X*p["axis_y"][i], cz-COM_X*p["axis_z"][i]
    ty, tz = cy+(LENGTH-COM_X)*p["axis_y"][i], cz+(LENGTH-COM_X)*p["axis_z"][i]
    ax_cross.plot([hy,ty], [hz,tz], color="#444444", lw=7, solid_capstyle="round")
    ax_cross.scatter([hy,ty],[hz,tz],c=[COL["head"],COL["tail"]],s=28,zorder=4)
    ax_cross.text(hy,hz,"H",ha="center",va="center",color="white",fontsize=5,weight="bold")
    ax_cross.text(ty,tz,"T",ha="center",va="center",color="white",fontsize=5,weight="bold")
    by,bz=p["B_y"][i],p["B_z"][i]; bn=max(np.hypot(by,bz),1e-12)
    ax_cross.arrow(0,0,.35*by/bn,.35*bz/bn,width=.012,color=COL["mag"],length_includes_head=True)
    if torque_prefix == "base":
        cyq=p["Tcontact_inferred_e1_Nmm"][i]; czq=p["Tcontact_inferred_e2_Nmm"][i]
    else:
        cyq=p["RM_e1_Nmm"][i]; czq=p["RM_e2_Nmm"][i]
    qn=max(np.hypot(cyq,czq),1e-12)
    ax_cross.arrow(cy,cz,.28*cyq/qn,.28*czq/qn,width=.012,color=COL["contact"],length_includes_head=True)
    bridge=p["head_gap"][i]<=20 and p["tail_gap"][i]<=20
    ax_cross.text(-.82,.84,(f"t={p['time'][i]*1e3:.2f} ms  bridge={'YES' if bridge else 'NO'}\n"
                              f"gap H/T={p['head_gap'][i]:.1f}/{p['tail_gap'][i]:.1f} um\n"
                              f"normal speed proxy={p['gap_rate'][i]:.1f} mm/s"),va="top",fontsize=6)
    ax_cross.set(xlim=(-.9,.9),ylim=(-.9,.9),aspect="equal",xlabel="e1 (mm)",ylabel="e2 (mm)")


def animations():
    base = pose(BASE_CASE, REPO/"calibration_analysis/StraightPipeControl/straight_exact_gap_timeseries.csv",
                HERE/"straight_bridge_torque_balance.csv", False)
    replay = pose(REPLAY_CASE, HERE/"prescribed_exact_gap_timeseries.csv",
                  HERE/"prescribed_reaction_vs_magnetic.csv", True)

    fig, axes = plt.subplots(1,2,figsize=(7.0,3.2))
    def update_base(i): draw_pose(axes[0],axes[1],base,i,"Free magnetic baseline","base"); return []
    FuncAnimation(fig,update_base,frames=len(base["time"]),interval=70).save(
        HERE/"Straight_Bridge_Torque_Balance.gif",writer=PillowWriter(fps=14),dpi=ANIMATION_RENDER_SCALE)
    plt.close(fig)

    fig, axes = plt.subplots(1,2,figsize=(7.0,3.2))
    def update_rep(i): draw_pose(axes[0],axes[1],replay,i,"Prescribed 30 deg replay","replay"); return []
    FuncAnimation(fig,update_rep,frames=len(replay["time"]),interval=70).save(
        HERE/"STRAIGHT_PRESCRIBED_WOBBLE_CONTACT_AUDIT_DualView.gif",writer=PillowWriter(fps=14),dpi=ANIMATION_RENDER_SCALE)
    plt.close(fig)

    fig, axes = plt.subplots(2,2,figsize=(7.0,6.0))
    def update_cmp(i):
        draw_pose(axes[0,0],axes[0,1],base,i,"Free magnetic baseline","base")
        draw_pose(axes[1,0],axes[1,1],replay,i,"Prescribed 30 deg replay","replay")
        return []
    FuncAnimation(fig,update_cmp,frames=min(len(base["time"]),len(replay["time"])),interval=70).save(
        HERE/"Baseline_vs_Prescribed_Contact.gif",writer=PillowWriter(fps=14),dpi=ANIMATION_RENDER_SCALE)
    plt.close(fig)


if __name__ == "__main__":
    torque_figures()
    animations()
