"""Render local-axis orbits, dual-view GIFs, and synchronized comparison."""
from pathlib import Path
import json
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import pandas as pd

from audit_common import AUDIT, CURVE_CSV, RP0


CASES = ("FRAME_CURRENT", "FRAME_LOCAL30", "FRAME_LOCAL40", "FRAME_LOCAL30_RAMP", "FRAME_LOCAL_BEST_CONTACT")
COLORS = {"FRAME_CURRENT": "#D55E00", "FRAME_LOCAL30": "#0072B2", "FRAME_LOCAL40": "#009E73",
          "FRAME_LOCAL30_RAMP": "#56B4E9", "FRAME_LOCAL_BEST_CONTACT": "#CC79A7"}
PIPE_RADIUS = .94
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 7.2, "axes.linewidth": .7, "axes.spines.right": False,
                     "axes.spines.top": False, "legend.frameon": False,
                     "svg.fonttype": "none", "pdf.fonttype": 42})


def save(fig, stem):
    fig.savefig(str(stem) + ".png", dpi=600, bbox_inches="tight")
    fig.savefig(str(stem) + ".svg", bbox_inches="tight")
    fig.savefig(str(stem) + ".pdf", bbox_inches="tight")
    plt.close(fig)


def existing_cases():
    return [case for case in CASES if (AUDIT / "cases" / case / "pose_diagnostic.csv").exists()]


def orbit_plot(case_id):
    folder = AUDIT / "cases" / case_id; pose = pd.read_csv(folder / "pose_diagnostic.csv")
    color = COLORS[case_id]; fig, ax = plt.subplots(figsize=(3.54, 3.25), constrained_layout=True)
    for angle, style in ((20, ":"), (30, "--"), (40, "-.")):
        radius = math.sin(math.radians(angle)); circle = plt.Circle((0,0), radius, fill=False, color="#B8BEC3", lw=.7, ls=style)
        ax.add_patch(circle)
        label_angle = math.radians(135.0)
        ax.text(radius*math.cos(label_angle), radius*math.sin(label_angle), "%d deg"%angle,
                color="#747C82", fontsize=6, ha="center", va="center")
    ax.plot(pose.q1, pose.q2, color=color, lw=1.2)
    contact = pose.contact_active.astype(bool)
    ax.scatter(pose.q1[contact], pose.q2[contact], s=5, color="#CC79A7", alpha=.38, label="Wall contact")
    ax.scatter(pose.q1.iloc[0], pose.q2.iloc[0], s=26, color="black", label="Start", zorder=4)
    ax.scatter(pose.q1.iloc[-1], pose.q2.iloc[-1], s=30, facecolors="none", edgecolors=color, label="End", zorder=4)
    arrow_indices = np.linspace(0, len(pose)-2, 10).astype(int)
    for i in arrow_indices:
        j = min(i + max(1, len(pose)//100), len(pose)-1)
        ax.annotate("", xy=(pose.q1.iloc[j],pose.q2.iloc[j]), xytext=(pose.q1.iloc[i],pose.q2.iloc[i]),
                    arrowprops=dict(arrowstyle="->", color=color, lw=.7))
    ax.axhline(0,color="#D3D7DA",lw=.6); ax.axvline(0,color="#D3D7DA",lw=.6)
    ax.set_aspect("equal"); ax.set(xlim=(-1.02,1.02),ylim=(-1.02,1.02),xlabel="a dot e1",ylabel="a dot e2",title=case_id+" robot-axis local orbit")
    ax.legend(fontsize=6, loc="lower left")
    save(fig, folder / "robot_axis_local_orbit")


def interpolated_pose(case_id, count=151):
    folder = AUDIT / "cases" / case_id; pose = pd.read_csv(folder / "pose_diagnostic.csv")
    proximity = pd.read_csv(folder / "proximity_diagnostic.csv")
    assert np.all(np.diff(pose.time_s) > 0), "pose time must increase"
    assert np.all(np.diff(proximity.time_s) > 0), "proximity time must increase"
    times = np.linspace(0, pose.time_s.iloc[-1], count); out = {"time_s": times}
    for column in pose.columns:
        if column != "time_s": out[column] = np.interp(times, pose.time_s, pose[column])
    sector = np.unwrap(np.radians(proximity.nearest_wall_sector_deg))
    out["sector_deg"] = np.degrees(np.interp(times, proximity.time_s, sector))
    return out


def draw_local(ax, data, i, title):
    ax.cla(); center=np.array([data[k][i] for k in ("center_x_mm","center_y_mm","center_z_mm")])
    e1=np.array([data[k][i] for k in ("e1_x","e1_y","e1_z")]); e2=np.array([data[k][i] for k in ("e2_x","e2_y","e2_z")])
    head=np.array([data[k][i] for k in ("head_x_mm","head_y_mm","head_z_mm")]); tail=np.array([data[k][i] for k in ("tail_x_mm","tail_y_mm","tail_z_mm")])
    hp=np.array([np.dot(head-center,e1),np.dot(head-center,e2)]); tp=np.array([np.dot(tail-center,e1),np.dot(tail-center,e2)])
    field=np.array([data[k][i] for k in ("Bx_T","By_T","Bz_T")]); bp=np.array([np.dot(field,e1),np.dot(field,e2)])*55
    ax.add_patch(plt.Circle((0,0),PIPE_RADIUS,facecolor="#DDE5EA",edgecolor="#69737B",alpha=.35,lw=1))
    ax.plot([hp[0],tp[0]],[hp[1],tp[1]],color="#0072B2",lw=8,solid_capstyle="round")
    ax.scatter(*hp,s=28,color="#D55E00",zorder=4); ax.scatter(*tp,s=20,color="#0072B2",zorder=4)
    ax.arrow(0,0,bp[0],bp[1],width=.012,head_width=.08,color="#009E73",length_includes_head=True)
    if data["contact_active"][i] > .5:
        angle=math.radians(data["sector_deg"][i]); ax.scatter(PIPE_RADIUS*math.cos(angle),PIPE_RADIUS*math.sin(angle),marker="x",s=35,color="#CC79A7")
    ax.text(.03,.97,"e2",transform=ax.transAxes,va="top",color="#747C82"); ax.text(.96,.51,"e1",transform=ax.transAxes,ha="right",color="#747C82")
    ax.axhline(0,color="#C8CED2",lw=.5); ax.axvline(0,color="#C8CED2",lw=.5); ax.set_aspect("equal")
    ax.set(xlim=(-1.08,1.08),ylim=(-1.08,1.08),xlabel="local e1 (mm)",ylabel="local e2 (mm)",title=title)


def dual_gif(case_id):
    folder=AUDIT/"cases"/case_id; data=interpolated_pose(case_id)
    curve=pd.read_csv(CURVE_CSV)[["x_mm","y_mm","z_mm"]].to_numpy(float); near=curve[np.linalg.norm(curve-RP0,axis=1)<7]
    lo=near.min(axis=0)-.5; hi=near.max(axis=0)+.5
    fig=plt.figure(figsize=(7.2,3.5),dpi=int("133")); ax=fig.add_subplot(1,2,1,projection="3d"); bx=fig.add_subplot(1,2,2)
    title_text = fig.suptitle("", fontsize=8.5)
    status_text = fig.text(.5,.01,"",ha="center",fontsize=6.2)
    def update(i):
        ax.cla(); ax.plot(near[:,0],near[:,1],near[:,2],color="#AEB5BA",lw=6,alpha=.35)
        head=np.array([data[k][i] for k in ("head_x_mm","head_y_mm","head_z_mm")]); tail=np.array([data[k][i] for k in ("tail_x_mm","tail_y_mm","tail_z_mm")]); rp=np.array([data[k][i] for k in ("rp_x_mm","rp_y_mm","rp_z_mm")])
        field=np.array([data[k][i] for k in ("Bx_T","By_T","Bz_T")]); field=field/max(np.linalg.norm(field),1e-15)
        ax.plot([head[0],tail[0]],[head[1],tail[1]],[head[2],tail[2]],color="#0072B2",lw=9,solid_capstyle="round"); ax.scatter(*head,s=22,color="#D55E00")
        ax.quiver(*rp,*field,length=.75,color="#009E73",linewidth=1.5); ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),zlim=(lo[2],hi[2])); ax.set_box_aspect(hi-lo); ax.view_init(elev=23,azim=-58); ax.set_axis_off(); ax.set_title("Fixed global camera")
        draw_local(bx,data,i,"COM-following tube-axis view")
        robot_w=(data["robot_phase_unwrapped_rad"][i]-data["robot_phase_unwrapped_rad"][0])/(2*math.pi); field_w=(data["field_phase_unwrapped_rad"][i]-data["field_phase_unwrapped_rad"][0])/(2*math.pi)
        title_text.set_text("%s | t=%.3f ms | robot winding=%+.3f | field winding=%+.3f"%(case_id,data["time_s"][i]*1e3,robot_w,field_w))
        status_text.set_text("tilt=%.1f deg | ds=%+.3f mm | contact=%s | robot phase=%+.1f deg | field phase=%+.1f deg"%(data["directed_tilt_deg"][i],data["delta_s_mm"][i],"ON" if data["contact_active"][i]>.5 else "OFF",math.degrees(data["robot_phase_unwrapped_rad"][i]),math.degrees(data["field_phase_unwrapped_rad"][i])))
        return [title_text, status_text]
    animation=FuncAnimation(fig,update,frames=len(data["time_s"]),interval=55,blit=False); animation.save(folder/(case_id+"_DualView_33p333.gif"),writer=PillowWriter(fps=18)); plt.close(fig)


def comparison_gif(case_ids):
    datasets=[interpolated_pose(case) for case in case_ids]; fig,axes=plt.subplots(1,len(case_ids),figsize=(9.6,3.25),dpi=int("120"))
    def update(i):
        for ax,case,data in zip(axes,case_ids,datasets):
            draw_local(ax,data,i,case)
            rw=(data["robot_phase_unwrapped_rad"][i]-data["robot_phase_unwrapped_rad"][0])/(2*math.pi)
            fw=(data["field_phase_unwrapped_rad"][i]-data["field_phase_unwrapped_rad"][0])/(2*math.pi)
            ax.text(.02,.02,"robot %+.2f | field %+.2f\ntilt %.1f deg | contact %s"%(rw,fw,data["directed_tilt_deg"][i],"ON" if data["contact_active"][i]>.5 else "OFF"),transform=ax.transAxes,fontsize=6,va="bottom")
        fig.suptitle("One-cycle magnetic-frame comparison | t=%.3f ms"%(datasets[0]["time_s"][i]*1e3),fontsize=9)
        return []
    animation=FuncAnimation(fig,update,frames=len(datasets[0]["time_s"]),interval=55,blit=False); animation.save(AUDIT/"Current_vs_Local30_vs_Local40_OneCycle.gif",writer=PillowWriter(fps=18)); plt.close(fig)


def winding_plot(case_ids):
    fig,axes=plt.subplots(2,1,figsize=(7.2,3.8),sharex=True,constrained_layout=True)
    for case in case_ids:
        pose=pd.read_csv(AUDIT/"cases"/case/"pose_diagnostic.csv"); t=pose.time_s*1e3
        rw=(pose.robot_phase_unwrapped_rad-pose.robot_phase_unwrapped_rad.iloc[0])/(2*math.pi); fw=(pose.field_phase_unwrapped_rad-pose.field_phase_unwrapped_rad.iloc[0])/(2*math.pi)
        axes[0].plot(t,rw,color=COLORS[case],label=case); axes[1].plot(t,fw,color=COLORS[case],label=case)
    axes[0].set(ylabel="Robot-axis winding"); axes[1].set(xlabel="Time (ms)",ylabel="Local B winding"); axes[0].legend(ncol=len(case_ids),fontsize=6)
    for letter,ax in zip("ab",axes): ax.text(-.08,1.05,letter,transform=ax.transAxes,fontweight="bold",fontsize=8)
    save(fig,AUDIT/"figures"/"dynamic_winding_comparison")


def main():
    cases=existing_cases()
    for case in cases: orbit_plot(case); dual_gif(case)
    mandatory=[case for case in ("FRAME_CURRENT","FRAME_LOCAL30","FRAME_LOCAL40") if case in cases]
    if len(mandatory)==3: comparison_gif(mandatory)
    if cases: winding_plot(cases)


if __name__=="__main__": main()
