"""Python-only publication figures for the cross-geometry causal audit."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
BLUE,ORANGE,RED,GREEN,PURPLE,GREY="#0072B2","#D55E00","#CC3311","#009E73","#882255","#5B6573"
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
 "font.size":7.2,"axes.titlesize":8,"axes.labelsize":7.2,"xtick.labelsize":7.2,"ytick.labelsize":7.2,
 "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.8,"legend.frameon":False,
 "pdf.fonttype":42,"svg.fonttype":"none"})


def save(fig,name):
    fig.savefig(HERE/f"{name}.png",dpi=600,bbox_inches="tight")
    fig.savefig(HERE/f"{name}.pdf",bbox_inches="tight")
    fig.savefig(HERE/f"{name}.svg",bbox_inches="tight")
    fig.savefig(HERE/f"{name}.tiff",dpi=600,bbox_inches="tight")
    plt.close(fig)


def draw_pose(ax,nodes,row):
    colors={"HEAD":ORANGE,"TAIL":GREEN,"BODY":GREY}
    for part in ("HEAD","TAIL","BODY"):
        subset=nodes[nodes.part==part]
        ax.scatter(np.radians(subset.wall_sector_deg),20-subset.gap_um,s=9,color=colors[part],alpha=.75,label=part)
    ax.set_theta_zero_location("E"); ax.set_theta_direction(1); ax.set_ylim(0,21)
    ax.set_thetagrids([0,45,90,135,180,225,315])
    ax.set_yticks([5,10,15,20]); ax.set_yticklabels(["15","10","5","0"])
    ax.set_title(f"{row.pose.replace('_',' ')} | {row.time_s*1e3:.3f} ms\ntilt {row.directed_tilt_deg:.1f} deg",pad=9)


def main():
    profile=pd.read_csv(HERE/"old_vs_freecad_head_profile.csv")
    summary=pd.read_csv(HERE/"cross_replay_bridge_summary.csv")
    poses=pd.read_csv(HERE/"bridge_pose_geometry.csv")
    nodes=pd.read_csv(HERE/"bridge_pose_nearwall_nodes.csv")

    fig,ax=plt.subplots(figsize=(3.54,2.45),constrained_layout=True)
    ax.plot(profile.axial_x_mm,profile.r_old_mm,color=GREY,lw=1,label="Old scaled mesh")
    ax.plot(profile.axial_x_mm,profile.r_new_mm,color=BLUE,lw=1.1,label="Exact FreeCAD mesh")
    ax.fill_between(profile.axial_x_mm,profile.r_old_mm,profile.r_new_mm,
                    where=(profile.axial_x_mm<-.55)&(profile.r_new_mm>profile.r_old_mm),color=ORANGE,alpha=.2,label="New HEAD outward")
    ax.set(xlabel="Body-fixed axial coordinate (mm)",ylabel="Outer radius (mm)",
           title="Exact CAD HEAD is fuller near the nose")
    ax.legend(loc="lower center");save(fig,"old_vs_freecad_head_profile")

    labels=[f"{g}-{t}" for g,t in zip(summary.geometry,summary.trajectory)]
    x=np.arange(len(summary)); width=.19
    fig,ax=plt.subplots(figsize=(3.54,2.55),constrained_layout=True)
    for j,(threshold,color) in enumerate(zip((20,10,5,0),(PURPLE,BLUE,ORANGE,GREY))):
        ax.bar(x+(j-1.5)*width,summary[f"longest{threshold}_ms"],width,color=color,label=f"{threshold} um")
    ax.axhline(.5,color=RED,lw=.75,ls="--"); ax.set_xticks(x,labels)
    ax.set(ylabel="Longest opposing bridge (ms)",title="Geometry substitution alone crosses the persistence screen")
    ax.legend(ncol=2);save(fig,"cross_replay_bridge_duration")

    fig,axes=plt.subplots(1,3,figsize=(7.2,2.65),subplot_kw={"projection":"polar"},constrained_layout=True)
    for ax,row in zip(axes,poses.itertuples()): draw_pose(ax,nodes[nodes.pose==row.pose],row)
    handles,labels=axes[0].get_legend_handles_labels(); fig.legend(handles,labels,loc="lower center",ncol=3)
    fig.suptitle("Observed longest bridge retains opposing near-wall sectors",fontsize=8)
    save(fig,"bridge_pose_geometry")
    for row in poses.itertuples():
        fig,ax=plt.subplots(figsize=(3.54,3.0),subplot_kw={"projection":"polar"},constrained_layout=True)
        draw_pose(ax,nodes[nodes.pose==row.pose],row); ax.legend(loc="lower left",bbox_to_anchor=(-.15,-.13),ncol=3)
        save(fig,f"bridge_pose_geometry_{row.pose}")


if __name__=="__main__":main()
