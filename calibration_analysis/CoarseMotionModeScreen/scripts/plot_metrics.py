"""Nature-style quantitative overview without an overall score."""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent; SCREEN=HERE.parent
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"svg.fonttype":"none","pdf.fonttype":42,"font.size":7.2,"axes.spines.right":False,"axes.spines.top":False,"axes.linewidth":.7,"legend.frameon":False})
COLORS={"geometry":"#0072B2","gradient":"#009E73","cone":"#D55E00","B":"#CC79A7","frequency":"#E69F00","friction":"#56B4E9","collision_damping":"#8C564B","hydro_translation":"#6F7F3F","hydro_rotation":"#9467BD","combo":"#222222"}

def main():
    df=pd.read_csv(SCREEN/"metrics/all_cases_metrics.csv"); plots=SCREEN/"plots"; plots.mkdir(exist_ok=True)
    fig=plt.figure(figsize=(7.2,5.2))
    gs=fig.add_gridspec(2,2,height_ratios=(1,1.18),hspace=.58,wspace=.38)
    axes=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,:])]
    for family,g in df.groupby("family",sort=False):
        kw=dict(s=20,color=COLORS[family],label=family.replace("_"," "),alpha=.85,edgecolor="white",linewidth=.35)
        axes[0].scatter(abs(g.delta_s_mm),abs(g.robot_phase_advance_deg),**kw)
        axes[1].scatter(g.max_directed_tilt_deg,g.KErot_KEtrans_ratio,**kw)
    axes[0].set(xlabel="|Axial displacement| (mm)",ylabel="|Robot phase advance| (deg)",title="Rotation versus translation")
    axes[1].set(xlabel="Maximum directed tilt (deg)",ylabel="Rotational/translational KE RMS",yscale="log",title="Tilt and energy partition")
    order=df.case_id.tolist(); x=np.arange(len(df)); bottom=np.zeros(len(df))
    for col,color,label in [("impact_fraction","#D55E00","Impact"),("sliding_fraction","#0072B2","Sliding"),("stuck_fraction","#666666","Stuck")]:
        axes[2].bar(x,df[col],bottom=bottom,color=color,width=.82,label=label); bottom+=df[col].to_numpy()
    axes[2].set(xticks=x,xticklabels=order,ylabel="Fraction of simulated time",title="Contact-state occupancy (no composite score)"); axes[2].tick_params(axis="x",labelrotation=90,labelsize=5.5,pad=1.5)
    for label in axes[2].get_xticklabels(): label.set_rotation_mode("anchor"); label.set_ha("right")
    axes[2].set_ylim(0,1)
    family_start=np.flatnonzero(df.family.ne(df.family.shift()).to_numpy())
    for boundary in family_start[1:]: axes[2].axvline(boundary-.5,color="#D0D0D0",linewidth=.55,zorder=0)
    axes[0].legend(fontsize=6,ncol=2,loc="best",handletextpad=.35,columnspacing=.7)
    axes[2].legend(fontsize=6,ncol=3,loc="upper right",handlelength=1.2,columnspacing=.8)
    for label,ax in zip("abc",axes): ax.text(-.12,1.10,label,transform=ax.transAxes,fontweight="bold",fontsize=8,va="top")
    fig.subplots_adjust(left=.10,right=.985,top=.94,bottom=.22)
    fig.savefig(plots/"motion_screen_quantitative_overview.png",dpi=600,bbox_inches="tight")
    fig.savefig(plots/"motion_screen_quantitative_overview.svg",bbox_inches="tight")
    fig.savefig(plots/"motion_screen_quantitative_overview.pdf",bbox_inches="tight")
    plt.close(fig)

if __name__=="__main__":main()
