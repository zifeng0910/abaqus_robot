"""Publication-grade preflight evidence for the stopped L2300 candidate."""
from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE=Path(__file__).resolve().parent
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
 "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.titlesize":8,
 "axes.spines.right":False,"axes.spines.top":False,"axes.linewidth":.8,
 "legend.frameon":False,"savefig.facecolor":"white"})


def save(fig,name):
    fig.savefig(HERE/f"{name}.png",dpi=600,bbox_inches="tight")
    fig.savefig(HERE/f"{name}.pdf",bbox_inches="tight")
    fig.savefig(HERE/f"{name}.svg",bbox_inches="tight")


def geometry_figure():
    p=pd.read_csv(HERE/"L2300_head_profile_fit.csv")
    g=json.loads((HERE/"Robot_parametric_L2p300_D0p815_WallWobble_geometry.json").read_text())
    x=p.nose_axial_x_mm.to_numpy(); old=p.old_L1800_convex_envelope_mm.to_numpy(); new=p.L2300_Bezier_head_mm.to_numpy()
    fig,ax=plt.subplots(figsize=(7.2,2.75),constrained_layout=True)
    ax.plot(x,old,color="#5B6573",lw=1.4,label="Old successful L1800 envelope")
    ax.plot(x,new,color="#0072B2",lw=1.6,label="L2300 smooth head")
    ax.plot([x[-1],1.8],[.4075,.4075],color="#5B6573",lw=1.4)
    ax.plot([x[-1],2.3],[.4075,.4075],color="#0072B2",lw=1.6)
    ax.vlines([1.8,2.3],0,.4075,colors=["#5B6573","#0072B2"],lw=[1.4,1.6])
    ax.fill_between([1.8,2.3],0,.4075,color="#009E73",alpha=.12)
    ax.text(2.05,.20,"+0.500 mm\nstraight cylinder",ha="center",va="center",color="#006D4F")
    ax.set(xlabel="Axial position from nose (mm)",ylabel="Radial envelope (mm)",xlim=(-.02,2.34),ylim=(0,.445),
           title="L2300 preserves the old head envelope and extends only the straight body")
    ax.legend(loc="lower center",ncol=2);ax.grid(color="#D9DDE2",lw=.45)
    ax.text(.99,.04,f"max outward {g['fit_source']['max_outward_excess_um']:.2f} um | RMS {g['fit_source']['profile_rmse_um']:.2f} um",
            transform=ax.transAxes,ha="right",fontsize=6)
    save(fig,"L1800_vs_L2300_geometry");plt.close(fig)


def mesh_figure():
    trials=pd.read_csv(HERE/"L2300_mesh_trial_summary.csv")
    fig,axes=plt.subplots(1,2,figsize=(7.2,2.75),constrained_layout=True)
    axes[0].plot(trials.trial,trials.HEAD_P95_deg,"o-",color="#D55E00",lw=1.4,label="HEAD P95")
    axes[0].plot(trials.trial,trials.overall_P95_deg,"o-",color="#0072B2",lw=1.2,label="Overall P95")
    axes[0].axhline(2,color="#CC3311",ls="--",lw=.9,label="Hard gate (2 deg)")
    axes[0].set(xlabel="Mesh audit trial",ylabel="Normal error P95 (deg)",xticks=trials.trial,
                title="HEAD faceting remains above the hard gate")
    axes[0].legend(fontsize=6);axes[0].grid(color="#D9DDE2",lw=.45)
    axes[1].bar(trials.trial,trials.elements/1000,color="#7A7F87",width=.65)
    axes[1].axhspan(25,35,color="#009E73",alpha=.13,label="Expected 25k-35k")
    axes[1].set(xlabel="Mesh audit trial",ylabel="C3D4 elements (thousands)",xticks=trials.trial,
                title="Curvature refinement cannot satisfy both targets")
    axes[1].legend(fontsize=6);axes[1].grid(axis="y",color="#D9DDE2",lw=.45)
    save(fig,"L2300_mesh_gate_failure");plt.close(fig)


if __name__=="__main__":geometry_figure();mesh_figure()
