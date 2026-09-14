"""Publication-grade mesh, profile, and offline replay evidence figures."""
from pathlib import Path
import json
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import pandas as pd


HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
AUDIT=REPO/"calibration_analysis"/"ReducedHydro_hidden_impact_audit"
sys.path.insert(0,str(AUDIT))
from audit_stage_a import mesh_properties

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
 "svg.fonttype":"none","pdf.fonttype":42,"font.size":7,"axes.titlesize":8,
 "axes.spines.right":False,"axes.spines.top":False,"axes.linewidth":.8,"legend.frameon":False})


def save(fig,name):
    fig.savefig(HERE/f"{name}.png",dpi=400,bbox_inches="tight")
    fig.savefig(HERE/f"{name}.pdf",bbox_inches="tight")
    fig.savefig(HERE/f"{name}.svg",bbox_inches="tight")
    fig.savefig(HERE/f"{name}.tiff",dpi=600,bbox_inches="tight")


def mesh_figure():
    dense_dir=REPO/"calibration_analysis"/"ReducedHydro_FreeCAD_L1800_validation"
    dense_job="Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_WallOn_Free_0083"
    simple_job="Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_Mesh060"
    meshes=[]
    for folder,job,tri_file in ((dense_dir,dense_job,dense_dir/"freecad_robot_surface_triangles_exact.csv"),
                                (HERE,simple_job,HERE/"current_mesh060_surface_triangles.csv")):
        mesh_properties_map,_,_=mesh_properties((folder/f"{job}.inp").read_text()); mesh=mesh_properties_map["Robot_SOLID"]
        tri=pd.read_csv(tri_file)[["n1","n2","n3"]].to_numpy(int)
        polys=np.asarray([[np.asarray(mesh["nodes"][int(n)])+mesh["shift"]-mesh["com"] for n in face] for face in tri])
        meshes.append(polys)
    dense=json.loads((dense_dir/"freecad_preflight_identity.json").read_text())
    simple=json.loads((HERE/"current_mesh060_audit.json").read_text())
    fig=plt.figure(figsize=(7.2,3.1),constrained_layout=True)
    labels=[f"Dense 0.035 mm\n{dense['abaqus_import']['node_count']:,} nodes | {dense['abaqus_import']['element_count']:,} C3D4\n"
            f"{dense['surface_triangle_count']:,} exterior triangles | P95 {dense['surface_normal_P95_deg']:.3f} deg",
            f"Simplified 0.060 mm\n{simple['node_count']:,} nodes | {simple['element_count']:,} C3D4\n"
            f"{simple['exterior_triangle_count']:,} exterior triangles | P95 {simple['surface_normal_P95_deg']:.3f} deg"]
    for i,(polys,title,color) in enumerate(zip(meshes,labels,("#7A7F87","#0072B2"))):
        ax=fig.add_subplot(1,2,i+1,projection="3d")
        ax.add_collection3d(Poly3DCollection(polys,facecolor=color,edgecolor=color,lw=.12,alpha=.10))
        ax.set(xlim=(-1.0,.85),ylim=(-.48,.48),zlim=(-.48,.48));ax.set_box_aspect((1.85,.96,.96));ax.view_init(elev=18,azim=-62)
        ax.set_title(title,pad=3);ax.set_axis_off()
    fig.suptitle("Rigid-body volume mesh reduction preserves the rounded HEAD contact surface",fontsize=9)
    save(fig,"Dense_vs_Simplified_RobotMesh");plt.close(fig)


def profile_figure():
    p=pd.read_csv(HERE/"head_clearance_profile_fit.csv");fit=json.loads((HERE/"head_clearance_fit.json").read_text())
    fig,ax=plt.subplots(figsize=(3.55,2.65),constrained_layout=True)
    ax.plot(p.nose_axial_x_mm,p.old_L1800_envelope_mm,color="#5B6573",lw=1.3,label="Old L1800 envelope")
    ax.plot(p.nose_axial_x_mm,p.current_FreeCAD_mm,color="#D55E00",lw=1.4,label="Current exact CAD")
    ax.plot(p.nose_axial_x_mm,p.HeadClearance_candidate_mm,color="#0072B2",lw=1.6,label="HeadClearance arc")
    ax.set(xlabel="Distance from nose tip (mm)",ylabel="Radial envelope (mm)",
           title="One-parameter circular-arc head fit")
    ax.legend(loc="lower right",fontsize=6);ax.grid(color="#D9DDE2",lw=.5)
    ax.text(.02,.98,f"Rhead = {fit['fitted_R_head_design_mm']:.6f} mm\nmax outward excess = {fit['max_outward_excess_um']:.3f} um",
            transform=ax.transAxes,va="top",fontsize=6)
    save(fig,"HeadClearance_profile_fit");plt.close(fig)


def replay_figure():
    mesh=pd.read_csv(HERE/"mesh060_replay_summary.csv")
    candidate=pd.read_csv(HERE/"head_clearance_replay_summary.csv")
    fig,axes=plt.subplots(1,2,figsize=(7.2,2.65),constrained_layout=True)
    x=np.arange(2);w=.19
    for j,(col,label,color) in enumerate((("longest20_ms","20 um","#0072B2"),("longest10_ms","10 um","#EE7733"),
                                         ("longest5_ms","5 um","#009E73"),("longest0_ms","0 um","#5B6573"))):
        axes[0].bar(x+(j-1.5)*w,mesh[col],w,label=label,color=color)
        axes[1].bar(x+(j-1.5)*w,candidate[col],w,label=label,color=color)
    axes[0].set_xticks(x,mesh.geometry);axes[1].set_xticks(x,candidate.trajectory.map({"OLD":"Old trajectory","NEW":"Current trajectory"}))
    axes[0].set_title("Mesh regression | current trajectory");axes[1].set_title("HeadClearance zero-Abaqus replay")
    for ax in axes:
        ax.axhline(.5,color="#CC3311",ls="--",lw=.8);ax.set_ylabel("Longest opposing bridge (ms)");ax.set_ylim(bottom=0)
    axes[1].legend(ncol=2,fontsize=6,loc="upper right")
    save(fig,"Offline_bridge_replay_gates");plt.close(fig)


def main():
    mesh_figure();profile_figure();replay_figure()


if __name__=="__main__":main()
