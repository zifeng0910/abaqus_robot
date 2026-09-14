"""Render synchronized fixed-camera L2300 and long/L1.8/L2.3 GIF evidence."""
from collections import Counter
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]; ROOT = REPO.parent
AUDIT = REPO/"calibration_analysis"/"ReducedHydro_hidden_impact_audit"
sys.path[:0] = [str(AUDIT), str(AUDIT/"normal_contact_damping_probe")]
from audit_stage_a import mesh_properties
from contact_probe_common import dense_rp
from exact_gap_audit import Wall

mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
                     "font.size":7,"axes.titlesize":8,"xtick.labelsize":5,"ytick.labelsize":5})
CENTER=np.array([-7.468174204284,-3.676918015967,-9.550745259298])


def exterior_faces(mesh):
    count=Counter()
    for _,a,b,c,d in mesh["elems"]:
        for face in ((a,b,c),(a,b,d),(a,c,d),(b,c,d)):count[tuple(sorted(face))]+=1
    return np.asarray([face for face,n in count.items() if n==1],int)


def load_case(folder,job,prefix,title,outcome,color,face_path=None,stride=4):
    meshes,rp,_=mesh_properties((folder/(job+".inp")).read_text());mesh=meshes["Robot_SOLID"]
    faces=(pd.read_csv(face_path)[["n1","n2","n3"]].to_numpy(int) if face_path else exterior_faces(mesh))[::stride]
    polygons=np.asarray([[np.asarray(mesh["nodes"][int(node)])+mesh["shift"] for node in face] for face in faces])
    t,data=dense_rp(folder/"candidate_private")
    return dict(folder=folder,job=job,prefix=prefix,title=title,outcome=outcome,color=color,mesh=mesh,
                rp=np.asarray(rp),polygons=polygons,t=t,data=data,
                axis=pd.read_csv(folder/(prefix+"_true_axis.csv")),
                phase=pd.read_csv(folder/(prefix+"_phase.csv" if prefix=="L2300" else prefix+"_local_phase.csv")),
                gap=pd.read_csv(folder/(prefix+"_exact_gap.csv")),
                bridge=pd.read_csv(folder/(prefix+"_bridge_timeline.csv")),
                events=pd.read_csv(folder/(prefix+"_contact_events.csv")),
                translation=pd.read_csv(folder/(prefix+"_translation.csv" if prefix=="L2300" else prefix+"_canonical_translation.csv")))


def moved(case,i):
    p=case["polygons"];return case["rp"]+case["data"]["U"][i]+Rotation.from_rotvec(case["data"]["UR"][i]).apply(
        (p-case["rp"]).reshape(-1,3)).reshape(p.shape)


def sample_state(case,time):
    i=int(np.searchsorted(case["t"],time).clip(0,len(case["t"])-1));axis=case["axis"];j=min(i,len(axis)-1)
    a=axis.loc[j,["axis_x","axis_y","axis_z"]].to_numpy(float)
    tangent=axis.loc[j,["tangent_x","tangent_y","tangent_z"]].to_numpy(float)
    head=axis.loc[j,["head_x_mm","head_y_mm","head_z_mm"]].to_numpy(float)
    tail=axis.loc[j,["tail_x_mm","tail_y_mm","tail_z_mm"]].to_numpy(float)
    e1=axis.loc[j,["e1_x","e1_y","e1_z"]].to_numpy(float);e2=axis.loc[j,["e2_x","e2_y","e2_z"]].to_numpy(float)
    phase_col="B_local_phase_rad" if "B_local_phase_rad" in case["phase"] else "phi_B_local_rad"
    pb=np.interp(time,case["phase"].time_s,case["phase"][phase_col]);b=np.cos(pb)*e1+np.sin(pb)*e2
    directed=np.degrees(np.arccos(np.clip(np.dot(a,tangent),-1,1)))
    gap_col="min_gap_um" if "min_gap_um" in case["gap"] else "gap_um"
    exact_gap=np.interp(time,case["gap"].time_s,case["gap"][gap_col])
    bridge_col="opposing_bridge_20um" if "opposing_bridge_20um" in case["bridge"] else "opposing_bridge"
    bridge=bool(np.interp(time,case["bridge"].time_s,case["bridge"][bridge_col])>=.5)
    event_ranges=case["events"][["start_s","end_s"]].to_numpy(float)
    contact=bool(np.any((time>=event_ranges[:,0])&(time<=event_ranges[:,1]))) if len(event_ranges) else False
    ds=np.interp(time,case["translation"].time_s,case["translation"].delta_s_mm)
    return i,head,tail,a,tangent,b,directed,exact_gap,contact,bridge,ds


def configure(ax,title):
    half=2.45;ax.set(xlim=(CENTER[0]-half,CENTER[0]+half),ylim=(CENTER[1]-half,CENTER[1]+half),zlim=(CENTER[2]-half,CENTER[2]+half))
    ax.set_box_aspect((1,1,1));ax.view_init(elev=19,azim=-57);ax.set_title(title,pad=2)
    ax.set_xlabel("x (mm)",labelpad=-3);ax.set_ylabel("y (mm)",labelpad=-3);ax.set_zlabel("z (mm)",labelpad=-4)


def draw(ax,case,time,wall_polys,curve):
    i,head,tail,axis,tangent,b,directed,gap,contact,bridge,ds=sample_state(case,time)
    ax.clear();ax.add_collection3d(Poly3DCollection(wall_polys,facecolor="#56B4E9",edgecolor="none",alpha=.09))
    ax.add_collection3d(Poly3DCollection(moved(case,i),facecolor=case["color"],edgecolor="#30343B",lw=.08,alpha=.87))
    ax.plot(curve[:,0],curve[:,1],curve[:,2],color="#667085",lw=.7,ls="--")
    mid=(head+tail)/2;ax.scatter(*head,s=19,color="#D55E00",depthshade=False);ax.scatter(*tail,s=19,color="#009E73",depthshade=False)
    for vector,color,scale in ((axis,"#0072B2",.9),(tangent,"#667085",.75),(b,"#CC3311",.7)):
        ax.quiver(*mid,*(vector*scale),color=color,linewidth=1.1,arrow_length_ratio=.16)
    ax.text2D(.02,.98,f"t={time*1e3:.4f} ms\ndirected tilt={directed:.1f} deg\nexact CAD gap={gap:+.2f} um\n"
              f"contact={'ON' if contact else 'off'} | bridge20={'YES' if bridge else 'no'}\nDelta s={ds:+.3f} mm\n{case['outcome']}",
              transform=ax.transAxes,va="top",fontsize=6)
    configure(ax,case["title"])


def main():
    decision=pd.read_csv(HERE/"L2300_runtime_identity.csv").iloc[0].classification
    current=load_case(HERE,"Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L2300_D0815_WallSupported_0083","L2300",
                      "L2.3 frozen CAD",decision,"#0072B2",HERE/"L2300_robot_surface_triangles.csv",3)
    long_dir=REPO/"calibration_analysis"/"ReducedHydro_zeta050_8p333_validation"
    long=load_case(long_dir,"Wobble_F30_G6L45_ReducedHydro_Zeta050_WallOn_Free_0083","zeta050_8p333",
                   "Long L~2.9","TOO LONG / JAM","#7A7F87",ROOT/"Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083_robot_surface_triangles_exact.csv",5)
    short_dir=REPO/"calibration_analysis"/"ReducedHydro_FreeCAD_L1800_validation"
    short=load_case(short_dir,"Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_WallOn_Free_0083","freecad_8p333",
                    "L1.8 frozen CAD","TOO SHORT / TUMBLE","#D55E00",short_dir/"freecad_robot_surface_triangles_exact.csv",4)
    wall=Wall();wall_polys=wall.tri[np.linalg.norm(wall.centers-CENTER,axis=1)<4.0][::3]
    curve=pd.read_csv(ROOT/"CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_TRUE_centerline_odb.csv")[["x_mm","y_mm","z_mm"]].to_numpy(float)
    curve=curve[np.linalg.norm(curve-CENTER,axis=1)<4.5]
    times=np.unique(np.r_[np.linspace(0,.008333,105),np.linspace(.0008,.0022,38),np.linspace(.006,.008333,38)])

    fig=plt.figure(figsize=(5.2,4.25),constrained_layout=True);ax=fig.add_subplot(111,projection="3d")
    def single_frame(k):
        draw(ax,current,times[k],wall_polys,curve);ax.set_title("L2.3 wall-supported wobble diagnostic | fixed camera",pad=2)
    animation=FuncAnimation(fig,single_frame,frames=len(times),interval=70,repeat=True)
    animation.save(HERE/"L2300_WallSupportedWobble_8p333.gif",writer=PillowWriter(fps=14));plt.close(fig)

    fig=plt.figure(figsize=(10.8,3.55),constrained_layout=True);axes=[fig.add_subplot(1,3,i+1,projection="3d") for i in range(3)]
    def compare_frame(k):
        for panel,case in zip(axes,(long,short,current)):draw(panel,case,times[k],wall_polys,curve)
        fig.suptitle("Length-regime comparison | synchronized time, camera, wall and vector scales",fontsize=9)
    animation=FuncAnimation(fig,compare_frame,frames=len(times),interval=70,repeat=True)
    animation.save(HERE/"Long_vs_L1800_vs_L2300_WallSupported_8p333.gif",writer=PillowWriter(fps=14));plt.close(fig)


if __name__=="__main__":main()
