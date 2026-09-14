"""Create synchronized fixed-camera evidence GIFs from existing trajectories only."""
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
REPO = HERE.parents[1]
ROOT = REPO.parent
AUDIT = REPO / "calibration_analysis" / "ReducedHydro_hidden_impact_audit"
sys.path[:0] = [str(AUDIT), str(AUDIT / "normal_contact_damping_probe")]
from audit_stage_a import mesh_properties
from contact_probe_common import dense_rp
from exact_gap_audit import Wall

mpl.rcParams.update({"font.family":"sans-serif", "font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
                     "font.size":7, "axes.titlesize":8, "xtick.labelsize":5, "ytick.labelsize":5,
                     "svg.fonttype":"none", "pdf.fonttype":42})


def exterior_faces(mesh):
    count = Counter()
    for _, a,b,c,d in mesh["elems"]:
        for face in ((a,b,c),(a,b,d),(a,c,d),(b,c,d)): count[tuple(sorted(face))] += 1
    return np.asarray([face for face,n in count.items() if n == 1], int)


def model_case(folder, job, prefix, title, color, face_file=None, stride=3):
    meshes, rp, _ = mesh_properties((folder/f"{job}.inp").read_text()); mesh=meshes["Robot_SOLID"]
    faces = (pd.read_csv(face_file)[["n1","n2","n3"]].to_numpy(int)
             if face_file else exterior_faces(mesh))
    faces = faces[::stride]
    base = np.asarray([[np.asarray(mesh["nodes"][int(n)])+mesh["shift"] for n in face] for face in faces])
    time, dense = dense_rp(folder/"candidate_private")
    return {"folder":folder,"mesh":mesh,"rp":np.asarray(rp),"base":base,"time":time,"dense":dense,
            "axis":pd.read_csv(folder/f"{prefix}_true_axis.csv"),
            "phase":pd.read_csv(folder/f"{prefix}_local_phase.csv"),
            "gap":pd.read_csv(folder/f"{prefix}_exact_gap.csv"),
            "bridge":pd.read_csv(folder/f"{prefix}_bridge_timeline.csv"),
            "events":pd.read_csv(folder/f"{prefix}_contact_events.csv"),
            "translation":pd.read_csv(folder/f"{prefix}_canonical_translation.csv"),
            "title":title,"color":color}


def moved(case, i, base=None):
    points = case["base"] if base is None else base
    return case["rp"] + case["dense"]["U"][i] + Rotation.from_rotvec(case["dense"]["UR"][i]).apply(
        (points-case["rp"]).reshape(-1,3)).reshape(points.shape)


def state(case, t):
    i=int(np.searchsorted(case["time"],t).clip(0,len(case["time"])-1)); axis=case["axis"]
    j=min(i,len(axis)-1)
    a=axis.loc[j,["axis_x","axis_y","axis_z"]].to_numpy(float)
    tangent=axis.loc[j,["tangent_x","tangent_y","tangent_z"]].to_numpy(float)
    head=axis.loc[j,["head_x_mm","head_y_mm","head_z_mm"]].to_numpy(float)
    tail=axis.loc[j,["tail_x_mm","tail_y_mm","tail_z_mm"]].to_numpy(float)
    e1=axis.loc[j,["e1_x","e1_y","e1_z"]].to_numpy(float); e2=axis.loc[j,["e2_x","e2_y","e2_z"]].to_numpy(float)
    assert np.all(np.diff(case["phase"].time_s)>0)
    phase=np.interp(t,case["phase"].time_s,case["phase"].phi_B_local_rad); b=np.cos(phase)*e1+np.sin(phase)*e2
    directed=np.degrees(np.arccos(np.clip(np.dot(a,tangent),-1,1))); folded=min(directed,180-directed)
    intervals=case["events"][["start_s","end_s"]].to_numpy(float)
    contact=bool(np.any((t>=intervals[:,0])&(t<=intervals[:,1]))) if len(intervals) else False
    bridge=bool(np.interp(t,case["bridge"].time_s,case["bridge"].opposing_bridge)>=.5)
    assert np.all(np.diff(case["gap"].time_s)>0) and np.all(np.diff(case["translation"].time_s)>0)
    gap=np.interp(t,case["gap"].time_s,case["gap"].gap_um)
    ds=np.interp(t,case["translation"].time_s,case["translation"].delta_s_mm)
    return i,head,tail,a,tangent,b,directed,folded,gap,contact,bridge,ds


def configure(ax, title):
    half=2.35; center=np.array([-7.468174204284,-3.676918015967,-9.550745259298])
    ax.set(xlim=(center[0]-half,center[0]+half),ylim=(center[1]-half,center[1]+half),zlim=(center[2]-half,center[2]+half))
    ax.set_box_aspect((1,1,1)); ax.view_init(elev=19,azim=-57); ax.set_title(title,pad=2)
    ax.set_xlabel("x (mm)",labelpad=-3);ax.set_ylabel("y (mm)",labelpad=-3);ax.set_zlabel("z (mm)",labelpad=-4)


def draw_case(ax,case,t,wall_polys,curve):
    i,head,tail,axis,tangent,b,directed,folded,gap,contact,bridge,ds=state(case,t)
    ax.clear(); robot=moved(case,i)
    ax.add_collection3d(Poly3DCollection(wall_polys,facecolor="#56B4E9",edgecolor="none",alpha=.09))
    ax.add_collection3d(Poly3DCollection(robot,facecolor=case["color"],edgecolor="#30343B",lw=.08,alpha=.86))
    ax.plot(curve[:,0],curve[:,1],curve[:,2],color="#5B6573",lw=.7,ls="--")
    mid=(head+tail)/2;ax.scatter(*head,s=18,color="#D55E00",depthshade=False);ax.scatter(*tail,s=18,color="#009E73",depthshade=False)
    for vector,color,scale in ((axis,"#0072B2",.9),(tangent,"#5B6573",.75),(b,"#CC3311",.7)):
        ax.quiver(*mid,*(vector*scale),color=color,linewidth=1.1,arrow_length_ratio=.16)
    ax.text2D(.02,.98,f"t={t*1e3:.4f} ms\ndirected/folded={directed:.1f}/{folded:.1f} deg\n"
              f"gap={gap:+.2f} um | contact={'ON' if contact else 'off'}\n"
              f"bridge20={'YES' if bridge else 'no'} | delta s={ds:+.3f} mm",
              transform=ax.transAxes,va="top",fontsize=6)
    configure(ax,case["title"])


def main():
    long_dir=REPO/"calibration_analysis"/"ReducedHydro_zeta050_8p333_validation"
    old_dir=REPO/"calibration_analysis"/"ReducedHydro_geometry_L1800_validation"
    new_dir=REPO/"calibration_analysis"/"ReducedHydro_FreeCAD_L1800_validation"
    long=model_case(long_dir,"Wobble_F30_G6L45_ReducedHydro_Zeta050_WallOn_Free_0083","zeta050_8p333",
                    "Long baseline | actual L~2.9 mm","#7A7F87",ROOT/"Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083_robot_surface_triangles_exact.csv",5)
    old=model_case(old_dir,"Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083","L1800_8p333",
                   "Old scaled mesh | L=1.800 mm","#009E73",None,2)
    new=model_case(new_dir,"Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_WallOn_Free_0083","freecad_8p333",
                   "Exact FreeCAD circular head | L=1.800 mm","#0072B2",new_dir/"freecad_robot_surface_triangles_exact.csv",4)
    wall=Wall(); center=np.array([-7.468174204284,-3.676918015967,-9.550745259298])
    wall_polys=wall.tri[np.linalg.norm(wall.centers-center,axis=1)<4.0][::3]
    curve=pd.read_csv(ROOT/"CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_TRUE_centerline_odb.csv")[["x_mm","y_mm","z_mm"]].to_numpy(float)
    curve=curve[np.linalg.norm(curve-center,axis=1)<4.5]
    times=np.unique(np.r_[np.linspace(0,.008333,90),np.linspace(0,.0015,42),
                          np.linspace(.0025,.0045,52),np.linspace(.0058,.0070,48)])
    fig=plt.figure(figsize=(10.8,3.55),constrained_layout=True);axes=[fig.add_subplot(1,3,i+1,projection="3d") for i in range(3)]
    def frame(k):
        for ax,case in zip(axes,(long,old,new)):draw_case(ax,case,times[k],wall_polys,curve)
        fig.suptitle("Shortening to L=1.800 mm suppresses persistent opposing-wall bridge | synchronized fixed camera",fontsize=9)
    frame(0)
    fig.savefig(HERE/"Long_vs_L1800Scaled_vs_L1800FreeCAD_keyframe.png",dpi=300,bbox_inches="tight")
    fig.savefig(HERE/"Long_vs_L1800Scaled_vs_L1800FreeCAD_keyframe.pdf",bbox_inches="tight")
    fig.savefig(HERE/"Long_vs_L1800Scaled_vs_L1800FreeCAD_keyframe.svg",bbox_inches="tight")
    if "--static-only" in sys.argv:
        plt.close(fig);return
    ani=FuncAnimation(fig,frame,frames=len(times),interval=65,repeat=True)
    ani.save(HERE/"Long_vs_L1800Scaled_vs_L1800FreeCAD_8p333.gif",writer=PillowWriter(fps=15));plt.close(fig)

    slow=np.linspace(.0058,.0070,97); fig=plt.figure(figsize=(7.2,3.55),constrained_layout=True)
    axes=[fig.add_subplot(1,2,i+1,projection="3d") for i in range(2)]
    full_refs=[]
    for case in (old,new):
        faces=exterior_faces(case["mesh"]); ids=np.unique(faces)
        points=np.asarray([np.asarray(case["mesh"]["nodes"][int(n)])+case["mesh"]["shift"] for n in ids])
        full_refs.append(points)
    near_cache=[]
    for case,base in zip((old,new),full_refs):
        rows=[]
        for t in slow:
            i=int(np.searchsorted(case["time"],t).clip(0,len(case["time"])-1)); moved_points=moved(case,i,base[:,None,:])[:,0,:]
            gaps,tris,closest=wall.query(moved_points); mask=gaps<=.020
            rows.append((moved_points[mask],np.unique(tris[mask]),float(gaps.min()*1e3)))
        near_cache.append(rows)
    def slow_frame(k):
        for panel,case,cache in zip(axes,(old,new),near_cache):
            draw_case(panel,case,slow[k],wall_polys,curve); points,sectors,gap=cache[k]
            if len(points): panel.scatter(points[:,0],points[:,1],points[:,2],s=8,color="#CC3311",depthshade=False)
            if len(sectors): panel.add_collection3d(Poly3DCollection(wall.tri[sectors],facecolor="#EE7733",edgecolor="none",alpha=.42))
        fig.suptitle("5.8-7.0 ms slow motion | red nodes and orange wall sectors are within 20 um",fontsize=8)
    slow_frame(0);ani=FuncAnimation(fig,slow_frame,frames=len(slow),interval=85,repeat=True)
    ani.save(HERE/"GeometryImprovement_SlowMotion_5p8_to_7p0ms.gif",writer=PillowWriter(fps=12));plt.close(fig)


if __name__=="__main__":main()
