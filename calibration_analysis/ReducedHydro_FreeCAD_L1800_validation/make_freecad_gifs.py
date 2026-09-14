"""Generate exact-CAD and old-vs-new fixed-camera validation animations."""
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
OLD = HERE.parent / "ReducedHydro_geometry_L1800_validation"
AUDIT = HERE.parent / "ReducedHydro_hidden_impact_audit"
sys.path[:0] = [str(AUDIT), str(AUDIT / "normal_contact_damping_probe")]
from audit_stage_a import mesh_properties
from exact_gap_audit import Wall
from contact_probe_common import dense_rp

NEW_JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_WallOn_Free_0083"
OLD_JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083"
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 7, "axes.titlesize": 8, "xtick.labelsize": 5, "ytick.labelsize": 5})


def exterior_faces(mesh):
    count = Counter()
    for _, a, b, c, d in mesh["elems"]:
        for face in ((a, b, c), (a, b, d), (a, c, d), (b, c, d)):
            count[tuple(sorted(face))] += 1
    return np.asarray([face for face, number in count.items() if number == 1], int)


def reference_polys(mesh, faces, stride):
    return np.asarray([[np.asarray(mesh["nodes"][int(node)]) + mesh["shift"] for node in face]
                       for face in faces[::stride]])


def move(polys, rp, displacement, rotvec):
    rotation = Rotation.from_rotvec(rotvec)
    return rp + displacement + rotation.apply((polys-rp).reshape(-1, 3)).reshape(polys.shape)


def configure(ax, center, title):
    half = 2.45
    ax.set(xlim=(center[0]-half, center[0]+half), ylim=(center[1]-half, center[1]+half),
           zlim=(center[2]-half, center[2]+half))
    ax.set_box_aspect((1, 1, 1)); ax.view_init(elev=19, azim=-57)
    ax.set_title(title, pad=2); ax.set_xlabel("x (mm)", labelpad=-3)
    ax.set_ylabel("y (mm)", labelpad=-3); ax.set_zlabel("z (mm)", labelpad=-4)


def draw_vectors(ax, axis, phase, i):
    head = axis.loc[i, ["head_x_mm", "head_y_mm", "head_z_mm"]].to_numpy(float)
    tail = axis.loc[i, ["tail_x_mm", "tail_y_mm", "tail_z_mm"]].to_numpy(float)
    center = (head+tail)/2
    tangent = axis.loc[i, ["tangent_x", "tangent_y", "tangent_z"]].to_numpy(float)
    robot_axis = axis.loc[i, ["axis_x", "axis_y", "axis_z"]].to_numpy(float)
    e1 = axis.loc[i, ["e1_x", "e1_y", "e1_z"]].to_numpy(float)
    e2 = axis.loc[i, ["e2_x", "e2_y", "e2_z"]].to_numpy(float)
    phi_b = phase.phi_B_local_rad.iloc[i]
    bvec = np.cos(phi_b)*e1 + np.sin(phi_b)*e2
    ax.scatter(*head, s=20, color="#D55E00", depthshade=False, label="HEAD")
    ax.scatter(*tail, s=20, color="#009E73", depthshade=False, label="TAIL")
    for vector, color, label, scale in ((robot_axis,"#0072B2","Directed axis",1.0),
                                         (tangent,"#5B6573","Local tangent",.85),
                                         (bvec,"#CC3311","B",.75)):
        ax.quiver(*center, *(vector*scale), color=color, linewidth=1.2,
                  arrow_length_ratio=.16, label=label)


def initial_polarity(axis):
    return np.sign(np.dot(axis.loc[0,["axis_x","axis_y","axis_z"]],
                          axis.loc[0,["tangent_x","tangent_y","tangent_z"]]))


def directed_state(axis, i):
    dot = np.dot(axis.loc[i,["axis_x","axis_y","axis_z"]],
                 axis.loc[i,["tangent_x","tangent_y","tangent_z"]])
    return np.degrees(np.arccos(np.clip(dot,-1,1))), ("INITIAL" if np.sign(dot)==initial_polarity(axis) else "REVERSED")


def main():
    new_meshes, new_rp, _ = mesh_properties((HERE/f"{NEW_JOB}.inp").read_text())
    old_meshes, old_rp, _ = mesh_properties((OLD/f"{OLD_JOB}.inp").read_text())
    new_mesh, old_mesh = new_meshes["Robot_SOLID"], old_meshes["Robot_SOLID"]
    new_faces = pd.read_csv(HERE/"freecad_robot_surface_triangles_exact.csv")[["n1","n2","n3"]].to_numpy(int)
    new_ref = reference_polys(new_mesh, new_faces, 4)
    old_ref = reference_polys(old_mesh, exterior_faces(old_mesh), 2)
    time, new_dense = dense_rp(HERE/"candidate_private")
    old_time, old_dense = dense_rp(OLD/"candidate_private")
    axis = pd.read_csv(HERE/"freecad_8p333_true_axis.csv")
    old_axis = pd.read_csv(OLD/"L1800_8p333_true_axis.csv")
    phase = pd.read_csv(HERE/"freecad_8p333_local_phase.csv")
    old_phase = pd.read_csv(OLD/"L1800_8p333_local_phase.csv")
    field = pd.read_csv(HERE/"freecad_8p333_field_orientation.csv")
    gap = pd.read_csv(HERE/"freecad_8p333_exact_gap.csv")
    bridge = pd.read_csv(HERE/"freecad_8p333_bridge_timeline.csv")
    contacts = pd.read_csv(HERE/"freecad_8p333_contact_events.csv")
    reversals = pd.read_csv(HERE/"freecad_head_tail_reversal_events.csv")
    wall = Wall(); center = new_ref.reshape(-1,3).mean(axis=0)
    wall_polys = wall.tri[np.linalg.norm(wall.centers-center,axis=1)<4.0][::3]
    curve = pd.read_csv(ROOT/"CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_TRUE_centerline_odb.csv")[["x_mm","y_mm","z_mm"]].to_numpy(float)
    curve = curve[np.linalg.norm(curve-center,axis=1)<4.5]
    special = np.r_[reversals.crossing_time_s,(contacts.start_s+contacts.end_s)/2,np.linspace(.0026,.0031,15)]
    frame_times = np.unique(np.r_[np.linspace(0,time[-1],120),special]); frame_times.sort()
    indices = np.searchsorted(time,frame_times).clip(0,len(time)-1)
    intervals = contacts[["start_s","end_s"]].to_numpy(float)

    fig=plt.figure(figsize=(5.2,4.25),constrained_layout=True); ax=fig.add_subplot(111,projection="3d")
    def draw_main(k):
        ax.clear(); i=int(indices[k]); t=time[i]
        polys=move(new_ref,new_rp,new_dense["U"][i],new_dense["UR"][i])
        ax.add_collection3d(Poly3DCollection(wall_polys,facecolor="#56B4E9",edgecolor="none",alpha=.10))
        ax.add_collection3d(Poly3DCollection(polys,facecolor="#0072B2",edgecolor="#30343B",lw=.10,alpha=.86))
        ax.plot(curve[:,0],curve[:,1],curve[:,2],color="#5B6573",lw=.8,ls="--"); draw_vectors(ax,axis,phase,i)
        tilt,polarity=directed_state(axis,i); active=bool(np.any((t>=intervals[:,0])&(t<=intervals[:,1])))
        ax.text2D(.02,.98,f"t = {t*1e3:.4f} ms\ndirected tilt = {tilt:.2f} deg\nfield tilt = {field.theta_B_deg.iloc[i]:.2f} deg\n"
                  f"polarity = {polarity}\nexact gap = {np.interp(t,gap.time_s,gap.gap_um):+.2f} um\n"
                  f"contact = {'ACTIVE' if active else 'off'} | bridge20 = {'YES' if bridge.opposing_bridge.iloc[i] else 'no'}",
                  transform=ax.transAxes,va="top"); configure(ax,center,"Exact FreeCAD L1800 | fixed camera")
    ani=FuncAnimation(fig,draw_main,frames=len(indices),interval=75,repeat=True)
    ani.save(HERE/"Wobble_F30_FreeCAD_L1800_D0815_ReducedHydro_8p333.gif",writer=PillowWriter(fps=13)); plt.close(fig)

    sync=np.unique(np.r_[np.linspace(0,min(time[-1],old_time[-1]),105),reversals.crossing_time_s,np.linspace(.0026,.0031,12)])
    ni=np.searchsorted(time,sync).clip(0,len(time)-1); oi=np.searchsorted(old_time,sync).clip(0,len(old_time)-1)
    fig=plt.figure(figsize=(7.2,3.55),constrained_layout=True); axes=[fig.add_subplot(1,2,j+1,projection="3d") for j in range(2)]
    def draw_compare(k):
        cases=((axes[0],old_ref,old_rp,old_dense,oi[k],old_axis,old_phase,"Old scaled L1800","#5B6573"),
               (axes[1],new_ref,new_rp,new_dense,ni[k],axis,phase,"Exact FreeCAD L1800","#0072B2"))
        for panel,ref,rp,dense,i,axial,ph,title,color in cases:
            panel.clear(); polys=move(ref,rp,dense["U"][i],dense["UR"][i])
            panel.add_collection3d(Poly3DCollection(wall_polys,facecolor="#56B4E9",edgecolor="none",alpha=.09))
            panel.add_collection3d(Poly3DCollection(polys,facecolor=color,edgecolor="#30343B",lw=.10,alpha=.86))
            panel.plot(curve[:,0],curve[:,1],curve[:,2],color="#5B6573",lw=.7,ls="--"); draw_vectors(panel,axial,ph,int(i))
            tilt,polarity=directed_state(axial,int(i)); panel.text2D(.02,.97,f"t = {sync[k]*1e3:.4f} ms\ndirected tilt = {tilt:.2f} deg\npolarity = {polarity}",transform=panel.transAxes,va="top")
            configure(panel,center,title)
        fig.suptitle("Old scaled versus exact FreeCAD geometry and trajectories | fixed camera",fontsize=8)
    ani=FuncAnimation(fig,draw_compare,frames=len(sync),interval=80,repeat=True)
    ani.save(HERE/"OldScaledL1800_vs_FreeCADL1800_8p333.gif",writer=PillowWriter(fps=12)); plt.close(fig)


if __name__=="__main__":
    main()
