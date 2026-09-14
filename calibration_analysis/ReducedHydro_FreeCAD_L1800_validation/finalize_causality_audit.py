"""Create pose-specific causal-audit tables from completed replay products."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import cross_geometry_trajectory_replay as replay


def main():
    summary=pd.read_csv(HERE/"cross_replay_bridge_summary.csv")
    timeline=pd.read_csv(HERE/"cross_replay_bridge_timeline.csv")
    d=summary[(summary.geometry=="NEW")&(summary.trajectory=="NEW")].iloc[0]
    segment=timeline[(timeline.geometry=="NEW")&(timeline.trajectory=="NEW")&
                     (timeline.time_s>=d.longest20_start_s)&(timeline.time_s<=d.longest20_end_s)]
    targets={"onset":d.longest20_start_s,"midpoint":(d.longest20_start_s+d.longest20_end_s)/2,
             "maximum_severity":segment.loc[segment.minimum_gap_um.idxmin(),"time_s"]}
    meshes,rp,_=replay.mesh_properties((HERE/f"{replay.NEW_JOB}.inp").read_text()); mesh=meshes["Robot_SOLID"]
    ids=np.unique(pd.read_csv(HERE/"freecad_robot_surface_triangles_exact.csv")[["n1","n2","n3"]].to_numpy())
    points0=np.asarray([mesh["nodes"][int(i)] for i in ids])+mesh["shift"]
    seed=np.array([.9647382600216,-.118874237214,.2348382536499]); axial=(points0-mesh["com"])@seed
    lo,hi=np.quantile(axial,[.10,.90]); parts=np.where(axial<=lo,"HEAD",np.where(axial>=hi,"TAIL","BODY"))
    time,data=replay.dense_rp(HERE/"candidate_private"); axis=pd.read_csv(HERE/"freecad_8p333_true_axis.csv")
    phase=pd.read_csv(HERE/"freecad_8p333_local_phase.csv"); wall=replay.Wall(); rows=[]; nodes=[]
    for label,target in targets.items():
        i=int(np.argmin(abs(time-target))); rotation=Rotation.from_rotvec(data["UR"][i])
        com=rp+data["U"][i]+rotation.apply(mesh["com"]-rp)
        points=com+rotation.apply(points0-mesh["com"]); gaps,tri,nearest=replay.exact_query(wall,points)
        e1=axis.loc[i,["e1_x","e1_y","e1_z"]].to_numpy(float); e2=axis.loc[i,["e2_x","e2_y","e2_z"]].to_numpy(float)
        selected=gaps*1e3<=20; angles=np.mod(np.arctan2(wall.normals[tri[selected]]@e2,wall.normals[tri[selected]]@e1),2*np.pi)
        centers=replay.clusters(angles); separation=max((np.degrees(np.arccos(np.clip(np.cos(a-b),-1,1))) for q,a in enumerate(centers) for b in centers[q+1:]),default=0.)
        for node,part,gap,angle in zip(ids[selected],parts[selected],gaps[selected]*1e3,np.degrees(angles)):
            nodes.append({"pose":label,"time_s":time[i],"node":int(node),"part":part,"gap_um":gap,"wall_sector_deg":angle})
        directed=float(axis.directed_tilt_deg.iloc[i]); azimuth=float(np.degrees(phase.phi_robot_rad.iloc[i]))
        row={"pose":label,"time_s":time[i],"directed_tilt_deg":directed,"local_azimuth_deg":azimuth,
             "head_gap_um":float((gaps[parts=="HEAD"]*1e3).min()),"tail_gap_um":float((gaps[parts=="TAIL"]*1e3).min()),
             "body_gap_um":float((gaps[parts=="BODY"]*1e3).min()),"wall_sector_A_deg":np.degrees(centers[0]) if centers else np.nan,
             "wall_sector_B_deg":np.degrees(centers[-1]) if len(centers)>1 else np.nan,"wall_sector_separation_deg":separation,
             "head_nodes_20um":int(np.sum(selected&(parts=="HEAD"))),"tail_nodes_20um":int(np.sum(selected&(parts=="TAIL"))),
             "body_nodes_20um":int(np.sum(selected&(parts=="BODY"))),"relative_to_first_flip":"after" if time[i]>.002830262 else "before"}
        rows.append(row)
    pd.DataFrame(rows).to_csv(HERE/"bridge_pose_geometry.csv",index=False)
    pd.DataFrame(nodes).to_csv(HERE/"bridge_pose_nearwall_nodes.csv",index=False)
    profile=pd.read_csv(HERE/"old_vs_freecad_head_profile.csv"); head=profile[profile.axial_x_mm<-.55]
    pd.DataFrame([{"head_max_outward_difference_um":head.new_minus_old_um.max(),
                   "head_max_difference_x_mm":head.loc[head.new_minus_old_um.idxmax(),"axial_x_mm"],
                   "head_positive_difference_start_x_mm":head.loc[head.new_minus_old_um>0,"axial_x_mm"].min(),
                   "head_positive_difference_end_x_mm":head.loc[head.new_minus_old_um>0,"axial_x_mm"].max()}]).to_csv(HERE/"head_profile_difference_summary.csv",index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__=="__main__":main()
