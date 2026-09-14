"""Four-way exact-wall replay of OLD/NEW geometry and OLD/NEW trajectories.

No Abaqus or ODB access occurs. A conservative triangle-centroid bound removes
points that cannot be within 20 um before the verified exact triangle query.
"""
from collections import Counter
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, cKDTree
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "ReducedHydro_geometry_L1800_validation"
AUDIT = HERE.parent / "ReducedHydro_hidden_impact_audit"
sys.path[:0] = [str(AUDIT), str(AUDIT / "normal_contact_damping_probe")]
from audit_stage_a import mesh_properties
from exact_gap_audit import Wall, closest
from contact_probe_common import dense_rp

NEW_JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_WallOn_Free_0083"
OLD_JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083"
THRESHOLDS = (20.0, 10.0, 5.0, 0.0)
DT = 1e-6


def exterior_ids(mesh):
    counts = Counter()
    for _, a, b, c, d in mesh["elems"]:
        for face in ((a,b,c),(a,b,d),(a,c,d),(b,c,d)):
            counts[tuple(sorted(face))] += 1
    return np.unique([node for face,count in counts.items() if count == 1 for node in face])


def unit(v):
    return v / np.maximum(np.linalg.norm(v,axis=-1,keepdims=True),1e-30)


def clusters(angles):
    if not len(angles): return []
    values=np.sort(np.mod(angles,2*np.pi)); gaps=np.diff(np.r_[values,values[0]+2*np.pi])
    cut=int(np.argmax(gaps)); ordered=np.r_[values[cut+1:],values[:cut+1]+2*np.pi]
    groups=np.split(ordered,np.where(np.diff(ordered)>np.radians(45))[0]+1)
    return [float(np.mod(np.angle(np.mean(np.exp(1j*g))),2*np.pi)) for g in groups if len(g)]


def intervals(times, flags):
    rows=[]; start=None; previous=None
    for time,flag in zip(times,np.asarray(flags,bool)):
        if flag and start is None: start=time
        if start is not None and (not flag or (previous is not None and time-previous>1.5*DT)):
            rows.append((start,previous,previous-start+DT)); start=time if flag else None
        previous=time
    if start is not None: rows.append((start,previous,previous-start+DT))
    return rows


def exact_query(wall, points):
    """Same certified centroid-radius bound as Wall.query, with k=8 start."""
    points=np.asarray(points,float); total=len(wall.tri); signed=np.empty(len(points))
    triangle=np.empty(len(points),int); nearest=np.empty((len(points),3)); remaining=np.arange(len(points)); k=min(8,total)
    while len(remaining):
        dc,indices=wall.tree.query(points[remaining],k=k)
        if k==1: dc,indices=dc[:,None],indices[:,None]
        tri=wall.tri[indices]; q,dist=closest(points[remaining,None,:],tri[:,:,0],tri[:,:,1],tri[:,:,2])
        row=np.arange(len(remaining)); best=np.argmin(dist,axis=1); dd=dist[row,best]; jj=indices[row,best]; qb=q[row,best]
        resolved=np.ones(len(remaining),bool) if k==total else dc[:,-1]-wall.radius>dd+1e-10
        out=remaining[resolved]; triangle[out]=jj[resolved]; nearest[out]=qb[resolved]
        signed[out]=dd[resolved]*np.sign(np.sum((points[out]-qb[resolved])*wall.normals[jj[resolved]],axis=1))
        remaining=remaining[~resolved]; k=min(k*2,total)
    return signed,triangle,nearest


def geometry(name, mesh, ids, axis_seed):
    points=np.asarray([mesh["nodes"][int(i)] for i in ids])+mesh["shift"]
    if name == "NEW":
        hull=ConvexHull(points).vertices; ids=np.asarray(ids)[hull]; points=points[hull]
    axial=(points-mesh["com"])@unit(np.asarray(axis_seed))
    lo,hi=np.quantile(axial,[.10,.90])
    labels=np.where(axial<=lo,"HEAD",np.where(axial>=hi,"TAIL","BODY"))
    return {"name":name,"mesh":mesh,"ids":np.asarray(ids,int),"points":points,
            "labels":labels,"com":np.asarray(mesh["com"])}


def trajectory(name, folder, mesh, rp, axis_file):
    time,data=dense_rp(folder/"candidate_private"); idx=np.arange(0,len(time),10,dtype=int)
    if idx[-1]!=len(time)-1: idx=np.r_[idx,len(time)-1]
    axis=pd.read_csv(axis_file).iloc[idx].reset_index(drop=True)
    return {"name":name,"time":time[idx],"U":data["U"][idx],"UR":data["UR"][idx],
            "mesh":mesh,"rp":np.asarray(rp),"axis":axis}


def replay(geom,traj,wall,center_tree):
    rotations=Rotation.from_rotvec(traj["UR"]).as_matrix()
    traj_com=traj["rp"]+traj["U"]+np.einsum("bij,j->bi",rotations,traj["mesh"]["com"]-traj["rp"])
    local=geom["points"]-geom["com"]
    e1_all=traj["axis"][["e1_x","e1_y","e1_z"]].to_numpy(float)
    e2_all=traj["axis"][["e2_x","e2_y","e2_z"]].to_numpy(float)
    timeline=[]
    for offset in range(0,len(traj["time"]),250):
        stop=min(len(traj["time"]),offset+250); rot=rotations[offset:stop]
        positions=traj_com[offset:stop,None,:]+np.einsum("bij,nj->bni",rot,local)
        flat=positions.reshape(-1,3)
        centroid_distance=wall.tree.query(flat,k=1)[0].reshape(stop-offset,len(local))
        centerline_distance=center_tree.query(flat,k=1)[0].reshape(stop-offset,len(local))
        # TRUE centerline samples are at most 0.270 mm apart; 0.70 mm is a
        # conservative inner cutoff for the nominal 1.005 mm lumen wall.
        possible=(centroid_distance <= wall.radius + .020001) & (centerline_distance >= .70)
        rows,cols=np.where(possible)
        exact_gap=np.full(possible.shape,.020001); exact_tri=np.full(possible.shape,-1,int)
        if len(rows):
            for chunk in range(0,len(rows),5000):
                part=slice(chunk,chunk+5000); points=positions[rows[part],cols[part]]
                gap,tri,_=exact_query(wall,points); exact_gap[rows[part],cols[part]]=gap; exact_tri[rows[part],cols[part]]=tri
        for local_i in range(stop-offset):
            global_i=offset+local_i; gaps=exact_gap[local_i]; tris=exact_tri[local_i]
            valid=tris>=0; minimum=float(gaps[valid].min()*1e3) if valid.any() else 20.001
            row={"geometry":geom["name"],"trajectory":traj["name"],"time_s":traj["time"][global_i],
                 "minimum_gap_um":minimum}
            e1=e1_all[global_i]; e2=e2_all[global_i]
            for threshold in THRESHOLDS:
                selected=valid & (gaps*1e3<=threshold); labels=geom["labels"][selected]
                angles=np.arctan2(wall.normals[tris[selected]]@e2,wall.normals[tris[selected]]@e1) if selected.any() else np.array([])
                centers=clusters(angles)
                separation=max((np.degrees(np.arccos(np.clip(np.cos(a-b),-1,1))) for q,a in enumerate(centers) for b in centers[q+1:]),default=0.0)
                tag=str(int(threshold)); row.update({
                    f"bridge_{tag}um":int(len(centers)>=2 and separation>=120),
                    f"head_count_{tag}um":int(np.sum(labels=="HEAD")),
                    f"tail_count_{tag}um":int(np.sum(labels=="TAIL")),
                    f"body_count_{tag}um":int(np.sum(labels=="BODY")),
                    f"sector_separation_{tag}um_deg":separation})
            timeline.append(row)
        print(f"{geom['name']} geometry + {traj['name']} trajectory: {stop}/{len(traj['time'])}",flush=True)
    return pd.DataFrame(timeline)


def profile(old_geom,new_geom):
    lo=max(((g["points"]-g["com"])@np.array([.9647382600216,-.118874237214,.2348382536499])).min() for g in (old_geom,new_geom))
    hi=min(((g["points"]-g["com"])@np.array([.9647382600216,-.118874237214,.2348382536499])).max() for g in (old_geom,new_geom))
    edges=np.linspace(lo,hi,181); centers=(edges[:-1]+edges[1:])/2; out={"axial_x_mm":centers}
    for key,g in (("old",old_geom),("new",new_geom)):
        q=g["points"]-g["com"]; axis=np.array([.9647382600216,-.118874237214,.2348382536499]); x=q@axis
        radius=np.linalg.norm(q-np.outer(x,axis),axis=1); values=np.full(len(centers),np.nan)
        bins=np.clip(np.digitize(x,edges)-1,0,len(centers)-1)
        for i in range(len(centers)):
            if np.any(bins==i): values[i]=radius[bins==i].max()
        values=pd.Series(values).interpolate(limit_direction="both").to_numpy(); out[f"r_{key}_mm"]=values
    frame=pd.DataFrame(out); frame["new_minus_old_um"]=(frame.r_new_mm-frame.r_old_mm)*1e3
    return frame


def summarize(timeline):
    rows=[]
    for (geometry_name,trajectory_name),group in timeline.groupby(["geometry","trajectory"],sort=False):
        row={"geometry":geometry_name,"trajectory":trajectory_name,"penetration_um":group.minimum_gap_um.min()}
        longest20=None
        for threshold in THRESHOLDS:
            tag=str(int(threshold)); events=intervals(group.time_s,group[f"bridge_{tag}um"])
            longest=max(events,key=lambda item:item[2]) if events else (np.nan,np.nan,0.0)
            row[f"longest{tag}_ms"]=longest[2]*1e3
            if threshold==20:
                longest20=longest; row["first_bridge_s"]=events[0][0] if events else np.nan
        if longest20[2]>0:
            segment=group[(group.time_s>=longest20[0])&(group.time_s<=longest20[1])&group.bridge_20um.astype(bool)]
            row.update({"HEAD":int(segment.head_count_20um.max()),"TAIL":int(segment.tail_count_20um.max()),
                        "BODY":int(segment.body_count_20um.max()),"longest20_start_s":longest20[0],"longest20_end_s":longest20[1]})
        else: row.update({"HEAD":0,"TAIL":0,"BODY":0,"longest20_start_s":np.nan,"longest20_end_s":np.nan})
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    new_meshes,new_rp,_=mesh_properties((HERE/f"{NEW_JOB}.inp").read_text()); new_mesh=new_meshes["Robot_SOLID"]
    old_meshes,old_rp,_=mesh_properties((OLD/f"{OLD_JOB}.inp").read_text()); old_mesh=old_meshes["Robot_SOLID"]
    seed=np.array([.9647382600216,-.118874237214,.2348382536499])
    new_ids=np.unique(pd.read_csv(HERE/"freecad_robot_surface_triangles_exact.csv")[["n1","n2","n3"]].to_numpy())
    old_geom=geometry("OLD",old_mesh,exterior_ids(old_mesh),seed); new_geom=geometry("NEW",new_mesh,new_ids,seed)
    old_traj=trajectory("OLD",OLD,old_mesh,old_rp,OLD/"L1800_8p333_true_axis.csv")
    new_traj=trajectory("NEW",HERE,new_mesh,new_rp,HERE/"freecad_8p333_true_axis.csv")
    wall=Wall()
    root=HERE.parents[2]
    centerline=pd.read_csv(root/"CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_TRUE_centerline_odb.csv")[["x_mm","y_mm","z_mm"]].to_numpy(float)
    center_tree=cKDTree(centerline); tables=[]
    for geom,traj in ((old_geom,old_traj),(new_geom,old_traj),(old_geom,new_traj),(new_geom,new_traj)):
        case_path=HERE/f"cross_replay_case_{geom['name'].lower()}_{traj['name'].lower()}.csv"
        if case_path.exists() and len(pd.read_csv(case_path,usecols=["time_s"]))==len(traj["time"]):
            table=pd.read_csv(case_path); print(f"reusing complete {case_path.name}")
        else:
            table=replay(geom,traj,wall,center_tree); table.to_csv(case_path,index=False)
        tables.append(table)
    timeline=pd.concat(tables,ignore_index=True); timeline.to_csv(HERE/"cross_replay_bridge_timeline.csv",index=False)
    summary=summarize(timeline); summary.to_csv(HERE/"cross_replay_bridge_summary.csv",index=False)
    head=profile(old_geom,new_geom); head.to_csv(HERE/"old_vs_freecad_head_profile.csv",index=False)
    a,b,c,d=[summary[(summary.geometry==g)&(summary.trajectory==t)].iloc[0] for g,t in (("OLD","OLD"),("NEW","OLD"),("OLD","NEW"),("NEW","NEW"))]
    if b.longest20_ms>.5 and c.longest20_ms<.5:
        causal="EXACT_HEAD_GEOMETRY_DOMINATES_BRIDGE_REINTRODUCTION"
    elif b.longest20_ms<.5 and c.longest20_ms>.5:
        causal="ROTATIONAL_TRAJECTORY_DOMINATES_BRIDGE_REINTRODUCTION"
    else:
        causal="GEOMETRY_TRAJECTORY_COUPLING_DRIVES_BRIDGE"
    identity={"wall":"SmoothWall114 exact triangles","pose_map":"COM(t)+R(t)*(x_ref-COM_ref)",
              "sample_interval_s":DT,"thresholds_um":list(THRESHOLDS),"centroid_bound":"d_centroid <= triangle_radius + 20.001 um",
              "centerline_prefilter":"nearest TRUE-centerline sample distance >= 0.70 mm; sample spacing <= 0.270 mm",
              "geometry_nodes":{"OLD":len(old_geom["ids"]),"NEW":len(new_geom["ids"])},
              "geometry_representation":"actual exterior mesh nodes; convex-hull extremal subset for convex NEW body",
              "old_old_regression_target_ms":.2588,"new_new_regression_target_ms":.547,"causal_classification":causal,
              "first_new_trajectory_flip_ms":2.830262}
    (HERE/"cross_replay_identity.json").write_text(json.dumps(identity,indent=2)+"\n")
    print(summary.to_string(index=False)); print(json.dumps(identity,indent=2))


if __name__=="__main__":
    main()
