"""Lightweight RP/contact/proximity metrics for every solved screening case."""

from pathlib import Path
import json
import math
import sys

import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.spatial.transform import Rotation

HERE=Path(__file__).resolve().parent; SCREEN=HERE.parent; REPO=SCREEN.parents[1]; ROOT=REPO.parent
AUDIT=REPO/"calibration_analysis/ReducedHydro_hidden_impact_audit"
VALIDATION=REPO/"calibration_analysis/ReducedHydro_L2300_frozenCAD_validation"
sys.path[:0]=[str(VALIDATION),str(AUDIT)]
from audit_frozen_cad_static import Wall, exact_query

RP0=np.array([-7.468174204284,-3.676918015967,-9.550745259298])
A0=np.array([0.9647382600216,-0.1188742372140,0.2348382536499]); A0/=np.linalg.norm(A0)
CURVE=pd.read_csv(ROOT/"CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_TRUE_centerline_odb.csv")


def get_array(z,name):
    keys=[k for k in z.files if k==name or k.startswith(name+" (Repeated:")]
    return max((z[k] for k in keys),key=len)


def dense(folder):
    z=np.load(folder/"private/rp_history_private.npz"); base=get_array(z,"V1"); t=base[:,0]
    data={p:np.column_stack([np.interp(t,get_array(z,p+str(i))[:,0],get_array(z,p+str(i))[:,1]) for i in (1,2,3)]) for p in ("U","UR","V","VR")}
    return t,data


def contact(folder,t):
    z=np.load(folder/"private/contact_history_private.npz")
    def one(prefix):
        cols=[]
        for i in (1,2,3):
            keys=[k for k in z.files if "|"+prefix+str(i)+" on surface " in k and "ASSEMBLY_ROBOT" in k]
            a=max((z[k] for k in keys),key=len); cols.append(np.interp(t,a[:,0],a[:,1]))
        return np.column_stack(cols)
    return one("CFN")+one("CFS")


def unit(a): return a/np.maximum(np.linalg.norm(a,axis=-1,keepdims=True),1e-30)


def projection(points):
    xyz=CURVE[["x_mm","y_mm","z_mm"]].to_numpy(float); arc=CURVE.arclength_mm.to_numpy(float)
    edge=np.diff(xyz,axis=0); edge2=np.einsum("ij,ij->i",edge,edge); tangent=unit(edge)
    rows=[]; previous=None
    for point in points:
        ids=np.arange(len(edge)) if previous is None else np.arange(max(0,previous-12),min(len(edge),previous+13))
        origin=xyz[ids]; frac=np.clip(np.einsum("ij,ij->i",point-origin,edge[ids])/edge2[ids],0,1)
        q=origin+frac[:,None]*edge[ids]; d=np.einsum("ij,ij->i",point-q,point-q); best=int(np.argmin(d))
        if d[best]>1:
            ids=np.arange(len(edge)); origin=xyz[:-1]; frac=np.clip(np.einsum("ij,ij->i",point-origin,edge)/edge2,0,1)
            q=origin+frac[:,None]*edge; d=np.einsum("ij,ij->i",point-q,point-q); best=int(np.argmin(d))
        segment=int(ids[best]); previous=segment
        rows.append((arc[segment]+frac[best]*(arc[segment+1]-arc[segment]),*tangent[segment]))
    return np.asarray(rows)


def frames(tangent):
    e1=np.empty_like(tangent)
    for i,a in enumerate(tangent):
        v=np.array([0.,0.,1.])-a[2]*a
        if np.linalg.norm(v)<1e-8: v=np.array([0.,1.,0.])-a[1]*a
        v/=np.linalg.norm(v)
        if i and np.dot(v,e1[i-1])<0:v*=-1
        e1[i]=v
    return e1,unit(np.cross(tangent,e1))


def intervals(flag):
    edge=np.diff(np.r_[False,flag,False].astype(int)); return list(zip(np.where(edge==1)[0],np.where(edge==-1)[0]-1))


def merge_ranges(ranges,max_gap):
    merged=[]
    for a,b in ranges:
        if merged and a-merged[-1][1]-1<=max_gap: merged[-1]=(merged[-1][0],b)
        else: merged.append((a,b))
    return merged


def proximity(item,folder,t,U,rot,indices,tangent,e1,e2):
    code="%04d"%round(item["length_mm"]*1000); mesh=SCREEN/"cases/_meshes"/("Robot_SCREENING_L%s_D0815_nodes.csv"%code)
    props=json.loads((REPO/"cad/freecad_parametric_robot/screening/coarse_motion_mode"/("Robot_SCREENING_L%s_D0815_geometry.json"%code)).read_text())
    nodes=pd.read_csv(mesh); local=nodes[["x_mm","y_mm","z_mm"]].to_numpy(float); x=local[:,0]
    surface=(x<=.26149478)|(x>=item["length_mm"]-1e-7)|(np.abs(np.linalg.norm(local[:,1:],axis=1)-.4075)<2e-4)
    local=local[surface][::2]; x=x[surface][::2]
    region=np.where(x<=.26149478,"HEAD",np.where(x>=item["length_mm"]-1e-7,"TAIL","BODY"))
    cad_com=np.asarray(props["center_of_mass_mm"]); initial=RP0+(local-cad_com)@rotation_x_to_axis(A0).T
    wall=Wall(); rows=[]
    for j,i in enumerate(indices):
        points=RP0+U[i]+rot[i].apply(initial-RP0); gaps,tri,_=exact_query(wall,points)
        angles=np.degrees(np.arctan2(wall.normals[tri]@e2[i],wall.normals[tri]@e1[i]))%360
        nearest=int(np.argmin(gaps))
        row={"time_s":t[i],"increment":int(i),"min_gap_um":gaps.min()*1e3,
             "HEAD_gap_um":gaps[region=="HEAD"].min()*1e3,"TAIL_gap_um":gaps[region=="TAIL"].min()*1e3,
             "nearest_wall_sector_deg":angles[nearest]}
        rows.append(row)
    result=pd.DataFrame(rows); result.to_csv(folder/"proximity_light.csv",index=False); return result


def rotation_x_to_axis(axis):
    x=np.array([1.,0.,0.]); v=np.cross(x,axis); c=np.dot(x,axis); k=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0.]])
    return np.eye(3)+k+k@k/(1+c)


def analyze(folder):
    item=json.loads((folder/"case_identity.json").read_text()); t,data=dense(folder); rot=Rotation.from_rotvec(data["UR"])
    axis=rot.apply(np.broadcast_to(A0,data["U"].shape)); com=RP0+data["U"]
    p=projection(com)
    forward_sign=1.0 if np.dot(axis[0],p[0,1:4])>=0 else -1.0
    s=forward_sign*p[:,0]; tangent=forward_sign*p[:,1:4]; e1,e2=frames(tangent)
    dot=np.einsum("ij,ij->i",axis,tangent); tilt=np.degrees(np.arccos(np.clip(dot,-1,1)))
    radial=axis-dot[:,None]*tangent
    phase=np.unwrap(np.arctan2(np.einsum("ij,ij->i",radial,e2),np.einsum("ij,ij->i",radial,e1)))
    vt=np.einsum("ij,ij->i",data["V"],tangent); force=contact(folder,t); active=np.linalg.norm(force,axis=1)>1e-8
    dt=float(np.median(np.diff(t))); event_ranges=merge_ranges(intervals(active),int(round(5e-6/dt))); states=np.full(len(t),"NONE",object)
    impacts=slides=stuck=0
    for a,b in event_ranges:
        duration=(b-a+1)*dt; separation=(b+round(20e-6/dt)<len(t) and not active[b+1:b+1+round(20e-6/dt)].any())
        if duration<1e-4 and separation: state="IMPACT"; impacts+=b-a+1
        elif np.mean(np.abs(vt[a:b+1]))>=50: state="SLIDING"; slides+=b-a+1
        else: state="STUCK"; stuck+=b-a+1
        states[a:b+1]=state
    sample=np.unique(np.r_[np.arange(0,len(t),100),len(t)-1]); prox_idx=np.unique(np.r_[np.arange(0,len(t),200),len(t)-1])
    prox=proximity(item,folder,t,data["U"],rot,prox_idx,tangent,e1,e2)
    support_head=prox.HEAD_gap_um<=20; support_tail=prox.TAIL_gap_um<=20
    sectors=np.unwrap(np.radians(prox.nearest_wall_sector_deg.to_numpy()))
    sector_transitions=int(np.sum(np.abs(np.diff(sectors))>math.radians(30)))
    props_code="%04d"%round(item["length_mm"]*1000)
    props=json.loads((REPO/"cad/freecad_parametric_robot/screening/coarse_motion_mode"/("Robot_SCREENING_L%s_D0815_geometry.json"%props_code)).read_text())
    ip=np.asarray(props["principal_inertia_tonne_mm2"]); ilocal=np.diag([ip.min(),ip.max(),ip.max()]); r0=rotation_x_to_axis(A0)
    initial_I=r0@ilocal@r0.T; rm=rot.as_matrix(); inertia=np.einsum("nij,jk,nlk->nil",rm,initial_I,rm)
    mass=item["mesh_mass_mg"]*1e-9; kt=.5*mass*np.einsum("ij,ij->i",data["V"],data["V"])*1e-3
    kr=.5*np.einsum("ni,nij,nj->n",data["VR"],inertia,data["VR"])*1e-3
    extrema=len(find_peaks(tilt,prominence=1.0)[0])+len(find_peaks(-tilt,prominence=1.0)[0])
    advance=float(np.degrees(phase[-1]-phase[0])); delta=float(s[-1]-s[0]); n=len(t)
    metrics=dict(item)
    metrics.update(status="VALID",max_directed_tilt_deg=float(tilt.max()),final_directed_tilt_deg=float(tilt[-1]),
        crossing_90deg=bool(np.any(dot<0)),polarity_reversal=bool(np.sign(dot[-1])!=np.sign(dot[0])),
        robot_phase_advance_deg=advance,phase_rate_sign=int(np.sign(advance)),delta_s_mm=delta,
        mean_abs_Vt_mm_s=float(np.mean(np.abs(vt))),contact_count=len(event_ranges),
        longest_contact_us=max(((b-a+1)*dt*1e6 for a,b in event_ranges),default=0.0),
        head_support_fraction=float(support_head.mean()),tail_support_fraction=float(support_tail.mean()),
        both_end_support_fraction=float((support_head&support_tail).mean()),wall_sector_span_deg=float(np.degrees(np.ptp(sectors))),
        wall_sector_transitions=sector_transitions,tilt_local_extrema_count=extrema,
        KErot_RMS_J=float(np.sqrt(np.mean(kr**2))),KEtrans_RMS_J=float(np.sqrt(np.mean(kt**2))),
        KErot_KEtrans_ratio=float(np.sqrt(np.mean(kr**2))/max(np.sqrt(np.mean(kt**2)),1e-30)),
        rotation_translation_index=float(abs(advance)/max(abs(delta),.05)),impact_fraction=impacts/n,
        sliding_fraction=slides/n,stuck_fraction=stuck/n)
    center=CURVE[["x_mm","y_mm","z_mm"]].to_numpy(float)[np.clip(np.searchsorted(CURVE.arclength_mm.to_numpy(float),forward_sign*s),0,len(CURVE)-1)]
    com_x=float(np.asarray(props["center_of_mass_mm"])[0]); head=com-com_x*axis; tail=com+(item["length_mm"]-com_x)*axis
    telemetry=pd.read_csv(next(folder.glob("*_telemetry.csv")))
    b=np.column_stack([np.interp(t,telemetry.t_s,telemetry[c]) for c in ("Bx_aba_T","By_aba_T","Bz_aba_T")])
    light=pd.DataFrame({"time_s":t[sample],"cycle_fraction":t[sample]*item["frequency_Hz"],
        "rp_x_mm":com[sample,0],"rp_y_mm":com[sample,1],"rp_z_mm":com[sample,2],
        "center_x_mm":center[sample,0],"center_y_mm":center[sample,1],"center_z_mm":center[sample,2],
        "head_x_mm":head[sample,0],"head_y_mm":head[sample,1],"head_z_mm":head[sample,2],
        "tail_x_mm":tail[sample,0],"tail_y_mm":tail[sample,1],"tail_z_mm":tail[sample,2],
        "B_x_T":b[sample,0],"B_y_T":b[sample,1],"B_z_T":b[sample,2],
        "axis_x":axis[sample,0],"axis_y":axis[sample,1],"axis_z":axis[sample,2],
        "tangent_x":tangent[sample,0],"tangent_y":tangent[sample,1],"tangent_z":tangent[sample,2],
        "e1_x":e1[sample,0],"e1_y":e1[sample,1],"e1_z":e1[sample,2],"e2_x":e2[sample,0],"e2_y":e2[sample,1],"e2_z":e2[sample,2],
        "directed_tilt_deg":tilt[sample],"robot_phase_deg":np.degrees(phase[sample]-phase[0]),
        "delta_s_mm":s[sample]-s[0],"Vt_mm_s":vt[sample],"contact_state":states[sample]})
    light.to_csv(folder/"pose_light.csv",index=False)
    (SCREEN/"metrics"/(item["case_id"]+"_metrics.json")).write_text(json.dumps(metrics,indent=2)+"\n")
    return metrics


def baseline_comparison(metrics):
    old=REPO/"calibration_analysis/ReducedHydro_L2300_frozenCAD_validation"
    tilt=pd.read_csv(old/"L2300_directed_tilt.csv"); phase=pd.read_csv(old/"L2300_phase.csv"); trans=pd.read_csv(old/"L2300_translation.csv")
    end=0.2/30; ti=tilt.time_s<=end; pi=phase.time_s<=end; si=trans.time_s<=end
    row={"quantity":["max_directed_tilt_deg","robot_phase_advance_deg","delta_s_mm"],
         "fine_0p060":[tilt.loc[ti,"directed_tilt_deg"].max(),phase.loc[pi,"robot_phase_advance_deg"].iloc[-1],trans.loc[si,"delta_s_mm"].iloc[-1]],
         "coarse_0p100":[metrics["max_directed_tilt_deg"],metrics["robot_phase_advance_deg"],metrics["delta_s_mm"]]}
    out=pd.DataFrame(row); out["same_sign_or_envelope"]=[True,np.sign(out.fine_0p060.iloc[1])==np.sign(out.coarse_0p100.iloc[1]),np.sign(out.fine_0p060.iloc[2])==np.sign(out.coarse_0p100.iloc[2])]
    out.to_csv(SCREEN/"metrics/coarse_baseline_comparison.csv",index=False)


def main():
    wanted=set(sys.argv[1:]); results=[]
    for identity in sorted((SCREEN/"cases").glob("*/case_identity.json")):
        folder=identity.parent; item=json.loads(identity.read_text())
        if wanted and item["case_id"] not in wanted: continue
        if item.get("status")!="SOLVED": continue
        print("ANALYZE",item["case_id"],flush=True); results.append(analyze(folder))
    existing=SCREEN/"metrics/all_cases_metrics.csv"
    old=pd.read_csv(existing) if existing.exists() else pd.DataFrame()
    combined=pd.concat([old,pd.DataFrame(results)],ignore_index=True) if len(old) else pd.DataFrame(results)
    if len(combined):
        combined=combined.drop_duplicates("case_id",keep="last")
        order={case_id:i for i,case_id in enumerate(pd.read_csv(SCREEN/"metrics/case_manifest.csv").case_id)}
        combined["_manifest_order"]=combined.case_id.map(order); combined=combined.sort_values("_manifest_order").drop(columns="_manifest_order")
        combined.to_csv(existing,index=False)
    if any(r["case_id"]=="GEO_230" for r in results): baseline_comparison(next(r for r in results if r["case_id"]=="GEO_230"))
    print(pd.DataFrame(results)[["case_id","max_directed_tilt_deg","robot_phase_advance_deg","delta_s_mm","contact_count","longest_contact_us"]].to_string(index=False))


if __name__=="__main__":main()
