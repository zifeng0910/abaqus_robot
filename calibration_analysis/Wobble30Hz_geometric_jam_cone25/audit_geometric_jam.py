"""Zero-cost geometry audit for the 30 Hz wall-on wobble case.

This is deliberately an audit only: it does not submit or modify an Abaqus job.
All distances are mm and use the validated Abaqus-frame SmoothWall114 triangles.
"""
from pathlib import Path
import json, math
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

ROOT = Path(r"J:\\abaqusfangzhen")
JOB = "Wobble_F30_G6L45_WallOn_Free_0083_WobbleSurvival"
OUT = ROOT / "abaqus_robot" / "calibration_analysis" / "Wobble30Hz_geometric_jam_cone25"
OUT.mkdir(parents=True, exist_ok=True)
frames_path = ROOT / f"{JOB}_bend_validation" / "robot_frames.csv"
axis_path = ROOT / f"{JOB}_true_axis_raw.csv"
tele_path = ROOT / f"{JOB}_telemetry.csv"
wall_nodes_path = ROOT / f"{JOB}_pipe_wall_nodes_exact.csv"
wall_el_path = ROOT / f"{JOB}_pipe_wall_triangles_exact.csv"

def write(df, name): df.to_csv(OUT/name, index=False)

fr = pd.read_csv(frames_path)
f0 = fr[fr.frame == fr.frame.min()].sort_values('node')
xyz0 = f0[['x_mm','y_mm','z_mm']].to_numpy(float)
rp0 = xyz0.mean(axis=0)
try:
    ax0 = pd.read_csv(axis_path).iloc[0]
    head = ax0[['head_x','head_y','head_z']].to_numpy(float)
    tail = ax0[['tail_x','tail_y','tail_z']].to_numpy(float)
except Exception:
    head, tail = xyz0[np.argmax(xyz0[:,0])], xyz0[np.argmin(xyz0[:,0])]
a0 = head-tail; a0 /= np.linalg.norm(a0)
# Use the directed body axis that follows increasing pipe arc (TAIL→HEAD).
# The mesh is geometrically almost symmetric, but this removes the 180°
# ambiguity when reporting cone tilt relative to the forward tangent.
a0 = -a0
q0 = xyz0-rp0
zq = q0 @ a0
rad = np.linalg.norm(q0-np.outer(zq,a0), axis=1)
L = float(zq.max()-zq.min())
diam_max, diam_p95, diam_med = 2*rad.max(), 2*np.percentile(rad,95), 2*np.median(rad)
head_q = q0[zq >= zq.max()-0.15*L]; tail_q = q0[zq <= zq.min()+0.15*L]
head_d = 2*np.max(np.linalg.norm(head_q-(head_q@a0)[:,None]*a0,axis=1))
tail_d = 2*np.max(np.linalg.norm(tail_q-(tail_q@a0)[:,None]*a0,axis=1))
geom = pd.DataFrame([dict(n_exterior_nodes=len(xyz0), axial_length_mm=L, max_diameter_mm=diam_max,
    p95_diameter_mm=diam_p95, median_diameter_mm=diam_med, head_diameter_mm=head_d,
    tail_diameter_mm=tail_d, rp_x_mm=rp0[0], rp_y_mm=rp0[1], rp_z_mm=rp0[2],
    head_x_mm=head[0], head_y_mm=head[1], head_z_mm=head[2], tail_x_mm=tail[0], tail_y_mm=tail[1], tail_z_mm=tail[2])])
write(geom,'actual_robot_geometry_dimensions.csv')

wn = pd.read_csv(wall_nodes_path).set_index('node')
we = pd.read_csv(wall_el_path)
tris = []
for _, r in we.iterrows():
    ids=[int(r.n1),int(r.n2),int(r.n3),int(r.n4)]
    tris += [(wn.loc[ids[0]].to_numpy(float),wn.loc[ids[1]].to_numpy(float),wn.loc[ids[2]].to_numpy(float)),
             (wn.loc[ids[0]].to_numpy(float),wn.loc[ids[2]].to_numpy(float),wn.loc[ids[3]].to_numpy(float))]
tri = np.asarray(tris,float); cent=tri.mean(axis=1)
ctree=cKDTree(cent)

tf = pd.read_csv(ROOT/'curvenew_CEL_xyrot56_exact.csv')
T = np.asarray(json.loads((ROOT/'abaqus_magpylib_frame_transform.json').read_text())['R_aba_to_mag'],float)
origin=np.asarray(json.loads((ROOT/'abaqus_magpylib_frame_transform.json').read_text())['origin_aba_mm'],float)
cv = tf[['x_mm','y_mm','z_mm']].to_numpy(float) @ T + origin
cv_tree=cKDTree(cv)
arc=tf.arclength_mm.to_numpy(float)
idx=int(np.argmin(np.linalg.norm(cv-rp0,axis=1)))
lo=max(0,idx-2); hi=min(len(cv)-1,idx+2)
tloc=cv[hi]-cv[lo]; tloc/=np.linalg.norm(tloc)
ref=np.array([0.,0.,1.]);
if abs(np.dot(ref,tloc))>0.9: ref=np.array([0.,1.,0.])
nloc=ref-np.dot(ref,tloc)*tloc; nloc/=np.linalg.norm(nloc); bloc=np.cross(tloc,nloc); bloc/=np.linalg.norm(bloc)
cpt=cv[idx]
# local wall section, with an explicitly stated slab so the estimate is reproducible
wp=wn[['x_mm','y_mm','z_mm']].to_numpy(float)
axproj=(wp-cpt)@tloc
sel=np.abs(axproj)<=0.45
uv=(wp[sel]-cpt)@np.column_stack([nloc,bloc])
local=pd.DataFrame([dict(station=idx, station_arc_mm=arc[idx], center_x_mm=cpt[0],center_y_mm=cpt[1],center_z_mm=cpt[2],
    tangent_x=tloc[0],tangent_y=tloc[1],tangent_z=tloc[2], n_x=nloc[0],n_y=nloc[1],n_z=nloc[2], b_x=bloc[0],b_y=bloc[1],b_z=bloc[2],
    slab_halfwidth_mm=.45, n_span_mm=np.ptp(uv[:,0]), b_span_mm=np.ptp(uv[:,1]),
    radial_median_mm=np.median(np.linalg.norm(uv,axis=1)), radial_min_mm=np.min(np.linalg.norm(uv,axis=1)),
    wall_nodes_in_slab=int(sel.sum()))])
write(local,'local_tube_geometry_dimensions.csv')

def rot_axis(axis):
    """Rotation taking a0 to axis."""
    v=np.cross(a0,axis); s=np.linalg.norm(v); c=np.dot(a0,axis)
    if s<1e-12: return np.eye(3) if c>0 else np.diag([1,-1,-1])
    K=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
    return np.eye(3)+K+K@K*((1-c)/(s*s))

def orient_normals():
    v1=tri[:,1]-tri[:,0]; v2=tri[:,2]-tri[:,0]
    nn=np.cross(v1,v2); nn/=np.maximum(np.linalg.norm(nn,axis=1)[:,None],1e-15)
    # orient towards the nearest centerline point (tube interior)
    nearest=cv[cv_tree.query(cent,k=1)[1]]
    flip=np.sum(nn*(nearest-cent),axis=1)<0
    nn[flip]*=-1
    return nn
tnorm=orient_normals()

def signed_min(points, k=24):
    ds, ids=ctree.query(points,k=min(k,len(tri)))
    if ds.ndim==1: ds=ds[:,None]; ids=ids[:,None]
    best=np.full(len(points),np.inf); bestsign=np.ones(len(points)); besttri=np.full(len(points),-1,int)
    for j in range(ds.shape[1]):
        tt=tri[ids[:,j]]; aa=tt[:,0]; bb=tt[:,1]; cc=tt[:,2]
        ab=bb-aa; ac=cc-aa; ap=points-aa
        d1=np.sum(ab*ap,axis=1); d2=np.sum(ac*ap,axis=1)
        # robust closest point, using barycentric regions
        out=aa.copy(); m=(d1<=0)&(d2<=0)
        bp=points-bb; d3=np.sum(ab*bp,axis=1); d4=np.sum(ac*bp,axis=1); mb=(d3>=0)&(d4<=d3)
        out[mb]=bb[mb]
        vc=d1*d4-d3*d2; mc=(d1>=0)&(d3<=0)&(vc<=0); v=d1/(d1-d3+1e-30); out[mc]=aa[mc]+v[mc,None]*ab[mc]
        cp=points-cc; d5=np.sum(ab*cp,axis=1); d6=np.sum(ac*cp,axis=1); mc2=(d6>=0)&(d5<=d6); out[mc2]=cc[mc2]
        vb=d5*d2-d1*d6; mb2=(d2>=0)&(d6<=0)&(vb<=0); w=d2/(d2-d6+1e-30); out[mb2]=aa[mb2]+w[mb2,None]*ac[mb2]
        va=d3*d6-d5*d4; mface=(va<=0)&((d4-d3)>=0)&((d5-d6)>=0); w2=(d4-d3)/(d4-d3+d5-d6+1e-30); out[mface]=bb[mface]+w2[mface,None]*(cc[mface]-bb[mface])
        rem=~(m|mb|mc|mc2|mb2|mface)
        den=va[rem]+vb[rem]+vc[rem]+1e-30; vv=vb[rem]/den; ww=vc[rem]/den; out[rem]=aa[rem]+vv[:,None]*ab[rem]+ww[:,None]*ac[rem]
        ud=np.linalg.norm(points-out,axis=1)
        sign=np.sign(np.sum((points-out)*tnorm[ids[:,j]],axis=1)); sign[sign==0]=1
        upd=ud<best; best[upd]=ud[upd]; bestsign[upd]=sign[upd]; besttri[upd]=ids[upd,j]
    return best*bestsign, besttri

# feasibility map: offsets allow COM to seek the best clearance, while theta/chi are orientation variables
thetas=np.arange(0,30.01,2.0); chis=np.arange(0,360,10.0); offs=np.linspace(-.24,.24,5)
rows=[]
for th in thetas:
  for ch in chis:
    ar=math.radians(th); cr=math.radians(ch)
    axis=math.cos(ar)*tloc+math.sin(ar)*(math.cos(cr)*nloc+math.sin(cr)*bloc)
    R=rot_axis(axis); best=(-1e9,0,0)
    for u in offs:
      for v in offs:
        pts=rp0 + u*nloc + v*bloc + (q0@R.T)
        g,_=signed_min(pts)
        gm=float(g.min())
        if gm>best[0]: best=(gm,u,v)
    rows.append(dict(theta_deg=th,chi_deg=ch,g_best_mm=best[0],u_best_mm=best[1],v_best_mm=best[2],
                     pass_geometric=best[0]>=0,pass_5um=best[0]>=.005,pass_20um=best[0]>=.020))
feas=pd.DataFrame(rows); write(feas,'orientation_feasibility_map.csv')
tm=[]
for ch,g in feas.groupby('chi_deg'):
    tm.append(dict(chi_deg=ch,theta_max_geometric_deg=g.loc[g.pass_geometric,'theta_deg'].max() if g.pass_geometric.any() else np.nan,
                   theta_max_5um_deg=g.loc[g.pass_5um,'theta_deg'].max() if g.pass_5um.any() else np.nan,
                   theta_max_20um_deg=g.loc[g.pass_20um,'theta_deg'].max() if g.pass_20um.any() else np.nan))
write(pd.DataFrame(tm),'theta_max_vs_azimuth.csv')

# dynamic overlay and exact multi-contact clusters
ax=pd.read_csv(axis_path); tel=pd.read_csv(tele_path)
mp=[]; cl=[]
for fi,g in fr.groupby('frame',sort=True):
    row=ax.iloc[min(int(fi),len(ax)-1)]; p=np.array([row.rp_x,row.rp_y,row.rp_z],float) if 'rp_x' in ax else np.array([row.rp_x_aba_mm,row.rp_y_aba_mm,row.rp_z_aba_mm],float)
    aa=-np.array([row.axis_x,row.axis_y,row.axis_z],float); aa/=max(np.linalg.norm(aa),1e-15)
    ci=int(np.argmin(np.linalg.norm(cv-p,axis=1))); l=max(0,ci-2); h=min(len(cv)-1,ci+2); tt=cv[h]-cv[l]; tt/=np.linalg.norm(tt)
    nn=nloc-np.dot(nloc,tt)*tt; nn/=max(np.linalg.norm(nn),1e-15); bb=np.cross(tt,nn); bb/=max(np.linalg.norm(bb),1e-15)
    th=math.degrees(math.acos(np.clip(np.dot(aa,tt),-1,1))); ch=(math.degrees(math.atan2(np.dot(aa,bb),np.dot(aa,nn)))%360)
    pts=g[['x_mm','y_mm','z_mm']].to_numpy(float); sg, ids=signed_min(pts); near=np.abs(sg)<.020
    normals=tnorm[ids[near]] if np.any(near) else np.empty((0,3)); side=np.sign(normals@nn) if len(normals) else np.array([])
    opposing=bool(np.any(side>0)&np.any(side<0)); pen=int(np.sum(sg<0)); near5=int(np.sum(sg<.005)); near20=int(np.sum(sg<.020))
    mp.append(dict(frame=int(fi),time_s=float(g.time_s.iloc[0]),theta_robot_deg=th,chi_robot_deg=ch,station=ci,actual_min_gap_mm=float(sg.min()),near5_nodes=near5,near20_nodes=near20,penetrating_nodes=pen,opposing_wall_clusters=opposing))
    cl.append(dict(frame=int(fi),time_s=float(g.time_s.iloc[0]),positive_normal_nodes=int(np.sum(side>0)) if len(side) else 0,negative_normal_nodes=int(np.sum(side<0)) if len(side) else 0,opposing_wall_clusters=opposing))
dyn=pd.DataFrame(mp); write(dyn,'dynamic_orientation_vs_feasibility.csv'); write(pd.DataFrame(cl),'multi_contact_cluster_history.csv')
first_bridge=dyn.loc[dyn.opposing_wall_clusters,'time_s'].min() if dyn.opposing_wall_clusters.any() else np.nan
persist=(dyn.actual_min_gap_mm<.005).rolling(10,min_periods=1).sum()>=8
classification='PERSISTENT_GEOMETRIC_JAM' if dyn.opposing_wall_clusters.any() and persist.any() else ('GEOMETRIC_WEDGE_LIMIT_REACHED' if dyn.opposing_wall_clusters.any() else 'NO_GEOMETRIC_WEDGE_EVIDENCE')
timeline=dyn[['frame','time_s','actual_min_gap_mm','near5_nodes','near20_nodes','penetrating_nodes','opposing_wall_clusters']].copy(); timeline['first_opposing_bridge_s']=first_bridge; timeline['classification']=classification; write(timeline,'jam_event_timeline.csv')
write(pd.DataFrame([dict(classification=classification,first_opposing_bridge_s=first_bridge,max_theta_deg=dyn.theta_robot_deg.max(),max_near20_nodes=dyn.near20_nodes.max(),max_penetrating_nodes=dyn.penetrating_nodes.max())]),'jam_classification.csv')

# figures
fig,axp=plt.subplots(figsize=(8,5)); sc=axp.scatter(feas.chi_deg,feas.theta_deg,c=feas.g_best_mm*1e3,s=12,cmap='coolwarm',vmin=-20,vmax=50); plt.colorbar(sc,ax=axp,label='best minimum gap (µm)'); axp.set(xlabel='azimuth χ (deg)',ylabel='tilt θ (deg)',title='SmoothWall114 orientation feasibility (coarse exact mesh)'); fig.tight_layout(); fig.savefig(OUT/'theta_feasibility_polar_map.png',dpi=180); plt.close(fig)
fig,axp=plt.subplots(figsize=(8,5)); sc=axp.scatter(feas.chi_deg,feas.theta_deg,c=feas.g_best_mm*1e3,s=10,cmap='coolwarm',vmin=-20,vmax=50); axp.plot(dyn.chi_robot_deg,dyn.theta_robot_deg,'k-',lw=1.5,label='dynamic path'); axp.legend(); axp.set(xlabel='χ (deg)',ylabel='θ (deg)',title='Dynamic trajectory over feasibility map'); fig.colorbar(sc,ax=axp,label='best minimum gap (µm)'); fig.tight_layout(); fig.savefig(OUT/'actual_trajectory_on_feasibility_map.png',dpi=180); plt.close(fig)
fig,axp=plt.subplots(2,1,sharex=True,figsize=(9,6)); axp[0].plot(dyn.time_s*1e3,dyn.actual_min_gap_mm*1e3); axp[0].axhline(0,color='k'); axp[0].set_ylabel('min gap (µm)'); axp[1].plot(dyn.time_s*1e3,dyn.theta_robot_deg); axp[1].set_ylabel('θ (deg)'); axp[1].set_xlabel('time (ms)'); fig.tight_layout(); fig.savefig(OUT/'theta_gap_phase_timeline.png',dpi=180); plt.close(fig)

report=f'''# Wobble30Hz geometric jam audit\n\nClassification: **{classification}**\n\n## Actual robot geometry\n```text\n{geom.to_string(index=False)}\n```\n\n## Local SmoothWall114 section\n```text\n{local.to_string(index=False)}\n```\n\nThe wall section uses the validated Abaqus-frame R3D4 triangles, split into 1,528 triangles, and a ±0.45 mm tangent slab. The orientation map uses the full 334-node exterior mesh, exact point-to-triangle nearest distances, and a ±0.24 mm transverse COM-offset grid. It is a coarse feasibility map (θ=2°, χ=10°); the geometric, 5 µm, and 20 µm thresholds are explicit in `orientation_feasibility_map.csv`.\n\n## Dynamic overlay\n- first opposing-wall bridge: `{first_bridge}` s\n- maximum tilt: {dyn.theta_robot_deg.max():.3f}°\n- maximum near-20-µm exterior nodes: {int(dyn.near20_nodes.max())}\n- maximum penetrating exterior nodes: {int(dyn.penetrating_nodes.max())}\n\nThe dynamic path is overlaid without altering the physical model. A cone-25 job is **not** submitted by this audit; it is permitted only if this geometry evidence is accepted as the cause of the tilt-and-jam.\n'''
(OUT/'Wobble30Hz_geometric_jam_and_cone25_report.md').write_text(report,encoding='utf-8')
print(report)
