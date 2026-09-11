"""Zero-cost Stage-A audit for phase-synchronised directional anti-impact design."""
from pathlib import Path
import math
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
JOB = "WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003"
GDIR = ROOT / "abaqus_robot" / "calibration_analysis" / JOB
def path(name):
    q=ROOT/name
    return q if q.exists() else GDIR/name

tele=pd.read_csv(path(JOB+"_telemetry.csv"))
vel=pd.read_csv(path(JOB+"_rp_velocity_history.csv"))
cp=pd.read_csv(GDIR/"cpress_peaks.csv")
gap=pd.read_csv(path(JOB+"_exact_wall_penetration.csv"))
rf=pd.read_csv(GDIR/"robot_frames.csv")
wn=pd.read_csv(path(JOB+"_pipe_wall_nodes_exact.csv"))
wt=pd.read_csv(path(JOB+"_pipe_wall_triangles_exact.csv"))
cl=pd.read_csv(ROOT/"curvenew_CEL_xyrot56_exact_abaqus.csv")
imp_path=ROOT/"WobbleCal_F30_Cone30_B10_WallOn_Free_003_wall_impulse_audit.csv"
if not imp_path.exists(): imp_path=ROOT/(JOB+"_wall_impulse_audit.csv")
imp=pd.read_csv(imp_path) if imp_path.exists() else pd.DataFrame()

def unit(v):
    v=np.asarray(v,float); n=np.linalg.norm(v)
    return v/n if n else np.zeros(3)
def tangent_at_s(s):
    a=cl.arclength_mm.values; xyz=cl[["x_mm","y_mm","z_mm"]].values
    i=int(np.clip(np.searchsorted(a,s)-1,0,len(a)-2)); return unit(xyz[i+1]-xyz[i])
def closest_tri(q,tri):
    # brute-force barycentric closest point, adequate for the 764 audit faces
    a,b,c=tri; ab=b-a; ac=c-a; ap=q-a; d1=ab@ap; d2=ac@ap
    if d1<=0 and d2<=0:return a
    bp=q-b; d3=ab@bp; d4=ac@bp
    if d3>=0 and d4<=d3:return b
    vc=d1*d4-d3*d2
    if vc<=0 and d1>=0 and d3<=0:return a+(d1/(d1-d3))*ab
    cp=q-c; d5=ab@cp; d6=ac@cp
    if d6>=0 and d5<=d6:return c
    vb=d5*d2-d1*d6
    if vb<=0 and d2>=0 and d6<=0:return a+(d2/(d2-d6))*ac
    va=d3*d6-d5*d4
    if va<=0 and d4>=d3 and d5>=d6:return b+((d4-d3)/((d4-d3)+(d5-d6)))*(c-b)
    den=1/(va+vb+vc); v=vb*den; w=vc*den; return a+ab*v+ac*w

node_xyz={int(r.node):np.array([r.x_mm,r.y_mm,r.z_mm],float) for r in wn.itertuples()}
tris=[]
for r in wt.itertuples():
    ids=[int(r.n1),int(r.n2),int(r.n3)]; tris.append((int(r.element),np.array([node_xyz[i] for i in ids])))
frames={int(k):g for k,g in rf.groupby("frame")}
def geometry(fr,node):
    g=frames.get(fr)
    if g is None:return None
    rr=g[g.node==node]; rr=rr.iloc[0] if len(rr) else g.iloc[0]
    q=rr[["x_mm","y_mm","z_mm"]].values.astype(float)
    best=min((np.linalg.norm(q-closest_tri(q,t)),eid,t,closest_tri(q,t)) for eid,t in tris)
    d,eid,tri,cpnt=best; cen=tri.mean(0); xyz=cl[["x_mm","y_mm","z_mm"]].values
    cpt=xyz[np.argmin(np.linalg.norm(xyz-cen,axis=1))]; n_away=unit(cpt-cen)
    return q,eid,cpnt,n_away,d

rows=[]
for r in cp.itertuples():
    if r.cpress_max_mpa<0.02 or not (0.0008 <= r.time_s <= 0.0020): continue
    fr=int(r.frame); g=geometry(fr,int(r.node))
    if g is None:continue
    q,eid,cpnt,naway,d=g; v=vel.iloc[min(fr,len(vel)-1)]; vcom=np.array([v.V1,v.V2,v.V3]); omg=np.array([v.VR1,v.VR2,v.VR3])
    xyz=frames[fr][["x_mm","y_mm","z_mm"]].values; com=xyz.mean(0); axis=unit(xyz[-1]-xyz[0]); vc=vcom+np.cross(omg,q-com)
    t_hat=tangent_at_s(float(tele.iloc[min(fr,len(tele)-1)].robot_arc_mm)); gr=gap.iloc[min(fr,len(gap)-1)]; tr=tele.iloc[min(fr,len(tele)-1)]
    rows.append(dict(frame=fr,time_s=r.time_s,cpress_max_mpa=r.cpress_max_mpa,contact_node=int(r.node),wall_element=eid,contact_x_mm=q[0],contact_y_mm=q[1],contact_z_mm=q[2],
        wall_normal_away_x=naway[0],wall_normal_away_y=naway[1],wall_normal_away_z=naway[2],vcom_n_mm_s=vcom@naway,vcontact_n_mm_s=vc@naway,vcontact_t_mm_s=vc@t_hat,
        t_center_x=t_hat[0],t_center_y=t_hat[1],t_center_z=t_hat[2],t_com_x=axis[0],t_com_y=axis[1],t_com_z=axis[2],n_dot_tcenter=naway@t_hat,n_dot_tcom=naway@axis,
        r_contact_x_mm=(q-com)[0],r_contact_y_mm=(q-com)[1],r_contact_z_mm=(q-com)[2],signed_gap_mm=gr.min_signed_gap_mm,nearest_wall_element=int(gr.nearest_wall_element),
        field_phase_rad=tr.instantaneous_phase_rad,field_phase_deg=np.degrees(tr.instantaneous_phase_rad)%360,geom_distance_mm=d))
audit=pd.DataFrame(rows); audit.to_csv(ROOT/"first_impact_vector_audit.csv",index=False); audit.to_csv(ROOT/"contact_normal_tangent_geometry.csv",index=False)

phase=[]
for i,r in tele.iterrows():
    if not (0.0008<=r.t_s<=0.0020):continue
    j=min(int(round(r.t_s/0.0001)),len(cp)-1); k=min(int(round(r.t_s/0.0001)),len(gap)-1)
    phase.append(dict(time_s=r.t_s,phase_rad=r.instantaneous_phase_rad,phase_deg=np.degrees(r.instantaneous_phase_rad)%360,cpress_mpa=cp.iloc[j].cpress_max_mpa,min_signed_gap_mm=gap.iloc[k].min_signed_gap_mm,Fmag_t_uN=r.force_tangent_N*1e6,Fmag_n_uN=r.force_normal_N*1e6,active_head=r.active_head,active_tail=r.active_tail))
pd.DataFrame(phase).to_csv(ROOT/"first_impact_phase_map.csv",index=False)

if len(imp):
    ii=imp.copy(); t=ii.time_s.values
    # The wall-impulse audit columns J* are already cumulative impulses (N*s),
    # not instantaneous forces; use their final values rather than integrating twice.
    integ=lambda c:float(ii[c].fillna(0).iloc[-1])
    th=tangent_at_s(float(tele.iloc[0].robot_arc_mm)); vf=vel.iloc[-1]; v0=vel.iloc[0]; dv=float(np.dot(np.array([vf.V1-v0.V1,vf.V2-v0.V2,vf.V3-v0.V3]),th)); mdv=10e-6*dv/1000.0
    jm=integ("Jmag_t_mN_s") if "Jmag_t_mN_s" in ii else 0.0; jc=integ("Jwall_on_robot_t_mN_s") if "Jwall_on_robot_t_mN_s" in ii else 0.0
    pd.DataFrame([dict(window_s=t[-1]-t[0],mass_kg=1e-5,delta_vt_mm_s=dv,m_delta_vt_mN_s=mdv,Jmag_t_mN_s=jm,Jcontact_t_mN_s=jc,resolved_sum_mN_s=jm+jc,residual_mN_s=mdv-jm-jc,relative_residual=(mdv-jm-jc)/max(abs(mdv),1e-12))]).to_csv(ROOT/"impulse_momentum_balance.csv",index=False)
    ii.to_csv(ROOT/"impulse_momentum_timeseries.csv",index=False)

if len(audit):
    fig,ax=plt.subplots(2,1,figsize=(8,6),sharex=True); ax[0].plot(audit.time_s*1e3,audit.vcontact_n_mm_s,'o-',label='Vcontact·n_away'); ax[0].axhline(0,color='k',lw=.7); ax[0].legend(); ax[1].plot(audit.time_s*1e3,audit.cpress_max_mpa,'o-',label='CPRESS max'); ax[1].legend(); ax[1].set_xlabel('time (ms)'); fig.tight_layout(); fig.savefig(ROOT/"first_impact_vector_geometry.png",dpi=180); plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4)); ax.plot(tele.t_s*1e3,np.degrees(tele.instantaneous_phase_rad)%360); ax.set_xlabel('time (ms)'); ax.set_ylabel('field phase (deg)'); fig.tight_layout(); fig.savefig(ROOT/"gap_cpress_Vt_vs_field_phase.png",dpi=180); plt.close(fig)

def min_norm_tensor(m,F):
    m=np.asarray(m,float); F=np.asarray(F,float)
    C=np.array([[m[0],0,0,m[1],m[2],0],[0,m[1],0,m[0],0,m[2]],[0,0,m[2],0,m[0],m[1]],[1,1,1,0,0,0]],float)
    x=np.linalg.lstsq(C,np.r_[F,0.0],rcond=None)[0]; return np.array([[x[0],x[3],x[4]],[x[3],x[1],x[5]],[x[4],x[5],x[2]]])
robot_axis=unit(cl[["x_mm","y_mm","z_mm"]].iloc[1].values-cl[["x_mm","y_mm","z_mm"]].iloc[0].values)
mvec=-1.168e-3*robot_axis
if len(audit):
    # Center just before the 1.3 ms CPRESS peak; fixed open-loop phase law.
    phase0=float(audit.loc[audit.cpress_max_mpa.idxmax(),"field_phase_deg"])-3.0
else: phase0=259.0
phase0%=360.0; halfwidth=7.0; ft_uN=400.0; fn_uN=3000.0
t0=tangent_at_s(float(tele.iloc[0].robot_arc_mm)); n0=unit(audit.iloc[0][["wall_normal_away_x","wall_normal_away_y","wall_normal_away_z"]].values) if len(audit) else np.array([0.,0.,1.]); Ftar=np.r_[ft_uN*1e-6*t0+fn_uN*1e-6*n0]; A=min_norm_tensor(mvec,Ftar); eig=np.linalg.eigvalsh(A)
pd.DataFrame([dict(phase_center_deg=phase0,ft_target_uN=ft_uN,fn_target_uN=fn_uN,m_x_Am2=mvec[0],m_y_Am2=mvec[1],m_z_Am2=mvec[2],t_hat_x=t0[0],t_hat_y=t0[1],t_hat_z=t0[2],n_away_x=n0[0],n_away_y=n0[1],n_away_z=n0[2],Axx_T_per_m=A[0,0],Ayy_T_per_m=A[1,1],Azz_T_per_m=A[2,2],Axy_T_per_m=A[0,1],Axz_T_per_m=A[0,2],Ayz_T_per_m=A[1,2],eig1_T_per_m=eig[0],eig2_T_per_m=eig[1],eig3_T_per_m=eig[2],fro_norm_T_per_m=np.linalg.norm(A),trace_T_per_m=np.trace(A),reconstruction_error_pct=100*np.linalg.norm(A@mvec-Ftar)/np.linalg.norm(Ftar),max_Baux_0p45mm_T=np.linalg.norm(A)*0.00045)]).to_csv(ROOT/"directional_tensor_design.csv",index=False)
pd.DataFrame([dict(profile="raised_cosine",phase_center_deg=phase0,half_width_deg=halfwidth,phase_start_deg=(phase0-halfwidth)%360,phase_end_deg=(phase0+halfwidth)%360,target_ft_uN=ft_uN,target_fn_uN=fn_uN,law="w=max(0,0.5*(1+cos(pi*dphi/halfwidth)))",open_loop=True)]).to_csv(ROOT/"phase_pulse_definition.csv",index=False)
R=np.asarray(json.load(open(ROOT/"abaqus_magpylib_frame_transform.json"))["R_aba_to_mag"],float)
nmag=unit(R@n0)
pd.DataFrame([dict(phase_center_deg=phase0,halfwidth_deg=halfwidth,ft_target_N=ft_uN*1e-6,fn_target_N=fn_uN*1e-6,normal_mag_x=nmag[0],normal_mag_y=nmag[1],normal_mag_z=nmag[2],coordinate_frame="MAGPYLIB_DXF",source_frame="ABAQUS",frame_transform="R_aba_to_mag")]).to_csv(ROOT/"directional_pulse_config.csv",index=False)

# Offline impulse estimate with the measured G6 phase cadence.
tt=np.linspace(0.0008,0.0020,1201); ph=np.interp(tt,tele.t_s,tele.instantaneous_phase_rad); ww=[]
for z in ph:
    dd=math.atan2(math.sin(z-math.radians(phase0)),math.cos(z-math.radians(phase0))); ww.append(.5*(1+math.cos(math.pi*dd/math.radians(halfwidth))) if abs(dd)<math.radians(halfwidth) else 0.)
ww=np.asarray(ww); impulse_t=float(np.trapz(ft_uN*1e-6*ww,tt)); impulse_n=float(np.trapz(fn_uN*1e-6*ww,tt)); vclose=float(abs(audit.vcontact_n_mm_s.iloc[audit.cpress_max_mpa.argmax()])*1e-3) if len(audit) else .5; dvn=impulse_n/1e-5
pd.DataFrame([dict(phase_center_deg=phase0,halfwidth_deg=halfwidth,window_start_s=tt[ww>0][0] if np.any(ww>0) else np.nan,window_end_s=tt[ww>0][-1] if np.any(ww>0) else np.nan,positive_tangential_impulse_Ns=impulse_t,away_normal_impulse_Ns=impulse_n,estimated_delta_v_away_m_s=dvn,impact_closing_speed_m_s=vclose,predicted_closing_reduction_pct=100*dvn/max(vclose,1e-12),gate_positive_tangent=impulse_t>0,gate_closing_reduction_30pct=(dvn/max(vclose,1e-12)>=.30))]).to_csv(ROOT/"phase_pulse_offline_impulse_estimate.csv",index=False)
# Field-direction and tensor-norm diagnostic curves for the selected one pulse.
fig,ax=plt.subplots(figsize=(8,4)); ax.plot(np.degrees(ph)%360,ww); ax.axvline(phase0,color='r',ls='--'); ax.set_xlabel('field phase (deg)'); ax.set_ylabel('pulse weight'); fig.tight_layout(); fig.savefig(ROOT/"target_force_direction.png",dpi=180); plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4)); norms=[]
for z in ph:
    norms.append(np.linalg.norm(min_norm_tensor(mvec,ft_uN*1e-6*t0+fn_uN*1e-6*n0)))
ax.plot(tt*1e3,norms); ax.set_xlabel('time (ms)'); ax.set_ylabel('||A|| (T/m)'); fig.tight_layout(); fig.savefig(ROOT/"gradient_tensor_norm_vs_phase.png",dpi=180); plt.close(fig)
print("Stage-A audit complete; contact rows",len(audit));
if len(audit):print(audit[["frame","time_s","cpress_max_mpa","contact_node","wall_element","vcontact_n_mm_s","n_dot_tcenter","field_phase_deg"]].to_string(index=False))
