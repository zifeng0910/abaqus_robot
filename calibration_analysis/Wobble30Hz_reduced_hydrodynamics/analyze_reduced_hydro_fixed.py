from pathlib import Path
import math, json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
JOB='Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083'; OUT=HERE
def unit(x):
    x=np.asarray(x,float); return x/max(np.linalg.norm(x),1e-30)
def rvmat(v):
    v=np.asarray(v,float); a=np.linalg.norm(v)
    if a<1e-14:return np.eye(3)
    k=v/a; K=np.array([[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[1],k[0],0.]])
    return np.eye(3)+math.sin(a)*K+(1-math.cos(a))*(K@K)
def logrot(R):
    a=math.acos(float(np.clip((np.trace(R)-1.)/2.,-1.,1.)))
    if a<1e-12:return np.zeros(3)
    return np.array([R[2,1]-R[1,2],R[0,2]-R[2,0],R[1,0]-R[0,1]])*a/(2*math.sin(a))
def interp(df,t,c):
    tc='time_s' if 'time_s' in df.columns else 't_s'
    return np.interp(t,df[tc].to_numpy(float),df[c].to_numpy(float))

rp=pd.read_csv(ROOT/(JOB+'_rp_fields.csv')); tel=pd.read_csv(ROOT/(JOB+'_telemetry.csv'))
cp=pd.read_csv(ROOT/(JOB+'_contact.csv')); gap=pd.read_csv(ROOT/(JOB+'_exact_wall_penetration.csv'))
hyd=pd.read_csv(ROOT/'reduced_hydro_runtime_load.csv',skiprows=2,skipinitialspace=True,
                names=['time_s','v1_mm_s','v2_mm_s','v3_mm_s','vr1','vr2','vr3',
                       'Fh1_N','Fh2_N','Fh3_N','Th1_Nmm','Th2_Nmm','Th3_Nmm'])
t=rp.time_s.to_numpy(float); ur=rp[['UR1','UR2','UR3']].to_numpy(float); rot=np.array([rvmat(x) for x in ur])
V=rp[['V1','V2','V3']].to_numpy(float); VR=rp[['VR1','VR2','VR3']].to_numpy(float)
om=np.zeros_like(ur)
for i in range(len(t)):
    j=1 if i==0 else i-1; k=i+1 if i<len(t)-1 else i; om[i]=logrot(rot[k]@rot[j].T)/max(t[k]-t[j],1e-30)
omega_diff=np.linalg.norm(om-VR,axis=1)

# True directed axis from the actual exterior node clouds (fixed end clusters).
rf=pd.read_csv(ROOT/'output'/(JOB+'_bend_validation')/'robot_frames.csv'); g0=rf[rf.frame==rf.frame.min()]
a_seed=unit([0.9647382600216,-0.1188742372140,0.2348382536499]); q=g0[['x_mm','y_mm','z_mm']].to_numpy(float)@a_seed; lo,hi=np.quantile(q,[.1,.9])
tail=set(g0.loc[q<=lo,'node']); head=set(g0.loc[q>=hi,'node']); heads=[]; tails=[]
for _,g in rf.groupby('frame',sort=True):
    heads.append(g[g.node.isin(head)][['x_mm','y_mm','z_mm']].to_numpy(float).mean(0)); tails.append(g[g.node.isin(tail)][['x_mm','y_mm','z_mm']].to_numpy(float).mean(0))
heads=np.asarray(heads); tails=np.asarray(tails); axis=np.array([unit(x-y) for x,y in zip(heads,tails)])

# Canonical global curve and continuous segment projection.
curve=pd.read_csv(ROOT/'CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_EXTSURF_centerline_odb.csv'); cv=curve[['x_mm','y_mm','z_mm']].to_numpy(float); S=curve.arclength_mm.to_numpy(float)
def project(p,prev):
    best=None
    for j in range(len(cv)-1):
        d=cv[j+1]-cv[j]; u=np.clip(np.dot(p-cv[j],d)/max(np.dot(d,d),1e-30),0,1); x=cv[j]+u*d; s=S[j]+u*(S[j+1]-S[j]); cost=np.dot(p-x,p-x)+(.02*(s-prev)**2 if prev is not None else 0)
        if best is None or cost<best[0]: best=(cost,s,x,unit(d))
    return best[1],best[2],best[3]
x0=np.array([-7.468174204284,-3.676918015967,-9.550745259298]); com=np.column_stack([x0[0]+rp.U1,x0[1]+rp.U2,x0[2]+rp.U3]); ss=[]; pp=[]; tang=[]; prev=None
for p in com:
    z=project(p,prev); ss.append(z[0]); pp.append(z[1]); tang.append(z[2]); prev=z[0]
ss=np.asarray(ss); pp=np.asarray(pp); tang=np.asarray(tang); e1=[]; e2=[]
for tt in tang:
    ee=np.array([0.,0.,1.])-tt[2]*tt
    if np.linalg.norm(ee)<1e-8: ee=np.array([0.,1.,0.])-tt[1]*tt
    ee=unit(ee); e1.append(ee); e2.append(unit(np.cross(tt,ee)))
e1=np.asarray(e1); e2=np.asarray(e2)
phi=np.unwrap(np.arctan2(np.einsum('ij,ij->i',axis,e2),np.einsum('ij,ij->i',axis,e1)))
B=np.column_stack([interp(tel,t,c) for c in ['Bx_aba_T','By_aba_T','Bz_aba_T']]); F=np.column_stack([interp(tel,t,c) for c in ['fx_aba_N','fy_aba_N','fz_aba_N']]); T=np.column_stack([interp(tel,t,c) for c in ['tx_aba_Nmm','ty_aba_Nmm','tz_aba_Nmm']]); phiB=np.unwrap(np.arctan2(np.einsum('ij,ij->i',B,e2),np.einsum('ij,ij->i',B,e1))); phicmd=interp(tel,t,'instantaneous_phase_rad')
def slope(tt,pp): return np.polyfit(tt-tt.mean(),pp,1)[0]/(2*np.pi) if len(tt)>=3 else np.nan
def wf(a,b):
    m=(t>=a)&(t<b); return m,slope(t[m],phi[m])
rates=np.full(len(t),np.nan)
for a,b in [(0,.001),(.001,.002),(.002,.004),(.004,.006333333),(.006333333,.0083331)]: rates[wf(a,b)[0]]=wf(a,b)[1]

# Hydro and rigid-body angular momentum/work in SI-consistent units.
fh=np.column_stack([interp(hyd,t,c) for c in ['Fh1_N','Fh2_N','Fh3_N']]); hh=np.column_stack([interp(hyd,t,c) for c in ['Th1_Nmm','Th2_Nmm','Th3_Nmm']]); hv=np.column_stack([interp(hyd,t,c) for c in ['v1_mm_s','v2_mm_s','v3_mm_s']])
m=1e-5; r=.61e-3; L=2.81e-3; Iax=.5*m*r*r; Itr=m*(3*r*r+L*L)/12.; xx=unit(np.cross([0,0,1.],a_seed)); yy=unit(np.cross(a_seed,xx)); Q=np.column_stack([xx,yy,a_seed]); Ibody=Q@np.diag([Itr,Itr,Iax])@Q.T; Ig=np.array([R@Ibody@R.T for R in rot]); H=np.einsum('nij,nj->ni',Ig,om); K=.5*np.einsum('ni,ni->n',om,H)
J=np.zeros_like(T); Jh=np.zeros_like(T); W=np.zeros(len(t)); Wh=np.zeros(len(t))
for i in range(1,len(t)):
    dt=t[i]-t[i-1]; J[i]=J[i-1]+.5*(T[i]+T[i-1])*dt; Jh[i]=Jh[i-1]+.5*(hh[i]+hh[i-1])*dt; W[i]=W[i-1]+.5*(T[i]@om[i]+T[i-1]@om[i-1])*dt*1e-3; Wh[i]=Wh[i-1]+.5*(hh[i]@om[i]+hh[i-1]@om[i-1])*dt*1e-3
Wnon=(K-K[0])-W; P_h=np.einsum('ij,ij->i',fh,hv); P_ht=np.einsum('ij,ij->i',hh,om); Vt=np.einsum('ij,ij->i',V,tang); Vs=np.gradient(ss,t)
data=pd.DataFrame({'time_s':t,'s_COM_mm':ss,'projection_residual_mm':np.linalg.norm(com-pp,axis=1),'Vt_direct_mm_s':Vt,'Vt_from_s_mm_s':Vs,'UR1':ur[:,0],'UR2':ur[:,1],'UR3':ur[:,2],'HEAD_x_mm':heads[:,0],'HEAD_y_mm':heads[:,1],'HEAD_z_mm':heads[:,2],'TAIL_x_mm':tails[:,0],'TAIL_y_mm':tails[:,1],'TAIL_z_mm':tails[:,2],'axis_x':axis[:,0],'axis_y':axis[:,1],'axis_z':axis[:,2],'phi_robot_rad':phi,'phi_B_local_rad':phiB,'phi_cmd_rad':phicmd,'f_robot_phase_equiv_Hz':rates,'omega_direct_norm_rad_s':np.linalg.norm(VR,axis=1),'omega_relative_norm_rad_s':np.linalg.norm(om,axis=1),'omega_wobble_norm_rad_s':np.linalg.norm(om-axis*np.einsum('ij,ij->i',om,axis)[:,None],axis=1),'omega_spin_rad_s':np.einsum('ij,ij->i',om,axis),'omega_reconstruction_diff_rad_s':omega_diff,'H_x_kgm2_s':H[:,0],'H_y_kgm2_s':H[:,1],'H_z_kgm2_s':H[:,2],'Krot_J':K,'Jmag_x_Nmm_s':J[:,0],'Jmag_y_Nmm_s':J[:,1],'Jmag_z_Nmm_s':J[:,2],'Jmag_norm_Nmm_s':np.linalg.norm(J,axis=1),'Jhydro_norm_Nmm_s':np.linalg.norm(Jh,axis=1),'Wmag_J':W,'Whydro_J':Wh,'Wnonmag_resolved_J':Wnon,'Fhydro_power_W':P_h,'Thydro_power_W':P_ht,'Tmag_norm_Nmm':np.linalg.norm(T,axis=1),'Fmag_norm_N':np.linalg.norm(F,axis=1)})
data.to_csv(OUT/'reduced_hydro_fixed_true_axis_phase.csv',index=False); data[['time_s','s_COM_mm','Vt_direct_mm_s','Vt_from_s_mm_s','projection_residual_mm']].to_csv(OUT/'reduced_hydro_fixed_canonical_translation.csv',index=False); data.to_csv(OUT/'reduced_hydro_fixed_angular_momentum.csv',index=False); data[['time_s','Jmag_x_Nmm_s','Jmag_y_Nmm_s','Jmag_z_Nmm_s','Jmag_norm_Nmm_s','Jhydro_norm_Nmm_s']].to_csv(OUT/'reduced_hydro_fixed_angular_impulse.csv',index=False); data[['time_s','Krot_J','Wmag_J','Whydro_J','Wnonmag_resolved_J','Fhydro_power_W','Thydro_power_W']].to_csv(OUT/'reduced_hydro_fixed_rotational_work.csv',index=False)
du=np.column_stack([np.gradient(rp[c].to_numpy(float),t) for c in ['U1','U2','U3']]); pd.DataFrame({'time_s':t,'V1':V[:,0],'V2':V[:,1],'V3':V[:,2],'dU1_dt':du[:,0],'dU2_dt':du[:,1],'dU3_dt':du[:,2],'V_vs_dU_abs_mm_s':np.linalg.norm(V-du,axis=1),'VR1_direct':VR[:,0],'VR2_direct':VR[:,1],'VR3_direct':VR[:,2]}).to_csv(OUT/'reduced_hydro_fixed_sensor_identity.csv',index=False)
pd.DataFrame([{'mass_target_mg':10.,'mass_target_kg':m,'robot_density_tonne_mm3':7.80906654321e-9,'inertia_source':'cylindrical envelope reconstruction; deck has no independent RP MASS card','Ibody_xx_kgm2':Ibody[0,0],'Ibody_yy_kgm2':Ibody[1,1],'Ibody_zz_kgm2':Ibody[2,2]}]).to_csv(OUT/'reduced_hydro_fixed_mass_inertia_audit.csv',index=False)

# Contact windows and corrected classification gate (n is row count, not len(Series)).
contact=pd.DataFrame({'time_s':gap.time_s,'exact_gap_mm':gap.min_signed_gap_mm,'nearest_node':gap.signed_node,'nearest_wall_element':gap.signed_wall_element,'CPRESS_max_MPa':cp.CPRESS_max_MPa,'CNORMF_max_N':cp.CNORMF_max_N,'CSHEARF_max_N':cp.CSHEARF_max_N}); contact.to_csv(OUT/'reduced_hydro_fixed_contact_timeline.csv',index=False)
g0t=float(contact.loc[contact.exact_gap_mm<=0,'time_s'].iloc[0]) if (contact.exact_gap_mm<=0).any() else np.nan; g5=float(contact.loc[contact.exact_gap_mm<=.005,'time_s'].iloc[0]) if (contact.exact_gap_mm<=.005).any() else np.nan; ctime=float(contact.loc[contact.CPRESS_max_MPa>.01,'time_s'].iloc[0]) if (contact.CPRESS_max_MPa>.01).any() else np.nan
jam=contact.assign(near_wall=contact.exact_gap_mm<=.005,geometric_contact=contact.exact_gap_mm<=0,solver_contact=contact.CPRESS_max_MPa>.01); jam.to_csv(OUT/'reduced_hydro_fixed_jam_timeline.csv',index=False)
summary=[]
for name,a,b in [('W0_precontact',0,g0t if np.isfinite(g0t) else .003),('W1_contact',g0t,g0t+.0005 if np.isfinite(g0t) else .0035),('W2_immediate_post',g0t+.0005 if np.isfinite(g0t) else .0035,g0t+.002 if np.isfinite(g0t) else .005),('W4_late',max(.006333333,g0t+.002 if np.isfinite(g0t) else .006333333),.0083331)]:
    if not np.isfinite(a): continue
    mask=(t>=a)&(t<b); summary.append({'window':name,'start_s':a,'end_s':b,'duration_ms':(b-a)*1000,'n':int(mask.sum()),'f_cmd_Hz':30.,'f_B_local_Hz':slope(t[mask],phiB[mask]),'f_robot_phase_equiv_Hz':slope(t[mask],phi[mask]),'R_phase_vs_30':abs(slope(t[mask],phi[mask]))/30 if mask.sum()>=3 else np.nan,'Tmag_RMS_Nmm':float(np.sqrt(np.mean(np.linalg.norm(T[mask],axis=1)**2))) if mask.any() else np.nan,'omega_wobble_RMS_rad_s':float(np.sqrt(np.mean(data.omega_wobble_norm_rad_s[mask]**2))) if mask.any() else np.nan,'gap_min_um':float(contact.loc[mask,'exact_gap_mm'].min()*1000) if mask.any() else np.nan,'CPRESS_max_MPa':float(contact.loc[mask,'CPRESS_max_MPa'].max()) if mask.any() else np.nan})
pd.DataFrame(summary).to_csv(OUT/'reduced_hydro_fixed_window_summary.csv',index=False)

# Figures required by the brief.
fig,ax=plt.subplots(2,1,figsize=(8,6),sharex=True); ax[0].plot(t,phi/(2*np.pi),label='robot local phase'); ax[0].plot(t,phiB/(2*np.pi),label='B local phase'); ax[0].legend(); ax[1].plot(t,data.omega_wobble_norm_rad_s,label='|omega_wobble|'); ax[1].set_xlabel('time (s)'); ax[1].legend(); fig.tight_layout(); fig.savefig(OUT/'reduced_hydro_fixed_phase_wobble.png',dpi=160); plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4)); ax.plot(t,K*1e6,label='Krot (uJ)'); ax.plot(t,W*1e6,label='Wmag (uJ)'); ax.plot(t,Wnon*1e6,label='Wnonmag resolved (uJ)'); ax.legend(); ax.set_xlabel('time (s)'); fig.tight_layout(); fig.savefig(OUT/'reduced_hydro_fixed_rotational_energy_work.png',dpi=160); plt.close(fig)
print('WROTE corrected reduced-hydro analyses; contact onset gap<=0:',g0t,'near-wall:',g5,'CPRESS:',ctime,'omega diff RMS:',float(np.sqrt(np.mean(omega_diff**2))))
