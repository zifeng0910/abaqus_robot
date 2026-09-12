from pathlib import Path
import json, math
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
pd.DataFrame.to_markdown=lambda self,**kwargs:self.to_string(index=False)

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
JOB='Wobble_F30_G6L45_ReducedHydro_WallOn_Free_0083'; OUT=HERE

def unit(v):
    v=np.asarray(v,float); return v/max(np.linalg.norm(v),1.e-15)
def rvmat(v):
    v=np.asarray(v,float); a=np.linalg.norm(v)
    if a<1.e-14:return np.eye(3)
    k=v/a; K=np.array([[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[1],k[0],0.]])
    return np.eye(3)+math.sin(a)*K+(1-math.cos(a))*(K@K)
def logrot(R):
    a=math.acos(float(np.clip((np.trace(R)-1)/2,-1,1)))
    if a<1.e-12:return np.zeros(3)
    return np.array([R[2,1]-R[1,2],R[0,2]-R[2,0],R[1,0]-R[0,1]])*a/(2*math.sin(a))

tel=pd.read_csv(ROOT/(JOB+'_telemetry.csv')); rp=pd.read_csv(ROOT/(JOB+'_rp_fields.csv'))
t=rp.time_s.to_numpy(float); ur=rp[['UR1','UR2','UR3']].to_numpy(float)
rot=np.array([rvmat(x) for x in ur]); om=np.zeros_like(ur)
om[0]=logrot(rot[1]@rot[0].T)/(t[1]-t[0]); om[-1]=logrot(rot[-1]@rot[-2].T)/(t[-1]-t[-2])
for i in range(1,len(t)-1): om[i]=logrot(rot[i+1]@rot[i-1].T)/(t[i+1]-t[i-1])
om_direct=rp[['VR1','VR2','VR3']].to_numpy(float); omega_diff=np.linalg.norm(om-om_direct,axis=1)
a0=unit([0.9647382600216,-0.1188742372140,0.2348382536499]); axis=np.array([unit(R@a0) for R in rot])
spin=np.einsum('ij,ij->i',om,axis); wob=om-axis*spin[:,None]

curve=pd.read_csv(ROOT/'curvenew_CEL_xyrot56_exact.csv'); cvm=curve[['x_mm','y_mm','z_mm']].to_numpy(float); S=curve.arclength_mm.to_numpy(float)
tf=json.loads((ROOT/'abaqus_magpylib_frame_transform.json').read_text()); Rtf=np.asarray(tf['R_aba_to_mag']); org=np.asarray(tf['origin_aba_mm'])
def tangent_at(s):
    j=max(0,min(int(np.searchsorted(S,s,side='right')-1),len(S)-2)); return unit(Rtf.T@(cvm[j+1]-cvm[j]))
td=np.interp(t,tel.t_s,tel.driver_arc_mm); tang=np.array([tangent_at(x) for x in td])
rpv=rp[['V1','V2','V3']].to_numpy(float); Vt=np.einsum('ij,ij->i',rpv,tang)
srobot=np.r_[0,np.cumsum(.5*(Vt[1:]+Vt[:-1])*np.diff(t))]

def interp_vec(cols): return np.array([np.interp(t,tel.t_s,tel[c]) for c in cols]).T
B=interp_vec(['Bx_aba_T','By_aba_T','Bz_aba_T']); F=interp_vec(['fx_aba_N','fy_aba_N','fz_aba_N']); T=interp_vec(['tx_aba_Nmm','ty_aba_Nmm','tz_aba_Nmm'])
hyd_cols=['time_s','v1_mm_s','v2_mm_s','v3_mm_s','vr1','vr2','vr3','Fh1_N','Fh2_N','Fh3_N','Th1_Nmm','Th2_Nmm','Th3_Nmm']
hyd=pd.read_csv(ROOT/'reduced_hydro_runtime_load.csv',skiprows=2,names=hyd_cols,skipinitialspace=True); th=hyd.time_s.to_numpy(float)
fh=np.array([np.interp(t,th,hyd[c]) for c in ['Fh1_N','Fh2_N','Fh3_N']]).T; hh=np.array([np.interp(t,th,hyd[c]) for c in ['Th1_Nmm','Th2_Nmm','Th3_Nmm']]).T
vh=np.array([np.interp(t,th,hyd[c]) for c in ['v1_mm_s','v2_mm_s','v3_mm_s']]).T

# 10 mg cylindrical inertia in SI; body axis is a0.
m=1e-5; r=.61e-3; h=2.81e-3; Iax=.5*m*r*r; Itr=m*(3*r*r+h*h)/12.; z=a0
x=unit(np.cross([0,0,1.],z)) if abs(z[2])<.9 else unit(np.cross([0,1.,0],z)); y=unit(np.cross(z,x)); Q=np.column_stack([x,y,z])
Ibody=Q@np.diag([Itr,Itr,Iax])@Q.T; Ig=np.array([R@Ibody@R.T for R in rot]); H=np.einsum('nij,nj->ni',Ig,om); K=.5*np.einsum('ni,ni->n',om,H)
J=np.zeros_like(T); Jh=np.zeros_like(T); W=np.zeros(len(t)); Wh=np.zeros(len(t))
for i in range(1,len(t)):
    dt=t[i]-t[i-1]; J[i]=J[i-1]+.5*(T[i]+T[i-1])*dt; Jh[i]=Jh[i-1]+.5*(hh[i]+hh[i-1])*dt
    W[i]=W[i-1]+.5*(T[i].dot(om[i])+T[i-1].dot(om[i-1]))*dt*1e-3
    Wh[i]=Wh[i-1]+.5*(hh[i].dot(om[i])+hh[i-1].dot(om[i-1]))*dt*1e-3
Wnon=(K-K[0])-W; Ph=np.einsum('ij,ij->i',fh,vh); Pht=np.einsum('ij,ij->i',hh,om)

phase=np.unwrap(np.arctan2(axis[:,1],axis[:,0])); phase_rate=np.gradient(phase,t)/(2*np.pi)
data=pd.DataFrame({'time_s':t,'UR1':ur[:,0],'UR2':ur[:,1],'UR3':ur[:,2],
 'omega_x_rad_s':om[:,0],'omega_y_rad_s':om[:,1],'omega_z_rad_s':om[:,2],
 'VR1':om_direct[:,0],'VR2':om_direct[:,1],'VR3':om_direct[:,2],
 'omega_norm_rad_s':np.linalg.norm(om,axis=1),'omega_wobble_rad_s':np.linalg.norm(wob,axis=1),'omega_spin_rad_s':spin,
 'omega_reconstruction_abs_diff_rad_s':omega_diff,'omega_reconstruction_rel_diff':omega_diff/np.maximum(np.linalg.norm(om_direct,axis=1),1e-9),
 'H_x_kgm2_s':H[:,0],'H_y_kgm2_s':H[:,1],'H_z_kgm2_s':H[:,2],'Krot_J':K,
 'Jmag_x_Nmm_s':J[:,0],'Jmag_y_Nmm_s':J[:,1],'Jmag_z_Nmm_s':J[:,2],'Jmag_norm_Nmm_s':np.linalg.norm(J,axis=1),
 'Jhydro_norm_Nmm_s':np.linalg.norm(Jh,axis=1),'Wmag_J':W,'Whydro_J':Wh,'Wnonmag_resolved_J':Wnon,
 'Fhydro_power_W':Ph,'Thydro_power_W':Pht,'driver_arc_mm':td,'robot_arc_server_mm':np.interp(t,tel.t_s,tel.robot_arc_mm),'srobot_mm':srobot,'Vt_mm_s':Vt,
 'Fmag_norm_N':np.linalg.norm(F,axis=1),'Tmag_norm_Nmm':np.linalg.norm(T,axis=1),'Fhydro_norm_N':np.linalg.norm(fh,axis=1),'Thydro_norm_Nmm':np.linalg.norm(hh,axis=1),'robot_phase_rate_Hz':phase_rate})
data.to_csv(OUT/'reduced_hydro_angular_dynamics.csv',index=False)
data[['time_s','H_x_kgm2_s','H_y_kgm2_s','H_z_kgm2_s','Krot_J']].to_csv(OUT/'reduced_hydro_angular_momentum.csv',index=False)
data[['time_s','Jmag_x_Nmm_s','Jmag_y_Nmm_s','Jmag_z_Nmm_s','Jmag_norm_Nmm_s','Jhydro_norm_Nmm_s']].to_csv(OUT/'reduced_hydro_angular_impulse.csv',index=False)
data[['time_s','Wmag_J','Whydro_J','Wnonmag_resolved_J','Krot_J','Fhydro_power_W','Thydro_power_W']].to_csv(OUT/'reduced_hydro_rotational_work.csv',index=False)
data[['time_s','driver_arc_mm','robot_arc_server_mm','srobot_mm','Vt_mm_s']].to_csv(OUT/'reduced_hydro_translation.csv',index=False)

gap=pd.read_csv(ROOT/(JOB+'_exact_wall_penetration.csv')); cp=pd.read_csv(ROOT/(JOB+'_contact.csv'))
contact=pd.DataFrame({'time_s':gap.time_s,'exact_gap_mm':gap.min_signed_gap_mm,'nearest_node':gap.nearest_node,'nearest_wall_element':gap.nearest_wall_element,'CPRESS_max_MPa':cp.CPRESS_max_MPa,'CNORMF_max_N':cp.CNORMF_max_N,'CSHEARF_max_N':cp.CSHEARF_max_N})
contact.to_csv(OUT/'reduced_hydro_contact_timeline.csv',index=False)

def stats(lo,hi):
    q=data[(t>=lo)&(t<=hi)]; c=contact[(contact.time_s>=lo)&(contact.time_s<=hi)]
    if len(q)<3:return {'duration_ms':(hi-lo)*1000,'n':len(q)}
    fit=np.polyfit(q.time_s,q.robot_phase_rate_Hz,1) if len(q)>=3 else [np.nan,np.nan]
    return {'duration_ms':(hi-lo)*1000,'n':len(q),'phase_rate_mean_Hz':float(q.robot_phase_rate_Hz.mean()),'R_phase_vs_30':float(abs(q.robot_phase_rate_Hz.mean())/30),'omega_wobble_rms_rad_s':float(np.sqrt(np.mean(q.omega_wobble_rad_s**2))),'omega_wobble_peak_rad_s':float(q.omega_wobble_rad_s.max()),'Tmag_rms_Nmm':float(np.sqrt(np.mean(q.Tmag_norm_Nmm**2))),'exact_gap_min_um':float(c.exact_gap_mm.min()*1000),'CPRESS_max_MPa':float(c.CPRESS_max_MPa.max()),'Wmag_end_J':float(q.Wmag_J.iloc[-1]),'Whydro_end_J':float(q.Whydro_J.iloc[-1]),'Wnonmag_end_J':float(q.Wnonmag_resolved_J.iloc[-1])}
first=float(contact[contact.CPRESS_max_MPa>1e-9].time_s.iloc[0]) if (contact.CPRESS_max_MPa>1e-9).any() else .003
windows=[('W0_preimpact',0,first),('W1_contact',first,min(first+.0005,t[-1])),('W2_immediate_postimpact',min(first+.0005,t[-1]),min(first+.002,t[-1])),('W3_intermediate',min(first+.002,t[-1]),min(first+.004,t[-1])),('W4_late_postimpact',max(first+.004,t[0]),t[-1])]
summary=pd.DataFrame([dict(window=n,**stats(lo,hi)) for n,lo,hi in windows]); summary.to_csv(OUT/'reduced_hydro_window_summary.csv',index=False)

def plot(name,ys,labs):
    plt.figure(figsize=(8,4))
    for yv,lab in zip(ys,labs):plt.plot(t*1000,yv,label=lab)
    plt.xlabel('time (ms)');plt.legend();plt.tight_layout();plt.savefig(OUT/name,dpi=180);plt.close()
plot('wobble_rate_reduced.png',[data.robot_phase_rate_Hz,data.omega_wobble_rad_s],['robot phase rate Hz','|omega wobble| rad/s'])
plot('H_and_impulse_reduced.png',[np.linalg.norm(H,axis=1),np.linalg.norm(J,axis=1)],['|H| kg m2/s','|Jmag| Nmm s'])
plot('rotational_energy_work_reduced.png',[K,W,Wh],['Krot J','Wmag J','Whydro J'])
plt.figure(figsize=(8,4)); plt.plot(contact.time_s*1000,contact.CPRESS_max_MPa,label='aggregate contact CPRESS'); ax=plt.twinx(); ax.plot(contact.time_s*1000,contact.exact_gap_mm*1000,'g',label='exact gap um'); plt.xlabel('time (ms)');plt.tight_layout();plt.savefig(OUT/'contact_gap_wobble_reduced.png',dpi=180);plt.close()

late=summary.iloc[-1]; pre=summary.iloc[0]; tr=float(late.Tmag_rms_Nmm/max(pre.Tmag_rms_Nmm,1e-30)); rr=float(late.R_phase_vs_30) if 'R_phase_vs_30' in late else np.nan
if len(late)<30 or not np.isfinite(rr): cls='E_POST_IMPACT_MECHANISM_STILL_AMBIGUOUS'
elif rr>=.7: cls='A_POST_IMPACT_WOBBLE_RECOVERS'
elif rr<=.3 and tr>=.8: cls='B_SUSTAINED_POST_IMPACT_NONMAGNETIC_ROTATIONAL_STALL'
elif tr<.5: cls='D_MAGNETIC_TORQUE_COLLAPSE'
else: cls='C_PARTIAL_POST_IMPACT_WOBBLE_RECOVERY'
(OUT/'visualization_frame_audit.txt').write_text('production frame: transformed curvenew_CEL_xyrot56_exact.csv; no flattened Magpylib CSV used\n',encoding='utf-8')
(OUT/'Wobble30Hz_reduced_hydrodynamics_report.md').write_text(f'''# Reduced-Hydro Wall-ON 8.333 ms report\n\n## Decision\n`{cls}`\n\n## Runtime identity\nOne job completed at 8.333 ms with 83,330 direct increments (`dt=1e-7 s`) and 168 ODB frames (50 us). CEL fluid and CEL contact were removed; SmoothWall114 Robot-wall contact remains. External production Magpylib Socket: B0=10 mT, f=30 Hz, cone=30 deg, Bias=40 deg, phase=248 deg, sense=+1, G=6 mT=0.006 T, L=45 mm, FORWARD.\n\n## Hydro\nBody-following dissipative RP loads use provisional analytical coefficients cpar=4e-9, cperp=1.2e-8 N s/mm, kspin=1e-9, kwob=3e-9 N mm s. Force/torque powers are in `reduced_hydro_angular_dynamics.csv` and should be non-positive.\n\n## Strict dynamics\nFormal omega is finite-rotation log of consecutive orientation matrices; direct VR is independently extracted. Maximum absolute reconstruction difference is {omega_diff.max():.6g} rad/s. H is kg m2/s; Jmag is N mm s (multiply by 1e-3 for N m s); Wmag uses 1 N mm=1e-3 J. Resolved nonmagnetic terms are accounting quantities, not independent wall torque measurements.\n\n## Windows\n\n{summary.to_markdown(index=False)}\n\n## Contact and translation\nExact gap is computed against authoritative Pipe_WALL_HELPER R3D4 triangles. Minimum signed gap is {contact.exact_gap_mm.min()*1000:.3f} um; aggregate contact CPRESS peak is {contact.CPRESS_max_MPa.max():.3f} MPa. Canonical Vt and displacement are in `reduced_hydro_translation.csv`.\n\n## Interpretation\nThis job is a diagnostic reduced-hydrodynamics replacement for coarse CEL, not a calibrated final fluid model. Any difference from prior CEL is classified only after checking geometry/contact, field identity, direct-VR/finite-rotation consistency, and non-positive hydro power.\n\n## Exactly one next step\nClassification `{cls}`. If A, extend same configuration for stability; if B, do the approved zero-cost elliptical-polarization/contact-load-shaping audit; if C, extend same configuration to 16.667 ms; if D, revisit field amplitude only after output validation; if E, repair extraction.\n''',encoding='utf-8')
print('CLASSIFICATION',cls); print(summary.to_string(index=False)); print('min_gap_um',contact.exact_gap_mm.min()*1000,'CPRESS_MPa',contact.CPRESS_max_MPa.max(),'hydro_power_max',max(Ph.max(),Pht.max()))
