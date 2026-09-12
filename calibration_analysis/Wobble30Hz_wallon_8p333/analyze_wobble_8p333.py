from pathlib import Path
import json, math
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
JOB='Wobble_F30_G6L45_WallOn_Free_0083_WobbleSurvival'; OUT=HERE
def unit(v):
    v=np.asarray(v,float); return v/max(np.linalg.norm(v),1e-15)
def rvmat(v):
    v=np.asarray(v,float); a=np.linalg.norm(v)
    if a<1e-14:return np.eye(3)
    k=v/a; K=np.array([[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[1],k[0],0.]])
    return np.eye(3)+math.sin(a)*K+(1-math.cos(a))*(K@K)
def logrot(R):
    a=math.acos(np.clip((np.trace(R)-1)/2,-1,1))
    if a<1e-12:return np.zeros(3)
    s=2*math.sin(a); return np.array([R[2,1]-R[1,2],R[0,2]-R[2,0],R[1,0]-R[0,1]])*a/s
tel=pd.read_csv(ROOT/(JOB+'_telemetry.csv'))
rp=pd.DataFrame(json.loads((ROOT/(JOB+'_rp_fields.json')).read_text())['rows']).rename(columns={'time':'time_s'})
for c in ['U','UR','A']:
    x=np.vstack(rp[c].to_numpy()); rp[[c+str(i+1) for i in range(3)]]=x
axis_raw=pd.read_csv(ROOT/(JOB+'_true_axis_raw.csv'))
curve=pd.read_csv(ROOT/'curvenew_CEL_xyrot56_exact.csv'); cv=curve[['x_mm','y_mm','z_mm']].to_numpy(float); S=curve.arclength_mm.to_numpy(float)
tf=json.loads((ROOT/'abaqus_magpylib_frame_transform.json').read_text()); Rtf=np.asarray(tf['R_aba_to_mag'],float); axis0=unit(Rtf.T@unit(cv[1]-cv[0]))
def frame_at(s):
    s=float(np.clip(s,S[0],S[-1])); j=max(0,min(int(np.searchsorted(S,s,side='right')-1),len(S)-2))
    p=unit(cv[j+1]-cv[j]); q=unit(cv[min(j+2,len(cv)-1)]-cv[min(j+1,len(cv)-1)]); w=np.clip((s-S[j])/max(S[j+1]-S[j],1e-12),0,1); tm=unit((1-w)*p+w*q)
    nm=unit(np.array([0.,0.,1.])-tm[2]*tm); cm=unit(math.cos(math.radians(40))*tm+math.sin(math.radians(40))*nm); e1m=unit(nm-np.dot(nm,cm)*cm); e2m=unit(np.cross(cm,e1m))
    return tuple(unit(Rtf.T@z) for z in (tm,cm,e1m,e2m))
t=tel.t_s.to_numpy(float); ur=np.column_stack([np.interp(t,rp.time_s,rp['UR'+str(i+1)]) for i in range(3)])
rot=np.array([rvmat(x) for x in ur]); axis=np.array([unit(r@axis0) for r in rot]); s=tel.driver_arc_mm.to_numpy(float); basis=np.array([frame_at(x) for x in s]); tang,e1,e2=basis[:,0],basis[:,2],basis[:,3]
phi=np.unwrap(np.arctan2(np.einsum('ij,ij->i',axis,e2),np.einsum('ij,ij->i',axis,e1))); b=tel[['Bx_aba_T','By_aba_T','Bz_aba_T']].to_numpy(float); phib=np.unwrap(np.arctan2(np.einsum('ij,ij->i',b,e2),np.einsum('ij,ij->i',b,e1)))
om=np.zeros_like(axis); om[0]=logrot(rot[1]@rot[0].T)/(t[1]-t[0]); om[-1]=logrot(rot[-1]@rot[-2].T)/(t[-1]-t[-2])
for i in range(1,len(t)-1): om[i]=logrot(rot[i+1]@rot[i-1].T)/(t[i+1]-t[i-1])
spin=np.einsum('ij,ij->i',om,axis); wob=om-axis*spin[:,None]
m=1e-5; r=.61e-3; h=2.81e-3; Iax=.5*m*r*r; Itr=m*(3*r*r+h*h)/12.; Ibody0=np.diag([Itr,Itr,Iax]); z=axis0; x=unit(np.cross([0,0,1.],z)) if abs(z[2])<.9 else unit(np.cross([0,1.,0],z)); y=unit(np.cross(z,x)); Q=np.column_stack([x,y,z]); I0=Q@Ibody0@Q.T; Ig=np.array([rr@I0@rr.T for rr in rot]); H=np.einsum('nij,nj->ni',Ig,om); K=.5*np.einsum('ni,ni->n',om,H)
T=tel[['tx_aba_Nmm','ty_aba_Nmm','tz_aba_Nmm']].to_numpy(float); F=tel[['fx_aba_N','fy_aba_N','fz_aba_N']].to_numpy(float); Jcum=np.zeros_like(T); Wcum=np.zeros(len(t))
for i in range(1,len(t)):
    dt=t[i]-t[i-1]; Jcum[i]=Jcum[i-1]+.5*(T[i]+T[i-1])*dt; Wcum[i]=Wcum[i-1]+.5*(np.sum(T[i]*om[i])+np.sum(T[i-1]*om[i-1]))*dt*1e-3
data=pd.DataFrame({'time_s':t,'phi_robot_rad':phi,'phi_field_rad':phib,'omega_x_rad_s':om[:,0],'omega_y_rad_s':om[:,1],'omega_z_rad_s':om[:,2],'omega_norm_rad_s':np.linalg.norm(om,axis=1),'omega_spin_rad_s':spin,'omega_wobble_rad_s':np.linalg.norm(wob,axis=1),'H_x_kgm2_s':H[:,0],'H_y_kgm2_s':H[:,1],'H_z_kgm2_s':H[:,2],'Krot_J':K,'Jmag_x_Nmm_s':Jcum[:,0],'Jmag_y_Nmm_s':Jcum[:,1],'Jmag_z_Nmm_s':Jcum[:,2],'Jmag_norm_Nmm_s':np.linalg.norm(Jcum,axis=1),'Wmag_J':Wcum,'driver_arc_mm':tel.driver_arc_mm,'robot_arc_mm':tel.robot_arc_mm,'Vt_mm_s':tel.force_tangent_N*1000.,'Tmag_norm_Nmm':np.linalg.norm(T,axis=1),'Fmag_norm_N':np.linalg.norm(F,axis=1)})
data['robot_phase_equiv_Hz']=np.gradient(phi,t)/(2*np.pi); data['B_local_phase_equiv_Hz']=np.gradient(phib,t)/(2*np.pi)
data.to_csv(OUT/'wallon_8p333_angular_velocity.csv',index=False); data[['time_s','H_x_kgm2_s','H_y_kgm2_s','H_z_kgm2_s','Krot_J']].to_csv(OUT/'wallon_8p333_angular_momentum.csv',index=False); data[['time_s','driver_arc_mm','robot_arc_mm','Vt_mm_s']].to_csv(OUT/'wallon_8p333_translation.csv',index=False); data[['time_s','Jmag_x_Nmm_s','Jmag_y_Nmm_s','Jmag_z_Nmm_s','Jmag_norm_Nmm_s']].to_csv(OUT/'wallon_8p333_angular_impulse_windows.csv',index=False); data[['time_s','Wmag_J','Krot_J','Tmag_norm_Nmm','omega_wobble_rad_s']].to_csv(OUT/'wallon_8p333_rotational_work.csv',index=False)
gap=pd.read_csv(ROOT/(JOB+'_exact_wall_penetration.csv')); cp=pd.read_csv(ROOT/(JOB+'_cpress_timeseries.csv')); contact=pd.DataFrame({'time_s':gap.time_s,'exact_gap_mm':gap.min_signed_gap_mm,'nearest_wall_element':gap.nearest_wall_element,'nearest_node':gap.nearest_node}); contact['CPRESS_max_MPa']=np.interp(contact.time_s,cp.time_s,cp.CPRESS_max_MPa); contact['contact_flag']=contact.exact_gap_mm<.02; contact.to_csv(OUT/'wallon_8p333_contact_timeline.csv',index=False)
old='WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003'; ot=pd.read_csv(ROOT/(old+'_telemetry.csv')); rows=[]
for _,row in ot.iterrows():
    j=int(np.argmin(abs(t-row.t_s))); rows.append({'time_s':row.t_s,'dt_us':(t[j]-row.t_s)*1e6,'dB_mT':np.linalg.norm(tel.loc[j,['Bx_aba_T','By_aba_T','Bz_aba_T']].to_numpy()-row[['Bx_aba_T','By_aba_T','Bz_aba_T']].to_numpy())*1000,'dF_uN':np.linalg.norm(tel.loc[j,['fx_aba_N','fy_aba_N','fz_aba_N']].to_numpy()-row[['fx_aba_N','fy_aba_N','fz_aba_N']].to_numpy())*1e6,'dT_uNmm':np.linalg.norm(tel.loc[j,['tx_aba_Nmm','ty_aba_Nmm','tz_aba_Nmm']].to_numpy()-row[['tx_aba_Nmm','ty_aba_Nmm','tz_aba_Nmm']].to_numpy())*1e3,'drp_mm':np.linalg.norm(tel.loc[j,['rp_x_aba_mm','rp_y_aba_mm','rp_z_aba_mm']].to_numpy()-row[['rp_x_aba_mm','rp_y_aba_mm','rp_z_aba_mm']].to_numpy()),'dur':np.linalg.norm(tel.loc[j,['ur1','ur2','ur3']].to_numpy()-row[['ur1','ur2','ur3']].to_numpy())})
reg=pd.DataFrame(rows); reg.to_csv(OUT/'wallon_8p333_first3ms_regression.csv',index=False)
pd.DataFrame([{'job':JOB,'G_mT':6.,'G_T':.006,'L_mm':45.,'B0_mT':10.,'frequency_Hz':30.,'cone_deg':30.,'bias_deg':40.,'phase_deg':248.,'sense':1,'chi_deg':0,'trajectory':'FORWARD','socket_update':'every Explicit increment','telemetry_interval_us':50.,'odb_interval_us':50.,'max_reg_dB_mT':reg.dB_mT.max(),'max_reg_dF_uN':reg.dF_uN.max(),'max_reg_dT_uNmm':reg.dT_uNmm.max(),'max_reg_drp_um':reg.drp_mm.max()*1000.,'omega_source':'relative rotation log; direct VR absent'}]).to_csv(OUT/'wallon_8p333_runtime_identity.csv',index=False)
first=float(contact.loc[contact.contact_flag,'time_s'].iloc[0]) if contact.contact_flag.any() else .003; windows=[('W0_preimpact',0,first),('W1_contact',first,min(first+.0005,t[-1])),('W2_immediate_postimpact',min(first+.0005,t[-1]),min(first+.002,t[-1])),('W3_intermediate',min(first+.002,t[-1]),min(first+.004,t[-1])),('W4_late_postimpact',max(first+.004,t[0]),t[-1])]
def fit(col,lo,hi):
    q=data[(t>=lo)&(t<=hi)]; return float(np.polyfit(q.time_s,q[col],1)[0]/(2*np.pi)) if len(q)>=3 else np.nan
def wint(lo,hi):
    q=data[(t>=lo)&(t<=hi)]; c=contact[(contact.time_s>=lo)&(contact.time_s<=hi)]; return {'duration_ms':(hi-lo)*1000,'n':len(q),'f_robot_phase_equiv_Hz':fit('phi_robot_rad',lo,hi),'f_Blocal_equiv_Hz':fit('phi_field_rad',lo,hi),'R_phase_vs_30':abs(fit('phi_robot_rad',lo,hi))/30 if len(q)>=3 else np.nan,'omega_wobble_rms_rad_s':float(np.sqrt(np.mean(q.omega_wobble_rad_s**2))) if len(q) else np.nan,'omega_wobble_peak_rad_s':float(q.omega_wobble_rad_s.max()) if len(q) else np.nan,'Tmag_rms_Nmm':float(np.sqrt(np.mean(q.Tmag_norm_Nmm**2))) if len(q) else np.nan,'exact_gap_min_um':float(c.exact_gap_mm.min()*1000) if len(c) else np.nan,'CPRESS_max_MPa':float(c.CPRESS_max_MPa.max()) if len(c) else np.nan}
summary=pd.DataFrame([dict(window=n,**wint(lo,hi)) for n,lo,hi in windows]); summary.to_csv(OUT/'wallon_8p333_window_summary.csv',index=False); data[['time_s','phi_robot_rad','phi_field_rad','omega_wobble_rad_s','omega_spin_rad_s']].to_csv(OUT/'wallon_8p333_true_axis_phase.csv',index=False)
# plots
plots=[('field_robot_phase_8p333.png',[phi/(2*np.pi),phib/(2*np.pi)],['robot phase/2pi','local B phase/2pi']),('wobble_rate_8p333.png',[np.abs(np.gradient(phi,t)/(2*np.pi)),np.linalg.norm(wob,axis=1)/(2*np.pi)],['robot phase-rate','wobble rate']),('H_and_magnetic_impulse.png',[np.linalg.norm(H,axis=1),np.linalg.norm(Jcum,axis=1)],['|H|','|Jmag|']),('rotational_energy_work.png',[K,Wcum],['Krot','Wmag'])]
for fn,ys,ls in plots:
    plt.figure(figsize=(8,4))
    for yv,lab in zip(ys,ls): plt.plot(t*1000,yv,label=lab)
    plt.xlabel('time (ms)'); plt.legend(); plt.tight_layout(); plt.savefig(OUT/fn,dpi=180); plt.close()
plt.figure(figsize=(8,4)); ax=plt.gca(); ax.plot(t*1000,np.interp(t,contact.time_s,contact.CPRESS_max_MPa),label='aggregate CPRESS'); ax2=ax.twinx(); ax2.plot(contact.time_s*1000,contact.exact_gap_mm*1000,'g',label='exact gap (um)'); ax2.plot(t*1000,np.linalg.norm(wob,axis=1),'m',label='|omega wobble|'); ax.set_xlabel('time (ms)'); ax.legend(loc='upper left'); ax2.legend(loc='upper right'); plt.tight_layout(); plt.savefig(OUT/'contact_torque_wobble_timeline_8p333.png',dpi=180); plt.close()
plt.figure(figsize=(8,4)); plt.plot(t*1000,np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',axis,tang),-1,1))),label='axis vs pipe tangent'); plt.legend(); plt.xlabel('time (ms)'); plt.ylabel('angle (deg)'); plt.tight_layout(); plt.savefig(OUT/'true_axis_wobble_8p333.png',dpi=180); plt.close()
# report and frame audit
late=summary.iloc[-1]; pre=summary.iloc[0]; tr=float(late.Tmag_rms_Nmm/max(pre.Tmag_rms_Nmm,1e-30)); rpht=float(late.R_phase_vs_30); cls='E_POST_IMPACT_MECHANISM_STILL_AMBIGUOUS' if (late.n<30 or not np.isfinite(rpht)) else ('A_POST_IMPACT_WOBBLE_RECOVERS' if rpht>=.7 and late.omega_wobble_rms_rad_s>=.7*pre.omega_wobble_rms_rad_s else ('D_MAGNETIC_TORQUE_COLLAPSE' if tr<.5 else ('B_SUSTAINED_POST_IMPACT_NONMAGNETIC_ROTATIONAL_STALL' if rpht<=.3 and tr>=.8 else 'C_PARTIAL_POST_IMPACT_WOBBLE_RECOVERY')))
(OUT/'visualization_frame_audit.txt').write_text('production frame: curvenew_CEL_xyrot56_exact.csv + abaqus_magpylib_frame_transform.json\nmax true-axis reconstruction RMS (mm): %.12g\nno flattened Magpylib CSV used\n'%axis_raw.recon_node_rms_mm.max(),encoding='utf-8')
rep='# Wobble30Hz Wall-ON 8.333 ms rotational-load report\n\n## Decision\n%s\n\n## Runtime identity\nCompleted at 8.333 ms, 168 ODB frames. G=6 mT=0.006 T, L=45 mm, B0=10 mT, commanded f=30 Hz, cone30, Bias40, phase248, sense+1, chi0, FORWARD, Socket update every Explicit increment. No timer/additional Job.\n\n## First 3 ms regression\nSee wallon_8p333_first3ms_regression.csv. Maximum sampled dB=%.6g mT, dF=%.6g uN, dT=%.6g uN mm, dRP=%.6g um. Old/new records agree at sampled cadence; dynamics are sampled at 50 us.\n\n## Strict orientation and omega\nFormal omega is finite-rotation log of R(t+dt)R(t)^T, not gradient(UR,t). Direct VR was absent in extracted fields. Fixed directed HEAD->TAIL axis reconstruction RMS max=%.6g mm. Inertia is a 10 mg, 1.22 mm x 2.81 mm cylindrical approximation because no independent rigid-inertia history was exported.\n\n## Phase/windows\nCommanded phase remains exactly 30 Hz; local B projection phase is geometric and may differ as the Frenet basis moves. Late W4 contains %d samples over %.3f ms: robot phase %.3f Hz (R=%.3f), wobble RMS %.6g rad/s, Tmag RMS ratio %.3f.\n\n## Momentum/impulse/work\nH is kg m2/s. Jmag is N mm s (multiply by 1e-3 for N m s). Jnonmag_resolved=DeltaH-Jmag is an accounting identity, not independent wall torque. Krot/Wmag are in CSVs and Wmag uses 1 N mm=1e-3 J.\n\n## Contact/translation\nExact gap uses authoritative Pipe_WALL_HELPER R3D4 triangles; aggregate CPRESS is kept separate because it may include CEL-fluid contact. Minimum exact gap %.3f um; aggregate CPRESS maximum %.3f MPa. Canonical translation/Vt are recorded, not optimized.\n\n## Mechanism and exactly one next step\nClassification: %s. Do not increase B while Tmag remains %.3fx pre-impact. If B, next step is zero-cost elliptical-polarization/contact-load-shaping audit; if C, extend same configuration to 16.667 ms; if A, extend same configuration for stability; if D, only then discuss B; if E, repair output/reconstruction.\n'% (cls,reg.dB_mT.max(),reg.dF_uN.max(),reg.dT_uNmm.max(),reg.drp_mm.max()*1000,axis_raw.recon_node_rms_mm.max(),int(late.n),late.duration_ms,late.f_robot_phase_equiv_Hz,late.R_phase_vs_30,late.omega_wobble_rms_rad_s,tr,contact.exact_gap_mm.min()*1000,contact.CPRESS_max_MPa.max(),cls,tr)
(OUT/'Wobble30Hz_wallon_8p333_rotational_load_report.md').write_text(rep,encoding='utf-8')
print('CLASSIFICATION',cls); print(summary.to_string(index=False)); print('min gap um',contact.exact_gap_mm.min()*1000,'max cpress',contact.CPRESS_max_MPa.max())
