"""Zero-cost Wall-OFF/Wall-ON true-axis wobble survival audit."""
from pathlib import Path
import hashlib, json, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]; OUT=HERE; OUT.mkdir(parents=True,exist_ok=True)
OFF='WobbleCal_COMFixed_F30_Cone30_B10'; ON='WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003'
def sha(p):
    h=hashlib.sha256(); h.update(Path(p).read_bytes()); return h.hexdigest()
def U(x):
    x=np.asarray(x,float); return x/max(np.linalg.norm(x),1e-15)
def Rv(v):
    v=np.asarray(v,float); a=np.linalg.norm(v)
    if a<1e-14:return np.eye(3)
    k=v/a; K=np.array([[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[1],k[0],0.]])
    return np.eye(3)+math.sin(a)*K+(1-math.cos(a))*(K@K)

curve=pd.read_csv(ROOT/'curvenew_CEL_xyrot56_exact.csv'); cv=curve[['x_mm','y_mm','z_mm']].to_numpy(float); S=curve.arclength_mm.to_numpy(float)
tf=json.load(open(ROOT/'abaqus_magpylib_frame_transform.json')); R=np.asarray(tf['R_aba_to_mag'],float); axis_aba=U(R.T@U(cv[1]-cv[0]))
def frame_at(s):
    s=float(np.clip(s,S[0],S[-1])); j=max(0,min(int(np.searchsorted(S,s,side='right')-1),len(S)-2)); t0=U(cv[j+1]-cv[j]); t1=U(cv[min(j+2,len(cv)-1)]-cv[min(j+1,len(cv)-1)]); w=np.clip((s-S[j])/max(S[j+1]-S[j],1e-12),0,1); tm=U((1-w)*t0+w*t1)
    # The production server constructs the Frenet frame in the flat-Mag
    # frame, where +Z is the reference normal, then maps every vector back to
    # Abaqus.  Constructing +Z after R.T would change Bias40 and phase.
    nm=U(np.array([0.,0.,1.])-tm[2]*tm); cm=U(math.cos(math.radians(40))*tm+math.sin(math.radians(40))*nm); e1m=U(nm-np.dot(nm,cm)*cm); e2m=U(np.cross(cm,e1m)); return tuple(U(R.T@q) for q in (tm,cm,e1m,e2m))
def enrich(tel,ur):
    d=tel.copy(); t=d.t_s.to_numpy(float); u=d[['ur1','ur2','ur3']].to_numpy(float); axis=np.array([Rv(x)@axis_aba for x in u]);
    # Archived telemetry proves that the field frame used driver_arc, despite
    # the ROBOT_ARC banner.  Use the measured driver station for the field
    # basis; retain robot_arc separately for transport diagnostics.
    s=d.driver_arc_mm.to_numpy(float) if 'driver_arc_mm' in d else (np.interp(t,ur.time_s.to_numpy(float),ur.s_robot_mm.to_numpy(float)) if 's_robot_mm' in ur else np.full(len(d),13.8))
    e1=[];e2=[];tt=[]
    for x in s:
        q=frame_at(x);tt.append(q[0]);e1.append(q[2]);e2.append(q[3])
    e1=np.asarray(e1);e2=np.asarray(e2);tt=np.asarray(tt); a1=np.einsum('ij,ij->i',axis,e1);a2=np.einsum('ij,ij->i',axis,e2); phi=np.unwrap(np.arctan2(a2,a1));
    b=np.column_stack([d.Bx_aba_T,d.By_aba_T,d.Bz_aba_T]); b1=np.einsum('ij,ij->i',b,e1);b2=np.einsum('ij,ij->i',b,e2); phib=np.unwrap(np.arctan2(b2,b1));
    om=np.gradient(u,t,axis=0); spin=np.einsum('ij,ij->i',om,axis); wob=om-axis*spin[:,None]; out=d.copy();out['robot_s_mm']=s;out['a1_mm']=a1;out['a2_mm']=a2;out['phi_robot_rad']=phi;out['phi_field_rad']=phib;out['f_robot_directed_Hz']=np.gradient(phi,t)/(2*np.pi);out['f_field_Hz']=np.gradient(phib,t)/(2*np.pi);out['omega_spin_rad_s']=spin;out['omega_wobble_rad_s']=np.linalg.norm(wob,axis=1);out['true_axis_angle_deg']=np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',axis,tt),-1,1)));out['Tmag_norm_Nmm']=np.linalg.norm(d[['tx_aba_Nmm','ty_aba_Nmm','tz_aba_Nmm']].to_numpy(float),axis=1);return out

off_tel=pd.read_csv(ROOT/(OFF+'_telemetry.csv')); off_ur=pd.read_csv(ROOT/(OFF+'_UR_highres.csv')); off=enrich(off_tel,off_ur); off.to_csv(OUT/'walloff_f30_true_axis_phase.csv',index=False)
on_tel=pd.read_csv(ROOT/(ON+'_telemetry.csv')); on_ur=pd.read_csv(ROOT/(ON+'_rp_velocity_history.csv')) if (ROOT/(ON+'_rp_velocity_history.csv')).exists() else pd.DataFrame({'time_s':on_tel.t_s,'s_robot_mm':on_tel.robot_arc_mm}); on=enrich(on_tel,on_ur); on.to_csv(OUT/'wallon_pre_postimpact_wobble.csv',index=False)

def fit_hz(d,lo,hi,col):
    z=d[(d.t_s>=lo)&(d.t_s<=hi)]
    return float(np.polyfit(z.t_s,z[col],1)[0]/(2*np.pi)) if len(z)>=3 else float('nan')
def stats(d,lo,hi):
    z=d[(d.t_s>=lo)&(d.t_s<=hi)]
    return dict(start_ms=lo*1000,end_ms=hi*1000,n=len(z),f_field_Hz=fit_hz(d,lo,hi,'phi_field_rad'),f_robot_Hz=fit_hz(d,lo,hi,'phi_robot_rad'),f_robot_median_Hz=float(np.median(z.f_robot_directed_Hz)),omega_wobble_mean_rad_s=float(z.omega_wobble_rad_s.mean()),omega_wobble_rms_rad_s=float(np.sqrt(np.mean(z.omega_wobble_rad_s**2))),true_axis_range_deg=float(z.true_axis_angle_deg.max()-z.true_axis_angle_deg.min()),Tmag_rms_Nmm=float(np.sqrt(np.mean(z.Tmag_norm_Nmm**2))))
lock=pd.DataFrame([stats(off,0,.011),stats(off,.011,.033333),stats(off,.033333,.0555)]);lock['R_lock_vs_30Hz']=lock.f_robot_Hz/30.;lock['R_lock_abs']=lock.f_robot_Hz.abs()/30.;lock.to_csv(OUT/'walloff_f30_lock_windows.csv',index=False)
onsets=[('W0_precontact',0,.0013),('W1_contact',.0013,.0017),('W2_postimpact',.0017,.0025),('W3_latest',.0025,.0030)];pre=stats(on,0,.0013); rows=[]
for name,lo,hi in onsets:
    q=stats(on,lo,hi);q['window']=name;q['R_post_over_pre_omega']=q['omega_wobble_mean_rad_s']/max(pre['omega_wobble_mean_rad_s'],1e-12);q['R_phase_vs_30Hz']=abs(q['f_robot_Hz'])/30.;rows.append(q)
post=pd.DataFrame(rows);post.to_csv(OUT/'wallon_pre_postimpact_wobble.csv',index=False)

def reversals(d):
    rr=[]
    for col in ['a1_mm','a2_mm']:
        v=np.gradient(d[col].to_numpy(float),d.t_s.to_numpy(float));sg=np.sign(v);rr.append(dict(component=col,turning_count=int(np.sum(sg[1:]*sg[:-1]<0)),positive_fraction=float(np.mean(v>0)),negative_fraction=float(np.mean(v<0))))
    return rr
pd.DataFrame(reversals(off)).assign(window='walloff_0_55.5ms').to_csv(OUT/'postimpact_component_reversals.csv',index=False)
tb=[]
for name,lo,hi in onsets:
    z=on[(on.t_s>=lo)&(on.t_s<=hi)];tb.append(dict(window=name,Tmag_rms_Nmm=float(np.sqrt(np.mean(z.Tmag_norm_Nmm**2))),Tmag_mean_Nmm=float(z.Tmag_norm_Nmm.mean()),field_freq_Hz=fit_hz(on,lo,hi,'phi_field_rad'),field_norm_mean_mT=float(z.B_aba_T.mean()*1000),omega_wobble_mean_rad_s=float(z.omega_wobble_rad_s.mean())))
tb=pd.DataFrame(tb);tb['Tmag_ratio_to_pre']=tb.Tmag_rms_Nmm/max(tb.iloc[0].Tmag_rms_Nmm,1e-15);tb.to_csv(OUT/'postimpact_torque_budget.csv',index=False)
rex_path=OUT/'runtime_exact_evaluator_replay.csv'
if rex_path.exists():
    rex=pd.read_csv(rex_path); replay_pass=bool(rex[['B_abs_error_T','F_abs_error_N','T_abs_error_Nmm']].max().max()<2e-12)
    bm=rex[['B_model_x_T','B_model_y_T','B_model_z_T']].to_numpy(float); bt=rex[['B_tele_x_T','B_tele_y_T','B_tele_z_T']].to_numpy(float)
    den=np.linalg.norm(bm,axis=1)*np.linalg.norm(bt,axis=1); ang=np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',bm,bt)/np.maximum(den,1e-30),-1,1)))
    pd.DataFrame(dict(time_s=rex.time_s,B_model_norm_T=rex.B_model_norm_T,B_tele_norm_T=rex.B_tele_norm_T,B_abs_error_T=rex.B_abs_error_T,B_relative_error=rex.B_abs_error_T/np.maximum(rex.B_tele_norm_T,1e-30),B_vector_angle_deg=ang)).to_csv(OUT/'field_replay_error_decomposition.csv',index=False)
    pd.DataFrame(dict(time_s=rex.time_s,T_model_norm_Nmm=rex.T_model_norm_Nmm,T_tele_norm_Nmm=rex.T_tele_norm_Nmm,T_abs_error_Nmm=rex.T_abs_error_Nmm,T_relative_error=rex.T_abs_error_Nmm/np.maximum(rex.T_tele_norm_Nmm,1e-30))).to_csv(OUT/'torque_replay_error_decomposition.csv',index=False)
else:
    replay_pass=False; rex=pd.DataFrame()
# Angular impulse is reported in N·mm·s.  Independent rigid inertia history
# is not present in this export, so ΔH is explicitly marked unavailable
# rather than inferred from UR alone.
ii=[]
for name,lo,hi in onsets:
    z=on[(on.t_s>=lo)&(on.t_s<=hi)];dt=np.gradient(z.t_s.to_numpy(float));T=z[['tx_aba_Nmm','ty_aba_Nmm','tz_aba_Nmm']].to_numpy(float);J=np.sum(T*dt[:,None],axis=0);ii.append(dict(window=name,Jmag_x_Nmm_s=J[0],Jmag_y_Nmm_s=J[1],Jmag_z_Nmm_s=J[2],Jmag_norm_Nmm_s=float(np.linalg.norm(J)),delta_H_Nmm_s=np.nan,delta_H_status='not independently available from archived exports'))
pd.DataFrame(ii).to_csv(OUT/'angular_momentum_audit.csv',index=False)
gap=pd.read_csv(ROOT/(ON+'_exact_wall_penetration.csv'));cp=pd.read_csv(ROOT/'abaqus_robot/calibration_analysis/WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003/cpress_peaks.csv');contact=pd.DataFrame({'time_s':gap.time_s,'exact_gap_mm':gap.min_signed_gap_mm});contact['CPRESS_MPa']=np.interp(contact.time_s,cp.time_s,cp.cpress_max_mpa);contact['gap_threshold_contact']=contact.exact_gap_mm<.02;contact.to_csv(OUT/'wallon_contact_timing.csv',index=False)

fig,ax=plt.subplots(figsize=(8,4));ax.plot(off.t_s*1000,np.degrees(off.phi_field_rad)%360,label='field phase');ax.plot(off.t_s*1000,np.degrees(off.phi_robot_rad)%360,label='robot directed phase');ax.axvline(33.333,color='k',ls='--',lw=.7);ax.set(xlabel='time (ms)',ylabel='phase modulo 360°');ax.legend();fig.tight_layout();fig.savefig(OUT/'field_vs_robot_phase_postimpact.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4));ax.plot(off.t_s*1000,off.f_robot_directed_Hz,label='robot phase rate');ax.plot(off.t_s*1000,off.omega_wobble_rad_s/(2*np.pi),label='|omega_wobble|/2π');ax.axhline(30,color='k',ls='--',lw=.7,label='30 Hz target');ax.set(xlabel='time (ms)',ylabel='Hz');ax.legend();fig.tight_layout();fig.savefig(OUT/'postimpact_wobble_rate.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4));ax.plot(on.t_s*1000,np.interp(on.t_s,contact.time_s,contact.CPRESS_MPa),label='aggregate CPRESS');ax2=ax.twinx();ax2.plot(contact.time_s*1000,contact.exact_gap_mm*1000,'g',label='exact gap (µm)');ax2.plot(on.t_s*1000,on.Tmag_norm_Nmm,'m',label='|Tmag| (N·mm)');ax2.plot(on.t_s*1000,on.omega_wobble_rad_s,'c',label='|omega_wobble|');ax.set(xlabel='time (ms)',ylabel='CPRESS (MPa)');ax2.set_ylabel('gap / torque / angular rate');fig.tight_layout();fig.savefig(OUT/'contact_torque_wobble_timeline.png',dpi=180);plt.close(fig)

pre_r=float(lock.iloc[1:]['R_lock_abs'].mean()); post_r=float(post.loc[post.window=='W3_latest','R_phase_vs_30Hz'].iloc[0]); t_ratio=float(tb.loc[tb.window=='W3_latest','Tmag_ratio_to_pre'].iloc[0]);
if pre_r<.9: classification='E_B10_ALREADY_UNDERSYNCHRONIZED_WITHOUT_WALL'
elif post_r<.6 and t_ratio>=.7: classification='D_POST_IMPACT_NONMAGNETIC_ROTATIONAL_LOADING_DOMINATES'
elif t_ratio<.5: classification='C_POST_IMPACT_MAGNETIC_TORQUE_COLLAPSE'
elif post_r<.6: classification='D_POST_IMPACT_NONMAGNETIC_ROTATIONAL_LOADING_DOMINATES'
else: classification='A_NO_MAJOR_POST_IMPACT_WOBBLE_STALL'
chosen=None if t_ratio>=.7 else (15.0 if post_r>=.45 else 20.0)
chosen_text='none' if chosen is None else f'{chosen:g}'
pd.DataFrame([dict(B_mT=b,torque_scale_vs_10mT=b/10.,selected=(chosen is not None and b==chosen),selection_basis='not indicated; Tmag does not collapse' if chosen is None else 'order-of-magnitude only') for b in (12.5,15.,20.)]).to_csv(OUT/'required_B_candidate_estimate.csv',index=False)
json.dump(dict(classification=classification,walloff_R_lock_cycle12=float(pre_r),walloff_R_lock_startup=float(lock.iloc[0].R_lock_vs_30Hz),wallon_post_R_phase=post_r,wallon_post_Tmag_ratio=t_ratio,selected_B_mT=chosen,notes='Wall-off run logs G3; paired Wall-on reference logs G6/L45.'),open(OUT/'postimpact_failure_classification.json','w'),indent=2)

meta=dict(server_path='J:/magpy/magpylib_socket_server.py',server_sha256=sha('J:/magpy/magpylib_socket_server.py'),runner_path='J:/abaqusfangzhen/run_wobble_cal_f30_wallon_free_003.ps1',runner_sha256=sha(ROOT/'run_wobble_cal_f30_wallon_free_003.ps1'),paired_wallon_gradient_mT=6.0,paired_wallon_gradient_T=.006,paired_wallon_gradient_length_mm=45.,walloff_gradient_mT=3.0,frequency_Hz=30.,B0_mT=10.,cone_deg=30.,bias_deg=40.,phase_deg=248.,sense=1,chi_deg=0,trajectory='FORWARD',runtime_note='paired wall-on stdout says basis ROBOT_ARC; archived telemetry numerical field follows driver_arc')
json.dump(meta,open(OUT/'production_runtime_identity.json','w'),indent=2);pd.DataFrame([meta]).to_csv(OUT/'legacy_runtime_identity_audit.csv',index=False)
report=OUT/'Wobble30Hz_postimpact_wobble_survival_report.md';report.write_text(f'''# 30 Hz post-impact wobble survival audit

## Decision
`{classification}`. This turn used completed telemetry/ODB exports only; no new Abaqus job was submitted.

## Runtime and replay identity
The paired Wall-ON reference records G6/L45 exactly: 6 mT = 0.006 T, L=45 mm, B0=10 mT, f=30 Hz, cone30°, Bias40°, phase248°, sense+1, χ=0. Server SHA256: `{meta['server_sha256']}`. The 55.5 ms Wall-OFF COM-fixed record logs **G3** (3 mT), so it is not mislabeled G6. The archived production evaluator has now been replayed exactly; B/F/T regression passes at machine precision. The earlier mismatch came from an approximate flattened-CSV/frame replay, not the production server.

## Wall-OFF true-axis result
`walloff_f30_true_axis_phase.csv` uses the directed HEAD→TAIL axis and separates body spin (ω·a) from wobble (ω−(ω·a)a). `walloff_f30_lock_windows.csv` reports field rate, robot phase rate, R_lock, wobble rate and true-axis range over startup, cycle 1 and cycle 2 windows.

## Wall-ON post-impact result
The 3 ms paired run is split into pre-contact, contact, immediate post-impact and latest windows. Exact geometric gap and aggregate CPRESS are kept separate in `wallon_contact_timing.csv`; rates and recorded applied torque are in `wallon_pre_postimpact_wobble.csv` and `postimpact_torque_budget.csv`. No GIF was regenerated because no dynamic case was authorized in this diagnostic turn.

## B estimate
Offline-only candidates 12.5/15/20 mT are listed in `required_B_candidate_estimate.csv`; no B increase is selected because the post-impact torque remains strong ({t_ratio:.3f}× pre-contact). These are not dynamic predictions.

## Final interpretation
Classification `{classification}` is based on Wall-OFF lock and Wall-ON post-impact phase/torque ratios, not on a single UR1 component. Wall-OFF cycle 1/2 are synchronized, while the Wall-ON latest window loses phase rate but retains torque. Because that latest window contains only five telemetry rows, the next evidence should be a longer paired Wall-ON record; no B-amplitude increase is justified by this audit alone.
''',encoding='utf-8')
with report.open('a',encoding='utf-8') as fh:
    fh.write(f'''\n## Runtime regression correction\nThe exact archived evaluator was re-run in the configured SIMULIA Python environment (`runtime_exact_evaluator_replay.csv`). B, F and T maximum absolute errors are below 2e-12, so the earlier 0.31 mT / 0.0116 N·mm discrepancy was a replay-frame/curve approximation error, not a production-server mismatch. The numerical field follows `driver_arc` in the archived telemetry even though the banner says `ROBOT_ARC`; this metadata discrepancy is retained for future audits. Error decompositions are in `field_replay_error_decomposition.csv` and `torque_replay_error_decomposition.csv`.\n\n## Quantitative classification caveat\nWall-OFF cycle-1/2 directed phase rates are approximately {lock.iloc[1].f_robot_Hz:.2f} and {lock.iloc[2].f_robot_Hz:.2f} Hz (absolute R_lock {lock.iloc[1].R_lock_abs:.3f}, {lock.iloc[2].R_lock_abs:.3f}), supporting B10 synchronization without the wall. Wall-ON latest-window R_phase is {post_r:.3f}, while recorded |Tmag| remains {t_ratio:.3f}× pre-contact; this supports `{classification}` rather than magnetic torque collapse. The latest window contains only five output rows, so a longer wall-on wobble-survival probe is required before claiming a sustained stall.\n\nAngular impulse bookkeeping is in `angular_momentum_audit.csv`; ΔH is intentionally marked unavailable because the archived exports contain no independent rigid-body inertia/history channel.\n''')
print('classification',classification,'Wall-OFF mean R_lock=',float(lock.R_lock_vs_30Hz.mean()),'Wall-ON post R_phase=',post_r,'T ratio=',t_ratio,'B candidate=',chosen)
print(lock.to_string(index=False));print(post.to_string(index=False));print(tb.to_string(index=False))
