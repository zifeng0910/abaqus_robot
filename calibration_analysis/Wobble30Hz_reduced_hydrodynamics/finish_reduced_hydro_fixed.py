from pathlib import Path
import numpy as np, pandas as pd
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]; JOB='Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083'
d=pd.read_csv(HERE/'reduced_hydro_fixed_true_axis_phase.csv'); tel=pd.read_csv(ROOT/(JOB+'_telemetry.csv')); hyd=pd.read_csv(ROOT/'reduced_hydro_runtime_load.csv',skiprows=2,names=['time_s','v1','v2','v3','vr1','vr2','vr3','Fh1','Fh2','Fh3','Th1','Th2','Th3'])
def it(df,t,c):
    tc='t_s' if 't_s' in df else 'time_s'; return np.interp(t,df[tc],df[c])
t=d.time_s.to_numpy(); F=np.column_stack([it(tel,t,c) for c in ['fx_aba_N','fy_aba_N','fz_aba_N']]); H=np.column_stack([it(hyd,t,c) for c in ['Fh1','Fh2','Fh3']]); V=d[['Vt_direct_mm_s']].to_numpy()[:,0]; tang=np.column_stack([np.gradient(d.s_COM_mm,t)*0+1,np.zeros(len(t)),np.zeros(len(t))])
# Legacy sparse-output diagnostic only. 6.15 ms is first CAPTURED contact,
# not first physical contact. Use ReducedHydro_hidden_impact_audit for conclusions.
tc=0.0061499997973442; mask=t<=tc; m=1e-5; Vxyz=np.column_stack([it(pd.read_csv(ROOT/(JOB+'_rp_fields.csv')),t,c) for c in ['V1','V2','V3']]); P=m*(Vxyz-Vxyz[0])*1e-3; J=np.zeros_like(P)
for i in range(1,len(t)): J[i]=J[i-1]+.5*(F[i]+H[i]+F[i-1]+H[i-1])*(t[i]-t[i-1])
rows=[]
for name,sel in [('before_first_sparse_contact_NOT_CONTACT_FREE',mask),('all',np.ones(len(t),bool))]:
    j=int(np.where(sel)[0][-1]); res=P[j]-J[j]; rows.append({'window':name,'end_time_s':t[j],'mDeltaV_x_kgm_s':P[j,0],'mDeltaV_y_kgm_s':P[j,1],'mDeltaV_z_kgm_s':P[j,2],'Jmag_plus_Jhyd_x_Ns':J[j,0],'Jmag_plus_Jhyd_y_Ns':J[j,1],'Jmag_plus_Jhyd_z_Ns':J[j,2],'residual_norm_kgm_s':np.linalg.norm(res),'relative_residual':np.linalg.norm(res)/max(np.linalg.norm(P[j]),1e-30)})
pd.DataFrame(rows).to_csv(HERE/'reduced_hydro_fixed_linear_momentum.csv',index=False)
gap=pd.read_csv(ROOT/(JOB+'_exact_wall_penetration.csv')); cp=pd.read_csv(ROOT/(JOB+'_contact.csv')); cons=gap.merge(cp,on=['time_s'],how='left'); cons['discrepancy_flag']=(cons.CPRESS_max_MPa>0.01)&(cons.min_signed_gap_mm>0); cons['contact_surface_name']='ROBOT_SOLID_SURF ↔ PIPE_WALL_HELPER_SURF'; cons['geometric_wall_source']='Pipe_WALL_HELPER R3D4'; cons.to_csv(HERE/'reduced_hydro_fixed_contact_gap_consistency.csv',index=False)
pd.DataFrame([{'case':'corrected Reduced-Hydro','job':JOB,'architecture':'Abaqus rigid robot + SmoothWall114 single-wall contact + local reduced hydro; CEL fluid removed','CPRESS_provenance':'aggregate wall contact output (not CEL)','min_exact_gap_um':gap.min_signed_gap_mm.min()*1000,'CPRESS_max_MPa':cp.CPRESS_max_MPa.max(),'first_exact_contact_ms':gap.loc[gap.min_signed_gap_mm<=0,'time_s'].iloc[0]*1000 if (gap.min_signed_gap_mm<=0).any() else np.nan,'first_solver_contact_ms':cp.loc[cp.CPRESS_max_MPa>.01,'time_s'].iloc[0]*1000 if (cp.CPRESS_max_MPa>.01).any() else np.nan}]).to_csv(HERE/'cel_vs_reduced_hydro_fixed_summary.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False)); print('WROTE momentum, gap consistency and comparison CSVs')
