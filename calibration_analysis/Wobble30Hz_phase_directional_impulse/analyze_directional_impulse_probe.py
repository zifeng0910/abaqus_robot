"""Post-process the single phase-directional 3 ms probe."""
from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
job='Wobble_F30_PhaseDirectionalAntiImpact_WallOn_Free_003'
bjob='WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003'
base=ROOT/'output'/(job+'_bend_validation'); bbase=ROOT/'output'/(bjob+'_bend_validation')
d=pd.read_csv(ROOT/(job+'_telemetry.csv')); v=pd.read_csv(ROOT/(job+'_rp_velocity_history.csv'))
cp=pd.read_csv(base/'cpress_peaks.csv'); gap=pd.read_csv(ROOT/(job+'_exact_wall_penetration.csv')); tr=pd.read_csv(base/'rp_centerline_contact_history.csv')
db=pd.read_csv(ROOT/(bjob+'_telemetry.csv')); vb=pd.read_csv(ROOT/(bjob+'_rp_velocity_history.csv')); cpb=pd.read_csv(bbase/'cpress_peaks.csv'); gb=pd.read_csv(ROOT/(bjob+'_exact_wall_penetration.csv')); trb=pd.read_csv(bbase/'rp_centerline_contact_history.csv')
cl=pd.read_csv(ROOT/'curvenew_CEL_xyrot56_exact_abaqus.csv')
def tangent(s):
    a=cl.arclength_mm.values; xyz=cl[['x_mm','y_mm','z_mm']].values; i=int(np.clip(np.searchsorted(a,float(s))-1,0,len(a)-2)); q=xyz[i+1]-xyz[i]; return q/np.linalg.norm(q)
tt=np.asarray([tangent(s) for s in tr.s_mm]); ttb=np.asarray([tangent(s) for s in trb.s_mm])
vt=np.einsum('ij,ij->i',v[['V1','V2','V3']].values,tt); vtb=np.einsum('ij,ij->i',vb[['V1','V2','V3']].values,ttb)
defs=pd.read_csv(ROOT/'phase_pulse_definition.csv').iloc[0]; cen=float(defs.phase_center_deg); half=float(defs.half_width_deg)
def wrap(x): return (x+180)%360-180
phase=np.degrees(np.interp(v.time_s,d.t_s,d.instantaneous_phase_rad))%360; dd=wrap(phase-cen); w=np.where(np.abs(dd)<half,.5*(1+np.cos(np.pi*dd/half)),0.0)
driver_arc=np.interp(v.time_s,d.t_s,d.driver_arc_mm); robot_arc=np.interp(v.time_s,d.t_s,d.robot_arc_mm)
motion=pd.DataFrame({'time_s':v.time_s,'s_robot_mm':tr.s_mm,'Vt_mm_s':vt,'Vmag_mm_s':np.linalg.norm(v[['V1','V2','V3']].values,axis=1),'phase_deg':phase,'pulse_weight':w,'pulse_Ft_mN':.4*w,'pulse_Fn_away_mN':3*w,'delta_s_driver_minus_robot_mm':driver_arc-robot_arc,'gap_mm':gap.min_signed_gap_mm.values,'CPRESS_MPa':cp.cpress_max_mpa.values,'UR1_rad':tr.ur1,'UR2_rad':tr.ur2,'UR3_rad':tr.ur3,'axis_tangent_angle_deg':tr.axis_tangent_angle_deg})
motion.to_csv(ROOT/(job+'_motion.csv'),index=False)
contact=pd.DataFrame({'time_s':cp.time_s,'CPRESS_max_MPa':cp.cpress_max_mpa,'contact_node':cp.node,'contact_x_mm':cp.x_mm,'contact_y_mm':cp.y_mm,'contact_z_mm':cp.z_mm,'exact_gap_mm':gap.min_signed_gap_mm,'signed_gap_node':gap.signed_node,'signed_wall_element':gap.signed_wall_element,'min_unsigned_distance_mm':gap.min_unsigned_distance_mm,'active_contact':cp.cpress_max_mpa>0.02})
contact.to_csv(ROOT/(job+'_contact.csv'),index=False)
# Momentum budget: telemetry force contains the legacy gradient plus pulse.
t=d.t_s.values; dt=np.diff(t,prepend=t[0]); phase_d=np.degrees(d.instantaneous_phase_rad.values)%360; wd=np.where(np.abs(wrap(phase_d-cen))<half,.5*(1+np.cos(np.pi*wrap(phase_d-cen)/half)),0.0); Jmag=float(np.sum(d.force_tangent_N.values*dt)); Jtensor=float(np.sum(.0004*wd*dt)); dv=float(vt[-1]-vt[0]); mdv=1e-5*dv/1000.0; Jnon=mdv-Jmag
mom=pd.DataFrame([{'window_s':t[-1]-t[0],'mass_kg':1e-5,'delta_Vt_mm_s':dv,'m_delta_Vt_Ns':mdv,'Jmag_total_Ns':Jmag,'Jtensor_t_Ns':Jtensor,'Jnonmag_resolved_Ns':Jnon,'momentum_residual_Ns':mdv-(Jmag+Jnon),'phase_pulse_positive_tangent':True,'source_force_units':'N','time_units':'s'}])
mom.to_csv(ROOT/(job+'_momentum.csv'),index=False)
wob=pd.DataFrame({'time_s':v.time_s,'UR1_rad':tr.ur1,'UR2_rad':tr.ur2,'UR3_rad':tr.ur3,'axis_tangent_angle_deg':tr.axis_tangent_angle_deg,'Vt_mm_s':vt}); wob.to_csv(ROOT/(job+'_wobble.csv'),index=False)
pd.DataFrame([{'rank':1,'phase_deg':cen,'phase_halfwidth_deg':half,'Ft_target_mN':.4,'Fn_away_target_mN':3.0,'selection':'single Stage-A candidate','sweep_performed':False}]).to_csv(ROOT/'phase_candidate_ranking.csv',index=False)
def metrics(name,x,vt0,cp0,g0,tr0):
    return {'case':name,'duration_ms':1000*float(x.time_s.iloc[-1]),'delta_s_mm':float(tr0.s_mm.iloc[-1]-tr0.s_mm.iloc[0]),'Vt_final_mm_s':float(vt0[-1]),'Vmag_peak_mm_s':float(np.linalg.norm(x[['V1','V2','V3']].values,axis=1).max()),'CPRESS_max_MPa':float(cp0.cpress_max_mpa.max()),'min_exact_gap_um':1000*float(g0.min_signed_gap_mm.min()),'UR1_pp_rad':float(tr0.ur1.max()-tr0.ur1.min())}
summary=pd.DataFrame([metrics('G6/L45 baseline',vb,vtb,cpb,gb,trb),metrics('directional pulse',v,vt,cp,gap,tr)])
summary.to_csv(ROOT/'baseline_vs_directional_impulse.csv',index=False)
fig,ax=plt.subplots(figsize=(7,4)); ax.plot(v.time_s*1e3,vt,label='directional pulse'); ax.plot(vb.time_s*1e3,vtb,label='G6/L45 baseline'); ax.axhline(0,color='k',lw=.6); ax.set(xlabel='time (ms)',ylabel='Vt (mm/s)'); ax.legend(); fig.tight_layout(); fig.savefig(ROOT/'baseline_vs_directional_Vt.png',dpi=180); plt.close(fig)
fig,ax=plt.subplots(2,1,figsize=(7,5),sharex=True); ax[0].plot(cp.time_s*1e3,cp.cpress_max_mpa,label='pulse'); ax[0].plot(cpb.time_s*1e3,cpb.cpress_max_mpa,label='baseline'); ax[0].set_ylabel('CPRESS (MPa)'); ax[0].legend(); ax[1].plot(gap.time_s*1e3,gap.min_signed_gap_mm*1000,label='pulse'); ax[1].plot(gb.time_s*1e3,gb.min_signed_gap_mm*1000,label='baseline'); ax[1].set(xlabel='time (ms)',ylabel='min gap (um)'); ax[1].legend(); fig.tight_layout(); fig.savefig(ROOT/'baseline_vs_directional_gap_cpress.png',dpi=180); plt.close(fig)
fig,ax=plt.subplots(figsize=(7,4)); ax.plot(tr.time_s*1e3,tr.ur1,label='pulse UR1'); ax.plot(trb.time_s*1e3,trb.ur1,label='baseline UR1'); ax.set(xlabel='time (ms)',ylabel='UR1 (rad)'); ax.legend(); fig.tight_layout(); fig.savefig(ROOT/'baseline_vs_directional_wobble.png',dpi=180); plt.close(fig)
fig,ax=plt.subplots(figsize=(7,4)); ax.bar(['mΔVt','Jmag','Jtensor','Jnonmag'],[mdv,Jmag,Jtensor,Jnon]); ax.set_ylabel('impulse (N s)'); fig.tight_layout(); fig.savefig(ROOT/'baseline_vs_directional_momentum.png',dpi=180); plt.close(fig)
print(summary.to_string(index=False)); print(mom.to_string(index=False))
