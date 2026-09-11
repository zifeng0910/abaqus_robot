"""Compare the chi=134 dynamic probe with the G6/L45 baseline."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
CL=pd.read_csv(ROOT/'curvenew_CEL_xyrot56_exact_abaqus.csv')
def tangent(s):
    a=CL.arclength_mm.to_numpy(float); x=CL[['x_mm','y_mm','z_mm']].to_numpy(float)
    i=int(np.clip(np.searchsorted(a,float(s))-1,0,len(a)-2)); v=x[i+1]-x[i]
    return v/max(np.linalg.norm(v),1e-15)
def load(job):
    d=pd.read_csv(ROOT/(job+'_telemetry.csv')); v=pd.read_csv(ROOT/(job+'_rp_velocity_history.csv'))
    b=ROOT/'output'/(job+'_bend_validation'); return d,v,pd.read_csv(b/'rp_centerline_contact_history.csv'),pd.read_csv(b/'cpress_peaks.csv'),pd.read_csv(ROOT/(job+'_exact_wall_penetration.csv'))
def metrics(name,job):
    d,v,tr,cp,g=load(job); s=tr.s_mm.to_numpy(float); tv=np.asarray([tangent(x) for x in s]); V=v[['V1','V2','V3']].to_numpy(float); vt=np.einsum('ij,ij->i',V,tv); hit=np.flatnonzero(cp.cpress_max_mpa.to_numpy(float)>.02); k=int(hit[0]) if len(hit) else 0; pre=max(0,k-2); post=min(len(vt)-1,k+2)
    return dict(case=name,job=job,duration_ms=1000*v.time_s.iloc[-1],delta_s_mm=s[-1]-s[0],Vt_final_mm_s=vt[-1],Vt_peak_mm_s=vt.max(),Vmag_peak_mm_s=np.linalg.norm(V,axis=1).max(),forward_fraction=np.mean(vt>0),first_contact_time_ms=1000*cp.time_s.iloc[k],first_impact_deltaVt_mm_s=vt[post]-vt[pre],CPRESS_max_MPa=cp.cpress_max_mpa.max(),exact_min_gap_um=1000*g.min_signed_gap_mm.min(),near_wall_frames_gap_lt20um=int(np.sum(g.min_signed_gap_mm<.020)),UR1_pp_rad=tr.ur1.max()-tr.ur1.min(),true_axis_angle_range_deg=tr.axis_tangent_angle_deg.max()-tr.axis_tangent_angle_deg.min(),true_axis_angle_max_deg=tr.axis_tangent_angle_deg.max(),contact_node_at_peak=int(cp.node.iloc[int(cp.cpress_max_mpa.argmax())]))
jobs=[('G6/L45 chi=0','WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003'),('chi=134 deg candidate','Wobble_F30_ConeFrameAzimuth_Chi134_WallOn_Free_003')]
out=pd.DataFrame([metrics(*x) for x in jobs]); out.to_csv(ROOT/'baseline_vs_chi_dynamic.csv',index=False)
for metric in ['delta_s_mm','Vt_final_mm_s','CPRESS_max_MPa','exact_min_gap_um','UR1_pp_rad','true_axis_angle_range_deg']:
    b=float(out.iloc[0][metric]); q=float(out.iloc[1][metric]); out.loc[1,metric+'_change_pct']=100*(q-b)/max(abs(b),1e-30)
out.to_csv(ROOT/'baseline_vs_chi_dynamic.csv',index=False)
for y,n,lab in [('Vt_final_mm_s','baseline_vs_chi_Vt.png','Vt (mm/s)'),('exact_min_gap_um','baseline_vs_chi_exact_gap.png','minimum exact gap (um)'),('true_axis_angle_range_deg','baseline_vs_chi_true_wobble.png','true axis angular range (deg)')]:
    plt.figure(figsize=(7,4))
    for _,job in jobs:
        d,v,tr,cp,g=load(job); s=tr.s_mm.to_numpy(float); tv=np.asarray([tangent(x) for x in s]); V=v[['V1','V2','V3']].to_numpy(float); vt=np.einsum('ij,ij->i',V,tv)
        if y=='Vt_final_mm_s': xx=v.time_s.to_numpy(float)*1000; z=vt
        elif y=='exact_min_gap_um': xx=g.time_s.to_numpy(float)*1000; z=g.min_signed_gap_mm.to_numpy(float)*1000
        else: xx=tr.time_s.to_numpy(float)*1000; z=tr.axis_tangent_angle_deg.to_numpy(float)
        plt.plot(xx,z,label=job)
    plt.xlabel('time (ms)'); plt.ylabel(lab); plt.grid(alpha=.3); plt.legend(fontsize=7); plt.tight_layout(); plt.savefig(ROOT/n,dpi=180); plt.close()
print(out.to_string(index=False))
