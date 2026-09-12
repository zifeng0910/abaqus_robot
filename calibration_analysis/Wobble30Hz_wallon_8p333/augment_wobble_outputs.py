from pathlib import Path
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
JOB='Wobble_F30_G6L45_WallOn_Free_0083_WobbleSurvival'
d=pd.read_csv(HERE/'wallon_8p333_angular_velocity.csv'); h=pd.read_csv(HERE/'wallon_8p333_angular_momentum.csv'); t=d.time_s.to_numpy()
tel=pd.read_csv(ROOT/(JOB+'_telemetry.csv')); T=tel[['tx_aba_Nmm','ty_aba_Nmm','tz_aba_Nmm']].to_numpy(float); om=d[['omega_x_rad_s','omega_y_rad_s','omega_z_rad_s']].to_numpy(float)
rows=[]
summary=pd.read_csv(HERE/'wallon_8p333_window_summary.csv')
for _,w in summary.iterrows():
    m=(t>=0)&(t<=t[-1]) if w.window=='W4_late_postimpact' else None
    # windows recovered from the summary boundaries by cumulative duration
    end=float(w.duration_ms/1000 + (0 if w.window=='W0_preimpact' else 0))
# use the exact data-driven boundaries reconstructed from first contact
gap=pd.read_csv(HERE.parent.parent.parent/(JOB+'_exact_wall_penetration.csv')); first=float(gap.loc[gap.min_signed_gap_mm<.02,'time_s'].iloc[0])
bounds=[('W0_preimpact',0,first),('W1_contact',first,min(first+.0005,t[-1])),('W2_immediate_postimpact',min(first+.0005,t[-1]),min(first+.002,t[-1])),('W3_intermediate',min(first+.002,t[-1]),min(first+.004,t[-1])),('W4_late_postimpact',max(first+.004,t[0]),t[-1])]
for name,lo,hi in bounds:
    sel=(t>=lo)&(t<=hi); ii=np.where(sel)[0]
    if len(ii)<2: continue
    J=np.trapezoid(T[ii],t[ii],axis=0) if hasattr(np,'trapezoid') else np.trapz(T[ii],t[ii],axis=0)
    dh=h.iloc[ii[-1]][['H_x_kgm2_s','H_y_kgm2_s','H_z_kgm2_s']].to_numpy(float)-h.iloc[ii[0]][['H_x_kgm2_s','H_y_kgm2_s','H_z_kgm2_s']].to_numpy(float)
    jw=dh-J*1e-3
    kw=float(h.iloc[ii[-1]].Krot_J-h.iloc[ii[0]].Krot_J); ww=float(d.iloc[ii[-1]].Wmag_J-d.iloc[ii[0]].Wmag_J)
    rows.append(dict(window=name,start_ms=lo*1000,end_ms=hi*1000,n=len(ii),DeltaH_x_kgm2_s=dh[0],DeltaH_y_kgm2_s=dh[1],DeltaH_z_kgm2_s=dh[2],DeltaH_norm_kgm2_s=float(np.linalg.norm(dh)),Jmag_x_Nmm_s=J[0],Jmag_y_Nmm_s=J[1],Jmag_z_Nmm_s=J[2],Jmag_norm_Nmm_s=float(np.linalg.norm(J)),Jnonmag_resolved_x_kgm2_s=jw[0],Jnonmag_resolved_y_kgm2_s=jw[1],Jnonmag_resolved_z_kgm2_s=jw[2],Jnonmag_resolved_norm_kgm2_s=float(np.linalg.norm(jw)),DeltaKrot_J=kw,Wmag_J=ww,Wnonmag_resolved_J=kw-ww))
pd.DataFrame(rows).to_csv(HERE/'wallon_8p333_angular_impulse_windows.csv',index=False)
# fixed-camera GIF
tri=pd.read_csv(ROOT/(JOB+'_pipe_wall_triangles_exact.csv')); wn=pd.read_csv(ROOT/(JOB+'_pipe_wall_nodes_exact.csv')); nd={int(r.node):[r.x_mm,r.y_mm,r.z_mm] for _,r in wn.iterrows()}; verts=[[nd[int(r[n])] for n in ('n1','n2','n3','n4')] for _,r in tri.iterrows()]; axraw=pd.read_csv(ROOT/(JOB+'_true_axis_raw.csv')); curve=pd.read_csv(ROOT/'curvenew_CEL_xyrot56_exact.csv'); Rtf=np.asarray(__import__('json').loads((ROOT/'abaqus_magpylib_frame_transform.json').read_text())['R_aba_to_mag'],float); origin=np.asarray(__import__('json').loads((ROOT/'abaqus_magpylib_frame_transform.json').read_text())['origin_aba_mm'],float); cv_aba=curve[['x_mm','y_mm','z_mm']].to_numpy(float)@Rtf+origin; ims=[]
for k in np.linspace(0,len(t)-1,min(84,len(t))).astype(int):
    fig=plt.figure(figsize=(7,5)); ax=fig.add_subplot(111,projection='3d'); ax.add_collection3d(Poly3DCollection(verts,alpha=.10,facecolor='gray',edgecolor='none')); ax.plot(cv_aba[:,0],cv_aba[:,1],cv_aba[:,2],'k-',lw=1,alpha=.7)
    rr=axraw.iloc[k]; hp=np.array([rr.head_x,rr.head_y,rr.head_z]); tp=np.array([rr.tail_x,rr.tail_y,rr.tail_z]); ax.plot([hp[0],tp[0]],[hp[1],tp[1]],[hp[2],tp[2]],'r-',lw=3); ax.scatter([hp[0]],[hp[1]],[hp[2]],c='r'); ax.scatter([tp[0]],[tp[1]],[tp[2]],c='b')
    allpts=np.vstack([wn[['x_mm','y_mm','z_mm']].to_numpy(float),cv_aba,axraw[['head_x','head_y','head_z','tail_x','tail_y','tail_z']].to_numpy(float).reshape(-1,3)]); lo=allpts.min(0)-.5; hi=allpts.max(0)+.5; ax.set_xlim(lo[0],hi[0]); ax.set_ylim(lo[1],hi[1]); ax.set_zlim(lo[2],hi[2]); ax.view_init(18,-60); ax.set_title('t=%.2f ms | phase-rate=%.1f Hz | |wobble|=%.0f rad/s'%(t[k]*1000,d.robot_phase_equiv_Hz.iloc[k],d.omega_wobble_rad_s.iloc[k])); fig.canvas.draw(); ims.append(Image.frombytes('RGB',fig.canvas.get_width_height(),fig.canvas.tostring_rgb())); plt.close(fig)
ims[0].save(HERE/'Wobble_F30_G6L45_WallOn_Free_8p333ms_WobbleSurvival.gif',save_all=True,append_images=ims[1:],duration=60,loop=0)
with (HERE/'Wobble30Hz_wallon_8p333_rotational_load_report.md').open('a',encoding='utf-8') as f:
    f.write('\\n## Integrated angular balance\\nThe window table now contains DeltaH, integrated magnetic angular impulse Jmag, resolved Jnonmag=DeltaH-Jmag, DeltaKrot, Wmag and Wnonmag. Jnonmag is an accounting identity, not an independent contact-torque measurement.\\n')
    for r in rows: f.write(' - %s: DeltaH=%.4g kg m2/s, |Jmag|=%.4g Nmm s, |Jnonmag|=%.4g kg m2/s, DeltaKrot=%.4g J, Wmag=%.4g J, Wnonmag=%.4g J.\\n'%(r['window'],r['DeltaH_norm_kgm2_s'],r['Jmag_norm_Nmm_s'],r['Jnonmag_resolved_norm_kgm2_s'],r['DeltaKrot_J'],r['Wmag_J'],r['Wnonmag_resolved_J']))
print('wrote impulse windows and GIF; late',rows[-1])
