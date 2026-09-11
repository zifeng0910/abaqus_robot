"""Zero-cost broad-to-local gradient timing audit on measured G6 trajectory."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
traj = pd.read_csv(ROOT/'gradient_length_scale_trajectory.csv')
cp = pd.read_csv(ROOT/'output/WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003_bend_validation/cpress_peaks.csv')
gap = pd.read_csv(ROOT/'WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003_exact_wall_penetration.csv')
t = traj.time_s.to_numpy(float)
broad = traj.Ft_L45p0_N.to_numpy(float)
local = traj.Ft_L3p75_N.to_numpy(float)
cp_i = np.interp(t, cp.time_s, cp.cpress_max_mpa)
gap_i = np.interp(t, gap.time_s, gap.min_signed_gap_mm)
quiet = (cp_i <= 0.10) & (gap_i >= 0)
first_onset = float(cp[cp.cpress_max_mpa>0.1].time_s.iloc[0])
clear = None
for i in range(int(np.searchsorted(t, first_onset)) , len(t)):
    j = np.searchsorted(t, t[i] + 0.0002)
    if j < len(t) and quiet[i:j].all():
        clear = float(t[i]); break
if clear is None:
    clear = 0.0017
safety = 0.0002
switch = clear + safety
transition = 0.0005
x = np.clip((t-switch)/transition, 0, 1)
w = x*x*(3-2*x)
gated = (1-w)*broad + w*local
df = traj[['time_s','s_robot_mm','s_driver_mm','delta_s_mm']].copy()
df['CPRESS_max_MPa'] = cp_i; df['exact_gap_um'] = gap_i*1000
df['Ft_broad_N'] = broad; df['Ft_local_N'] = local; df['Ft_gated_N'] = gated
df['weight'] = w; df['dFt_gated_N_per_s'] = np.gradient(gated,t)
df.to_csv(ROOT/'broad_to_local_force_continuity.csv',index=False)
imax = int(np.argmax(cp_i)); ig = int(np.argmin(gap_i))
events = [('t0',0.0),('first_contact_onset',first_onset),
          ('CPRESS_peak',float(t[imax])),('min_gap',float(t[ig])),('contact_clear',clear),
          ('switch_start',switch),('switch_end',switch+transition),('t_end',float(t[-1]))]
pd.DataFrame(events,columns=['event','time_s']).to_csv(ROOT/'broad_to_local_event_timeline.csv',index=False)
bud=[]
for name,a in [('post_switch',switch),('post_transition',switch+transition)]:
    q=(t>=a); bud.append({'window':name,'start_ms':a*1000,'end_ms':t[-1]*1000,
      'J_forward_estimated_Ns':float(np.trapz(gated[q],t[q])),'mean_Ft_uN':float(gated[q].mean()*1e6),
      'positive_fraction':float(np.mean(gated[q]>0))})
pd.DataFrame(bud).to_csv(ROOT/'broad_to_local_impulse_estimate.csv',index=False)
cont={'t_clear_ms':clear*1000,'switch_start_ms':switch*1000,'switch_end_ms':(switch+transition)*1000,
      'max_abs_dFt_uN':float(np.max(np.abs(np.diff(gated)))*1e6),
      'max_abs_dFt_per_s':float(np.max(np.abs(np.gradient(gated,t)))),
      'mean_Ft_gated_uN':float(gated.mean()*1e6),'min_Ft_gated_uN':float(gated.min()*1e6),
      'positive_fraction':float(np.mean(gated>0))}
pd.DataFrame([cont]).to_csv(ROOT/'broad_to_local_continuity_summary.csv',index=False)
fig,ax=plt.subplots(4,1,figsize=(9,8),sharex=True)
ax[0].plot(t*1000,traj.Vt_global_mm_s,color='k'); ax[0].set_ylabel('Vt (mm/s)')
ax[1].plot(t*1000,broad*1e6,label='broad L45'); ax[1].plot(t*1000,local*1e6,label='local L3.75'); ax[1].plot(t*1000,gated*1e6,label='gated'); ax[1].set_ylabel('Ft (µN)'); ax[1].legend()
ax[2].plot(cp.time_s*1000,cp.cpress_max_mpa,color='tab:purple'); ax[2].set_ylabel('CPRESS (MPa)')
ax[3].plot(gap.time_s*1000,gap.min_signed_gap_mm*1000,color='tab:green'); ax[3].set_ylabel('gap (µm)'); ax[3].set_xlabel('time (ms)')
for a in ax: a.axvline(switch*1000,color='k',ls='--',lw=.8); a.axvline((switch+transition)*1000,color='k',ls=':',lw=.8)
fig.tight_layout(); fig.savefig(ROOT/'event_timeline_force_contact_velocity.png',dpi=180); plt.close(fig)
plt.figure(figsize=(7,4)); plt.plot(t*1000,traj.s_robot_mm-traj.s_robot_mm.iloc[0],color='k'); plt.xlabel('time (ms)'); plt.ylabel('robot Δs (mm)'); plt.tight_layout(); plt.savefig(ROOT/'canonical_s_vs_time.png',dpi=180); plt.close()
plt.figure(figsize=(7,4)); plt.plot(t*1000,broad*1e6,label='broad L45'); plt.plot(t*1000,local*1e6,label='local L3.75'); plt.plot(t*1000,gated*1e6,label='gated'); plt.xlabel('time (ms)'); plt.ylabel('Ft (µN)'); plt.legend(); plt.tight_layout(); plt.savefig(ROOT/'Ft_broad_local_gated.png',dpi=180); plt.close()
print('clear=%.3f ms switch=%.3f ms end=%.3f ms' % (clear*1000,switch*1000,(switch+transition)*1000))
print(pd.DataFrame([cont]).to_string(index=False)); print(pd.DataFrame(bud).to_string(index=False))
