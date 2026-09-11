"""Post-process local-gradient probe against G6 baseline (no solver run)."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'output'
CASES = {
    'G6_L45': 'WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003',
    'Local_G6_L3p75': 'Wobble_F30_LocalForwardGradient_WallOn_Free_0083',
}
MASS = 1e-5
cl = pd.read_csv(ROOT / 'curvenew_CEL_xyrot56_exact_abaqus.csv')
cs = cl.arclength_mm.to_numpy(float)
xyz = cl[['x_mm', 'y_mm', 'z_mm']].to_numpy(float)
tv = np.gradient(xyz, cs, axis=0)
tv /= np.maximum(np.linalg.norm(tv, axis=1)[:, None], 1e-15)

def tangent(s):
    q = float(np.clip(s, cs[0], cs[-1]))
    v = np.array([np.interp(q, cs, tv[:, j]) for j in range(3)])
    return v / max(np.linalg.norm(v), 1e-15)

def nearest(df, t):
    a = df.time_s.to_numpy(float)
    i = int(np.searchsorted(a, t)); i = min(max(i, 0), len(df)-1)
    if i and abs(a[i-1]-t) < abs(a[i]-t): i -= 1
    return df.iloc[i]

summaries = []
histories, telems, contacts = {}, {}, {}
for label, job in CASES.items():
    h = pd.read_csv(OUT / f'{job}_bend_validation' / 'rp_centerline_contact_history.csv')
    vh = pd.read_csv(ROOT / f'{job}_rp_velocity_history.csv')
    vt = np.array([np.dot(nearest(vh, r.time_s)[['V1','V2','V3']].to_numpy(float), tangent(r.s_mm)) for _, r in h.iterrows()])
    h['Vt_canonical_mm_s'] = vt
    h['dsdt_mm_s'] = np.gradient(h.s_mm.to_numpy(float), h.time_s.to_numpy(float))
    tel = pd.read_csv(ROOT / f'{job}_telemetry.csv')
    cp = pd.read_csv(OUT / f'{job}_bend_validation' / 'cpress_peaks.csv')
    ex = pd.read_csv(ROOT / f'{job}_exact_wall_penetration.csv')
    histories[label], telems[label], contacts[label] = h, tel, (cp, ex)
    onset_rows = cp[cp.cpress_max_mpa > 0.1]
    onset = float(onset_rows.time_s.iloc[0]) if len(onset_rows) else np.nan
    imax = cp.cpress_max_mpa.idxmax(); igap = ex.min_signed_gap_mm.idxmin()
    ft, tm = tel.force_tangent_N.to_numpy(float), tel.t_s.to_numpy(float)
    summaries.append({
        'case': label, 'job': job, 'duration_ms': h.time_s.iloc[-1]*1000,
        's_start_mm': h.s_mm.iloc[0], 's_end_mm': h.s_mm.iloc[-1], 'delta_s_mm': h.s_mm.iloc[-1]-h.s_mm.iloc[0],
        'Vt_mean_mm_s': vt.mean(), 'Vt_final_mm_s': vt[-1], 'Vt_min_mm_s': vt.min(), 'Vt_max_mm_s': vt.max(),
        'positive_Vt_fraction': np.mean(vt > 0), 'precontact_mean_Vt_mm_s': vt[h.time_s.to_numpy(float) < onset].mean() if np.isfinite(onset) and np.any(h.time_s < onset) else np.nan,
        'Ft_mean_uN': ft.mean()*1e6, 'Ft_min_uN': ft.min()*1e6, 'Ft_max_uN': ft.max()*1e6, 'Jmag_t_Ns': np.trapz(ft, tm),
        'contact_onset_ms': onset*1000, 'CPRESS_max_MPa': cp.cpress_max_mpa.max(), 'CPRESS_peak_time_ms': cp.loc[imax, 'time_s']*1000,
        'min_gap_um': ex.min_signed_gap_mm.min()*1000, 'gap_time_ms': ex.loc[igap, 'time_s']*1000,
        'UR1_pp_rad': h.ur1.max()-h.ur1.min(), 'UR2_pp_rad': h.ur2.max()-h.ur2.min(), 'UR3_pp_rad': h.ur3.max()-h.ur3.min(),
        'max_abs_axis_angle_deg': np.abs(h.axis_tangent_angle_deg).max(),
    })
summary = pd.DataFrame(summaries)
summary.to_csv(ROOT / 'local_gradient_8p333ms_canonical_motion.csv', index=False)

windows = []
for label, h in histories.items():
    t, v = h.time_s.to_numpy(float), h.Vt_canonical_mm_s.to_numpy(float)
    onset = float(summary.loc[summary.case == label, 'contact_onset_ms'].iloc[0]) / 1000
    stages = [('precontact', 0, onset), ('impact', max(0, onset-.0002), min(t[-1], onset+.0002)),
              ('postimpact', min(t[-1], onset+.0002), min(t[-1], onset+.0015)), ('late', min(t[-1], onset+.0015), t[-1]+1e-12)]
    for name, a, b in stages:
        q = (t >= a) & (t < b)
        windows.append({'case': label, 'window': name, 'start_ms': a*1000, 'end_ms': b*1000,
                        'mean_Vt_mm_s': v[q].mean() if q.any() else np.nan, 'median_Vt_mm_s': np.median(v[q]) if q.any() else np.nan,
                        'positive_fraction': np.mean(v[q] > 0) if q.any() else np.nan})
win = pd.DataFrame(windows)
win.to_csv(ROOT / 'local_gradient_8p333ms_velocity_windows.csv', index=False)

budget = []
for _, r in summary.iterrows():
    h = histories[r.case]; dp = MASS * (h.Vt_canonical_mm_s.iloc[-1] - h.Vt_canonical_mm_s.iloc[0]) * 1e-3
    budget.append({'case': r.case, 'delta_p_t_Ns': dp, 'Jmag_t_Ns': r.Jmag_t_Ns,
                   'magnetic_only_residual_Ns': dp-r.Jmag_t_Ns, 'residual_over_delta_p': (dp-r.Jmag_t_Ns)/(abs(dp)+1e-30)})
pd.DataFrame(budget).to_csv(ROOT / 'local_gradient_8p333ms_momentum_budget.csv', index=False)

events = []
for label, (cp, _) in contacts.items():
    for _, r in cp.loc[cp.cpress_max_mpa.nlargest(8).index].iterrows():
        events.append({'case': label, 'time_ms': r.time_s*1000, 'CPRESS_MPa': r.cpress_max_mpa,
                       'instance': r.instance, 'node': r.node, 'x_mm': r.x_mm, 'y_mm': r.y_mm, 'z_mm': r.z_mm})
pd.DataFrame(events).to_csv(ROOT / 'local_gradient_8p333ms_contact.csv', index=False)
summary[['case','UR1_pp_rad','UR2_pp_rad','UR3_pp_rad','max_abs_axis_angle_deg']].to_csv(ROOT / 'local_gradient_8p333ms_wobble.csv', index=False)

colors = {'G6_L45':'#377eb8', 'Local_G6_L3p75':'#e41a1c'}
plt.figure(figsize=(7,4));
for k,h in histories.items(): plt.plot(h.time_s*1000, h.s_mm-h.s_mm.iloc[0], label=k, color=colors[k])
plt.axhline(0,color='k',lw=.6); plt.xlabel('time (ms)'); plt.ylabel('canonical Δs (mm)'); plt.legend(); plt.tight_layout(); plt.savefig(ROOT/'canonical_s_vs_time_G6_vs_localgradient.png',dpi=180); plt.close()
plt.figure(figsize=(7,4));
for k,h in histories.items(): plt.plot(h.time_s*1000, h.Vt_canonical_mm_s, label=k, color=colors[k])
plt.axhline(0,color='k',lw=.6); plt.xlabel('time (ms)'); plt.ylabel('canonical Vt (mm/s)'); plt.legend(); plt.tight_layout(); plt.savefig(ROOT/'Vt_vs_time_G6_vs_localgradient.png',dpi=180); plt.close()
fig, ax = plt.subplots(1,3, figsize=(11,3.5))
for j,c in enumerate(['UR1_pp_rad','UR2_pp_rad','UR3_pp_rad']): ax[j].bar(summary.case,summary[c],color=[colors[x] for x in summary.case]); ax[j].set_title(c); ax[j].tick_params(axis='x',rotation=25)
fig.tight_layout(); fig.savefig(ROOT/'wobble_G6_vs_localgradient.png',dpi=180); plt.close(fig)
fig, ax = plt.subplots(2,1, figsize=(7,5), sharex=True)
for k,(cp,ex) in contacts.items(): ax[0].plot(cp.time_s*1000,cp.cpress_max_mpa,label=k,color=colors[k]); ax[1].plot(ex.time_s*1000,ex.min_signed_gap_mm*1000,label=k,color=colors[k])
ax[0].set_ylabel('CPRESS max (MPa)'); ax[1].set_ylabel('min signed gap (µm)'); ax[1].set_xlabel('time (ms)'); ax[0].legend(); ax[0].axhline(10,color='k',ls='--',lw=.6); ax[1].axhline(0,color='k',lw=.6); fig.tight_layout(); fig.savefig(ROOT/'contact_G6_vs_localgradient.png',dpi=180); plt.close(fig)

report_path = ROOT/'Wobble30Hz_local_gradient_lengthscale_report.md'
old = report_path.read_text(encoding='utf-8')
base = old.split('\n## 6.', 1)[0]
with open(report_path,'w',encoding='utf-8') as f:
    r = summary.set_index('case')
    f.write(base.rstrip()+'\n\n')
    def md_table(df):
        cols=list(df.columns); out=['|'+'|'.join(str(c) for c in cols)+'|','|'+'|'.join(['---']*len(cols))+'|']
        for _,row in df.iterrows(): out.append('|'+'|'.join(str(row[c]) for c in cols)+'|')
        return '\n'.join(out)
    f.write('\n\n## 6. Dynamic 8.333 ms validation (G=6 mT, L=3.75 mm)\n\n')
    f.write('The production Gaussian gradient was tested in the Wall-ON CEL model with no other physical parameter change.\n\n')
    f.write(md_table(summary)); f.write('\n\n')
    f.write('The local profile ends at **Δs = %.4f mm** and **Vt = %.1f mm/s**, while the G6/L45 baseline ends at **Δs = %.4f mm** and **Vt = %.1f mm/s**. The candidate CPRESS peak is %.2f MPa and exact minimum signed gap %.3f µm.\n\n' % (r.loc['Local_G6_L3p75','delta_s_mm'],r.loc['Local_G6_L3p75','Vt_final_mm_s'],r.loc['G6_L45','delta_s_mm'],r.loc['G6_L45','Vt_final_mm_s'],r.loc['Local_G6_L3p75','CPRESS_max_MPa'],r.loc['Local_G6_L3p75','min_gap_um']))
    f.write('### Stage-wise canonical velocity\n\n'+md_table(win)+'\n\n')
    f.write('The offline force increase is real, but concentrating it near the measured positive Δs (~5.1 mm) injects a large force during the contact-sensitive initial state. In the coupled solve the first-contact CPRESS and penetration worsen, and no post-impact forward recovery appears. UR1 peak-to-peak remains %.3f rad, so this is not loss of magnetic torque; it is a translation/contact failure.\n\n' % r.loc['Local_G6_L3p75','UR1_pp_rad'])
    f.write('**Decision: `LOCAL_GRADIENT_PROFILE_DOES_NOT_RECOVER_FORWARD_MOTION`.** Do not submit a 16.667 ms continuation or increase G / further shrink L in this scalar Gaussian architecture. Preserve G6 as the softer-contact wobble reference and move to impulse-cancellation or another spatial/temporal field architecture.\n')
print(summary.to_string(index=False))
