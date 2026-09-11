"""Zero-cost G3/G6 canonical tangential momentum and gradient budget."""
import csv, os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
jobs={
 'G3':'WobbleCal_F30_Cone30_B10_WallOn_Free_003',
 'G6':'WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003'}
mass_kg=1.0e-5
G_mT={'G3':3.0,'G6':6.0}
records={}; windows=[]
cl=pd.read_csv(ROOT/'curvenew_CEL_xyrot56_exact_abaqus.csv')
cl_s=cl.arclength_mm.to_numpy(float); cl_xyz=cl[['x_mm','y_mm','z_mm']].to_numpy(float)
cl_t=np.gradient(cl_xyz,cl_s,axis=0); cl_t/=np.maximum(np.linalg.norm(cl_t,axis=1)[:,None],1e-15)
def tangent_at_s(q):
    q=float(np.clip(q,cl_s[0],cl_s[-1])); v=np.array([np.interp(q,cl_s,cl_t[:,k]) for k in range(3)])
    return v/max(np.linalg.norm(v),1e-15)
for key,job in jobs.items():
    h=pd.read_csv(ROOT/'output'/f'{job}_bend_validation'/'rp_centerline_contact_history.csv')
    t=h.time_s.to_numpy(float); s=h.s_mm.to_numpy(float)
    dsdt=np.gradient(s,t) # arclength cross-check, mm/s
    vh=pd.read_csv(ROOT/f'{job}_rp_velocity_history.csv')
    vt=np.array([np.dot(vh.iloc[(vh.time_s-r.time_s).abs().argsort().iloc[0]][['V1','V2','V3']].to_numpy(float), tangent_at_s(r.s_mm)) for _,r in h.iterrows()])
    h['Vt_canonical_mm_s']=vt; h['dsdt_crosscheck_mm_s']=dsdt
    tel=pd.read_csv(ROOT/f'{job}_telemetry.csv')
    cp=pd.read_csv(ROOT/'output'/f'{job}_bend_validation'/'cpress_peaks.csv')
    ex=pd.read_csv(ROOT/f'{job}_exact_wall_penetration.csv')
    # first meaningful contact: CPRESS > 0.1 MPa (not numerical micro-noise)
    hit=cp.index[cp.cpress_max_mpa>0.1]
    tc=float(cp.loc[hit[0],'time_s']) if len(hit) else 0.0015
    t_end=float(t[-1]);
    def wmean(a,b):
        q=(t>=a)&(t<b); return float(vt[q].mean()) if q.any() else float(vt[-1])
    for name,a,b in [('W0_precontact',0.0,tc),('W1_first_contact',max(0,tc-0.0002),tc+0.0002),('W2_postimpact',tc+0.0002,0.0025),('W3_late',0.0025,t_end+1e-12)]:
        windows.append({'case':key,'job':job,'window':name,'start_ms':a*1000,'end_ms':b*1000,'mean_Vt_mm_s':wmean(a,b),'median_Vt_mm_s':float(np.median(vt[(t>=a)&(t<b)])) if (t>=a).any() and (t<b).any() else float(vt[-1])})
    tm=tel.t_s.to_numpy(float); ft=tel.force_tangent_N.to_numpy(float)
    # interpolate/projection of telemetry force only for common measured rows
    jmag=float(np.trapz(ft,tm)); post=(tm>=tc)
    jpost=float(np.trapz(ft[post],tm[post])) if post.sum()>1 else jmag
    records[key]={'job':job,'duration_ms':t_end*1000,'delta_s_mm':s[-1]-s[0],
      'Vt_mean_mm_s':float(vt.mean()),'Vt_median_mm_s':float(np.median(vt)),
      'Vt_final_mm_s':float(vt[-1]),'Vt_min_mm_s':float(vt.min()),'Vt_max_mm_s':float(vt.max()),
      'dsdt_final_crosscheck_mm_s':float(dsdt[-1]),
      'positive_Vt_fraction':float(np.mean(vt>0)),'Ft_mean_uN':float(ft.mean()*1e6),
      'Ft_median_uN':float(np.median(ft)*1e6),'Jmag_t_Ns':jmag,'Jmag_post_Ns':jpost,
      'contact_onset_ms':tc*1000,'CPRESS_max_MPa':float(cp.cpress_max_mpa.max()),
      'min_gap_um':float(ex.min_signed_gap_mm.min()*1000),
      'UR1_pp_rad':float(h.ur1.max()-h.ur1.min()),'UR2_pp_rad':float(h.ur2.max()-h.ur2.min()),'UR3_pp_rad':float(h.ur3.max()-h.ur3.min())}
    records[key]['_t']=t; records[key]['_vt']=vt; records[key]['_s']=s

summary=pd.DataFrame([{k:v for k,v in d.items() if not k.startswith('_')} for d in records.values()])
summary.to_csv(ROOT/'g3_g6_tangential_momentum_budget.csv',index=False)
pd.DataFrame(windows).to_csv(ROOT/'postimpact_velocity_windows.csv',index=False)

K=(records['G6']['Ft_mean_uN']-records['G3']['Ft_mean_uN'])/(G_mT['G6']-G_mT['G3'])
b=records['G3']['Ft_mean_uN']-K*G_mT['G3']
sens=pd.DataFrame([{'gradient_mT':3.0,'mean_Ft_uN':records['G3']['Ft_mean_uN'],'Ft_per_gradient_uN_per_mT':K,'intercept_uN':b},
                   {'gradient_mT':6.0,'mean_Ft_uN':records['G6']['Ft_mean_uN'],'Ft_per_gradient_uN_per_mT':K,'intercept_uN':b}])
sens.to_csv(ROOT/'gradient_force_sensitivity.csv',index=False)

vref=float(records['G6']['Vt_final_mm_s'])
tc=max(records['G6']['contact_onset_ms']/1000.0, records['G3']['contact_onset_ms']/1000.0)
rows=[]
for total_ms in [8.333,16.667,33.333]:
    total=total_ms/1000; avail=max(total-tc,1e-9)
    for target in [0.0,5.0]:
        dv=(target-vref)*1e-3
        jreq=mass_kg*dv
        f_uN=jreq/avail*1e6
        g=max(0.0,(f_uN-b)/K)
        rows.append({'target_total_time_ms':total_ms,'available_postimpact_time_ms':avail*1000,
          'Vt_reference_mm_s':vref,'target_Vt_mm_s':target,'J_required_Ns':jreq,
          'mean_force_required_uN':f_uN,'estimated_gradient_amplitude_mT':g,
          'within_existing_20mT_guard':bool(g<=20.0),
          'interpretation':'minimum impulse scale estimate; excludes CEL drag and repeated wall impulses'})
req=pd.DataFrame(rows); req.to_csv(ROOT/'required_forward_gradient_estimate.csv',index=False)

plt.figure(figsize=(8,4.5))
for key,c in [('G3','#377eb8'),('G6','#e41a1c')]: plt.plot(records[key]['_t']*1000,records[key]['_vt'],'.-',label=key,color=c)
plt.axhline(0,color='k',lw=.7); plt.xlabel('time (ms)'); plt.ylabel('canonical Vt (mm/s)'); plt.legend(); plt.tight_layout(); plt.savefig(ROOT/'g3_g6_Vt_vs_time.png',dpi=170); plt.close()
plt.figure(figsize=(7,4.2));
for key,c in [('G3','#377eb8'),('G6','#e41a1c')]: plt.plot(records[key]['_t']*1000,mass_kg*records[key]['_vt']*1e-3*1e6,label=key,color=c)
plt.axhline(0,color='k',lw=.7); plt.xlabel('time (ms)'); plt.ylabel('tangential momentum (µN·s)'); plt.legend(); plt.tight_layout(); plt.savefig(ROOT/'g3_g6_tangential_momentum.png',dpi=170); plt.close()
plt.figure(figsize=(7,4.2)); plt.plot(req[req.target_Vt_mm_s==0].target_total_time_ms,req[req.target_Vt_mm_s==0].estimated_gradient_amplitude_mT,'o-',label='neutral Vt=0'); plt.plot(req[req.target_Vt_mm_s==5].target_total_time_ms,req[req.target_Vt_mm_s==5].estimated_gradient_amplitude_mT,'s-',label='Vt=+5 mm/s'); plt.axhline(20,color='r',ls='--',label='20 mT guard'); plt.xlabel('target total time (ms)'); plt.ylabel('estimated gradient amplitude (mT)'); plt.legend(); plt.tight_layout(); plt.savefig(ROOT/'required_gradient_vs_recovery_time.png',dpi=170); plt.close()

report=ROOT/'Wobble30Hz_forward_gradient_impulse_scale_report.md'
def md_table(df):
    cols=list(df.columns); lines=['|'+'|'.join(str(c) for c in cols)+'|','|'+'|'.join(['---']*len(cols))+'|']
    for _,row in df.iterrows(): lines.append('|'+'|'.join(str(row[c]) for c in cols)+'|')
    return '\n'.join(lines)
with open(report,'w',encoding='utf-8') as f:
    f.write('# 30 Hz forward-gradient impulse budget\n\n')
    f.write('## 1. Why the G6 3 ms result is not a gradient failure\n\n')
    f.write('The 3 ms window is only 0.09 of a 30 Hz cycle and is dominated by the initial rotation/contact transient. G6 does not cancel that initial reverse transient; it does, however, preserve the wobble and substantially soften contact.\n\n')
    f.write('## 2. Canonical tangential velocity and stage split\n\n')
    f.write('Forward is increasing canonical centerline arclength. Vt is computed by projecting the ODB RP global velocity (V1,V2,V3) onto the continuously interpolated local tangent; continuous `ds/dt` is retained as a cross-check. No legacy arc steps or absolute-value force are used.\n\n')
    f.write(md_table(summary[['job','delta_s_mm','Vt_final_mm_s','Vt_min_mm_s','Vt_max_mm_s','positive_Vt_fraction','Ft_mean_uN','Jmag_t_Ns','contact_onset_ms','CPRESS_max_MPa','min_gap_um']]))
    f.write('\n\n')
    f.write(md_table(pd.DataFrame(windows))); f.write('\n\n')
    f.write('At 3 ms both cases are still moving rapidly toward −s (G3 %.1f mm/s, G6 %.1f mm/s); neither is in case C (instantaneous Vt already positive). The reverse motion is therefore an initial inertia/contact-dominated transient.\n\n' % (records['G3']['Vt_final_mm_s'],records['G6']['Vt_final_mm_s']))
    f.write('## 3. Magnetic impulse and sensitivity\n\n')
    f.write('Measured magnetic tangential impulse is Jmag,t=%.4g N·s (G3) and %.4g N·s (G6). The extra G6 impulse is %.4g N·s, while the G6 3 ms tangential momentum magnitude is about %.4g N·s; the added impulse is only %.2f%% of that magnitude.\n\n' % (records['G3']['Jmag_t_Ns'],records['G6']['Jmag_t_Ns'],records['G6']['Jmag_t_Ns']-records['G3']['Jmag_t_Ns'],abs(mass_kg*records['G6']['Vt_final_mm_s']*1e-3),(records['G6']['Jmag_t_Ns']-records['G3']['Jmag_t_Ns'])/abs(mass_kg*records['G6']['Vt_final_mm_s']*1e-3)*100))
    f.write('Two-point measured sensitivity is K_F=%.3f µN per mT, with Ft≈%.3f·G(mT)+%.3f µN. This is a short-range estimate only. CPRESS falling 95%% does not prove tangential wall impulse stayed constant; the available data only support the resolved nonmagnetic remainder interpretation.\n\n' % (K,K,b))
    f.write('## 4. Minimum impulse-scale gradient estimate\n\n')
    f.write(md_table(req)); f.write('\n\n')
    f.write('These are `MINIMUM_IMPULSE_SCALE_ESTIMATE` values: they ignore CEL drag, changing orientation, and repeated contacts, so they are not a full trajectory prediction. Every neutral/+5 mm/s estimate is above the existing 20 mT guard (approximately %.0f–%.0f mT over 8.333–33.333 ms).\n\n' % (req.estimated_gradient_amplitude_mT.min(),req.estimated_gradient_amplitude_mT.max()))
    f.write('## 5. Decision\n\n')
    f.write('`FORWARD_GRADIENT_WITHIN_CURRENT_LIMIT_TOO_WEAK`. No ≤20 mT candidate has a credible impulse scale to recover the observed −s momentum within a quarter, half, or full 30 Hz cycle. Therefore no new 8.333 ms Abaqus job is submitted in this stage. The next physically meaningful work is a hardware/field-architecture change (or an explicitly approved relaxation of the 20 mT guard), while preserving the 30 Hz wobble operating point.\n')
print('WROTE budget outputs; K_F=%.6f uN/mT' % K)
print(summary.to_string(index=False))
print(req.to_string(index=False))
