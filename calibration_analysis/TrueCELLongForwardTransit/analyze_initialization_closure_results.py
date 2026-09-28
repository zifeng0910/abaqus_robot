import csv,json,math
from pathlib import Path
import numpy as np
ROOT=Path('J:/abaqusfangzhen/abaqus_robot/calibration_analysis/TrueCELLongForwardTransit'); OUT=ROOT/'initialization_closure'; P=ROOT/'axial_boundary_closure/facefix'; C=OUT/'child'
AX=np.array([.9762799602464296,-.0618320247446977,.2074951564186524]); AX/=np.linalg.norm(AX); N=np.array([.07770950949965555,.9945677258930866,-.06925511349479567]); N/=np.linalg.norm(N); ROCK=np.cross(AX,N); ROCK/=np.linalg.norm(ROCK); AX0=AX
def load(folder):
 z=np.load(folder/'rp_history.npz');t=z['U1'][:,0]
 def vec(p):return np.column_stack([np.interp(t,z[f'{p}{i}'][:,0],z[f'{p}{i}'][:,1]) for i in (1,2,3)])
 u,ur,v,vr=vec('U'),vec('UR'),vec('V'),vec('VR'); th=np.linalg.norm(ur,axis=1);k=np.divide(ur,th[:,None],out=np.zeros_like(ur),where=th[:,None]>1e-14); body=AX0*np.cos(th[:,None])+np.cross(k,AX0)*np.sin(th[:,None])+k*(k@AX0)[:,None]*(1-np.cos(th[:,None])); alpha=np.degrees(np.arctan2(body@N,body@AX)); return {'t':t,'s':u@AX,'v':v@AX,'alpha':alpha,'omega':vr@ROCK}
def val(m,k,t):return np.interp(t,m['t'],m[k])
def zeros(m,a,b):
 t=m['t'];v=m['v'];q=(t>a+1e-10)&(t<b-1e-10);tt=t[q];vv=v[q];o=[]
 for i in range(len(tt)-1):
  if vv[i]*vv[i+1]<0:o.append(float(tt[i]-vv[i]*(tt[i+1]-tt[i])/(vv[i+1]-vv[i])))
 return o
def cycle(name,m,a):
 b=a+.01;t=np.r_[a,m['t'][(m['t']>a)&(m['t']<b)],b];s=val(m,'s',t);v=val(m,'v',t);neg=v<0
 return {'case':name,'cycle':int(a/.01)+1,'delta_s_mm':float(s[-1]-s[0]),'mean_v_s_mm_s':float((s[-1]-s[0])/.01),'end_v_s_mm_s':float(v[-1]),'minimum_v_s_mm_s':float(v.min()),'maximum_v_s_mm_s':float(v.max()),'minimum_s_mm':float(s.min()),'maximum_s_mm':float(s.max()),'velocity_crosses_zero':bool(zeros(m,a,b)),'zero_crossing_ms':'|'.join(format(x*1000,'.6f') for x in zeros(m,a,b)),'negative_duration_ms':float(np.sum(np.diff(t)[neg[:-1]])*1000),'reverse_displacement_mm':float(np.sum(np.maximum(0,-np.diff(s))))}
def write(name,rows):
 with (OUT/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def energy(folder,name):
 z=np.load(folder/'energy_history.npz');o=[]
 for k in sorted(z.files):
  q=z[k][z[k][:,0]<=.02000001];o.append({'case':name,'quantity':k,'at_0ms':float(q[0,1]),'at_10ms':float(np.interp(.01,q[:,0],q[:,1])),'at_20ms':float(q[-1,1]),'min_0_20ms':float(q[:,1].min()),'max_0_20ms':float(q[:,1].max())})
 return o
p=load(P);c=load(C);tp,tc=p['t'],c['t'];j=np.searchsorted(tc,tp);common=[]
for i,x in enumerate(tp):
 for jj in (j[i]-1,j[i]):
  if 0<=jj<len(tc) and abs(tc[jj]-x)<=5e-8:common.append({'time_s':float(x),'parent_s_mm':float(p['s'][i]),'child_s_mm':float(c['s'][jj]),'delta_s_mm':float(c['s'][jj]-p['s'][i]),'parent_v_s_mm_s':float(p['v'][i]),'child_v_s_mm_s':float(c['v'][jj]),'delta_v_s_mm_s':float(c['v'][jj]-p['v'][i]),'parent_rocking_deg':float(p['alpha'][i]),'child_rocking_deg':float(c['alpha'][jj]),'parent_omega_rad_s':float(p['omega'][i]),'child_omega_rad_s':float(c['omega'][jj])});break
write('common_timestamp_FACEFIX20_vs_INITFIX20.csv',common); cycles=[]
for name,m in [('FACEFIX20',p),('INITFIX20',c)]:cycles += [cycle(name,m,0),cycle(name,m,.01)]
write('cycle_metrics_FACEFIX20_vs_INITFIX20.csv',cycles)
events=[]
for name,m in [('FACEFIX20',p),('INITFIX20',c)]:
 t=m['t'];v=m['v']; zc=zeros(m,.01,.02); idx=np.flatnonzero((t>.0085)&(t<.01)&(v<0)); onset=float(t[idx[0]]) if len(idx) else math.nan; ii=np.flatnonzero((t>.0085)&(t<.01)&(v<-.3)); events.append({'case':name,'native_negative_v_onset_ms':onset*1000 if np.isfinite(onset) else 'NOT RESOLVABLE','zero_crossings_ms':'|'.join(format(x*1000,'.6f') for x in zc),'max_reverse_v_mm_s':float(v.min()),'max_reverse_v_time_ms':float(t[v.argmin()]*1000),'v_negative_duration_cycle2_ms':cycles[1 if name=='FACEFIX20' else 3]['negative_duration_ms'],'extra_deceleration_onset_ms':float(t[ii[0]]*1000) if len(ii) else 'NOT RESOLVABLE'})
write('event_timings_FACEFIX20_vs_INITFIX20.csv',events); er=energy(P,'FACEFIX20')+energy(C,'INITFIX20');write('energy_work_FACEFIX20_vs_INITFIX20.csv',er)
# Native field/EVF summaries, including startup 0-2 ms extrema and final volume.
field=[]
for name,folder in [('FACEFIX20',ROOT/'axial_boundary_closure/facefix'),('INITFIX20',C)]:
 fp=folder/'fluid_summary_native.csv'
 if fp.exists():
  rows=list(csv.DictReader(fp.open()));
  for r in rows:
   if float(r['time_s'])<=.0200001:r['case']=name;field.append(r)
write('fluid_volume_EVF_history_FACEFIX20_vs_INITFIX20.csv',field)
# Earliest material divergence: absolute displacement > 0.01 mm sustained 25 us.
arr=np.array([[r['time_s'],r['delta_s_mm']] for r in common]); mask=np.abs(arr[:,1])>.01; earliest=float(arr[np.flatnonzero(mask)[0],0]) if mask.any() else math.nan
summary={'parent':'F100_G2P20_ZEROPRESSURE_FACEFIX20','child':'F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX20','common_timestamp_count':len(common),'common_timestamp_tolerance_s':5e-8,'max_abs_displacement_difference_mm':float(np.max(np.abs(arr[:,1]))),'rms_displacement_difference_mm':float(np.sqrt(np.mean(arr[:,1]**2))),'max_abs_velocity_difference_mm_s':float(np.max(np.abs([r['delta_v_s_mm_s'] for r in common]))),'rms_velocity_difference_mm_s':float(np.sqrt(np.mean([r['delta_v_s_mm_s']**2 for r in common]))),'earliest_abs_displacement_difference_gt_0p01mm_ms':earliest*1000 if np.isfinite(earliest) else 'NOT REACHED','cycles':cycles,'events':events,'energy':er,'fluid_history_rows':len(field),'startup_window':'0–2 ms native field rows retained; interpretation diagnostic only'}
(OUT/'initialization_closure_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
