import csv, json, math
from pathlib import Path
import numpy as np

ROOT=Path(r'J:/abaqusfangzhen/abaqus_robot/calibration_analysis/TrueCELLongForwardTransit')
OUT=ROOT/'initialization_closure'; P=ROOT/'initialization_closure'/'child'; C=OUT/'fluidfric0'
AX=np.array([.9762799602464296,-.0618320247446977,.2074951564186524]); AX/=np.linalg.norm(AX)
N=np.array([.07770950949965555,.9945677258930866,-.06925511349479567]); N/=np.linalg.norm(N)
ROCK=np.cross(AX,N); ROCK/=np.linalg.norm(ROCK)

def load(folder):
 z=np.load(folder/'rp_history.npz'); t=z['U1'][:,0]
 def vec(prefix): return np.column_stack([np.interp(t,z[f'{prefix}{i}'][:,0],z[f'{prefix}{i}'][:,1]) for i in (1,2,3)])
 u,ur,v,vr=vec('U'),vec('UR'),vec('V'),vec('VR')
 th=np.linalg.norm(ur,axis=1); k=np.divide(ur,th[:,None],out=np.zeros_like(ur),where=th[:,None]>1e-14)
 body=AX*np.cos(th[:,None])+np.cross(k,AX)*np.sin(th[:,None])+k*(k@AX)[:,None]*(1-np.cos(th[:,None]))
 return {'t':t,'s':u@AX,'v':v@AX,'alpha':np.degrees(np.arctan2(body@N,body@AX)),'omega':vr@ROCK}

def zeros(t,v,a,b):
 out=[]
 for i in range(len(t)-1):
  if t[i]<a or t[i+1]>b: continue
  if v[i]*v[i+1]<0: out.append(float(t[i]-v[i]*(t[i+1]-t[i])/(v[i+1]-v[i])))
 return out

def interp(m,key,t): return np.interp(t,m['t'],m[key])

def negative_duration(t,v,a,b):
 ts=[a]+[x for x in t if a<x<b]+[b]
 # direct interpolation avoids assumptions about native timestamp spacing
 vs=np.interp(ts,t,v); total=0.0
 for x0,x1,y0,y1 in zip(ts[:-1],ts[1:],vs[:-1],vs[1:]):
  if y0<0 and y1<0: total+=x1-x0
  elif y0*y1<0:
   xc=x0-y0*(x1-x0)/(y1-y0)
   total += (xc-x0) if y0<0 else (x1-xc)
 return total

def cycle(name,m,a):
 b=a+.01; ts=np.r_[a,m['t'][(m['t']>a)&(m['t']<b)],b]; s=interp(m,'s',ts); v=interp(m,'v',ts); zz=zeros(m['t'],m['v'],a,b)
 reverse=0.0
 for x0,x1,y0,y1 in zip(ts[:-1],ts[1:],s[:-1],s[1:]):
  if y1<y0: reverse += y0-y1
 return {'case':name,'cycle':int(a/.01)+1,'delta_s_mm':float(s[-1]-s[0]),'mean_v_s_mm_s':float((s[-1]-s[0])/.01),'end_v_s_mm_s':float(v[-1]),'minimum_v_s_mm_s':float(v.min()),'maximum_v_s_mm_s':float(v.max()),'minimum_s_mm':float(s.min()),'maximum_s_mm':float(s.max()),'zero_crossing_ms':'|'.join(f'{z*1000:.6f}' for z in zz),'negative_v_duration_ms':negative_duration(m['t'],m['v'],a,b)*1000,'reverse_displacement_mm':float(reverse)}

def energy(folder,name):
 z=np.load(folder/'energy_history.npz'); out=[]
 for k in sorted(z.files):
  q=z[k]; out.append({'case':name,'quantity':k,'at_0ms':float(q[0,1]),'at_10ms':float(np.interp(.01,q[:,0],q[:,1])),'at_20ms':float(q[-1,1]),'min_0_20ms':float(q[:,1].min()),'max_0_20ms':float(q[:,1].max())})
 return out

parent=load(P); child=load(C); cycles=[cycle('INITFIX20',parent,0),cycle('INITFIX20',parent,.01),cycle('FLUIDFRIC0',child,0),cycle('FLUIDFRIC0',child,.01)]
with (OUT/'cycle_metrics_INITFIX20_vs_FLUIDFRIC0.csv').open('w',newline='') as f: w=csv.DictWriter(f,fieldnames=cycles[0].keys()); w.writeheader(); w.writerows(cycles)
common=[]
for i,t in enumerate(parent['t']):
 j=int(np.argmin(abs(child['t']-t))); common.append({'time_s':float(t),'delta_s_mm':float(child['s'][j]-parent['s'][i]),'delta_v_s_mm_s':float(child['v'][j]-parent['v'][i]),'parent_s_mm':float(parent['s'][i]),'child_s_mm':float(child['s'][j]),'parent_v_s_mm_s':float(parent['v'][i]),'child_v_s_mm_s':float(child['v'][j]),'delta_alpha_deg':float(child['alpha'][j]-parent['alpha'][i]),'delta_omega_rad_s':float(child['omega'][j]-parent['omega'][i]),'time_delta_s':float(child['t'][j]-t)})
with (OUT/'common_timestamp_INITFIX20_vs_FLUIDFRIC0.csv').open('w',newline='') as f: w=csv.DictWriter(f,fieldnames=common[0].keys()); w.writeheader(); w.writerows(common)
arr=np.array([r['delta_s_mm'] for r in common]); vel=np.array([r['delta_v_s_mm_s'] for r in common]); early=np.flatnonzero((np.abs(arr)>.01)&(np.array([r['time_s'] for r in common])<.008));
events=[]
for name,m in [('INITFIX20',parent),('FLUIDFRIC0',child)]:
 z1=zeros(m['t'],m['v'],0,.01); z2=zeros(m['t'],m['v'],.01,.02); events.append({'case':name,'cycle1_zero_crossing_ms':'|'.join(f'{z*1000:.6f}' for z in z1),'cycle2_zero_crossing_ms':'|'.join(f'{z*1000:.6f}' for z in z2),'min_v_time_ms':float(m['t'][np.argmin(m['v'])]*1000),'min_v_mm_s':float(m['v'].min())})
with (OUT/'event_timings_INITFIX20_vs_FLUIDFRIC0.csv').open('w',newline='') as f: w=csv.DictWriter(f,fieldnames=events[0].keys()); w.writeheader(); w.writerows(events)
energies=energy(P,'INITFIX20')+energy(C,'FLUIDFRIC0')
with (OUT/'energy_work_INITFIX20_vs_FLUIDFRIC0.csv').open('w',newline='') as f: w=csv.DictWriter(f,fieldnames=energies[0].keys()); w.writeheader(); w.writerows(energies)
summary={'parent':'F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX20','child':'F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX_FLUIDFRIC0_20','common_timestamp_count':len(common),'max_abs_displacement_difference_mm':float(abs(arr).max()),'rms_displacement_difference_mm':float(np.sqrt(np.mean(arr**2))),'max_abs_velocity_difference_mm_s':float(abs(vel).max()),'rms_velocity_difference_mm_s':float(np.sqrt(np.mean(vel**2))),'earliest_abs_displacement_gt_0p01mm_ms':float(common[early[0]]['time_s']*1000) if len(early) else None,'cycles':cycles,'events':events,'energy':energies,'startup_0_8ms_max_abs_displacement_difference_mm':float(abs(arr[np.array([r['time_s'] for r in common])<=.008]).max())}
(OUT/'fluidfric0_closure_summary.json').write_text(json.dumps(summary,indent=2)+'\n'); print(json.dumps(summary,indent=2))
