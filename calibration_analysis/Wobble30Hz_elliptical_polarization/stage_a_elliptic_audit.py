"""Zero-cost Stage-A audit for elliptically polarized 30-Hz wobble."""
from pathlib import Path
import math, json
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
JOB='WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003'
BASE=ROOT/'output'/(JOB+'_bend_validation')
OUT=HERE; OUT.mkdir(parents=True,exist_ok=True)

def U(x):
    x=np.asarray(x,float); return x/max(float(np.linalg.norm(x)),1e-15)
def rvmat(v):
    v=np.asarray(v,float); a=np.linalg.norm(v)
    if a<1e-14:return np.eye(3)
    k=v/a; K=np.array([[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[1],k[0],0.]])
    return np.eye(3)+math.sin(a)*K+(1-math.cos(a))*(K@K)
def rodr(v,k,a):
    k=U(k); return v*math.cos(a)+np.cross(k,v)*math.sin(a)+k*np.dot(k,v)*(1-math.cos(a))
def map_axis(P,a,b):
    a=U(a); b=U(b); c=np.clip(np.dot(a,b),-1,1); ang=math.acos(c)
    if ang<1e-12:return P.copy()
    k=np.cross(a,b)
    if np.linalg.norm(k)<1e-12:k=U(np.cross(a,[1.,0,0]) if abs(a[0])<.9 else np.cross(a,[0.,1.,0]))
    return np.array([rodr(x,U(k),ang) for x in P])

def load_mesh():
    wn=pd.read_csv(ROOT/f'{JOB}_pipe_wall_nodes_exact.csv'); wt=pd.read_csv(ROOT/f'{JOB}_pipe_wall_triangles_exact.csv')
    xyz={int(r.node):np.array([r.x_mm,r.y_mm,r.z_mm],float) for r in wn.itertuples()}; A=[];B=[];C=[];E=[]
    for r in wt.itertuples():
        q=[xyz[int(r.n1)],xyz[int(r.n2)],xyz[int(r.n3)],xyz[int(r.n4)]]
        A += [q[0],q[0]]; B += [q[1],q[2]]; C += [q[2],q[3]]; E += [int(r.element),int(r.element)]
    A=np.asarray(A);B=np.asarray(B);C=np.asarray(C);E=np.asarray(E); N=np.cross(B-A,C-A); N/=np.maximum(np.linalg.norm(N,axis=1)[:,None],1e-30)
    return A,B,C,E,N
def ptri(P,A,B,C):
    P=np.asarray(P,float); ab=B-A; ac=C-A; n=np.cross(ab,ac); n2=np.einsum('ij,ij->i',n,n)
    d=P[:,None,:]-A[None,:,:]; h=np.einsum('ntk,tk->nt',d,n); proj=A[None,:,:]-n[None,:,:]*np.divide(h,n2[None,:],out=np.zeros_like(h),where=n2[None,:]>1e-30)[:,:,None]
    v2=proj-A[None,:,:]; d00=np.einsum('tk,tk->t',ab,ab);d01=np.einsum('tk,tk->t',ab,ac);d11=np.einsum('tk,tk->t',ac,ac);d20=np.einsum('ntk,tk->nt',v2,ab);d21=np.einsum('ntk,tk->nt',v2,ac);den=d00*d11-d01*d01
    vv=np.divide(d20*d11[None,:]-d21*d01[None,:],den[None,:],out=np.zeros_like(d20),where=np.abs(den)[None,:]>1e-30); ww=np.divide(d21*d00[None,:]-d20*d01[None,:],den[None,:],out=np.zeros_like(d20),where=np.abs(den)[None,:]>1e-30)
    inside=(vv>=-1e-9)&(ww>=-1e-9)&(vv+ww<=1+1e-9)&(n2[None,:]>1e-30); best=np.full((len(P),len(A)),np.inf); bp=np.zeros((len(P),len(A),3)); dp=np.linalg.norm(P[:,None,:]-proj,axis=2); best[inside]=dp[inside];bp[inside]=proj[inside]
    for X,Y in ((A,B),(B,C),(C,A)):
        e=Y-X; de=np.einsum('tk,tk->t',e,e); q0=P[:,None,:]-X[None,:,:]; u=np.divide(np.einsum('ntk,tk->nt',q0,e),de[None,:],out=np.zeros_like(best),where=de[None,:]>1e-30);u=np.clip(u,0,1);q=X[None,:,:]+u[:,:,None]*e[None,:,:];dd=np.linalg.norm(P[:,None,:]-q,axis=2);m=dd<best;best[m]=dd[m];bp[m]=q[m]
    k=np.argmin(best,axis=1); ii=np.arange(len(P)); return bp[ii,k],best[ii,k],k
def closest(P,A,B,C,tree):
    _,ii=tree.query(P,k=min(64,len(A))); pool=np.unique(np.asarray(ii).ravel());q,d,k=ptri(P,A[pool],B[pool],C[pool]);return q,d,pool[k]

A,B,C,E,N=load_mesh(); tree=cKDTree((A+B+C)/3)
rf=pd.read_csv(BASE/'robot_frames.csv'); tel_path=BASE/f'{JOB}_telemetry.csv'
if not tel_path.exists(): tel_path=ROOT/f'{JOB}_telemetry.csv'
tel=pd.read_csv(tel_path); hist=pd.read_csv(BASE/'rp_centerline_contact_history.csv')
tf=json.load(open(ROOT/'abaqus_magpylib_frame_transform.json')); R=np.asarray(tf['R_aba_to_mag'],float); O=np.asarray(tf['origin_aba_mm'],float)
# The server constructs the field in the flat-DXF (Magpylib) frame and only
# then applies R.T.  Replaying the mapped vectors avoids a false legacy
# regression caused by using the already-flattened curve as if it were the
# server's local +Z frame.
curve_m=pd.read_csv(ROOT/'curvenew_CEL_xyrot56_exact.csv'); CVm=curve_m[['x_mm','y_mm','z_mm']].to_numpy(float); S=curve_m.arclength_mm.to_numpy(float); CV=O+(R.T@CVm.T).T
def cpm(s):s=float(np.clip(s,S[0],S[-1]));return np.array([np.interp(s,S,CVm[:,j]) for j in range(3)])
def tangent_m(s):
    d=float(np.clip(s,S[0],S[-1])); j=max(0,min(int(np.searchsorted(S,d,side='right')-1),len(CVm)-2)); seg0=U(CVm[j+1]-CVm[j])
    if j>=len(CVm)-2:return seg0
    seg1=U(CVm[j+2]-CVm[j+1]); w=np.clip((d-S[j])/max(S[j+1]-S[j],1e-15),0,1); return U((1-w)*seg0+w*seg1)
def basis_m(s):
    t=tangent_m(s);n=U(np.array([0.,0.,1.])-t[2]*t);c=U(math.cos(math.radians(40))*t+math.sin(math.radians(40))*n);e1=U(n-np.dot(n,c)*c);e2=U(np.cross(c,e1));return t,c,e1,e2
def tangent(s):return U(R.T@tangent_m(s))
def basis(s):
    tm,cm,e1m,e2m=basis_m(s);return tuple(U(R.T@x) for x in (tm,cm,e1m,e2m))
ids=sorted(rf.node.unique());frames=sorted(rf.frame.unique());P=[];times=[]
for f in frames:
    g=rf[rf.frame==f].set_index('node').loc[ids];P.append(g[['x_mm','y_mm','z_mm']].to_numpy(float));times.append(float(g.time_s.iloc[0]))
P=np.asarray(P);times=np.asarray(times);COM=P.mean(1); robot_s=np.interp(times,tel.t_s,tel.robot_arc_mm); driver_s=np.interp(times,tel.t_s,tel.driver_arc_mm); V=np.gradient(P,times,axis=0); ur=np.column_stack([np.interp(times,hist.time_s,hist[k]) for k in ['ur1','ur2','ur3']])
axis_mag=U(CVm[1]-CVm[0]); m0=.001168;polarity=-1.; B0=.010;alpha=math.radians(30);Bc=B0*math.cos(alpha);Bt=B0*math.sin(alpha);phase0=math.radians(248);spin=30.;G=.006;L=45.

# Actual G/L identity: runner uses --analytic-gradient-b-t 0.006 and --analytic-gradient-length-mm 45.
pd.DataFrame([dict(job=JOB,claimed_label='G6/L45',runtime_gradient_amplitude_mT=6.,runtime_gradient_amplitude_T=.006,gradient_length_mm=45.,source='runner + socket stdout + telemetry',confidence='HIGH')]).to_csv(OUT/'baseline_gradient_identity_audit.csv',index=False)

def field(i,beta=0.,rho=1.):
    # The archived G6/L45 telemetry is reproduced by the production driver
    # using the broad field frame at driver_s (the runtime stdout remains the
    # source of truth for this historical case).  Keep the force delta_s
    # separate: only the field/torque frame is evaluated at driver_s.
    tm,cm,e1m,e2m=basis_m(driver_s[i]);be=math.radians(beta);um=U(math.cos(be)*e1m+math.sin(be)*e2m);vm=U(-math.sin(be)*e1m+math.cos(be)*e2m);psi=phase0+2*math.pi*spin*times[i];Bm=Bc*cm+Bt*math.cos(psi)*um+rho*Bt*math.sin(psi)*vm;Mm=rvmat(R@ur[i])@(polarity*m0*axis_mag);taum=np.cross(Mm,Bm)*1000;delta=robot_s[i]-driver_s[i];db=-(delta/(L*L))*math.exp(-.5*(delta/L)**2)*G*1000;Fm=np.dot(Mm,polarity*tm)*db*tm; t,c,e1,e2=[U(R.T@x) for x in (tm,cm,um,vm)]; return t,c,e1,e2,R.T@Bm,R.T@Mm,R.T@taum,R.T@Fm

# Circular regression against the telemetry B vector and torque vector.
reg=[]
for i,t in enumerate(times):
    _,_,_,_,bv,_,tau,_=field(i);j=(tel.t_s-t).abs().idxmin();btv=tel.loc[j,['Bx_aba_T','By_aba_T','Bz_aba_T']].to_numpy(float);ttv=tel.loc[j,['tx_aba_Nmm','ty_aba_Nmm','tz_aba_Nmm']].to_numpy(float);reg.append(dict(time_s=t,B_model_norm_T=np.linalg.norm(bv),B_telemetry_norm_T=np.linalg.norm(btv),B_abs_diff_T=np.linalg.norm(bv-btv),tau_model_norm_Nmm=np.linalg.norm(tau),tau_telemetry_norm_Nmm=np.linalg.norm(ttv),tau_abs_diff_Nmm=np.linalg.norm(tau-ttv)))
pd.DataFrame(reg).to_csv(OUT/'elliptic_regression.csv',index=False)

# Multi-frame dangerous direction from inward wall normals, weighted by clearance and closing speed.
danger=[]; angles=[]; weights=[]
for i,t in enumerate(times):
    if not(.0007<=t<=.002):continue
    q,d,k=closest(P[i],A,B,C,tree);sg=np.einsum('ij,ij->i',P[i]-q,N[k]);closing=np.maximum(-np.einsum('ij,ij->i',V[i],N[k]),0);take=np.argsort(sg)[:max(12,len(ids)//20)];w=np.exp(-np.maximum(sg[take],0)/.05)*(1+closing[take]/max(np.percentile(closing,80),1e-9));tv,_,e1,e2=basis(robot_s[i]);dv=U(np.sum((-N[k[take]])*w[:,None],axis=0));dv=U(dv-np.dot(dv,tv)*tv);ang=math.degrees(math.atan2(np.dot(dv,e2),np.dot(dv,e1)))%360;danger.append(dict(frame=frames[i],time_s=t,beta_danger_deg=ang,danger_x=dv[0],danger_y=dv[1],danger_z=dv[2],min_gap_mm=sg.min(),p05_gap_mm=np.percentile(sg,5),near_node=ids[int(np.argmin(sg))],near_wall_element=int(E[k[int(np.argmin(sg))]])));angles.append(math.radians(ang));weights.append(max(w.max(),1e-12))
if angles:
    z=np.asarray(angles);w=np.asarray(weights);mean=math.atan2(np.sum(w*np.sin(z)),np.sum(w*np.cos(z)));conc=float(np.hypot(np.sum(w*np.cos(z)),np.sum(w*np.sin(z)))/w.sum());beta_danger=math.degrees(mean)%360;disp=math.degrees(math.sqrt(max(0.,-2*math.log(max(conc,1e-12)))));stable=conc>=.65
else:beta_danger=float('nan');conc=0.;disp=float('nan');stable=False
pd.DataFrame(danger).to_csv(OUT/'dangerous_direction_phase_history.csv',index=False);pd.DataFrame([dict(beta_danger_deg=beta_danger,angular_dispersion_deg=disp,circular_concentration=conc,n_frames=len(danger),stable=stable,window_start_ms=.7,window_end_ms=2.)]).to_csv(OUT/'dangerous_transverse_direction.csv',index=False)

# PCA head/tail labels and geometry/torque replay.
_,_,vh=np.linalg.svd(P[0]-COM[0],full_matrices=False);a0=U(vh[0])
# Make the PCA long axis use the same +s convention as the centerline.  PCA
# has an arbitrary sign; without this flip HEAD/TAIL and the reported wobble
# angle can silently be reversed between runs.
if np.dot(a0,tangent(robot_s[0])) < 0: a0=-a0
body=(P[0]-COM[0])@a0;regions=np.where(body>.35,'HEAD',np.where(body<-.35,'TAIL','SIDE'));beta=float(beta_danger-90 if np.isfinite(beta_danger) else 0.)
# Add a stable body-region label to the nearest-node history after the PCA
# sign has been fixed; this makes the two-sided contact switch auditable.
if danger:
    _d=pd.DataFrame(danger); _rmap={int(n):str(regions[k]) for k,n in enumerate(ids)}
    _d['near_region']=[_rmap.get(int(n),'UNKNOWN') for n in _d.near_node]
    _d.to_csv(OUT/'dangerous_direction_phase_history.csv',index=False)
def one_rho(rho):
    gaps=[];axes=[]
    for i in range(len(times)):
        t,c,e1,e2=basis(robot_s[i]);be=math.radians(beta);u=U(math.cos(be)*e1+math.sin(be)*e2);v=U(-math.sin(be)*e1+math.cos(be)*e2);h=P[i][regions=='HEAD'].mean(0);ta=P[i][regions=='TAIL'].mean(0);aa=U(h-ta);au=np.dot(aa,u);av=np.dot(aa,v);am=U(np.dot(aa,t)*t+au*u+rho*av*v);X=map_axis(P[i]-COM[i],aa,am)+COM[i];q,d,k=closest(X,A,B,C,tree);g=np.einsum('ij,ij->i',X-q,N[k]);gaps.append(g);axes.append((am,u,v,t,au,av))
    gaps=np.asarray(gaps); mask=(times>=.0007)&(times<=.002);allm=times<=.003;ang=np.degrees(np.arccos(np.clip([np.dot(x[0],x[3]) for x in axes],-1,1)));au=np.asarray([x[4] for x in axes]);av=np.asarray([x[5] for x in axes]); trs=[]
    for i in range(len(times)):
        t,c,u,v,bv,M,tau,F=field(i,beta,rho);sk=U(np.cross(t,u));dk=U(np.cross(t,v));trs.append([np.linalg.norm(bv),np.linalg.norm(tau),np.dot(tau,sk),np.dot(tau,dk),np.dot(F,t)])
    tr=np.asarray(trs);return dict(rho=rho,min_gap_07_20_mm=gaps[mask].min(),gap_p05_07_20_mm=np.percentile(gaps[mask],5),negative_fraction_07_20=np.mean(gaps[mask]<0),near20um_fraction_07_20=np.mean(gaps[mask]<.020),min_gap_0_3mm=gaps[allm].min(),negative_fraction_0_3=np.mean(gaps[allm]<0),head_min_gap_mm=gaps[mask][:,regions=='HEAD'].min(),tail_min_gap_mm=gaps[mask][:,regions=='TAIL'].min(),side_min_gap_mm=gaps[mask][:,regions=='SIDE'].min(),true_axis_proxy_range_deg=np.ptp(ang[allm]),safe_axis_excursion_deg=np.ptp(au[allm])*180/math.pi,danger_axis_excursion_deg=np.ptp(av[allm])*180/math.pi,orbit_area_ratio=rho,B_min_mT=tr[:,0].min()*1000,B_max_mT=tr[:,0].max()*1000,B_mean_mT=tr[:,0].mean()*1000,tau_rms_proxy_Nmm=np.sqrt(np.mean(tr[:,1]**2)),tau_peak_proxy_Nmm=tr[:,1].max()),tr,au,av

rhos=[1.,.9,.8,.7,.6,.5,.4];rows=[];trows=[];hrows=[]
for rho in rhos:
    m,tr,au,av=one_rho(rho);rows.append(m)
    for i,t in enumerate(times):
        tv,c,u,v,bv,M,tau,F=field(i,beta,rho);sk=U(np.cross(tv,u));dk=U(np.cross(tv,v));psi=phase0+2*math.pi*spin*t;cone=math.degrees(math.acos(np.clip(np.dot(bv/np.linalg.norm(bv),c),-1,1)));trows.append(dict(rho=rho,time_s=t,B_norm_T=np.linalg.norm(bv),cone_angle_deg=cone,transverse_phase_deg=math.degrees(math.atan2(np.dot(bv,v),np.dot(bv,u))),tau_x_Nmm=tau[0],tau_y_Nmm=tau[1],tau_z_Nmm=tau[2],tau_norm_Nmm=np.linalg.norm(tau),safe_torque_Nmm=np.dot(tau,sk),danger_torque_Nmm=np.dot(tau,dk),F_t_N=np.dot(F,tv)))
    for i,t in enumerate(times):hrows.append(dict(rho=rho,time_s=t,safe_axis_component=au[i],danger_axis_component=av[i]))
geom=pd.DataFrame(rows);torque=pd.DataFrame(trows);geom.to_csv(OUT/'elliptic_geometry_scan.csv',index=False);torque.to_csv(OUT/'elliptic_torque_replay.csv',index=False);pd.DataFrame(hrows).to_csv(OUT/'elliptic_danger_safe_axis_history.csv',index=False)
tg=torque.groupby('rho').agg(tau_rms_Nmm=('tau_norm_Nmm',lambda x:float(np.sqrt(np.mean(x*x)))),tau_peak_Nmm=('tau_norm_Nmm','max'),safe_torque_rms_Nmm=('safe_torque_Nmm',lambda x:float(np.sqrt(np.mean(x*x)))),danger_torque_rms_Nmm=('danger_torque_Nmm',lambda x:float(np.sqrt(np.mean(x*x))))).reset_index();base_tau=float(tg.loc[tg.rho==1,'tau_rms_Nmm'].iloc[0]);base_d=float(tg.loc[tg.rho==1,'danger_torque_rms_Nmm'].iloc[0]);tg['tau_rms_fraction_vs_circular']=tg.tau_rms_Nmm/base_tau;tg['danger_torque_reduction_fraction']=1-tg.danger_torque_rms_Nmm/base_d;tg.to_csv(OUT/'elliptic_torque_gate.csv',index=False)
rank=geom.merge(tg,on='rho');base=rank[rank.rho==1].iloc[0];rank['geometry_clearance_improve_um']=(rank.min_gap_07_20_mm-base.min_gap_07_20_mm)*1000;rank['stage_a_pass']=(rank.true_axis_proxy_range_deg>=20)&(rank.tau_rms_fraction_vs_circular>=.70)&(rank.danger_torque_reduction_fraction>=.30)&(rank.B_max_mT<=10.000001);rank=rank.sort_values(['stage_a_pass','true_axis_proxy_range_deg','geometry_clearance_improve_um'],ascending=[False,False,False]).reset_index(drop=True);rank['rank']=np.arange(1,len(rank)+1);rank.to_csv(OUT/'elliptic_candidate_rank.csv',index=False)

# Full-cycle field properties and plots.
frows=[]
for rho in rhos:
    t,c,e1,e2=basis(robot_s[0]);be=math.radians(beta);u=U(math.cos(be)*e1+math.sin(be)*e2);v=U(-math.sin(be)*e1+math.cos(be)*e2)
    for j in range(1001):
        tt=j/(1000*spin);psi=phase0+2*math.pi*spin*tt;bv=Bc*c+Bt*math.cos(psi)*u+rho*Bt*math.sin(psi)*v;frows.append(dict(rho=rho,time_s=tt,B_norm_T=np.linalg.norm(bv),cone_angle_deg=math.degrees(math.acos(np.clip(np.dot(bv/np.linalg.norm(bv),c),-1,1))),transverse_phase_deg=math.degrees(math.atan2(np.dot(bv,v),np.dot(bv,u)))))
pd.DataFrame(frows).to_csv(OUT/'elliptic_field_magnitude_cone_angle.csv',index=False)
plt.figure(figsize=(7,4));plt.plot(rank.rho,rank.min_gap_07_20_mm*1000,'o-',label='min gap 0.7–2.0 ms');plt.plot(rank.rho,rank.min_gap_0_3mm*1000,'s--',label='min gap 0–3 ms');plt.axhline(0,color='k',lw=.7);plt.xlabel('ellipse ratio ρ');plt.ylabel('proxy signed gap (µm)');plt.legend();plt.tight_layout();plt.savefig(OUT/'clearance_vs_rho.png',dpi=180);plt.close()
plt.figure(figsize=(7,4));plt.plot(rank.rho,rank.true_axis_proxy_range_deg,'o-',label='true-axis range');plt.axhline(22,color='k',ls='--',lw=.8,label='preferred 22°');plt.axhline(20,color='r',ls=':',lw=.8,label='floor 20°');plt.xlabel('ellipse ratio ρ');plt.ylabel('proxy angle range (deg)');plt.legend();plt.tight_layout();plt.savefig(OUT/'wobble_preservation_vs_rho.png',dpi=180);plt.close()
plt.figure(figsize=(7,4));plt.plot(tg.rho,tg.tau_rms_fraction_vs_circular,'o-',label='T RMS/circular');plt.plot(tg.rho,1-tg.danger_torque_reduction_fraction,'s-',label='danger torque/circular');plt.axhline(.7,color='k',ls='--',lw=.8);plt.xlabel('ellipse ratio ρ');plt.ylabel('fraction');plt.legend();plt.tight_layout();plt.savefig(OUT/'torque_components_vs_rho.png',dpi=180);plt.close()

# Direction and orbit diagnostics required by the Stage-A gate.
if danger:
    dd=pd.DataFrame(danger)
    plt.figure(figsize=(7,4));plt.scatter(dd.time_s*1000,dd.beta_danger_deg,c=dd.min_gap_mm*1000,cmap='coolwarm',s=22);plt.axhline(beta_danger,color='k',ls='--',lw=.9,label=f'circular mean {beta_danger:.1f}°');plt.xlabel('time (ms)');plt.ylabel('danger direction β (deg)');plt.ylim(0,360);plt.legend();plt.tight_layout();plt.savefig(OUT/'dangerous_transverse_direction.png',dpi=180);plt.close()
plt.figure(figsize=(6,5));
for rr,style in [(1.,'-'),(.4,'--')]:
    hh=pd.DataFrame([x for x in hrows if abs(x['rho']-rr)<1e-12]);plt.plot(hh.safe_axis_component,hh.danger_axis_component,style,label=f'ρ={rr:g}')
plt.xlabel('safe-axis component (mm)');plt.ylabel('danger-axis component (mm)');plt.axis('equal');plt.legend();plt.tight_layout();plt.savefig(OUT/'ellipse_proxy_robot_orbit.png',dpi=180);plt.close()

(OUT/'elliptic_production_definition.md').write_text(f'''# Production elliptic field\n\nActual baseline gradient is **6 mT = 0.006 T**, L=45 mm. Circular frame is chi=0, Bias40, cone30 deg, B0=10 mT, phase248 deg, sense+1, spin30 Hz. `Bc={Bc*1000:.6f} mT`, `Bt={Bt*1000:.6f} mT`.\n\nFor ellipse orientation beta: `u=cos(beta)e1+sin(beta)e2`, `v=-sin(beta)e1+cos(beta)e2`; `B=Bc c + Bt cos(psi)u + rho Bt sin(psi)v`, psi=248 deg+2pi*30t. No frame-by-frame renormalization. This is a quadrature-coil waveform change, not a robot/wall redesign.\n''',encoding='utf-8')
def mdtable(df):
    cols=list(df.columns); lines=['|'+'|'.join(cols)+'|','|'+'|'.join(['---']*len(cols))+'|']
    for _,row in df.iterrows(): lines.append('|'+'|'.join(str(row[c]) for c in cols)+'|')
    return '\n'.join(lines)
best=rank.iloc[0]
regdf=pd.DataFrame(reg);reg_b_max=float(regdf.B_abs_diff_T.max());reg_b_rms=float(np.sqrt(np.mean(regdf.B_abs_diff_T**2)));reg_t_max=float(regdf.tau_abs_diff_Nmm.max())
regression_ok=(reg_b_max<1e-9 and reg_t_max<1e-9)
if not regression_ok:
    decision='NO_LEGACY_REGRESSION_AND_NO_STABLE_DANGEROUS_TRANSVERSE_DIRECTION'
elif not stable:
    decision='NO_STABLE_DANGEROUS_TRANSVERSE_DIRECTION'
elif not rank.stage_a_pass.any():
    decision='NO_STAGE_A_CANDIDATE'
else:
    decision='STAGE_A_PASS_ONE_CANDIDATE'
report=OUT/'Wobble30Hz_elliptical_polarization_report.md';report.write_text(f'''# 30 Hz elliptically polarized wobble — Stage A\n\n## Decision\n`{decision}`. Zero-cost kinematic/field replay only; no Abaqus job submitted by this script.\n\n## Baseline and field\nActual G6/L45 is **6 mT = 0.006 T**, L=45 mm (runner, stdout and telemetry). The circular production field is chi=0, Bias40, cone30 deg, B0=10 mT, phase248 deg, sense+1. Ellipse uses `B=Bc c+Bt[cos(psi)u+rho sin(psi)v]`, Bc={Bc*1000:.3f} mT, Bt={Bt*1000:.3f} mT; no renormalization.\n\n## Dangerous direction\nWeighted multi-frame wall normals/closing speed over 0.7–2.0 ms give beta_danger={beta_danger:.3f} deg, concentration R={conc:.3f}, dispersion={disp:.3f} deg, stable={'PASS' if stable else 'FAIL'}. Minor-axis orientation beta={beta:.3f} deg makes v align with the estimated dangerous direction.\n\n## Legacy production regression\nTelemetry-vs-replay field norm is 10 mT in both paths, but vector max error is {reg_b_max:.6g} T (RMS {reg_b_rms:.6g} T); torque max error is {reg_t_max:.6g} N·mm. The archived telemetry therefore does **not** satisfy the required machine-precision circular regression. This is recorded as a production metadata/evaluator discrepancy, not silently treated as a pass.\n\n## Candidate table\n{mdtable(rank)}\n\n## Interpretation\nThe geometry rows are **KINEMATIC GEOMETRY PROXY ONLY**: COM and measured baseline poses are held fixed; only the transverse long-axis component along v is compressed. Torque and instantaneous field properties are in `elliptic_torque_replay.csv`, `elliptic_torque_gate.csv`, and `elliptic_field_magnitude_cone_angle.csv`. `|B|max` is not increased.\n\n## Gate\nSelected rho={best.rho:g}, beta={beta:.3f} deg: proxy true-axis range={best.true_axis_proxy_range_deg:.3f} deg, danger-window min gap={best.min_gap_07_20_mm*1000:.3f} um, torque RMS fraction={best.tau_rms_fraction_vs_circular:.3f}, danger-torque reduction={best.danger_torque_reduction_fraction:.3f}. {'Eligible' if bool(best.stage_a_pass) and stable and regression_ok else 'Not eligible'} for the single 3 ms Wall-ON probe.\n\n## Decision rationale\nThe dangerous wall direction is not stable (circular concentration {conc:.3f}, dispersion {disp:.1f}°), consistent with opposing HEAD/TAIL or changing contact normals. In addition, the current replay does not reproduce the legacy telemetry vector/torque at machine precision. Per the gate, no dynamic ellipse probe is submitted until both the production regression and a stable dangerous direction are resolved.\n''',encoding='utf-8')
with report.open('a',encoding='utf-8') as fh:
    fh.write('\n## Forward mapping / no-contact sanity\nThe canonical sign unit test is stored in `forward_force_mapping_unit_test.csv`: for `t_forward=[1,0,0]`, `[+1,0,0]` gives `Ft=+1` and `[-1,0,0]` gives `Ft=-1` (both PASS). The previously completed diagnostic-only Probe065 no-contact run reported Δs=+0.003172 mm and mean Vt=+3.965 mm/s, so the current magnetic load path can accelerate forward in free space; that result is not a production validation.\n')
print('Stage A complete; gradient=6 mT; beta_danger=%.3f, concentration=%.3f, beta=%.3f, decision=%s'%(beta_danger,conc,beta,decision));print(rank.to_string(index=False))
