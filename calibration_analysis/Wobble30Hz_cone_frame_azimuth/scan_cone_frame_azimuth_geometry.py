"""Zero-cost rigid cone-frame azimuth geometry proxy on G6/L45."""
from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import sys

ROOT=Path(__file__).resolve().parent
JD=ROOT/'abaqus_robot/calibration_analysis/WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003'
rf=pd.read_csv(JD/'robot_frames.csv'); tel=pd.read_csv(JD/'WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003_telemetry.csv')
wn=pd.read_csv(ROOT/'Wobble_F30_PhaseDirectionalAntiImpact_WallOn_Free_003_pipe_wall_nodes_exact.csv')
wm={int(r.node):np.array([r.x_mm,r.y_mm,r.z_mm],float) for r in wn.itertuples()}
wt=pd.read_csv(ROOT/'Wobble_F30_PhaseDirectionalAntiImpact_WallOn_Free_003_pipe_wall_triangles_exact.csv')
tri=[]
for r in wt.itertuples():
    q=[int(r.n1),int(r.n2),int(r.n3),int(r.n4)]
    tri += [[wm[q[0]],wm[q[1]],wm[q[2]]],[wm[q[0]],wm[q[2]],wm[q[3]]]]
tri=np.asarray(tri,float); tree=cKDTree(tri.mean(axis=1)); tn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]); tn/=np.maximum(np.linalg.norm(tn,axis=1)[:,None],1e-30)
curve=pd.read_csv(ROOT/'curvenew_CEL_xyrot56_exact_abaqus.csv'); C=curve[['x_mm','y_mm','z_mm']].to_numpy(float); S=curve.arclength_mm.to_numpy(float)
tc=tri.mean(axis=1)
for j,x in enumerate(tc):
    q=C[np.argmin(np.sum((C-x)**2,axis=1))]
    if np.dot(tn[j],q-x)<0: tn[j]*=-1.0
def U(v): return np.asarray(v,float)/max(float(np.linalg.norm(v)),1e-15)
def cp(s): s=float(np.clip(s,S[0],S[-1])); return np.array([np.interp(s,S,C[:,j]) for j in range(3)])
def tt(s): return U(cp(min(S[-1],s+.01))-cp(max(S[0],s-.01)))
def R(v,k,a): return v*math.cos(a)+np.cross(k,v)*math.sin(a)+k*np.dot(k,v)*(1-math.cos(a))
def ptri(p,T):
    a,b,c=T[:,0],T[:,1],T[:,2]; ab=b-a; ac=c-a; ap=p-a; d1=np.einsum('ij,ij->i',ab,ap); d2=np.einsum('ij,ij->i',ac,ap); out=np.empty_like(a); done=np.zeros(len(T),bool)
    m=(d1<=0)&(d2<=0); out[m]=a[m]; done|=m; bp=p-b; d3=np.einsum('ij,ij->i',ab,bp); d4=np.einsum('ij,ij->i',ac,bp); m=(~done)&(d3>=0)&(d4<=d3); out[m]=b[m]; done|=m
    vc=d1*d4-d3*d2; m=(~done)&(vc<=0)&(d1>=0)&(d3<=0); z=d1/(d1-d3+1e-30); out[m]=a[m]+z[m,None]*ab[m]; done|=m; cpv=p-c; d5=np.einsum('ij,ij->i',ab,cpv); d6=np.einsum('ij,ij->i',ac,cpv); m=(~done)&(d6>=0)&(d5<=d6); out[m]=c[m]; done|=m
    vb=d5*d2-d1*d6; m=(~done)&(vb<=0)&(d2>=0)&(d6<=0); z=d2/(d2-d6+1e-30); out[m]=a[m]+z[m,None]*ac[m]; done|=m; va=d3*d6-d5*d4; m=(~done)&(va<=0)&((d4-d3)>=0)&((d5-d6)>=0); z=(d4-d3)/((d4-d3)+(d5-d6)+1e-30); out[m]=b[m]+z[m,None]*(c[m]-b[m]); done|=m
    m=~done; den=va[m]+vb[m]+vc[m]+1e-30; v=vb[m]/den; z=vc[m]/den; out[m]=a[m]+ab[m]*v[:,None]+ac[m]*z[:,None]; return out
def gaps(P):
    ids=np.atleast_2d(tree.query(P,k=min(48,len(tri)))[1]); out=[]
    for i,p in enumerate(P):
        ii=np.atleast_1d(ids[i]); q=ptri(p,tri[ii]); d=np.linalg.norm(q-p,axis=1); j=int(np.argmin(d)); out.append((q[j],d[j],int(ii[j])))
    return out
ids=sorted(rf.node.unique()); f0=rf[rf.frame==rf.frame.min()].set_index('node').loc[ids][['x_mm','y_mm','z_mm']].to_numpy(float); c0=f0.mean(axis=0); ax0=tt(0.0); body=(f0-c0)@ax0; regions=np.where(body>.35,'HEAD',np.where(body<-.35,'TAIL','SIDE'))
frs=sorted(rf.frame.unique()); P=np.asarray([rf[rf.frame==f].set_index('node').loc[ids][['x_mm','y_mm','z_mm']].to_numpy(float) for f in frs]); times=np.asarray([float(rf[rf.frame==f].time_s.iloc[0]) for f in frs]); com=P.mean(axis=1); sr=np.interp(times,tel.t_s,tel.robot_arc_mm)
def run(chi):
    X=[]
    for i in range(len(P)):
        k=tt(sr[i]); X.append(com[i]+np.asarray([R(v-com[i],k,math.radians(chi)) for v in P[i]]))
    X=np.asarray(X); V=np.gradient(X,times,axis=0); rows=[]
    for i in range(len(X)):
        gg=gaps(X[i]); sg=np.asarray([float(np.dot(X[i,j]-gg[j][0],tn[gg[j][2]])) for j in range(len(ids))]); close=np.maximum(-np.asarray([np.dot(V[i,j],tn[gg[j][2]]) for j in range(len(ids))]),0); rows.append(dict(frame=frs[i],time_s=times[i],chi_deg=chi,min_gap_mm=sg.min(),p05_gap_mm=np.percentile(sg,5),low10_mean_gap_mm=np.mean(np.sort(sg)[:max(1,len(sg)//10)]),frac_gap_lt0=np.mean(sg<0),frac_gap_lt20um=np.mean(sg<.020),max_closing_mm_s=close.max(),head_min_gap_mm=sg[regions=='HEAD'].min(),tail_min_gap_mm=sg[regions=='TAIL'].min(),side_min_gap_mm=sg[regions=='SIDE'].min(),head_max_closing_mm_s=close[regions=='HEAD'].max(),tail_max_closing_mm_s=close[regions=='TAIL'].max(),side_max_closing_mm_s=close[regions=='SIDE'].max(),worst_node=ids[int(np.argmin(sg))],worst_tri=gg[int(np.argmin(sg))][2]))
    return pd.DataFrame(rows)
chi_values=([0]+list(range(int(sys.argv[1]),int(sys.argv[2])+1,int(sys.argv[3]))) if len(sys.argv)>=4 else list(range(0,360,5)))
allm=pd.concat([run(x) for x in chi_values],ignore_index=True); danger=allm[(allm.time_s>=.0007)&(allm.time_s<=.0020)]
def ag(g):
    return pd.Series(dict(chi_deg=g.chi_deg.iloc[0],min_gap_mm=g.min_gap_mm.min(),p05_gap_mm=g.p05_gap_mm.min(),low10_mean_gap_mm=g.low10_mean_gap_mm.min(),frac_frames_gap_lt0=np.mean(g.min_gap_mm<0),frac_frames_gap_lt20um=np.mean(g.min_gap_mm<.020),max_closing_mm_s=g.max_closing_mm_s.max(),head_min_gap_mm=g.head_min_gap_mm.min(),tail_min_gap_mm=g.tail_min_gap_mm.min(),side_min_gap_mm=g.side_min_gap_mm.min(),head_max_closing_mm_s=g.head_max_closing_mm_s.max(),tail_max_closing_mm_s=g.tail_max_closing_mm_s.max(),side_max_closing_mm_s=g.side_max_closing_mm_s.max(),danger_frames=len(g)))
scan=danger.groupby('chi_deg',sort=False).apply(ag).reset_index(drop=True); b=scan[scan.chi_deg==0].iloc[0]; scan['min_gap_improve_um']=(scan.min_gap_mm-b.min_gap_mm)*1000; scan['negative_exposure_reduction_pct']=100*(b.frac_frames_gap_lt0-scan.frac_frames_gap_lt0)/max(b.frac_frames_gap_lt0,1e-12); scan['closing_reduction_pct']=100*(b.max_closing_mm_s-scan.max_closing_mm_s)/max(b.max_closing_mm_s,1e-12); scan['gate_pass']=((scan.min_gap_improve_um>=50)|(scan.negative_exposure_reduction_pct>=50)|(scan.closing_reduction_pct>=30))&(scan.head_min_gap_mm>=b.head_min_gap_mm-.001)&(scan.tail_min_gap_mm>=b.tail_min_gap_mm-.001); scan.to_csv(ROOT/'cone_frame_azimuth_geometry_scan.csv',index=False)
rank=scan.sort_values(['gate_pass','min_gap_mm','p05_gap_mm','low10_mean_gap_mm','max_closing_mm_s'],ascending=[False,False,False,False,True]).reset_index(drop=True); rank['rank']=np.arange(1,len(rank)+1); rank.to_csv(ROOT/'cone_frame_azimuth_rank.csv',index=False); best=float(rank.iloc[0].chi_deg); allm[allm.chi_deg.isin([0,best])].to_csv(ROOT/'chi_candidate_exact_gap_history.csv',index=False); rank[rank.chi_deg.isin([0,best])][['chi_deg','head_min_gap_mm','tail_min_gap_mm','side_min_gap_mm','head_max_closing_mm_s','tail_max_closing_mm_s','side_max_closing_mm_s']].to_csv(ROOT/'chi_candidate_head_tail_clearance.csv',index=False); rank[['chi_deg','max_closing_mm_s','head_max_closing_mm_s','tail_max_closing_mm_s','side_max_closing_mm_s']].to_csv(ROOT/'chi_candidate_closing_speed_proxy.csv',index=False)
T=tt(sr[0]); N=U(np.array([0.,0.,1.])-T[2]*T); B=U(np.cross(T,N)); beta=math.radians(40); cc=U(math.cos(beta)*T+math.sin(beta)*N); e1=U(N-np.dot(N,cc)*cc); e2=U(np.cross(cc,e1)); defs=[]
for x in [0,best]:
    c=R(cc,T,math.radians(x)); q=R(e1,T,math.radians(x)); z=R(e2,T,math.radians(x)); defs.append(dict(chi_deg=x,t_x=T[0],t_y=T[1],t_z=T[2],c_x=c[0],c_y=c[1],c_z=c[2],e1_x=q[0],e1_y=q[1],e1_z=q[2],e2_x=z[0],e2_y=z[1],e2_z=z[2],c_dot_t=np.dot(c,T),c_dot_e1=np.dot(c,q),c_dot_e2=np.dot(c,z),e1_dot_e2=np.dot(q,z),bias_deg=math.degrees(math.acos(np.clip(np.dot(c,T),-1,1))),handedness=np.dot(np.cross(q,z),c)))
pd.DataFrame(defs).to_csv(ROOT/'cone_frame_current_definition.csv',index=False)
for y,n,lab in [('min_gap_mm','clearance_vs_chi.png','minimum signed gap (mm)'),('max_closing_mm_s','closing_speed_vs_chi.png','max outward closing proxy (mm/s)')]:
    plt.figure(figsize=(8,4)); plt.plot(scan.chi_deg,scan[y],'-o',ms=2); plt.axvline(best,color='r',ls='--',label=f'best χ={best:.0f}°'); plt.xlabel('χ (deg)'); plt.ylabel(lab); plt.grid(alpha=.3); plt.legend(); plt.tight_layout(); plt.savefig(ROOT/n,dpi=160); plt.close()
fig,ax=plt.subplots(figsize=(6,6)); p=f0-c0
for x,col,lab in [(0,'k','χ=0'),(best,'r',f'χ={best:.0f}°')]:
    q=np.asarray([R(v,T,math.radians(x)) for v in p]); ax.scatter(q[:,1],q[:,2],s=8,alpha=.5,c=col,label=lab)
ax.set_xlabel('global y relative COM (mm)'); ax.set_ylabel('global z relative COM (mm)'); ax.set_aspect('equal'); ax.grid(alpha=.3); ax.legend(); fig.tight_layout(); fig.savefig(ROOT/'baseline_vs_candidate_cross_section_orbit.png',dpi=180); plt.close(fig)
print('baseline danger:',b.to_dict()); print('best:',rank.iloc[0].to_dict()); print(rank.head(8)[['rank','chi_deg','min_gap_mm','min_gap_improve_um','closing_reduction_pct','gate_pass']].to_string(index=False))
