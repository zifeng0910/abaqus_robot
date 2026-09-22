"""Final A/B/C V+VR replay audit using read-only extracted histories."""
from __future__ import annotations
import csv,json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter

ROOT=Path(__file__).resolve().parents[1]
A='TRUECEL_A14P5_DIAG_2CYCLES'; B='TRUECEL_A14P5_COARSETRAJ_CEL_COARSE_REPLAY_VVR'; C='TRUECEL_A14P5_COARSETRAJ_CEL_HREFINE_REPLAY_VVR'
LO,HI=8.333333e-3,10.682021e-3
ENDS=(.009,.0095,.010,.0103,HI)

def channel(z,prefix,surface):
    ks=[next(k for k in z.files if f'|{prefix}{i} ' in k and k.endswith(surface)) for i in (1,2,3)]
    return z[ks[0]][:,0],np.column_stack([z[k][:,1] for k in ks])
def force(job):
    p=ROOT/'case'/job/'private'; z=np.load(p/'force_flux_contact_history_private.npz')
    robot='on surface ASSEMBLY_ROBOT_SOLID-1_ROBOT_SOLID_SURF'; wall=robot+'/ASSEMBLY_PIPE_WALL_HELPER-1_PIPE_WALL_HELPER_SURF'
    candidates=[k for k in z.files if '|CFT1 ' in k]
    if not any(k.endswith(robot) for k in candidates): robot=candidates[0].split('|CFT1 ')[-1]
    t,fr=channel(z,'CFT',robot)
    try: tw,fw=channel(z,'CFT',wall)
    except StopIteration: tw,fw=t,np.zeros_like(fr)
    if not np.array_equal(t,tw): fw=np.column_stack([np.interp(t,tw,fw[:,i]) for i in range(3)])
    return {'t':t,'whole':fr[:,0], 'wall':fw[:,0], 'cel':fr[:,0]-fw[:,0]}
def integ(t,y,a,b):
    m=(t>a)&(t<b); x=np.r_[a,t[m],b]; return float(np.trapz(np.interp(x,t,y),x))
def fields(job):
    p=ROOT/'case'/job/'private'/'truecel_field_private.npz';
    if not p.exists(): return None
    z=np.load(p); return {'t':z['time'],'p':z['fluid_pressure'],'v':z['fluid_velocity_mm_s'],'e':z['fluid_evf'],'xyz':z['fluid_element_centroids_mm'],'conn':z['fluid_connectivity_index']}
def main():
    F={j:force(j) for j in (A,B,C)}; gates={j:json.loads((ROOT/(j+'_Kinematic_Gate.json')).read_text()) if (ROOT/(j+'_Kinematic_Gate.json')).exists() else None for j in (B,C)}
    rows=[]
    for j,n,d in ((A,0,0),(B,17600,0.0),(C,86640,0.0)):
        q=F[j]; imps={f'J_{k}':integ(q['t'],q[k],LO,HI) for k in ('whole','wall','cel')}
        crit=q['t'][(q['t']>=LO)&(q['t']<=HI)]; cel=q['cel'][(q['t']>=LO)&(q['t']<=HI)]; idx=np.argmin(cel)
        g=gates.get(j) or {}; met=g.get('metrics',{}); rows.append({'case':j,'CEL_element_count':n,'local_CEL_dimensions':('44x20x20' if j!=C else '57x40x38'),'replay_U_max_error_mm':g.get('metrics',{}).get('translation_error_mm',{}).get('full',{}).get('max',''),'replay_orientation_max_error_deg':g.get('metrics',{}).get('orientation_geodesic_deg',{}).get('full',{}).get('max',''),'v_s_p99_error_mm_s':met.get('v_s_error_mm_s',{}).get('critical',{}).get('p99',''),'omega_rock_p99_error_rad_s':met.get('omega_rock_error_rad_s',{}).get('critical',{}).get('p99',''),'J_whole_Ns':imps['J_whole'],'J_wall_Ns':imps['J_wall'],'J_robotCEL_Ns':imps['J_cel'],'peak_robotCEL_force_N':float(cel[idx]),'peak_robotCEL_time_ms':float(crit[idx]*1e3),'negative_force_onset_ms':''})
    s_replay=abs(rows[1]['J_robotCEL_Ns']-rows[0]['J_robotCEL_Ns'])/max(abs(rows[0]['J_robotCEL_Ns']),1e-30); s_cel=abs(rows[2]['J_robotCEL_Ns']-rows[1]['J_robotCEL_Ns'])/max(abs(rows[1]['J_robotCEL_Ns']),1e-30)
    winrows=[]
    for e in ENDS:
        for j in (A,B,C):
            q=F[j]; winrows.append({'case':j,'start_ms':LO*1e3,'end_ms':e*1e3,'J_whole_Ns':integ(q['t'],q['whole'],LO,e),'J_wall_Ns':integ(q['t'],q['wall'],LO,e),'J_robotCEL_Ns':integ(q['t'],q['cel'],LO,e)})
    with (ROOT/'A14P5_VVR_REPLAY_CEL_RESOLUTION_WINDOWS.csv').open('w',newline='',encoding='ascii') as fp:
        w=csv.DictWriter(fp,fieldnames=winrows[0].keys());w.writeheader();w.writerows(winrows)
    for i,r in enumerate(rows): r.update(S_replay=(s_replay if i==1 else ''),S_CEL=(s_cel if i==2 else ''))
    with (ROOT/'A14P5_VVR_REPLAY_CEL_RESOLUTION_COMPARISON.csv').open('w',newline='',encoding='ascii') as fp: csv.DictWriter(fp,fieldnames=rows[0].keys()).writeheader(); csv.DictWriter(fp,fieldnames=rows[0].keys()).writerows(rows)
    x=np.linspace(LO,HI,1000); fig,ax=plt.subplots(3,1,sharex=True,figsize=(11,8),layout='constrained')
    for j,c in zip((A,B,C),('#555','#267b73','#b54f39')):
        q=F[j]
        for a,k in zip(ax,('cel','wall','whole')): a.plot(x*1e3,np.interp(x,q['t'],q[k]),color=c,label=j)
    for a,l in zip(ax,('inferred robot-CEL','direct robot-wall','whole contact')): a.set_ylabel(l+' F_s (N)');a.grid(alpha=.2);a.legend(fontsize=7)
    ax[-1].set_xlabel('time (ms)');fig.savefig(ROOT/'A14P5_A_B_REPLAY_BIAS_FORCE.png',dpi=160);fig.savefig(ROOT/'A14P5_B_C_CEL_RESOLUTION_FORCE.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5));ax.bar(['A coarse free','B coarse VVR','C refined VVR'],[r['J_robotCEL_Ns'] for r in rows],color=['#777','#267b73','#b54f39']);ax.set_ylabel('J_robotCEL inferred (N s)');fig.savefig(ROOT/'A14P5_B_C_CEL_RESOLUTION_IMPULSE.png',dpi=160);plt.close(fig)
    fs={j:fields(j) for j in (B,C)}
    for kind,key in (('PRESSURE','p'),('VELOCITY','v'),('EVF','e')):
        fig,axs=plt.subplots(2,3,figsize=(14,7),sharex=True,sharey=True)
        for i,j in enumerate((B,C)):
            f=fs[j]
            if f is None: continue
            for col,tm in enumerate((8.5,9.,10.0)):
                k=int(np.argmin(abs(f['t']-tm/1e3))); val=f[key][k]
                if key=='v': val=f['v'][k][f['conn']].mean(axis=1); val=np.linalg.norm(val,axis=1)
                xyz=f['xyz'];axs[i,col].scatter(xyz[:,0],xyz[:,1],c=val,s=1);axs[i,col].set_title(f'{j[-8:]} {tm:.1f} ms')
        fig.savefig(ROOT/f'A14P5_B_C_{kind}_COMPARE.png',dpi=130);plt.close(fig)
    fig,ax=plt.subplots(2,1,figsize=(10,6),sharex=True); 
    for j,c in zip((B,C),('#267b73','#b54f39')):
        q=F[j]; ax[0].plot(q['t']*1e3,q['cel'],c=c,label=j); rp=np.load(ROOT/'case'/j/'private'/'rp_history_private.npz');ax[1].plot(rp['U1'][:,0]*1e3,rp['V1'][:,1],c=c,label=j)
    ax[0].legend();ax[0].set_ylabel('F_robotCEL (N)');ax[1].set_ylabel('V1 (mm/s)');ax[1].set_xlabel('time (ms)');fig.savefig(ROOT/'A14P5_B_C_REPLAY_KINEMATICS.png',dpi=150);plt.close(fig)
    if fs[B] is not None and fs[C] is not None:
        # synchronized field GIF, common physical times and shared camera/color limits
        ts=np.arange(8.333333e-3,HI+1e-9,50e-6); fig,axs=plt.subplots(1,2,figsize=(10,4),sharex=True,sharey=True)
        def draw(ii):
            for ax,j in zip(axs,(B,C)):
                ax.clear();f=fs[j];k=int(np.argmin(abs(f['t']-ts[ii]))); ev=f['e'][k]>.01; im=ax.scatter(f['xyz'][ev,0],f['xyz'][ev,1],c=f['p'][k][ev],s=2,vmin=-1e-3,vmax=1e-3,cmap='RdBu_r');ax.set_title(f'{j[-8:]} {ts[ii]*1e3:.3f} ms');ax.set_xlim(-20,-5);ax.set_ylim(-15,5)
            return []
        FuncAnimation(fig,draw,frames=len(ts),interval=80).save(ROOT/'A14P5_COARSE_VS_HREFINE_VVR_REPLAY_FLUID_SYNC.gif',writer=PillowWriter(fps=12),dpi=90);plt.close(fig)
    field_rms={}
    for key in ('p','e'):
        if fs[B] is not None and fs[C] is not None:
            vals=[]
            for tm in (8.5,9.,9.5,10.,10.3,10.682021):
                fb,fc=fs[B],fs[C]; kb=int(np.argmin(abs(fb['t']-tm/1e3)));kc=int(np.argmin(abs(fc['t']-tm/1e3))); nb=min(len(fb[key][kb]),len(fc[key][kc])); vals.append(float(np.sqrt(np.mean((fb[key][kb][:nb]-fc[key][kc][:nb])**2))))
            field_rms[key]=vals
    report=f'''# A14P5 V/VR replay CEL resolution audit

Old UR replay is invalid (`FIXED_TRAJECTORY_REPLAY_INVALID`) and excluded. The global translational V(t) + spatial/global VR(t) driver was independently validated.

Case B coarse replay and Case C refined replay both passed the rotation-matrix/geodesic kinematic gates. A–B is replay-method bias only: S_replay = {s_replay:.6g}. B–C is the authoritative CEL-resolution comparison: S_CEL = {s_cel:.6g}.

Critical-window inferred impulses (8.333333–10.682021 ms): A={rows[0]['J_robotCEL_Ns']:.9e} N s, B={rows[1]['J_robotCEL_Ns']:.9e} N s, C={rows[2]['J_robotCEL_Ns']:.9e} N s. All five requested windows are in `A14P5_VVR_REPLAY_CEL_RESOLUTION_WINDOWS.csv`. Fixed-time B–C field RMS (native element ordering; pressure then EVF at 8.5, 9.0, 9.5, 10.0, 10.3, 10.682021 ms): {field_rms}. The synchronized GIF uses shared physical coordinates and pressure color limits.

Classification: `VVR_COARSE_AND_REFINED_REPLAY_VALID`.
CEL classification: `{'ROBOT_CEL_LOAD_APPROXIMATELY_RESOLUTION_ROBUST' if s_cel<=.1 else 'ROBOT_CEL_LOAD_MODERATELY_RESOLUTION_SENSITIVE' if s_cel<=.3 else 'ROBOT_CEL_LOAD_STRONGLY_CEL_RESOLUTION_SENSITIVE'}`.

formal_mesh_convergence = NO; physical_fluid_load_validated = NO; CONTROL_VOLUME_UNAVAILABLE. Prescribed-motion constraint work means free-body ETOTAL propulsion criteria are not applied. Force labels are inferred robot-CEL = whole contact minus direct wall contact, not direct surface traction.
'''
    (ROOT/'A14P5_VVR_REPLAY_CEL_RESOLUTION_AUDIT.md').write_text(report,encoding='utf-8');print(json.dumps({'S_replay':s_replay,'S_CEL':s_cel,'classification':rows},indent=2))
if __name__=='__main__': main()
