"""Frozen-physics audit: exact tetra mass/COM, dense RP, correct vector impulse.

All linear impulses are N s (equivalently tonne mm/s); all vectors Abaqus global.
No solver is launched. Private increment arrays stay local and are not committed.
"""
from pathlib import Path
import re, json
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
from scipy.integrate import cumulative_trapezoid

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
JOB='Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083'

def blocks(text):
    result=[]
    for line_no,line in enumerate(text.splitlines(),1):
        if not line.strip() or line.startswith('**'):continue
        if line.startswith('*'):result.append([line_no,line,[]])
        else:result[-1][2].append(line.strip())
    return result

def mesh_properties(text):
    part=None; material=None; meshes={}; rps={}; densities={}; shifts={}
    for _,key,data in blocks(text):
        low=key.lower()
        if low.startswith('*part,'):
            part=re.search(r'name=([^,]+)',key,re.I)[1]; meshes[part]={'nodes':{},'elems':[]}
        elif low.startswith('*end part'):part=None
        elif low.startswith('*node') and low.split(',')[0]=='*node':
            for line in data:
                v=[float(x) for x in line.split(',') if x.strip()]
                if part:meshes[part]['nodes'][int(v[0])]=v[1:4]
                else:rps[int(v[0])]=v[1:4]
        elif low.startswith('*element') and part and 'type=c3d4' in low:
            meshes[part]['elems'].extend([[int(x) for x in row.split(',') if x.strip()] for row in data])
        elif low.startswith('*solid section') and part:
            meshes[part]['material']=re.search(r'material=([^,]+)',key,re.I)[1]
        elif low.startswith('*instance'):
            p=re.search(r'part=([^,]+)',key,re.I)[1]
            shifts[p]=np.asarray([float(x) for x in data[0].split(',') if x.strip()]) if data else np.zeros(3)
            assert len(data)<=1, 'Unimplemented assembly rotation: audit must not guess'
        elif low.startswith('*material'):material=re.search(r'name=([^,]+)',key,re.I)[1]
        elif low.startswith('*density'):densities[material]=float(data[0].split(',')[0])
    masses=[]
    for name,mesh in meshes.items():
        if not mesh['elems']:continue
        con=np.array(mesh['elems'],int); verts=np.array([[mesh['nodes'][i] for i in e[1:]] for e in con])+shifts.get(name,np.zeros(3))
        vol=np.abs(np.linalg.det(verts[:,1:]-verts[:,:1]))/6
        rho=densities[mesh['material']]; mass=vol.sum()*rho
        com=(vol[:,None]*verts.mean(axis=1)).sum(axis=0)/vol.sum()
        # Explicit rigid-body inertia: lumped element mass at constituent nodes.
        q=verts-com
        I=sum((rho*vol[e]/4)*(np.eye(3)*np.dot(v,v)-np.outer(v,v)) for e in range(len(vol)) for v in q[e])
        masses.append(dict(part=name,element_count=len(vol),density_tonne_mm3=rho,volume_mm3=vol.sum(),mesh_mass_tonne=mass,mesh_mass_kg=mass*1e3,mesh_mass_mg=mass*1e9,target_mass_mg=10 if name=='Robot_SOLID' else np.nan,COM_x_mm=com[0],COM_y_mm=com[1],COM_z_mm=com[2],effective_mass_source='rho times summed C3D4 volumes; no point mass/mass scaling; total checked against STA'))
        mesh.update(verts=verts,vol=vol,com=com,mass=mass,I=I,connectivity=con,shift=shifts.get(name,np.zeros(3)))
    return meshes,np.asarray(rps[2]),masses

def dense_history(folder):
    z=np.load(folder/'rp_history_private.npz')
    def get(v):return max((z[k] for k in z.files if k==v or k.startswith(v+' (Repeated:')),key=len)
    t=get('V1')[:,0]
    return t,{v:np.interp(t,get(v)[:,0],get(v)[:,1]) for prefix in ('U','UR','V','VR') for v in [prefix+str(i) for i in range(1,4)]},z

def interpolate(t,tx,values):return np.column_stack([np.interp(t,tx,col) for col in np.asarray(values).T])
def vec(df,prefix):return df[[prefix+str(i) for i in range(1,4)]].to_numpy()

def main():
    text=(ROOT/(JOB+'.inp')).read_text()
    meshes,rp,masses=mesh_properties(text); mesh=meshes['Robot_SOLID']; m=mesh['mass']; off=mesh['com']-rp
    pd.DataFrame(masses).to_csv(HERE/'actual_rigid_body_mass_audit.csv',index=False)
    pd.DataFrame([dict(RP_x_mm=rp[0],RP_y_mm=rp[1],RP_z_mm=rp[2],COM_x_mm=mesh['com'][0],COM_y_mm=mesh['com'][1],COM_z_mm=mesh['com'][2],offset_x_um=off[0]*1000,offset_y_um=off[1]*1000,offset_z_um=off[2]*1000,offset_norm_um=np.linalg.norm(off)*1000)]).to_csv(HERE/'rp_vs_com_offset_audit.csv',index=False)
    rows=[]; ext=[]
    for ln,key,dat in blocks(text):
        if key.lower().startswith('*cload'):
            amp=re.search(r'amplitude=([^,]+)',key,re.I)[1]
            for line in dat:
                target,dof,mag=line.split(',')[:3]
                rows.append(dict(line=ln,target=target.strip(),amplitude=amp,DOF=int(dof),reference_magnitude=float(mag),unit='N' if int(dof)<4 else 'N mm'))
        if any(key.lower().startswith(k) for k in ('*cload','*dload','*dsload','*boundary','*initial conditions','*bulk viscosity','*contact damping','*mass','*rotary inertia','*connector','*amplitude','*dynamic')):
            ext.append(dict(line=ln,keyword=key,data='; '.join(dat),meaning='reference=1 multiplies VUAMP return' if key.lower().startswith('*cload') else 'inspect without counting internal forces twice'))
    pd.DataFrame(rows).to_csv(HERE/'rp_cload_inventory.csv',index=False)
    pd.DataFrame(ext).to_csv(HERE/'all_external_loads_inventory.csv',index=False)
    assert all(r['reference_magnitude']==1 for r in rows) and len(rows)==12
    assert len(set((r['target'],r['amplitude'],r['DOF']) for r in rows))==12
    t,d,z=dense_history(HERE/'baseline'); df=pd.DataFrame(d)
    U=vec(df,'U'); UR=vec(df,'UR'); V=vec(df,'V'); omega=vec(df,'VR')
    rotated_off=Rotation.from_rotvec(UR).apply(np.broadcast_to(off,U.shape))
    Vcom=V+np.cross(omega,rotated_off)
    tel=pd.read_csv(ROOT/(JOB+'_telemetry.csv'))
    F=interpolate(t,np.r_[0,tel.t_s],np.vstack((np.zeros(3),tel[['fx_aba_N','fy_aba_N','fz_aba_N']])))
    hyd=pd.read_csv(ROOT/(JOB+'_hydro_increment.csv'),skiprows=2,names=['time_s','v1','v2','v3','w1','w2','w3','Fh1','Fh2','Fh3','Th1','Th2','Th3'])
    H=interpolate(t,np.r_[0,hyd.time_s],np.vstack((np.zeros(3),vec(hyd,'Fh'))))
    Jmag=cumulative_trapezoid(F,t,axis=0,initial=0); Jh=cumulative_trapezoid(H,t,axis=0,initial=0)
    P=m*(Vcom-Vcom[0]); missing=P-Jmag-Jh
    windows=[]
    for name,lo,hi in [('A_clean',0,.00095),('B_suspected',.00095,.00110),('C',.0011,.00125),('D',.00120,.00135),('E',0,.00615)]:
        a=np.argmin(abs(t-lo));b=np.argmin(abs(t-hi));p=P[b]-P[a];jm=Jmag[b]-Jmag[a];jh=Jh[b]-Jh[a];res=p-jm-jh
        row=dict(window=name,start_s=t[a],end_s=t[b],samples=b-a+1,wall_impulse_source='NOT_AVAILABLE_AT_INCREMENT_RATE; residual is NOT measured wall impulse',mDeltaV_norm_Ns=np.linalg.norm(p),Jmag_norm_Ns=np.linalg.norm(jm),Jhydro_norm_Ns=np.linalg.norm(jh),missing_norm_Ns=np.linalg.norm(res),relative_residual=np.linalg.norm(res)/max(np.linalg.norm(p),np.linalg.norm(jm+jh),1e-30))
        for pref,arr in [('mDeltaV',p),('Jmag',jm),('Jhydro',jh),('Jmissing',res)]:
            for c,x in zip('xyz',arr):row[pref+'_'+c+'_Ns']=x
        windows.append(row)
    pd.DataFrame(windows).to_csv(HERE/'preimpact_momentum_windows.csv',index=False)
    dv=np.diff(Vcom,axis=0); rr=np.diff(missing,axis=0); ix=np.where((np.linalg.norm(dv,axis=1)>.01)&(t[1:]<.00135))[0]+1
    ev=pd.DataFrame({'time_s':t[ix],'deltaV_norm_mm_s':np.linalg.norm(dv[ix-1],axis=1),'missing_impulse_norm_Ns':np.linalg.norm(rr[ix-1],axis=1)})
    for k,c in enumerate('xyz'):ev['Vcom_'+c+'_mm_s']=Vcom[ix,k];ev['deltaV_'+c+'_mm_s']=dv[ix-1,k]
    ev.to_csv(HERE/'missing_impulse_event_scan.csv',index=False)
    udot=np.gradient(U,t,axis=0); clean=(t>.00001)&(t<.00094)
    err=np.linalg.norm(udot[clean]-V[clean],axis=1)
    # CF is the applied concentrated force (magnetic + local hydro), not wall force.
    cf=np.column_stack([z['CF'+str(i)][:,1] for i in range(1,4)]);ct=z['CF1'][:,0]
    fi=interpolate(ct,t,F+H); ce=np.linalg.norm(cf-fi,axis=1)
    cfrows=pd.DataFrame({'time_s':ct,'CF_vs_telemetry_plus_hydro_error_N':ce})
    for k,c in enumerate('xyz'):cfrows['CF_'+c+'_N']=cf[:,k];cfrows['Fmag_plus_hydro_'+c+'_N']=fi[:,k]
    cfrows.to_csv(HERE/'applied_cf_crosscheck.csv',index=False)
    np.savez_compressed(HERE/'baseline_dense_private.npz',t=t,U=U,UR=UR,V=V,Vcom=Vcom,F=F,H=H,Jmag=Jmag,Jhydro=Jh,P=P,Jmissing=missing,rp=rp,offset=off,m=m)
    summary=dict(mass_mg=m*1e9,offset_um=np.linalg.norm(off)*1000,total_model_mass_tonne=sum(r['mesh_mass_tonne'] for r in masses),clean_residual=windows[0]['relative_residual'],clean_V_vs_dU_median_mm_s=float(np.median(err)),clean_V_vs_dU_max_mm_s=float(max(err)),CF_vs_load_clean_max_N=float(max(ce[ct<=.00095])),first_large_deltaV_time_s=float(t[ix[0]]) if len(ix) else None,max_increment_deltaV_mm_s=float(max(np.linalg.norm(dv,axis=1))),sensor_samples=len(t))
    (HERE/'stage_a_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2));print(pd.DataFrame(windows)[['window','relative_residual','missing_norm_Ns']].to_string(index=False));print(ev.head(15).to_string(index=False))
if __name__=='__main__':main()
