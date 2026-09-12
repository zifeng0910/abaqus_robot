"""Independent surface-resultant contact impulse vs robot linear momentum."""
import json,sys
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
from scipy.integrate import cumulative_trapezoid
from audit_stage_a import HERE,ROOT,JOB,dense_history,vec,interpolate,mesh_properties
from exact_gap_audit import Wall,tests

NEW='Wobble_F30_G6L45_ReducedHydroFixed_ContactAudit_0012'

def main():
    t0,d,z=dense_history(HERE/'diagnostic');t=np.round(t0/1e-7)*1e-7
    assert np.all(np.diff(t)>0)
    dd=pd.DataFrame(d);U=vec(dd,'U');UR=vec(dd,'UR');V=vec(dd,'V');w=vec(dd,'VR')
    meshes,rp,masses=mesh_properties((ROOT/(NEW+'.inp')).read_text());mesh=meshes['Robot_SOLID'];m=mesh['mass'];off=mesh['com']-rp
    Vcom=V+np.cross(w,Rotation.from_rotvec(UR).apply(np.broadcast_to(off,V.shape)))
    base=np.load(HERE/'baseline_dense_private.npz')
    # Same increment index, not interpolation across a discontinuity.
    identity=[]
    for name,val in [('U',U),('UR',UR),('V',V),('Vcom',Vcom)]:
        err=np.linalg.norm(val-base[name][:len(t)],axis=1)
        identity.append(dict(quantity=name,max_error=float(max(err)),RMS_error=float(np.sqrt(np.mean(err**2))),samples=len(t),matching='same explicit increment index'))
    pd.DataFrame(identity).to_csv(HERE/'instrumentation_dynamics_regression.csv',index=False)
    ch=np.load(HERE/'diagnostic'/'contact_history_private.npz')
    def force(prefix,surface):
        a=[]
        for i in (1,2,3):
            key=next(k for k in ch.files if '|'+prefix+str(i)+' on surface ' in k and surface in k)
            ar=ch[key];assert len(ar)==len(t)
            a.append(ar[:,1])
        return np.column_stack(a)
    N=force('CFN','ASSEMBLY_ROBOT');S=force('CFS','ASSEMBLY_ROBOT');T=force('CFT','ASSEMBLY_ROBOT');wallT=force('CFT','ASSEMBLY_PIPE')
    Fwall=N+S
    Ftotalhist=np.column_stack([z['CF'+str(i)][:,1] for i in (1,2,3)])
    tel=pd.read_csv(ROOT/(NEW+'_telemetry.csv'));hyd=pd.read_csv(ROOT/(NEW+'_hydro_increment.csv'),skiprows=2,names=['time_s','v1','v2','v3','w1','w2','w3','Fh1','Fh2','Fh3','Th1','Th2','Th3'])
    F=interpolate(t,np.r_[0,tel.t_s],np.vstack((np.zeros(3),tel[['fx_aba_N','fy_aba_N','fz_aba_N']])))
    H=interpolate(t,np.r_[0,hyd.time_s],np.vstack((np.zeros(3),vec(hyd,'Fh'))))
    jf=cumulative_trapezoid(F,t,axis=0,initial=0);jh=cumulative_trapezoid(H,t,axis=0,initial=0);jw=cumulative_trapezoid(Fwall,t,axis=0,initial=0)
    jn=cumulative_trapezoid(N,t,axis=0,initial=0);js=cumulative_trapezoid(S,t,axis=0,initial=0)
    P=m*(Vcom-Vcom[0]);res=P-jf-jh-jw
    active=np.linalg.norm(Fwall,axis=1)>1e-8
    edges=np.diff(np.r_[False,active,False].astype(int));starts=np.where(edges==1)[0];ends=np.where(edges==-1)[0]-1
    events=[]
    for a,b in zip(starts,ends):
        events.append(dict(start_s=t[a],end_s=t[b],active_duration_us=(b-a+1)*.1,samples=b-a+1,peak_force_N=float(max(np.linalg.norm(Fwall[a:b+1],axis=1)))))
    pd.DataFrame(events).to_csv(HERE/'contact_event_catalog.csv',index=False)
    windows=[]
    for name,lo,hi in [('A_clean',0,.00095),('B_first_impact',.00095,.00110),('C_after_first',.00110,.00125),('D_second_impact_clipped',.0012,.0013),('all_1p3ms',0,.0013),('first_impact_detail',.00098,.00108)]:
        a=int(round(lo/1e-7));b=int(round(hi/1e-7));p=P[b]-P[a];fm=jf[b]-jf[a];fh=jh[b]-jh[a];fw=jw[b]-jw[a];rr=p-fm-fh-fw
        row=dict(window=name,start_s=t[a],end_s=t[b],mDeltaV_norm_Ns=np.linalg.norm(p),Jmag_norm_Ns=np.linalg.norm(fm),Jhydro_norm_Ns=np.linalg.norm(fh),Jwall_norm_Ns=np.linalg.norm(fw),residual_norm_Ns=np.linalg.norm(rr),relative_residual=np.linalg.norm(rr)/max(np.linalg.norm(p),np.linalg.norm(fm+fh+fw),1e-30),relative_without_wall=np.linalg.norm(p-fm-fh)/max(np.linalg.norm(p),np.linalg.norm(fm+fh),1e-30))
        for pref,vv in [('mDeltaV',p),('Jmag',fm),('Jhydro',fh),('Jwall',fw),('residual',rr)]:
            for c,x in zip('xyz',vv):row[pref+'_'+c+'_Ns']=x
        windows.append(row)
    pd.DataFrame(windows).to_csv(HERE/'first_impact_momentum_closure.csv',index=False)
    # Field resultant sum comparison is a separate summation/sign check.
    fields=pd.read_csv(HERE/'diagnostic'/'contact_fields.csv');ft=fields.time_s.to_numpy();fi=np.rint(ft/1e-7).astype(int)
    fsum=fields[['Fn_x_N','Fn_y_N','Fn_z_N']].to_numpy()+fields[['Fs_x_N','Fs_y_N','Fs_z_N']].to_numpy()
    ferr=np.linalg.norm(fsum-Fwall[fi],axis=1)
    wall=Wall();tests();faces=pd.read_csv(ROOT/(JOB+'_robot_surface_triangles_exact.csv'));ids=np.unique(faces[['n1','n2','n3']].to_numpy());nodes=np.array([mesh['nodes'][i] for i in ids])+mesh['shift']
    # Exact distance from every exterior node on every increment in the event.
    ind=np.arange(9800,10801);gaps=[]
    if '--reuse-gap' in sys.argv:
        gaps=pd.read_csv(HERE/'first_impact_gap_history.csv')
        assert np.allclose(gaps.time_s,t[ind],rtol=0,atol=1e-12)
    else:
        for i in ind:
            pos=rp+U[i]+Rotation.from_rotvec(UR[i]).apply(nodes-rp)
            gap,j,q=wall.query(pos);k=np.argmin(gap)
            gaps.append(dict(time_s=t[i],gap_um=gap[k]*1000,contact_node=int(ids[k]),wall_element=int(wall.elements[j[k]]),triangle_index=int(j[k])))
        gaps=pd.DataFrame(gaps);gaps.to_csv(HERE/'first_impact_gap_history.csv',index=False)
    # CPRESS is only at 0.5 us field samples; do NOT invent an increment CPRESS.
    history=gaps.copy();history['speed_mm_s']=np.linalg.norm(Vcom[ind],axis=1)
    fieldmap={int(round(row.time_s/1e-7)):row for row in fields.itertuples()}
    history['CPRESS_MPa']=[fieldmap[i].CPRESS_MPa if i in fieldmap else np.nan for i in ind]
    history['CPRESS_node']=[fieldmap[i].CPRESS_node if i in fieldmap else np.nan for i in ind]
    for prefix,val in [('Vcom',Vcom),('Fmag',F),('Fhydro',H),('Fwall',Fwall),('Jwall_cumulative',jw),('mDeltaV',P),('residual',res)]:
        for k,c in enumerate('xyz'):history[prefix+'_'+c]=val[ind,k]
    history['force_unit']='N';history['velocity_unit']='mm/s';history['impulse_unit']='N s';history['gap_source']='all 217 exterior nodes vs true wall triangles, increment orientation'
    history.to_csv(HERE/'first_impact_increment_history.csv',index=False)
    forcehist=pd.DataFrame({'time_s':t,'contact_active':active.astype(int)})
    for prefix,val in [('Fn_robot',N),('Fs_robot',S),('Ftotal_robot',T),('Ftotal_wall',wallT)]:
        for k,c in enumerate('xyz'):forcehist[prefix+'_'+c+'_N']=val[:,k]
    forcehist.to_csv(HERE/'first_impact_contact_force_history.csv',index=False)
    a,b=int(round(.001/1e-7)),int(round(.0010501/1e-7));dv=Vcom[b]-Vcom[a]
    im=max(range(len(t)),key=lambda i:np.linalg.norm(Fwall[i]));ig=gaps.gap_um.idxmin()
    sm=dict(events=events,mass_mg=m*1e9,RP_COM_offset_um=np.linalg.norm(off)*1000,dt_quantization_max_s=float(max(abs(t-t0))),physics_regression_max_V_mm_s=identity[2]['max_error'],action_reaction_max_N=float(max(np.linalg.norm(T+wallT,axis=1))),CFT_vs_CFN_plus_CFS_max_N=float(max(np.linalg.norm(T-Fwall,axis=1))),field_sum_vs_history_max_N=float(max(ferr)),applied_CF_vs_logged_force_max_N=float(max(np.linalg.norm(Ftotalhist-F-H,axis=1))),clean_contact_max_N=float(max(np.linalg.norm(Fwall[t<=.00095],axis=1))),clean_residual=windows[0]['relative_residual'],first_impact_residual=windows[1]['relative_residual'],whole_residual=windows[4]['relative_residual'],first_impact_Jwall_Ns=[windows[1]['Jwall_'+c+'_Ns'] for c in 'xyz'],first_impact_Jwall_norm_Ns=windows[1]['Jwall_norm_Ns'],first_impact_Jnormal_Ns=list(jn[11000]-jn[9500]),first_impact_Jshear_Ns=list(js[11000]-js[9500]),deltaV_1_to_1p0501_xyz_mm_s=dv.tolist(),deltaV_1_to_1p0501_norm_mm_s=float(np.linalg.norm(dv)),peak_force_N=float(np.linalg.norm(Fwall[im])),peak_force_time_ms=t[im]*1000,min_event_gap_um=float(gaps.gap_um.min()),min_event_gap_time_ms=float(gaps.iloc[ig].time_s*1000),min_event_gap_node=int(gaps.iloc[ig].contact_node),min_event_gap_wall_element=int(gaps.iloc[ig].wall_element),peak_CPRESS_MPa=float(fields.CPRESS_MPa.max()),peak_CPRESS_time_ms=float(fields.loc[fields.CPRESS_MPa.idxmax()].time_s*1000))
    sm['decision']='FORCE_TRANSFER_VALID_HIDDEN_WALL_IMPACT_CONFIRMED' if sm['first_impact_residual']<.1 and sm['clean_residual']<.1 and sm['physics_regression_max_V_mm_s']<1e-3 else 'MOMENTUM_AUDIT_STILL_AMBIGUOUS'
    (HERE/'diagnostic_summary.json').write_text(json.dumps(sm,indent=2,default=lambda x:x.item()))
    np.savez_compressed(HERE/'diagnostic_dense_private.npz',t=t,U=U,UR=UR,V=V,Vcom=Vcom,F=F,H=H,Fwall=Fwall,Jmag=jf,Jhydro=jh,Jwall=jw,P=P,residual=res)
    # Public plotting source, not huge raw telemetry.
    public=pd.DataFrame({'time_s':t,'speed_mm_s':np.linalg.norm(Vcom,axis=1),'Fwall_norm_N':np.linalg.norm(Fwall,axis=1),'residual_norm_Ns':np.linalg.norm(res,axis=1)})
    for prefix,val in [('Vcom',Vcom),('P',P),('Jmag',jf),('Jhydro',jh),('Jwall',jw)]:
        for k,c in enumerate('xyz'):public[prefix+'_'+c]=val[:,k]
    public.to_csv(HERE/'momentum_plot_source.csv',index=False)
    print(json.dumps(sm,indent=2,default=lambda x:x.item()));print(pd.DataFrame(windows)[['window','relative_residual','relative_without_wall']].to_string(index=False))
if __name__=='__main__':main()
