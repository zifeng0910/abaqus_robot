"""Read-only ODB history/field inventory and robot-side contact export."""
from odbAccess import openOdb
import csv, json, sys
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation

job, dest = sys.argv[1:3]
out = Path(dest); out.mkdir(parents=True, exist_ok=True)
o = openOdb(job + '.odb', readOnly=True)
st = list(o.steps.values())[-1]
inventory=[]; arrays={}; contact_history={}; assembly_history={}
for rn, region in st.historyRegions.items():
    for vn, ho in region.historyOutputs.items():
        a=np.asarray(ho.data, dtype=float)
        dt=np.diff(a[:,0]) if len(a)>1 else np.array([float('nan')])
        inventory.append(dict(region=rn,variable=vn,samples=len(a),dt_min_s=float(np.nanmin(dt)),dt_median_s=float(np.nanmedian(dt)),dt_max_s=float(np.nanmax(dt))))
        if 'Node ASSEMBLY.2' == rn:
            arrays[vn]=a
        if rn == 'Assembly ASSEMBLY':
            assembly_history[vn]=a
        if vn.startswith(('CFN','CFS','CFT')):
            contact_history[rn+'|'+vn]=a
with (out/'existing_odb_history_inventory.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(inventory[0]));w.writeheader();w.writerows(inventory)
np.savez_compressed(out/'rp_history_private.npz',**arrays)
np.savez_compressed(out/'contact_history_private.npz',**contact_history)
np.savez_compressed(out/'assembly_history_private.npz',**assembly_history)
freq={}; contacts=[]; poses=[]; node_checks=[]
rp=o.rootAssembly.nodeSets['RP_ROBOT']
robot=o.rootAssembly.instances['ROBOT_SOLID-1']
initial={n.label:np.array(n.coordinates,dtype=float) for n in robot.nodes}
rp0=np.array(list(rp.nodes[0])[0].coordinates,dtype=float)
for frame in st.frames:
    r=dict(time_s=float(frame.frameValue))
    for vn in ('U','UR','V','VR','A','AR','CF','RF','RM'):
        if vn in frame.fieldOutputs:
            v=frame.fieldOutputs[vn].getSubset(region=rp).values
            if v:
                try: a=v[0].dataDouble
                except: a=v[0].data
                for j,x in enumerate(a):r[vn+str(j+1)]=float(x)
    poses.append(r)
    if len(poses)%20==1 and 'U' in frame.fieldOutputs and 'UR1' in r:
        rot=Rotation.from_rotvec([r['UR'+str(i)] for i in (1,2,3)])
        trans=np.array([r['U'+str(i)] for i in (1,2,3)])
        err=[]
        for v in frame.fieldOutputs['U'].getSubset(region=robot).values:
            try: a=np.array(v.dataDouble)
            except: a=np.array(v.data)
            pred=rp0+trans+rot.apply(initial[v.nodeLabel]-rp0)
            err.append(np.linalg.norm(pred-initial[v.nodeLabel]-a))
        node_checks.append(dict(time_s=r['time_s'],max_reconstruction_error_mm=max(err),nodes=len(err)))
    c=dict(time_s=float(frame.frameValue),CPRESS_MPa=0.,CPRESS_node=-1,CPRESS_instance='',COPEN_min_mm=float('nan'),Fn_x_N=0.,Fn_y_N=0.,Fn_z_N=0.,Fs_x_N=0.,Fs_y_N=0.,Fs_z_N=0.,wall_Fn_x_N=0.,wall_Fn_y_N=0.,wall_Fn_z_N=0.,wall_Fs_x_N=0.,wall_Fs_y_N=0.,wall_Fs_z_N=0.)
    for key,fo in frame.fieldOutputs.items():
        if key.startswith(('CPRESS','CNORMF','CSHEARF','COPEN')):
            freq.setdefault(key,[]).append(float(frame.frameValue))
            for v in fo.values:
                inst=v.instance.name if v.instance else 'ASSEMBLY'
                try: a=v.dataDouble
                except: a=v.data
                if key.startswith('COPEN'):
                    val=float(a if np.isscalar(a) else a[0])
                    if abs(val)<1 and (np.isnan(c['COPEN_min_mm']) or val<c['COPEN_min_mm']):c['COPEN_min_mm']=val
                elif key.startswith('CPRESS'):
                    val=float(a if np.isscalar(a) else a[0])
                    if val>c['CPRESS_MPa']:c.update(CPRESS_MPa=val,CPRESS_node=v.nodeLabel,CPRESS_instance=inst)
                elif inst in ('ROBOT_SOLID-1','PIPE_WALL_HELPER-1'):
                    prefix=('wall_' if inst=='PIPE_WALL_HELPER-1' else '')+('Fn_' if key.startswith('CNORMF') else 'Fs_')
                    for k,x in zip('xyz',a):c[prefix+k+'_N']+=float(x)
    contacts.append(c)
for filename,rows in [('rp_fields.csv',poses),('contact_fields.csv',contacts),('orientation_node_crosscheck.csv',node_checks)]:
    keys=sorted(set(k for row in rows for k in row));keys=['time_s']+[k for k in keys if k!='time_s']
    with (out/filename).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
fr=[]
for k,times in freq.items():
    dt=np.diff(times)
    fr.append(dict(variable=k,source='field',samples=len(times),dt_median_s=float(np.median(dt)),dt_min_s=float(min(dt)),dt_max_s=float(max(dt))))
with (out/'contact_output_frequency_audit.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(fr[0]));w.writeheader();w.writerows(fr)
print('HISTORY REGIONS',list(st.historyRegions.keys()));print('RP variables',[(k,len(v)) for k,v in arrays.items()]);print('CONTACT',fr)
o.close()
