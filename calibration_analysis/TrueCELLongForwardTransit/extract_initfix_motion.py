from __future__ import print_function
import sys,os
import numpy as np
from odbAccess import openOdb
odb_path,out=sys.argv[1:3]; os.makedirs(out,exist_ok=True)
axis=np.array([.9762799602464296,-.0618320247446977,.2074951564186524]); axis/=np.linalg.norm(axis)
n=np.array([.07770950949965555,.9945677258930866,-.06925511349479567]); n/=np.linalg.norm(n)
origin=np.array([-10.591104434430079,-3.464137020895571,-10.215768586540108])+3.2*axis-.015*n
odb=openOdb(odb_path,readOnly=True)
try:
 step=list(odb.steps.values())[-1]; insts=odb.rootAssembly.instances
 rname=next(k for k in insts.keys() if 'ROBOT' in k.upper()); inst=insts[rname]
 labels=np.array([x.label for x in inst.nodes],int); coords=np.array([x.coordinates for x in inst.nodes],float)
 label_to_index={int(v):i for i,v in enumerate(labels)}
 times=[]; frames=[]; stride=max(1,len(step.frames)//81)
 for idx,fr in enumerate(step.frames):
  t=float(fr.frameValue)
  if t>.0200001 or idx%stride: continue
  if 'U' not in fr.fieldOutputs: continue
  u=np.zeros_like(coords)
  for v in fr.fieldOutputs['U'].getSubset(region=inst).values:
   j=label_to_index.get(int(v.nodeLabel))
   if j is not None: u[j]=v.data[:3]
  frames.append(coords+u); times.append(t)
 np.savez_compressed(os.path.join(out,'robot_motion_frames.npz'),time=np.array(times),coords=np.array(frames),initial=coords,axis=axis,n=n,origin=origin)
 print('EXTRACTED',rname,len(times),coords.shape)
finally: odb.close()
