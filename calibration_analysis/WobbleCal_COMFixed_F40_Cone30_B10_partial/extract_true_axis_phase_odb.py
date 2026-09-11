"""Extract fixed geometric HEAD/TAIL vectors from a rigid-body ODB.

Run with Abaqus Python:
  abaqus python extract_true_axis_phase_odb.py JOB OUT.csv

The ODB has no explicit HEAD/TAIL node sets, so fixed sets are defined once
from the 5% extreme projections of the reference robot mesh along its first
principal axis.  They are then carried by the RP finite rotation vector.
"""
from __future__ import print_function
from odbAccess import openOdb
import csv, math, sys

def dot(a,b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def sub(a,b): return (a[0]-b[0],a[1]-b[1],a[2]-b[2])
def add(a,b): return (a[0]+b[0],a[1]+b[1],a[2]+b[2])
def mul(a,q): return (a[0]*q,a[1]*q,a[2]*q)
def norm(a): return math.sqrt(max(0.0,dot(a,a)))
def unit(a):
    n=norm(a); return mul(a,1.0/n) if n>1e-15 else (1.0,0.0,0.0)
def rot(v,rv):
    ang=norm(rv)
    if ang<1e-14: return v
    k=mul(rv,1.0/ang); c=math.cos(ang); s=math.sin(ang)
    return add(add(mul(v,c),mul((k[1]*v[2]-k[2]*v[1],k[2]*v[0]-k[0]*v[2],k[0]*v[1]-k[1]*v[0]),s)),mul(k,dot(k,v)*(1-c)))
def vec(v): return (float(v[0]),float(v[1]),float(v[2]))
def field_vec(frame,key,region):
    if key not in frame.fieldOutputs: return (0.0,0.0,0.0)
    vs=frame.fieldOutputs[key].getSubset(region=region).values
    return vec(vs[0].data) if vs else (0.0,0.0,0.0)
def flatten_nodes(obj):
    out=[]
    try:
        for x in obj:
            if hasattr(x,'coordinates'): out.append(x)
            else: out.extend(flatten_nodes(x))
    except TypeError:
        pass
    return out

job=sys.argv[1]; out=sys.argv[2]
odb=openOdb(job+'.odb',readOnly=True)
root=odb.rootAssembly; step=list(odb.steps.values())[-1]
rp_nodes=flatten_nodes(root.nodeSets['RP_ROBOT'].nodes)
if len(rp_nodes)==0: raise RuntimeError('RP_ROBOT empty')
rp=rp_nodes[0]; rp0=vec(rp.coordinates)
robot=root.instances['ROBOT_SOLID-1']; nodes=list(robot.nodes)
pts=[vec(n.coordinates) for n in nodes]
c=tuple(sum(p[k] for p in pts)/len(pts) for k in range(3))
# Power iteration for first principal axis.
cov=[[0.0]*3 for _ in range(3)]
for p in pts:
 q=sub(p,c)
 for i in range(3):
  for j in range(3): cov[i][j]+=q[i]*q[j]
axis=(1.0,0.0,0.0)
for _ in range(60):
 w=tuple(sum(cov[i][j]*axis[j] for j in range(3)) for i in range(3)); axis=unit(w)
proj=[dot(sub(p,c),axis) for p in pts]
lo=sorted(proj)[max(0,int(0.05*len(proj)))]
hi=sorted(proj)[min(len(proj)-1,int(0.95*len(proj)))]
tail=[(n,p) for n,p,q in zip(nodes,pts,proj) if q<=lo]
head=[(n,p) for n,p,q in zip(nodes,pts,proj) if q>=hi]
if not head or not tail: raise RuntimeError('failed to define fixed head/tail groups')
step_region=root.nodeSets['RP_ROBOT']
with open(out,'w') as f:
 w=csv.writer(f)
 w.writerow(['frame','time_s','ur1','ur2','ur3','rp_x','rp_y','rp_z','head_x','head_y','head_z','tail_x','tail_y','tail_z','axis_x','axis_y','axis_z','recon_node_rms_mm','head_count','tail_count'])
 for fi,fr in enumerate(step.frames):
  ur=field_vec(fr,'UR',step_region); u=field_vec(fr,'U',step_region); rpcur=add(rp0,u)
  h=[add(rpcur,rot(sub(p,rp0),ur)) for n,p in head]
  t=[add(rpcur,rot(sub(p,rp0),ur)) for n,p in tail]
  hc=tuple(sum(p[k] for p in h)/len(h) for k in range(3)); tc=tuple(sum(p[k] for p in t)/len(t) for k in range(3)); ah=unit(sub(hc,tc))
  # Compare reconstructed rigid positions to nodal U where available.
  rms=0.0; ncheck=0
  if 'U' in fr.fieldOutputs:
   uv={v.nodeLabel:vec(v.data) for v in fr.fieldOutputs['U'].getSubset(region=robot).values}
   for n,p in (head[:10]+tail[:10]):
    if n.label in uv:
     pr=add(p,uv[n.label]); pp=add(rpcur,rot(sub(p,rp0),ur)); d=sub(pr,pp); rms+=dot(d,d); ncheck+=1
  rms=math.sqrt(rms/ncheck) if ncheck else float('nan')
  w.writerow([fi,fr.frameValue,ur[0],ur[1],ur[2],rpcur[0],rpcur[1],rpcur[2],hc[0],hc[1],hc[2],tc[0],tc[1],tc[2],ah[0],ah[1],ah[2],rms,len(head),len(tail)])
print('wrote',out,'frames',len(step.frames),'head',len(head),'tail',len(tail),'axis0',axis)
odb.close()
