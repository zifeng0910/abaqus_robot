from pathlib import Path
import re, json, math, hashlib, shutil, csv
import numpy as np
from scipy.spatial import ConvexHull
from scipy.ndimage import label
ROOT=Path(__file__).resolve().parent; BASE=ROOT/'case'/'F100_G2P20_ZEROPRESSURE_FACEFIX20'; JOB='F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX20'; DST=ROOT/'case'/JOB
if DST.exists(): raise RuntimeError('refuse overwrite existing child')
deck=(BASE/(BASE.name+'.inp')).read_text(encoding='ascii')
def part(name):
 m=re.search(r'\*Part, name='+re.escape(name)+r'\n(.*?)\*End Part',deck,re.S|re.I)
 if not m: raise RuntimeError(name)
 return m.group(1)
def parse_part(name,etype):
 p=part(name); nm=re.search(r'\*Node\n(.*?)(?=^\*|\Z)',p,re.S|re.M|re.I); nodes={}
 for ln in nm.group(1).splitlines():
  if ln.strip(): a=[x.strip() for x in ln.split(',')]; nodes[int(a[0])]=np.array([float(x) for x in a[1:4]])
 em=re.search(r'\*Element, type='+etype+r'\n(.*?)(?=^\*|\Z)',p,re.S|re.M|re.I); els={}
 for ln in em.group(1).splitlines():
  if ln.strip(): a=[x.strip() for x in ln.split(',')]; els[int(a[0])]=[int(x) for x in a[1:]]
 return nodes,els
fnodes,fels=parse_part('FLUID_EULERIAN','EC3D8R'); rnodes,rels=parse_part('ROBOT_SOLID','C3D4')
# Fluid cells are numbered in the generated EC3D8R order: i fastest, then j, then k.
labels=sorted(fels); centers=np.array([np.mean([fnodes[n] for n in fels[e]],axis=0) for e in labels]);
# exact structured spacing recovered from first cell
v0=np.array([fnodes[n] for n in fels[1]]); ds=np.linalg.norm(v0[1]-v0[0]); dq=np.linalg.norm(v0[3]-v0[0]); cellvol=ds*dq*dq
hull=ConvexHull(np.array(list(rnodes.values()))); eqn=hull.equations[:,:3]; eqd=hull.equations[:,3]
# tube axis and origin from audited closure constants
axis=np.array([.9762799602464296,-.0618320247446977,.2074951564186524]); axis/=np.linalg.norm(axis); n=np.array([.07770950949965555,.9945677258930866,-.06925511349479567]); n/=np.linalg.norm(n); b=np.cross(axis,n); b/=np.linalg.norm(b)
pipe_origin=np.array([-10.591104434430079,-3.464137020895571,-10.215768586540108])+3.2*axis-.015*n; R=.667345
# parent initialized set
m0=re.search(r'\*Elset, elset=FLUID_CEL_INIT_ALL, instance=Fluid_EULERIAN-1\n(.*?)\n\*',deck,re.S|re.I); parent=set()
for x in re.findall(r'\d+',m0.group(1)): parent.add(int(x))
# Deterministic 8^3 subcell quadrature, exact geometry predicates.
q=(np.arange(8)+.5)/8-.5; uu,vv,ww=np.meshgrid(q,q,q,indexing='ij'); local=np.column_stack([uu.ravel()*ds,vv.ravel()*dq,ww.ravel()*dq]);
fractions=[]; rows=[]; grid=np.zeros((44,20,20),float)
for idx,e in enumerate(labels):
 c=centers[idx]; pts=c+local[:,0,None]*axis+local[:,1,None]*n+local[:,2,None]*b; rel=pts-pipe_origin; ss=rel@axis; qn=rel@n; qb=rel@b
 inside_tube=(np.abs(ss)<=6.0)&(qn*qn+qb*qb<=R*R)
 inside_robot=np.all(pts@eqn.T+eqd<=1e-10,axis=1)
 ok=inside_tube & ~inside_robot; frac=float(ok.mean()); fractions.append(frac)
 i=idx%44; j=(idx//44)%20; k=idx//(44*20); grid[i,j,k]=frac
 radial=np.sqrt(qn[len(qn)//2]**2+qb[len(qb)//2]**2); robot_any=bool(np.any(inside_robot)); tube_edge=bool(abs(radial-R)<=math.sqrt(2)*dq/2 or robot_any)
 rows.append({'element_id':e,'i':i,'j':j,'k':k,'cell_volume_mm3':cellvol,'parent_evf':1.0 if e in parent else 0.0,'child_evf':frac,'inside_tube_fraction':float(inside_tube.mean()),'robot_occupied_fraction':float(inside_robot.mean()),'robot_adjacent':robot_any,'tube_wall_adjacent':bool(abs(radial-R)<=math.sqrt(2)*dq/2),'axial_end_cell':bool(i in (0,43)),'physically_accessible':frac>0,'whole_cell_parent_excluded_but_accessible':bool(frac>0 and e not in parent),'legitimate_solid_exclusion':bool(frac==0 and (inside_robot.mean()>0.99 or inside_tube.mean()<0.01))})
# 6-neighbor components among occupied cells
def comps(g):
 lab,n=label(g>1e-12,structure=np.array([[[0,0,0],[0,1,0],[0,0,0]],[[0,1,0],[1,1,1],[0,1,0]],[[0,0,0],[0,1,0],[0,0,0]]]))
 sizes=np.bincount(lab.ravel())[1:]; return int(n),float(sizes.max()/max((g>1e-12).sum(),1)),lab
np0=grid.copy(); np1=grid.copy(); parentgrid=np.zeros_like(grid)
for e in parent: ii=labels.index(e); i=ii%44;j=(ii//44)%20;k=ii//(44*20);parentgrid[i,j,k]=1
npar,fracpar,lpar=comps(parentgrid); nchild,fracchild,lchild=comps(np1)
vol_parent=len(parent)*cellvol; vol_child=float(np.sum(grid)*cellvol); domain=17600*cellvol; tube_nom=math.pi*R*R*12; robot_vol=hull.volume; accessible_nom=tube_nom-robot_vol; deficit_parent=100*(accessible_nom-vol_parent)/accessible_nom; deficit_child=100*(accessible_nom-vol_child)/accessible_nom
# group EVFs rounded to 6 decimals; full cells grouped separately.
groups={}
for e,f in zip(labels,fractions):
 if f>1e-12: groups.setdefault(round(f,6),[]).append(e)
# modify assembly init elset block and initial conditions
oldblock=m0.group(0); elset_lines=[]
for qf,es in sorted(groups.items(),reverse=True):
 name='INITFIX_EVF_'+format(qf,'.6f').replace('.','P')
 elset_lines.append(f'*Elset, elset={name}, instance=Fluid_EULERIAN-1')
 for z in range(0,len(es),16): elset_lines.append(', '.join(str(x) for x in es[z:z+16]))
newblock='\n'.join(elset_lines)
deck2=deck.replace(oldblock,newblock,1)
oldic='*Initial Conditions, type=VOLUME FRACTION\nFLUID_CEL_INIT_ALL, Fluid_EULERIAN-1_WATER, 1.0'
ic='*Initial Conditions, type=VOLUME FRACTION\n'+'\n'.join(f'INITFIX_EVF_{format(qf,".6f").replace(".","P")}, Fluid_EULERIAN-1_WATER, {qf:.6f}' for qf in sorted(groups,reverse=True))
if oldic not in deck2: raise RuntimeError('initial condition block not found')
deck2=deck2.replace(oldic,ic,1)
# preserve all physics and only add a provenance header
header='** INITFIX20: deterministic 8x8x8 subcell EVF initialization from FACEFIX20. Mesh/physics/boundaries unchanged.\n'
deck2=header+deck2
DST.mkdir(parents=True)
(DST/(JOB+'.inp')).write_text(deck2,encoding='ascii'); shutil.copy2(BASE/'vuamp_precomputed_truecel.f90',DST/'vuamp_precomputed_truecel.f90'); shutil.copy2(BASE/'magnetic_field_gradient_table_B0P11_A14P5.dat',DST/'magnetic_field_gradient_table_B0P11_A14P5.dat')
for r in rows: r['parent_component_id']=int(lpar[r['i'],r['j'],r['k']]); r['child_component_id']=int(lchild[r['i'],r['j'],r['k']])
with (ROOT/'initialization_cell_map.csv').open('w',newline='') as f: w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
aud={'parent_job':BASE.name,'child_job':JOB,'mesh_elements':len(labels),'cell_volume_mm3':cellvol,'domain_volume_mm3':domain,'tube_nominal_volume_mm3':tube_nom,'robot_convex_hull_volume_mm3':robot_vol,'nominal_accessible_liquid_volume_mm3':accessible_nom,'parent_initialized_liquid_volume_mm3':vol_parent,'child_initialized_liquid_volume_mm3':vol_child,'parent_void_volume_mm3':domain-vol_parent,'child_void_volume_mm3':domain-vol_child,'parent_deficit_percent':deficit_parent,'child_deficit_percent':deficit_child,'parent_initialized_elements':len(parent),'child_positive_evf_elements':int(sum(f>0 for f in fractions)),'child_partial_evf_elements':int(sum(1e-12<f<1-1e-12 for f in fractions)),'parent_connected_components':npar,'child_connected_components':nchild,'parent_largest_component_fraction':fracpar,'child_largest_component_fraction':fracchild,'parent_artificial_accessible_empty_volume_mm3':float(sum(r['cell_volume_mm3'] for r in rows if r['whole_cell_parent_excluded_but_accessible'])),'child_artificial_accessible_empty_volume_mm3':0.0,'subcell_quadrature':'8x8x8 midpoint samples; fluid inside axial cylinder and outside robot convex hull; conservative no-solid-overlap predicate','limitations':'Convex-hull robot exclusion is conservative for concave details; subcell EVF is an approximation, not exact CAD Boolean.'}
(ROOT/'initialization_volume_audit.json').write_text(json.dumps(aud,indent=2)+'\n'); (DST/'initialization_manifest.json').write_text(json.dumps({'job':JOB,'parent':BASE.name,'status':'PREPARED','physical_change':'Eulerian initial volume fractions only','mesh_unchanged':True,'boundary_sets_unchanged':True,'contact_unchanged':True,'quadrature':'8x8x8','volume_audit':aud},indent=2)+'\n')
print(json.dumps(aud,indent=2))
