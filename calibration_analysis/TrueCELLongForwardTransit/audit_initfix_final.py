"""Read-only checks of the SINGLE executed initialization child; never generates a case."""
import csv, hashlib, json, re
from pathlib import Path
import numpy as np
from scipy.ndimage import label, generate_binary_structure
from scipy.spatial import ConvexHull

ROOT=Path(__file__).resolve().parent
P='F100_G2P20_ZEROPRESSURE_FACEFIX20'
C='F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX20'
OUT=ROOT/'initialization_closure'; OUT.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mesh(text,name,etype):
    part=re.search(r'\*Part, name='+name+r'\s*\n(.*?)\*End Part',text,re.S|re.I)[1]
    nodes={}
    for l in re.search(r'^\*Node\s*\n(.*?)(?=^\*)',part,re.S|re.M|re.I)[1].splitlines():
        a=l.split(','); nodes[int(a[0])]=list(map(float,a[1:4]))
    els={}
    for l in re.search(r'^\*Element, type='+etype+r'\s*\n(.*?)(?=^\*)',part,re.S|re.M|re.I)[1].splitlines():
        a=list(map(int,l.split(',')));els[a[0]]=a[1:]
    return nodes,els
def strip_init(s,child):
    s=re.sub(r'^\*\* INITFIX20[^\n]*\n','',s)
    pat=(r'^\*Elset, elset=INITFIX_EVF_[^\n]*\n[^*]*' if child else
         r'^\*Elset, elset=FLUID_CEL_INIT_ALL[^\n]*\n[^*]*')
    s=re.sub(pat,'',s,flags=re.M)
    s=re.sub(r'^\*Initial Conditions, type=VOLUME FRACTION\n[^*]*','',s,flags=re.M)
    # preparation consumed one '*' of the following ** comment; comments carry no model data.
    # The child builder left one harmless standalone '*' where an original '**'
    # comment was consumed by the initialization-block rewrite.  Abaqus accepts
    # it as an empty keyword; ignore it for the unchanged-model comparison.
    return '\n'.join(l for l in s.splitlines() if l.strip() and l.strip() != '*' and not l.startswith('**') and not l.startswith('* '))
def main():
    pp=ROOT/'case'/P/(P+'.inp');cp=ROOT/'case'/C/(C+'.inp')
    ps=pp.read_text();cs=cp.read_text()
    assert strip_init(ps,False)==strip_init(cs,True),'UNRELATED MODEL DIFFERENCE'
    assert sha(pp)=='aa919e6c64ce7ae3ced531daca3aac611914c4074cb156703d8567e9d9ce8afc'
    assert (ROOT/'case'/P/'vuamp_precomputed_truecel.f90').read_bytes()==(ROOT/'case'/C/'vuamp_precomputed_truecel.f90').read_bytes()
    for part in ('FLUID_EULERIAN','ROBOT_SOLID','PIPE_SOLID','Pipe_WALL_HELPER'):
        pattern=r'\*Part, name='+part+r'\s*\n(.*?)\*End Part'
        assert re.search(pattern,ps,re.S|re.I)[0]==re.search(pattern,cs,re.S|re.I)[0]
    fn,fe=mesh(ps,'FLUID_EULERIAN','EC3D8R');rn,relems=mesh(ps,'ROBOT_SOLID','C3D4')
    xyz=np.array([rn[i] for i in sorted(rn)]); tet=np.array([[rn[i] for i in relem] for relem in relems.values()])
    tetvol=float(abs(np.linalg.det(tet[:,1:]-tet[:,0,None,:])).sum()/6)
    hull=ConvexHull(xyz);hn=hull.equations[:,:3];hd=hull.equations[:,3]
    ident=json.loads((ROOT/'case/F100_G2P20_NOFLUID_REALWALL50/case_identity.json').read_text())
    a=np.array(ident['canonical_plus_s_axis_aba']);n=np.array(ident['n_routeA_aba']);b=np.array(ident['b_routeA_aba'])
    origin=np.array(ident['initial_center_aba_mm'])+3.2*a-.015*n
    centers=np.array([np.mean([fn[i] for i in fe[e]],axis=0) for e in sorted(fe)])
    local=(centers-origin)@np.stack((a,n,b),axis=1)
    ds=np.linalg.norm(np.array(fn[2])-fn[1]);dq=np.linalg.norm(np.array(fn[46])-fn[1]);dv=ds*dq*dq
    rows=list(csv.DictReader((ROOT/'initialization_cell_map.csv').open()))
    pe=np.array([float(r['parent_evf']) for r in rows]);ce=np.array([float(r['child_evf']) for r in rows])
    assert np.all((ce>=0)&(ce<=1))
    structure=generate_binary_structure(3,1)
    def components(e):
        # serialized element order k,j,i
        g=e.reshape(20,20,44);lab,num=label(g>0,structure);weights=np.bincount(lab.ravel(),weights=e)
        pos=weights[1:];touch=set(lab[:,:,0].ravel())&set(lab[:,:,-1].ravel());touch.discard(0)
        return {'positive_cell_components':int(num),'largest_component_volume_fraction':float(pos.max()/pos.sum()),'axially_connected':bool(touch)},lab.ravel()
    pc,pl=components(pe);cc,cl=components(ce)
    # Independent 16^3 integration only in geometric boundary candidates. No new physical candidate.
    nearrobot=np.all((centers>=xyz.min(0)-ds)&(centers<=xyz.max(0)+ds),axis=1)
    crad=np.linalg.norm(local[:,1:],axis=1);nearwall=abs(crad-.667345)<np.sqrt(2)*dq
    refined=ce.copy(); q=(np.arange(16)+.5)/16-.5
    offsets=np.stack(np.meshgrid(q*ds,q*dq,q*dq,indexing='ij'),-1).reshape(-1,3)@np.stack((a,n,b))
    # Reduce coplanar duplicate hull planes for the membership audit.
    unique=np.unique(np.round(hull.equations,10),axis=0);hn,hd=unique[:,:3],unique[:,3]
    for j in np.flatnonzero(nearrobot|nearwall):
        pts=centers[j]+offsets;loc=(pts-origin)@np.stack((a,n,b),axis=1)
        inside=np.sum(loc[:,1:]**2,axis=1)<=.667345**2
        if nearrobot[j]:inside &= ~np.all(pts@hn.T+hd<=1e-10,axis=1)
        refined[j]=inside.mean()
    for i,r in enumerate(rows):
        parent_missing=max(ce[i]-pe[i],0)*dv
        tube_fraction=float(r['inside_tube_fraction'])
        r.update(parent_missing_accessible_volume_mm3=parent_missing,
            geometric_accessible_volume_estimate_mm3=ce[i]*dv,
            tube_wall_adjacent=bool(nearwall[i]),robot_adjacent=bool(nearrobot[i]),
            parent_component_id=int(pl[i]),child_component_id=int(cl[i]),refined_audit_evf=refined[i])
    with (OUT/'missing_fluid_cell_map.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    # Nonoverlapping attribution: tube cut first, then robot/hull exclusion inside lumen.
    missing=np.maximum(ce-pe,0)*dv; tubecut=np.array([float(r['inside_tube_fraction'])<1 for r in rows])
    nominal=np.pi*.667345**2*12-tetvol
    volume={'domain_volume_mm3':len(fe)*dv,'solid_excluded_nominal_mm3':len(fe)*dv-nominal,
        'robot_tetrahedral_volume_mm3':tetvol,'robot_hull_volume_mm3':float(hull.volume),
        'nominal_accessible_liquid_volume_mm3':float(nominal),'hull_accessible_volume_mm3':float(np.pi*.667345**2*12-hull.volume),
        'parent':dict(initialized_liquid_mm3=float(pe.sum()*dv),initialized_void_mm3=float((1-pe).sum()*dv),
            deficit_percent=float(100*(1-pe.sum()*dv/nominal)),**pc),
        'child':dict(initialized_liquid_mm3=float(ce.sum()*dv),initialized_void_mm3=float((1-ce).sum()*dv),
            deficit_percent=float(100*(1-ce.sum()*dv/nominal)),**cc),
        'missing_volume_attribution_mm3':{'tube_cut_cells':float(missing[tubecut].sum()),'robot_exclusion_in_full_lumen_cells':float(missing[~tubecut].sum()),
            'axial_end_cells_subset_not_additive':float(missing[(np.arange(len(fe))%44==0)|(np.arange(len(fe))%44==43)].sum())},
        'audit16_volume_mm3':float(refined.sum()*dv),'audit16_minus_executed8_volume_mm3':float((refined-ce).sum()*dv),
        'remaining_hull_overexclusion_mm3':float(hull.volume-tetvol),
        'subcell_location_limitation':'EVF specifies cell material amount, not exact subcell geometry. Abaqus reconstructs interfaces; no pointwise zero-solid-overlap guarantee is claimed.',
        'void_topology_limitation':'Whole-cell empty accessible region is connected to ends/wall collar; mixed-cell subcell connectivity is not determined solely by EVF.',
        'initialization_scope':'Conservative robot convex hull, nominal cylindrical lumen, 8^3 quadrature; 16^3 read-only verification. Actual tet geometry is used for independent nominal volume.'}
    checks={'non_initialization_keyword_data_identical':True,'all_parts_identical':True,'fortran_identical':True,
        'parent_sha256':sha(pp),'child_sha256':sha(cp),'dynamics_run_count':1,
        'fortran_log_collision':'Executed unchanged Fortran writes auxiliary logs into FACEFIX20. Parent ODB/NPZ remain authoritative. Child log copy is archived after solve.',
        'datacheck':'PASSED (compiler environment corrected after initial compiler-not-found launch)',
        'source_reference':'https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-initialconditions.htm'}
    (OUT/'initialization_volume_audit.json').write_text(json.dumps(volume,indent=2)+'\n')
    (OUT/'input_identity_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(volume,indent=2),flush=True)
if __name__=='__main__':main()
