"""Verified nearest point to actual wall triangles; no centreline distance proxy.

Positive means lumen side. Centreline is used ONLY to orient triangle normals.
All triangle candidates are bounded by centroid distance minus maximal radius.
Rigid-node reconstruction is checked against independently exported ODB nodes.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from audit_stage_a import HERE,ROOT,JOB,mesh_properties

def closest(P,A,B,C):
    # Broadcasting supports P(N,1,3), triangles(N,K,3).
    e=B-A;f=C-A;n=np.cross(e,f);n2=np.sum(n*n,axis=-1)
    h=np.sum((P-A)*n,axis=-1)
    q=P-n*(h/n2)[...,None]
    ee=np.sum(e*e,axis=-1);ef=np.sum(e*f,axis=-1);ff=np.sum(f*f,axis=-1)
    w=q-A;we=np.sum(w*e,axis=-1);wf=np.sum(w*f,axis=-1);den=ee*ff-ef*ef
    u=(ff*we-ef*wf)/den;v=(ee*wf-ef*we)/den
    inside=(u>=-1e-13)&(v>=-1e-13)&(u+v<=1+1e-13)
    d=np.sum((P-q)**2,axis=-1);d=np.where(inside,d,np.inf)
    for aa,bb in ((A,B),(B,C),(C,A)):
        edge=bb-aa;lam=np.clip(np.sum((P-aa)*edge,axis=-1)/np.sum(edge*edge,axis=-1),0,1)
        qq=aa+lam[...,None]*edge;dd=np.sum((P-qq)**2,axis=-1);use=dd<d
        q=np.where(use[...,None],qq,q);d=np.minimum(dd,d)
    return q,np.sqrt(d)

def tests():
    A=np.array([[0.,0,0]]); B=np.array([[1.,0,0]]);C=np.array([[0.,1,0]])
    for p,want in [([.2,.3,2],[.2,.3,0]),([2,0,0],[1,0,0]),([.7,.7,1],[.5,.5,0]),([-.2,-.3,-1],[0,0,0])]:
        q,d=closest(np.array([p]),A,B,C);assert np.allclose(q,[want]),(p,q,want)
    return 'PASS interior, edges, vertices and both sides'

class Wall:
    def __init__(self):
        df=pd.read_csv(ROOT/(JOB+'_pipe_wall_nodes_exact.csv')).set_index('node');el=pd.read_csv(ROOT/(JOB+'_pipe_wall_triangles_exact.csv'))
        tri=[];labels=[]
        for row in el.itertuples():
            for ids in ((row.n1,row.n2,row.n3),(row.n1,row.n3,row.n4)):
                tri.append(df.loc[list(ids)].to_numpy());labels.append(row.element)
        self.tri=np.array(tri);self.elements=np.array(labels);mid=self.tri.mean(axis=1)
        self.centers=mid;self.tree=cKDTree(mid);self.radius=np.linalg.norm(self.tri-mid[:,None,:],axis=2).max()
        curve=pd.read_csv(ROOT/'CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_EXTSURF_centerline_odb.csv')[['x_mm','y_mm','z_mm']].to_numpy()
        a=curve[:-1];e=np.diff(curve,axis=0);lam=np.clip(np.sum((mid[:,None,:]-a)*e,axis=2)/np.sum(e*e,axis=1),0,1)
        q=a[None,:,:]+lam[:,:,None]*e;nearest=q[np.arange(len(mid)),np.argmin(np.linalg.norm(mid[:,None,:]-q,axis=2),axis=1)]
        normal=np.cross(self.tri[:,1]-self.tri[:,0],self.tri[:,2]-self.tri[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
        normal*=np.where(np.sum(normal*(nearest-mid),axis=1)>0,1,-1)[:,None];self.normals=normal
    def query(self,p,all_faces=False):
        k=len(self.tri) if all_faces else 64
        while True:
            dc,ix=self.tree.query(p,k=k)
            tr=self.tri[ix];q,dist=closest(p[:,None,:],tr[:,:,0,:],tr[:,:,1,:],tr[:,:,2,:]);idx=np.argmin(dist,axis=1)
            dd=dist[np.arange(len(p)),idx]
            if k==len(self.tri) or np.all(dc[:,-1]-self.radius>dd+1e-10):break
            k=min(k*2,len(self.tri))
        j=ix[np.arange(len(p)),idx];qq=q[np.arange(len(p)),idx]
        signed=dd*np.sign(np.sum((p-qq)*self.normals[j],axis=1))
        return signed,j,qq

def main():
    tests();wall=Wall();meshes,rp,_=mesh_properties((ROOT/(JOB+'.inp')).read_text());mesh=meshes['Robot_SOLID']
    faces=pd.read_csv(ROOT/(JOB+'_robot_surface_triangles_exact.csv'));ids=np.unique(faces[['n1','n2','n3']].to_numpy())
    nodes=np.array([mesh['nodes'][i] for i in ids])+mesh['shift'];d=np.load(HERE/'baseline_dense_private.npz'); t=d['t']
    orig=pd.read_csv(ROOT/'output'/ (JOB+'_lub_baseline')/'robot_frames.csv')
    regress=[];sparse=[]
    for ti,df in orig.groupby('time_s'):
        i=np.argmin(abs(t-ti));P=rp+d['U'][i]+Rotation.from_rotvec(d['UR'][i]).apply(nodes-rp)
        actual=df.set_index('node').loc[ids][['x_mm','y_mm','z_mm']].to_numpy()
        err=np.linalg.norm(P-actual,axis=1).max();regress.append(err)
        g,j,q=wall.query(actual); k=np.argmin(g)
        sparse.append(dict(time_s=ti,gap_um=g[k]*1000,node=int(ids[k]),wall_element=int(wall.elements[j[k]]),triangle_index=int(j[k])))
    pd.DataFrame(sparse).to_csv(HERE/'verified_sparse_wall_gap.csv',index=False)
    event_idx=np.flatnonzero((t>=.00095-1e-10)&(t<=.00110+1e-10))[::5]
    # Dense 0.5-us reconstruction is not independent contact force evidence.
    rows=[]
    for i in event_idx:
        P=rp+d['U'][i]+Rotation.from_rotvec(d['UR'][i]).apply(nodes-rp)
        g,j,q=wall.query(P);k=np.argmin(g)
        rows.append(dict(time_s=t[i],gap_um=g[k]*1000,node=int(ids[k]),wall_element=int(wall.elements[j[k]]),triangle_index=int(j[k]),source='existing increment orientation and U; not new contact output'))
    pd.DataFrame(rows).to_csv(HERE/'existing_dense_geometry_event.csv',index=False)
    summary=dict(unit_tests=tests(),surface_nodes=len(ids),surface_faces=len(faces),wall_quads=len(wall.tri)//2,reconstruction_max_error_um=max(regress)*1000,sparse_min_gap_um=min(x['gap_um'] for x in sparse),event_min_gap_um=min(x['gap_um'] for x in rows),projection_fix='q=P-n*h/n2, legacy used A-n*h/n2',sign='positive lumen side; closest actual wall triangle')
    (HERE/'geometry_audit_summary.json').write_text(json.dumps(summary,indent=2));print(summary);print(pd.DataFrame(sparse).iloc[19:23].to_string(index=False))
if __name__=='__main__':main()
