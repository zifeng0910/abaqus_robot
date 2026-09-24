"""Cheap phase-matched rocking check of the existing 30-ms F100 result."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / 'case' / 'TRUECEL_B0P11_G2P20_A14P5_F100_FAST'
KEYS = ('UR1','UR2','UR3','VR1','VR2','VR3','V1','V2','V3','U1','U2','U3')

def unit(x):
    x=np.asarray(x,float); return x/np.linalg.norm(x)

def main():
    ident=json.loads((CASE/'case_identity.json').read_text())
    a=np.load(CASE/'private/rp_history_private.npz'); b=np.load(CASE/'F100_CYCLE3_RESTART_private/rp_history_private.npz')
    h={}
    for k in KEYS:
        first=a[k].astype(float); second=b[k].astype(float); second[:,0]+=first[-1,0]
        h[k]=np.vstack((first,second[second[:,0]>first[-1,0]+1e-12]))
    a.close();b.close()
    t=h['UR1'][:,0]; ur=np.column_stack([h[f'UR{i}'][:,1] for i in (1,2,3)])
    vr=np.column_stack([h[f'VR{i}'][:,1] for i in (1,2,3)])
    v=np.column_stack([h[f'V{i}'][:,1] for i in (1,2,3)])
    rock=unit(ident['b_routeA_aba']); n=unit(ident['n_routeA_aba']); c=unit(ident['canonical_plus_s_axis_aba'])
    axis0=unit(ident['head_tail_axis_aba'])
    axis=Rotation.from_rotvec(ur).apply(np.broadcast_to(axis0,ur.shape))
    alpha=np.rad2deg(np.arctan2(axis@n,axis@c))
    omega=vr@rock; vs=v@c
    phase=np.linspace(0,1,2001); local=phase*.01
    aa=np.array([np.interp((k-1)*.01+local,t,alpha) for k in (1,2,3)])
    ww=np.array([np.interp((k-1)*.01+local,t,omega) for k in (1,2,3)])
    base=aa[0]-aa[0].mean()
    metrics=[]
    for k in range(3):
        y=aa[k]; yc=y-y.mean()
        coeff=np.sum(yc*np.exp(-2j*np.pi*phase))/len(phase)
        coeff0=np.sum(base*np.exp(-2j*np.pi*phase))/len(phase)
        lag=((np.angle(coeff)-np.angle(coeff0)+np.pi)%(2*np.pi)-np.pi)*180/np.pi
        metrics.append({'cycle':k+1,'rocking_peak_to_peak_deg':float(np.ptp(y)),
                        'fundamental_amplitude_deg':float(2*abs(coeff)),
                        'phase_lag_vs_C1_deg':float(lag),
                        'rms_phase_matched_alpha_error_deg':float(np.sqrt(np.mean((y-aa[0])**2))),
                        'p95_abs_omega_difference_rad_s':float(np.percentile(abs(ww[k]-ww[0]),95))})
    mask=(t>=.02)&(vs<0); onset=float(t[np.flatnonzero(mask)[0]])
    s=np.column_stack([h[f'U{i}'][:,1] for i in (1,2,3)])@c
    back=np.maximum.accumulate(s)-s; imax=int(np.argmax(back));
    out={'metrics':metrics,'cycle3_negative_velocity_onset_ms':onset*1000,
         'max_backtrack_time_ms':float(t[imax]*1000),
         'classification':('ROCKING_REMAINS_PHASE_COHERENT_THROUGH_CYCLE3'
             if metrics[2]['fundamental_amplitude_deg']>=.7*metrics[0]['fundamental_amplitude_deg']
             and abs(metrics[2]['phase_lag_vs_C1_deg'])<15
             and metrics[2]['rms_phase_matched_alpha_error_deg']<5
             else 'ROCKING_CONTROL_DEGRADES_IN_CYCLE3')}
    # The saved contact history contains only whole-surface General Contact
    # resultants.  It cannot be interpreted as a robot-wall pair force.
    contact=np.load(CASE/'F100_CYCLE3_RESTART_private/contact_history_private.npz')
    out['direct_robot_wall_resultant']={
        'available':bool(any('/ASSEMBLY_PIPE_WALL_HELPER' in k and
                             'on surface ASSEMBLY_ROBOT_SOLID' in k for k in contact.files)),
        'interpretation':'Pair-isolated history absent; whole General Contact includes robot-fluid and is not a wall resultant.'}
    contact.close()
    import re
    deck=(CASE/f'{CASE.name}.inp').read_text(encoding='latin1')
    def nodes(name):
        block=re.search(r'(?ms)^\*Part, name='+re.escape(name)+r'\s*$.*?^\*End Part\s*$',deck).group(0)
        body=re.search(r'(?ms)^\*Node\s*$\n(.*?)(?=^\*)',block).group(1)
        return np.array([[float(x) for x in line.split(',')[1:4]] for line in body.splitlines() if line.strip()])
    robot=nodes('Robot_SOLID'); wall=nodes('Pipe_WALL_HELPER')
    rp0=np.asarray(ident['initial_center_aba_mm']); pipe=rp0-float(ident['initial_axial_shift_mm'])*c-float(ident['radial_offset_n_mm'])*n
    ss=(wall-pipe)@c; ring=wall[abs(ss-ss.min())<1e-7]
    xy=np.column_stack(((ring-pipe)@n,(ring-pipe)@rock)); xy=xy[np.argsort(np.arctan2(xy[:,1],xy[:,0]))]
    normals=[]; radii=[]
    for p,q in zip(xy,np.roll(xy,-1,axis=0)):
        edge=q-p; nv=np.array([edge[1],-edge[0]]);nv=unit(nv)
        if nv@(p+q)<0:nv=-nv
        normals.append(nv[0]*n+nv[1]*rock);radii.append(p@nv)
    normals=np.asarray(normals); radius=float(np.mean(radii))
    material=(robot-rp0)@c; regions={'TAIL':robot[material<=material.min()+.25],
                                     'HEAD':robot[material>=material.max()-.25]}
    events=[]
    for label,ti in [('negative_velocity_onset',onset),('maximum_backtrack',float(t[imax]))]:
        u=np.array([np.interp(ti,h[f'U{i}'][:,0],h[f'U{i}'][:,1]) for i in (1,2,3)])
        urt=np.array([np.interp(ti,h[f'UR{i}'][:,0],h[f'UR{i}'][:,1]) for i in (1,2,3)])
        gaps={}
        for name,pts0 in regions.items():
            pts=rp0+u+Rotation.from_rotvec(urt).apply(pts0-rp0)
            gaps[name]=float(radius-np.max((pts-pipe)@normals.T))
        events.append({'event':label,'time_ms':ti*1000,'signed_gap_mm':gaps})
    out['wall_gap_events']=events
    (ROOT/'TRUECEL_F100_G2P20_ROCKING_CHECK.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
