import csv, math, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
telemetry = sys.argv[1] if len(sys.argv) > 1 else 'WobbleCal_F30_Cone30_B10_WallOn_Free_003_telemetry.csv'
centerline = sys.argv[2] if len(sys.argv) > 2 else 'curvenew_CEL_xyrot56_exact_abaqus.csv'
out_csv = sys.argv[3] if len(sys.argv) > 3 else 'forward_gradient_force_replay.csv'

def tangent_table(path):
    rows = list(csv.DictReader(open(path, newline='')))
    s = np.array([float(r['arclength_mm']) for r in rows])
    p = np.array([[float(r['x_mm']), float(r['y_mm']), float(r['z_mm'])] for r in rows])
    # Increasing canonical arclength is the sole forward convention.
    d = np.gradient(p, s, axis=0, edge_order=1)
    d /= np.maximum(np.linalg.norm(d, axis=1)[:, None], 1e-15)
    return s, p, d

def tangent_at(sq, s, tang):
    q = min(max(float(sq), float(s[0])), float(s[-1]))
    out = np.array([np.interp(q, s, tang[:, k]) for k in range(3)])
    return out / max(np.linalg.norm(out), 1e-15)

def rotvec(ur):
    a = float(np.linalg.norm(ur))
    if a < 1e-14: return np.eye(3)
    u = ur/a
    K = np.array([[0.,-u[2],u[1]],[u[2],0.,-u[0]],[-u[1],u[0],0.]])
    return np.eye(3) + math.sin(a)*K + (1-math.cos(a))*(K@K)

s, p, tang = tangent_table(centerline)
rows = list(csv.DictReader(open(telemetry, newline='')))
axis0 = tang[0]              # --robot-axis-tangent in production server
polarity = -1.0              # baseline WallOn_Free_003
m = 0.001168                # A m^2, fixed production value
L = 45.0
Gvals = [(-0.006, 'G6_minus'), (0.006, 'G6_plus'), (-0.003, 'G3_minus'), (0.003, 'G3_plus')]
out=[]
for r in rows:
    t=float(r['t_s']); ur=np.array([float(r['ur1']),float(r['ur2']),float(r['ur3'])])
    dr=float(r['driver_arc_mm']); rr=float(r['robot_arc_mm'])
    td=tangent_at(dr,s,tang); tr=tangent_at(rr,s,tang)
    delta=rr-dr; env=math.exp(-0.5*(delta/L)**2)
    R=rotvec(ur); moment=R@(polarity*m*axis0)
    # Production analytic-gradient branch: grad_axis=polarity*t_robot,
    # gdir=t_robot; force=(m dot grad_axis)*dB/ds*gdir.  A signed G is
    # used here only for mathematical replay; production CLI accepts G>=0.
    vals={}
    for G,label in Gvals:
        dbds=-(delta/(L*L))*env*G*1000.0
        grad_axis=polarity*tr
        f=(moment@grad_axis)*dbds*tr*float(r['drive_scale'])
        ft=float(f@td)
        vals[label]=(f,ft)
    rec={'t_s':t,'driver_arc_mm':dr,'robot_arc_mm':rr,'delta_s_mm':delta,
         'tangent_driver_x':td[0],'tangent_driver_y':td[1],'tangent_driver_z':td[2],
         'tangent_robot_x':tr[0],'tangent_robot_y':tr[1],'tangent_robot_z':tr[2],
         'moment_x_Am2':moment[0],'moment_y_Am2':moment[1],'moment_z_Am2':moment[2],
         'telemetry_ft_N':float(r.get('force_tangent_N','nan'))}
    for label,(f,ft) in vals.items():
        rec[label+'_ft_N']=ft; rec[label+'_fnorm_N']=float(np.linalg.norm(f)); rec[label+'_fx_N']=f[0]; rec[label+'_fy_N']=f[1]; rec[label+'_fz_N']=f[2]
    out.append(rec)
fields=list(out[0])
with open(out_csv,'w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=fields); w.writeheader(); w.writerows(out)

def stats(label):
    a=np.array([x[label+'_ft_N'] for x in out]); t=np.array([x['t_s'] for x in out])
    imp=float(np.trapz(a,t)); return dict(mean=float(a.mean()),median=float(np.median(a)),min=float(a.min()),max=float(a.max()),positive_fraction=float(np.mean(a>0)),impulse_Ns=imp)
print('Wrote', out_csv)
for label in ['G3_plus','G3_minus','G6_plus','G6_minus']:
    print(label, stats(label))
print('Baseline telemetry Ft stats:', stats('telemetry') if 'telemetry_ft_N' in out[0] else '')
