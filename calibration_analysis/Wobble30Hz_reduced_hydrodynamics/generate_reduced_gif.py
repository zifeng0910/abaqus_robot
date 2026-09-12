from pathlib import Path
import json, numpy as np, pandas as pd, matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]; JOB='Wobble_F30_G6L45_ReducedHydro_WallOn_Free_0083'
rf=pd.read_csv(ROOT/'output'/(JOB+'_bend_validation')/'robot_frames.csv')
dyn=pd.read_csv(HERE/'reduced_hydro_angular_dynamics.csv')
wn=pd.read_csv(ROOT/(JOB+'_pipe_wall_nodes_exact.csv')); tri=pd.read_csv(ROOT/(JOB+'_pipe_wall_triangles_exact.csv'))
nd={int(r.node):np.array([float(r.x_mm),float(r.y_mm),float(r.z_mm)]) for _,r in wn.iterrows()}
verts=[[nd[int(r[n])] for n in ('n1','n2','n3','n4')] for _,r in tri.iterrows()]
curve=pd.read_csv(ROOT/'curvenew_CEL_xyrot56_exact.csv'); tf=json.loads((ROOT/'abaqus_magpylib_frame_transform.json').read_text())
cv=curve[['x_mm','y_mm','z_mm']].to_numpy(float)@np.asarray(tf['R_aba_to_mag'])+np.asarray(tf['origin_aba_mm'])
axis=np.array([.9647382600216,-.1188742372140,.2348382536499]); axis/=np.linalg.norm(axis)
allpts=np.vstack([wn[['x_mm','y_mm','z_mm']].to_numpy(float),cv]); lo=allpts.min(0)-.4; hi=allpts.max(0)+.4; ims=[]
for k in np.linspace(0,int(rf.frame.max()),84).astype(int):
    g=rf[rf.frame==k]; P=g[['x_mm','y_mm','z_mm']].to_numpy(float); c=P.mean(0); q=P-c; hp=c+q[np.argmax(q@axis)]; tp=c+q[np.argmin(q@axis)]; row=dyn.iloc[np.argmin(abs(dyn.time_s-g.time_s.iloc[0]))]
    fig=plt.figure(figsize=(8,6)); ax=fig.add_subplot(111,projection='3d'); ax.add_collection3d(Poly3DCollection(verts,alpha=.10,facecolor='gray',edgecolor='none')); ax.plot(cv[:,0],cv[:,1],cv[:,2],'k-',lw=1); ax.scatter([hp[0]],[hp[1]],[hp[2]],c='r',s=25); ax.scatter([tp[0]],[tp[1]],[tp[2]],c='b',s=25); ax.plot([tp[0],hp[0]],[tp[1],hp[1]],[tp[2],hp[2]],'r-',lw=3)
    ax.set_xlim(lo[0],hi[0]); ax.set_ylim(lo[1],hi[1]); ax.set_zlim(lo[2],hi[2]); ax.view_init(18,-60); ax.set_xlabel('X (mm)'); ax.set_ylabel('Y (mm)'); ax.set_zlabel('Z (mm)'); ax.set_title('Reduced-Hydro | t=%.3f ms | phase-rate=%.2f Hz | wobble=%.1f rad/s | Vt=%.1f mm/s'%(row.time_s*1000,row.robot_phase_rate_Hz,row.omega_wobble_rad_s,row.Vt_mm_s)); fig.canvas.draw(); ims.append(Image.frombytes('RGB',fig.canvas.get_width_height(),fig.canvas.tostring_rgb())); plt.close(fig)
out=HERE/'Wobble_F30_G6L45_ReducedHydro_WallOn_Free_8p333ms.gif'; ims[0].save(out,save_all=True,append_images=ims[1:],duration=60,loop=0); print('WROTE',out,'frames',len(ims))
