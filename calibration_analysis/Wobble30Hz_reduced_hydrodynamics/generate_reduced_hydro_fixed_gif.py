from pathlib import Path
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]; JOB='Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083'
rf=pd.read_csv(ROOT/'output'/(JOB+'_bend_validation')/'robot_frames.csv'); dyn=pd.read_csv(HERE/'reduced_hydro_fixed_true_axis_phase.csv'); wn=pd.read_csv(ROOT/(JOB+'_pipe_wall_nodes_exact.csv')); tri=pd.read_csv(ROOT/(JOB+'_pipe_wall_triangles_exact.csv'))
nd={int(r.node):np.array([r.x_mm,r.y_mm,r.z_mm],float) for _,r in wn.iterrows()}; verts=[[nd[int(r[n])] for n in ('n1','n2','n3','n4')] for _,r in tri.iterrows()]
cv=pd.read_csv(ROOT/'CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_EXTSURF_centerline_odb.csv')[['x_mm','y_mm','z_mm']].to_numpy(float)
allpts=np.vstack([wn[['x_mm','y_mm','z_mm']].to_numpy(float),cv]); lo=allpts.min(0)-.5; hi=allpts.max(0)+.5; ims=[]
for k in np.linspace(0,int(rf.frame.max()),84).astype(int):
    g=rf[rf.frame==k]; P=g[['x_mm','y_mm','z_mm']].to_numpy(float); row=dyn.iloc[np.argmin(abs(dyn.time_s-g.time_s.iloc[0]))]; hp=np.array([row.HEAD_x_mm,row.HEAD_y_mm,row.HEAD_z_mm]); tp=np.array([row.TAIL_x_mm,row.TAIL_y_mm,row.TAIL_z_mm])
    fig=plt.figure(figsize=(8,6)); ax=fig.add_subplot(111,projection='3d'); ax.add_collection3d(Poly3DCollection(verts,alpha=.10,facecolor='gray',edgecolor='none')); ax.plot(cv[:,0],cv[:,1],cv[:,2],'k-',lw=1); ax.plot([tp[0],hp[0]],[tp[1],hp[1]],[tp[2],hp[2]],'r-',lw=3); ax.scatter(*hp,c='r',s=28,label='HEAD'); ax.scatter(*tp,c='b',s=28,label='TAIL'); ax.legend(loc='upper left')
    ax.set_xlim(lo[0],hi[0]); ax.set_ylim(lo[1],hi[1]); ax.set_zlim(lo[2],hi[2]); ax.view_init(18,-60); ax.set_xlabel('X (mm)'); ax.set_ylabel('Y (mm)'); ax.set_zlabel('Z (mm)'); ax.set_title('Reduced-Hydro fixed | t=%.3f ms | local phase-rate=%.2f Hz | |omega_wobble|=%.1f rad/s | Vt=%.1f mm/s | gap=%.2f um'%(row.time_s*1000,row.f_robot_phase_equiv_Hz,row.omega_wobble_norm_rad_s,row.Vt_direct_mm_s,row.projection_residual_mm*1000)); fig.canvas.draw(); ims.append(Image.frombytes('RGB',fig.canvas.get_width_height(),fig.canvas.tostring_rgb())); plt.close(fig)
out=HERE/'Wobble_F30_ReducedHydroFixed_WallOn_8p333.gif'; ims[0].save(out,save_all=True,append_images=ims[1:],duration=60,loop=0); print('WROTE',out,'frames',len(ims))
