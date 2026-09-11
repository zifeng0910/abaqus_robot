"""Render a fixed-camera four-panel CEL-FSI bend validation GIF."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image

job = sys.argv[1]
centerline_file = Path(sys.argv[2])
bend_start = float(sys.argv[3]) if len(sys.argv) > 3 else 13.49
bend_end = float(sys.argv[4]) if len(sys.argv) > 4 else 18.56
stride = max(int(sys.argv[5]) if len(sys.argv) > 5 else 1, 1)
cone_frame_azimuth_deg = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0
base = Path('output') / (job + '_bend_validation')
dyn = pd.read_csv(base / 'rp_centerline_contact_history.csv')
robot = pd.read_csv(base / 'robot_frames.csv')
pipe = pd.read_csv(base / 'pipe_wall_nodes.csv')
evf = pd.read_csv(base / 'evf_frames.csv')
cp = pd.read_csv(base / 'cpress_peaks.csv')
cl = pd.read_csv(centerline_file)
telemetry_path = Path(job + '_telemetry.csv')
telemetry = pd.read_csv(telemetry_path) if telemetry_path.exists() else pd.DataFrame()

cl_s = cl.arclength_mm.to_numpy(float)
cl_xyz = cl[['x_mm','y_mm','z_mm']].to_numpy(float)
cl_tan = np.gradient(cl_xyz, cl_s, axis=0, edge_order=1)
cl_tan /= np.maximum(np.linalg.norm(cl_tan, axis=1)[:, None], 1e-15)

def tangent_at_s(sq):
    q = np.clip(float(sq), cl_s[0], cl_s[-1])
    v = np.array([np.interp(q, cl_s, cl_tan[:, k]) for k in range(3)])
    return v / max(np.linalg.norm(v), 1e-15)

frames = sorted(dyn.frame.unique())[::stride]
all_xyz = np.vstack([
    pipe[['x_mm','y_mm','z_mm']].values,
    robot[['x_mm','y_mm','z_mm']].values,
    evf[['x_mm','y_mm','z_mm']].values,
])
mins, maxs = all_xyz.min(axis=0), all_xyz.max(axis=0)
span = np.maximum(maxs-mins, 1.0)
mins -= 0.05*span
maxs += 0.05*span

bcl = cl[(cl.arclength_mm >= bend_start) & (cl.arclength_mm <= bend_end)]
if bcl.empty:
    bcl = cl
bmins = bcl[['x_mm','y_mm','z_mm']].values.min(axis=0)-2.0
bmaxs = bcl[['x_mm','y_mm','z_mm']].values.max(axis=0)+2.0

r0 = robot[robot.frame == frames[0]]
r0_min = r0[['x_mm','y_mm']].values.min(axis=0)
r0_max = r0[['x_mm','y_mm']].values.max(axis=0)
traj = dyn[['x_mm','y_mm','z_mm']].values
evf_sum = evf.groupby('frame').evf.sum().reindex(frames).fillna(0.0)
evf_rel = evf_sum / max(evf_sum.iloc[0], 1.0e-30)
dyn_render = dyn[dyn.frame.isin(frames)]

def bbox(ax, lo, hi, **kw):
    x0,y0=lo; x1,y1=hi
    ax.plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],**kw)

imgs=[]
pngs=[]
for fi in frames:
    d = dyn[dyn.frame == fi].iloc[0]
    rb = robot[robot.frame == fi]
    ef = evf[evf.frame == fi]
    c = cp[cp.frame == fi].iloc[0]
    fig = plt.figure(figsize=(14,10), dpi=105)
    ax1=fig.add_subplot(2,2,1)
    ax2=fig.add_subplot(2,2,2,projection='3d')
    ax3=fig.add_subplot(2,2,3)
    ax4=fig.add_subplot(2,2,4)

    # Front / overall. All limits remain fixed for every frame.
    ax1.scatter(pipe.x_mm,pipe.y_mm,s=1.0,c='#3182bd',alpha=.30,label='Pipe wall')
    if not ef.empty:
        ax1.scatter(ef.x_mm,ef.y_mm,s=3,c=ef.evf,cmap='YlOrBr',vmin=0,vmax=1,alpha=.42,label='CEL EVF')
    ax1.scatter(rb.x_mm,rb.y_mm,s=6,c='#d62728',alpha=.85,label='Robot')
    ax1.plot(cl.x_mm,cl.y_mm,'k--',lw=1,label='Centerline')
    ax1.plot(traj[:fi+1,0],traj[:fi+1,1],c='#31a354',lw=2,label='RP trajectory')
    bbox(ax1,r0_min,r0_max,color='k',ls='--',lw=1)
    ax1.set(xlim=(mins[0],maxs[0]),ylim=(mins[1],maxs[1]),aspect='equal',title='Front / overall (fixed world camera)')
    ax1.legend(fontsize=7,loc='best')

    # Fixed 3-D world view.
    ax2.scatter(pipe.x_mm,pipe.y_mm,pipe.z_mm,s=.6,c='#6baed6',alpha=.18)
    if not ef.empty:
        ax2.scatter(ef.x_mm,ef.y_mm,ef.z_mm,s=2,c=ef.evf,cmap='YlOrBr',vmin=0,vmax=1,alpha=.35)
    ax2.scatter(rb.x_mm,rb.y_mm,rb.z_mm,s=5,c='#d62728',alpha=.8)
    ax2.plot(cl.x_mm,cl.y_mm,cl.z_mm,'k--',lw=1)
    ax2.plot(traj[:fi+1,0],traj[:fi+1,1],traj[:fi+1,2],c='#31a354',lw=2)
    # Magnetic diagnostics: B vector from the socket telemetry and the
    # authoritative +s gradient direction, both in the Abaqus frame.
    if not telemetry.empty:
        ti = int(np.argmin(np.abs(telemetry.t_s.to_numpy(float) - float(d.time_s))))
        tr = telemetry.iloc[ti]
        p0 = np.array([float(d.x_mm), float(d.y_mm), float(d.z_mm)])
        bvec = np.array([float(tr.get('Bx_aba_T', 0.0)),
                         float(tr.get('By_aba_T', 0.0)),
                         float(tr.get('Bz_aba_T', 0.0))])
        if np.linalg.norm(bvec) > 1e-15:
            bdir = bvec / np.linalg.norm(bvec)
            ax2.quiver(*p0, *(3.0*bdir), color='#d62728', linewidth=1.6,
                       arrow_length_ratio=0.18, label='B vector (Abaqus)')
        gdir = tangent_at_s(float(tr.get('driver_arc_mm', d.s_mm)))
        ax2.quiver(*p0, *(3.0*gdir), color='#2ca02c', linewidth=1.6,
                   arrow_length_ratio=0.18, label='+s gradient')
        if abs(cone_frame_azimuth_deg) > 1.0e-12:
            zref = np.array([0.0, 0.0, 1.0])
            nref = zref - np.dot(zref, gdir)*gdir
            nref /= max(np.linalg.norm(nref), 1.0e-15)
            beta = np.deg2rad(40.0)
            caxis = np.cos(beta)*gdir + np.sin(beta)*nref
            ang = np.deg2rad(cone_frame_azimuth_deg)
            caxis = (caxis*np.cos(ang) + np.cross(gdir,caxis)*np.sin(ang) +
                     gdir*np.dot(gdir,caxis)*(1.0-np.cos(ang)))
            caxis /= max(np.linalg.norm(caxis), 1.0e-15)
            ax2.quiver(*p0, *(3.0*caxis), color='#6a3d9a', linewidth=1.8,
                       arrow_length_ratio=0.18, label=f'cχ (χ={cone_frame_azimuth_deg:.0f}°)')
        ax2.text2D(0.02, 0.96, f'red=B, green=+s, purple=cχ (χ={cone_frame_azimuth_deg:.0f}°)', transform=ax2.transAxes,
                    fontsize=7)
    ax2.set(xlim=(mins[0],maxs[0]),ylim=(mins[1],maxs[1]),zlim=(mins[2],maxs[2]),title='3-D CEL-FSI trajectory')
    ax2.view_init(elev=24,azim=-58)

    # Bend detail, also fixed.
    ax3.scatter(pipe.x_mm,pipe.y_mm,s=1.2,c='#3182bd',alpha=.28)
    if not ef.empty:
        ax3.scatter(ef.x_mm,ef.y_mm,s=4,c=ef.evf,cmap='YlOrBr',vmin=0,vmax=1,alpha=.45)
    ax3.scatter(rb.x_mm,rb.y_mm,s=8,c='#d62728',alpha=.9)
    ax3.plot(bcl.x_mm,bcl.y_mm,'k--',lw=1.4)
    ax3.plot(traj[:fi+1,0],traj[:fi+1,1],c='#31a354',lw=2)
    ax3.set(xlim=(bmins[0],bmaxs[0]),ylim=(bmins[1],bmaxs[1]),aspect='equal',title=f'Bend detail: s={d.s_mm:.3f} mm')

    # Quantitative timeline.
    ax4.plot(dyn_render.time_s,dyn_render.s_mm,label='centerline s (mm)',c='#31a354')
    ax4.plot(dyn_render.time_s,dyn_render.axis_tangent_angle_deg,label='axis/tangent angle (deg)',c='#756bb1')
    ax4.plot(cp.time_s,cp.cpress_max_mpa,label='max CPRESS (MPa)',c='#de2d26',alpha=.8)
    ax4.plot(dyn_render.time_s,evf_rel.values,label='EVF total / initial',c='#e6550d')
    ax4.axhline(bend_start,color='gray',ls='--',lw=.8)
    ax4.axhline(bend_end,color='gray',ls=':',lw=.8)
    ax4.axhline(15,color='#756bb1',ls='--',lw=.8)
    ax4.axvline(d.time_s,color='k',lw=1)
    ax4.set(xlim=(dyn.time_s.min(),dyn.time_s.max()),title='Motion / contact / CEL diagnostics',xlabel='time (s)')
    ax4.legend(fontsize=7,loc='best')

    delta_s = float(d.s_mm - dyn.s_mm.iloc[0])
    fig.suptitle(f'CEL-FSI route | {job} | t={d.time_s:.6f} s | CPRESS={c.cpress_max_mpa:.4g} MPa | Δs={delta_s:+.4f} mm')
    fig.tight_layout()
    p=base/(job+f'_bend_frame_{int(fi):05d}.png')
    fig.savefig(p)
    plt.close(fig)
    pngs.append(p)
    imgs.append(Image.open(p).convert('RGB'))

gif = Path(job + '_CEL_FSI_bend_validation.gif')
if imgs:
    imgs[0].save(gif,save_all=True,append_images=imgs[1:],duration=70,loop=0,optimize=False)
    final_png=Path(job + '_CEL_FSI_bend_validation_final.png')
    imgs[-1].save(final_png)
    print('GIF',gif.resolve(),'frames',len(imgs))
    print('PNG',final_png.resolve())
