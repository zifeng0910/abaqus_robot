"""Python-only publication figures for the exact FreeCAD validation."""
from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'ReducedHydro_geometry_L1800_validation'
BLUE,ORANGE,RED,GREEN,PURPLE,GREY='#0072B2','#D55E00','#CC3311','#009E73','#882255','#5B6573'
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
 'font.size':7.2,'axes.titlesize':8,'axes.labelsize':7.2,'xtick.labelsize':7.2,'ytick.labelsize':7.2,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.8,'legend.frameon':False,
 'pdf.fonttype':42,'svg.fonttype':'none'})


def save(fig,name):
    fig.savefig(HERE/f'{name}.png',dpi=600,bbox_inches='tight')
    fig.savefig(HERE/f'{name}.pdf',bbox_inches='tight')
    fig.savefig(HERE/f'{name}.svg',bbox_inches='tight')
    fig.savefig(HERE/f'{name}.tiff',dpi=600,bbox_inches='tight')
    plt.close(fig)


def profile(points,axis,rp):
    q=points-rp; x=q@axis; r=np.linalg.norm(q-np.outer(x,axis),axis=1)
    return x,r


def main():
    ident=json.loads((HERE/'freecad_preflight_identity.json').read_text())
    s=json.loads((HERE/'freecad_8p333_summary.json').read_text())
    old_s=json.loads((OLD/'L1800_8p333_summary.json').read_text())
    d=pd.read_csv(HERE/'freecad_8p333_directed_tilt.csv');od=pd.read_csv(OLD/'L1800_8p333_true_axis.csv')
    op=pd.read_csv(OLD/'L1800_8p333_local_phase.csv');p=pd.read_csv(HERE/'freecad_8p333_local_phase.csv')
    f=pd.read_csv(HERE/'freecad_8p333_field_orientation.csv');g=pd.read_csv(HERE/'freecad_8p333_exact_gap.csv')
    c=pd.read_csv(HERE/'freecad_8p333_contact_events.csv');b=pd.read_csv(HERE/'freecad_8p333_bridge_timeline.csv')
    q=pd.read_csv(HERE/'freecad_8p333_torque.csv');tr=pd.read_csv(HERE/'freecad_8p333_canonical_translation.csv')

    geom=pd.read_csv(HERE/'freecad_vs_abaqus_geometry_audit.csv')
    fig,ax=plt.subplots(figsize=(3.54,2.45),constrained_layout=True)
    x=np.arange(len(geom));ax.bar(x-.18,geom.CAD,.36,color=GREY,label='FreeCAD STEP');ax.bar(x+.18,geom.Abaqus_mesh,.36,color=BLUE,label='Abaqus mesh')
    ax.set(xticks=x,xticklabels=geom.metric,ylabel='Value (mm or mm³)',title='Abaqus mesh preserves the authoritative CAD envelope');ax.legend();save(fig,'freecad_vs_abaqus_geometry')

    ns=pd.read_csv(HERE/'freecad_robot_surface_nodes_global.csv');newp=ns[['x_mm','y_mm','z_mm']].to_numpy();axis=np.array([.9647382600216,-.1188742372140,.2348382536499]);rp=np.array(ident['RP_mm']);xn,rn=profile(newp,axis,rp)
    from analyze_freecad_8p333 import mesh_properties
    om=mesh_properties((OLD/'Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083.inp').read_text())[0]['Robot_SOLID'];on=np.array(list(om['nodes'].values()))+om['shift'];xo,ro=profile(on,axis,rp)
    fig,ax=plt.subplots(figsize=(3.54,2.2),constrained_layout=True);ax.scatter(xo,ro,s=.8,color=GREY,alpha=.35,label='Old scaled mesh');ax.scatter(xn,rn,s=.8,color=BLUE,alpha=.35,label='Exact CAD mesh');ax.set(xlabel='HEAD-to-TAIL coordinate (mm)',ylabel='Radius (mm)',title='Rounded-head CAD changes volume, not nominal L/D');ax.legend();save(fig,'old_scaled_vs_freecad_robot_shape')

    olddot=np.einsum('ij,ij->i',od[['axis_x','axis_y','axis_z']],od[['tangent_x','tangent_y','tangent_z']]);old_direct=np.degrees(np.arccos(np.clip(olddot,-1,1)))
    fig,ax=plt.subplots(figsize=(3.54,2.35),constrained_layout=True);ax.plot(od.time_s*1e3,old_direct,color=GREY,lw=.9,label='Old scaled L1800');ax.plot(d.time_s*1e3,d.directed_tilt_deg,color=BLUE,lw=1.1,label='Exact CAD L1800');ax.axhline(90,color=RED,lw=.7,ls='--');ax.set(xlabel='Time (ms)',ylabel='Directed tilt (deg)',title='Directed tilt exposes HEAD-TAIL reversal');ax.legend();save(fig,'old_scaled_vs_freecad_directed_tilt')

    fig,ax=plt.subplots(figsize=(3.54,2.35),constrained_layout=True);ax.plot(op.time_s*1e3,op.phi_robot_advance_deg,color=GREY,lw=.9,label='Old scaled L1800');ax.plot(p.time_s*1e3,p.phi_robot_advance_deg,color=BLUE,lw=1.1,label='Exact CAD L1800');ax.set(xlabel='Time (ms)',ylabel='Robot phase advance (deg)',title='Exact head tests whether late phase reversal persists');ax.legend();save(fig,'old_scaled_vs_freecad_phase')

    fig,ax=plt.subplots(figsize=(3.54,2.4),constrained_layout=True);ax.plot(f.time_s*1e3,d.directed_tilt_deg,color=BLUE,lw=1,label='Robot directed tilt');ax.plot(f.time_s*1e3,f.theta_B_deg,color=ORANGE,lw=1,label='Field tilt');ax.plot(f.time_s*1e3,f.robot_B_misalignment_deg,color=PURPLE,lw=.8,label='Robot-field misalignment');ax.set(xlabel='Time (ms)',ylabel='Angle (deg)',title='Robot overshoot is separated from field orientation');ax.legend();save(fig,'robot_vs_field_tilt_freecad')

    fig,axs=plt.subplots(3,1,figsize=(3.54,4.0),sharex=True,constrained_layout=True,gridspec_kw={'height_ratios':[2,1,1]});axs[0].plot(g.time_s*1e3,g.gap_um,color=GREY,lw=.65);axs[0].axhline(0,color='black',lw=.6);axs[0].set(ylabel='Exact gap (µm)',title='Exact gap distinguishes impacts from opposing-wall bridge');axs[1].vlines((c.start_s+c.end_s)*.5e3,0,c.peak_force_N,color=RED,lw=.8);axs[1].set(ylabel='Peak force (N)');axs[2].fill_between(b.time_s*1e3,0,b.opposing_bridge,step='post',color=PURPLE,alpha=.75);axs[2].set(xlabel='Time (ms)',ylabel='Bridge',ylim=(0,1.05));save(fig,'gap_contact_bridge_freecad')

    fig,ax=plt.subplots(figsize=(3.54,2.35),constrained_layout=True);ax.plot(q.time_s*1e3,q.Tmag_norm_Nmm,color=ORANGE,lw=.9,label='Magnetic');ax.plot(q.time_s*1e3,q.Thydro_norm_Nmm,color=GREEN,lw=.9,label='Hydrodynamic');ax.set_yscale('symlog',linthresh=1e-6);ax.set(xlabel='Time (ms)',ylabel='Torque (N mm)',title='Torque audit tests drive collapse and dissipation');ax.legend();save(fig,'magnetic_vs_hydro_torque_freecad')

    fig,axs=plt.subplots(2,1,figsize=(3.54,3.1),sharex=True,constrained_layout=True);axs[0].plot(tr.time_s*1e3,tr.delta_s_mm,color=BLUE,lw=1);axs[0].set(ylabel='Δs (mm)',title='Canonical translation remains a diagnostic outcome');axs[1].plot(tr.time_s*1e3,tr.Vt_mm_s,color=GREY,lw=.9);axs[1].set(xlabel='Time (ms)',ylabel='Tangent velocity (mm/s)');save(fig,'canonical_translation_freecad')

    reversals=pd.read_csv(HERE/'freecad_head_tail_reversal_events.csv')
    fig,ax=plt.subplots(figsize=(3.54,2.55),constrained_layout=True)
    ax.plot(d.time_s*1e3,d.directed_tilt_deg,color=BLUE,lw=1.05)
    ax.axhline(90,color=GREY,lw=.7,ls='--')
    for row in reversals.itertuples():
        ax.axvline(row.crossing_time_ms,color=RED,lw=.65,alpha=.85)
        ax.text(row.crossing_time_ms,94,str(row.event_number),ha='center',va='bottom',fontsize=6,color=RED)
    ax.fill_between(d.time_s*1e3,90,d.directed_tilt_deg,where=d.directed_tilt_deg<90,
                    color=ORANGE,alpha=.12,label='Polarity reversed')
    ax.set(xlabel='Time (ms)',ylabel='Directed tilt (deg)',ylim=(0,180),
           title='Four transverse crossings reverse HEAD-TAIL polarity')
    ax.legend(loc='lower left');save(fig,'directed_tilt_and_head_tail_reversals')

    sensitivity=pd.read_csv(HERE/'freecad_bridge_threshold_sensitivity.csv')
    fig,ax=plt.subplots(figsize=(3.54,2.4),constrained_layout=True)
    colors=[PURPLE,BLUE,GREY,GREY]
    bars=ax.bar(sensitivity.threshold_um.astype(str),sensitivity.longest_bridge_ms,
                color=colors,width=.62)
    ax.axhline(.5,color=RED,lw=.75,ls='--',label='Persistent screen (0.5 ms)')
    for bar,value in zip(bars,sensitivity.longest_bridge_ms):
        ax.text(bar.get_x()+bar.get_width()/2,value+.018,f'{value:.3f}',ha='center',va='bottom',fontsize=6)
    ax.set(xlabel='Near-wall threshold (µm)',ylabel='Longest opposing bridge (ms)',ylim=(0,.64),
           title='Persistent bridge is specific to the 20 µm criterion')
    ax.legend(loc='upper right');save(fig,'bridge_threshold_sensitivity')


if __name__=='__main__':main()
