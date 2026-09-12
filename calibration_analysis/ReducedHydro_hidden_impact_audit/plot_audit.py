"""Quantitative-grid evidence: sparse frames missed a resolved contact pulse.

Contract: Python; 183 mm report figures; raw deterministic traces, no statistical
replicates or smoothing; original and diagnostic runs are explicitly labelled.
Blue=momentum, orange=contact, grey=magnetic/hydro. PNG/PDF editable text outputs.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,'axes.labelsize':9,'axes.titlesize':10,'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':300})
BLUE='#21618C';ORANGE='#C76728';GRAY='#616A6B'
def save(fig,name):
    fig.savefig(HERE/(name+'.png'),dpi=300)
    fig.savefig(HERE/(name+'.pdf'))
    fig.savefig(HERE/(name+'.svg'))
    plt.close(fig)
def main():
    d=pd.read_csv(HERE/'momentum_plot_source.csv');g=pd.read_csv(HERE/'first_impact_gap_history.csv');sp=pd.read_csv(HERE/'verified_sparse_wall_gap.csv');cp=pd.read_csv(HERE/'diagnostic'/'contact_fields.csv')
    t=d.time_s.to_numpy()*1e3; ev=g.time_s.to_numpy()*1e3
    assert np.all(np.diff(t)>0) and np.all(np.diff(ev)>0)
    fig,ax=plt.subplots(2,1,figsize=(7.2,5.2),sharex=True,layout='constrained')
    ax[0].plot(t,d.speed_mm_s,color=BLUE,label='Direct COM speed (0.1-us history)')
    ax[0].set_ylabel('COM speed (mm/s)');ax[0].legend(loc='lower right')
    ax[1].plot(ev,g.gap_um,color=ORANGE,label='Corrected wall gap (rigid-node reconstruction)')
    sel=sp[(sp.time_s>=.00098)&(sp.time_s<=.00108)]
    ax[1].scatter(sel.time_s*1e3,sel.gap_um,color=GRAY,s=24,zorder=3,label='Original 50-us field frames')
    ax[1].axhline(0,color='black',lw=.6);ax[1].set_ylabel('Signed gap (µm)');ax[1].set_xlabel('Time (ms)');ax[1].legend(loc='upper left');ax[1].set_xlim(.98,1.08)
    fig.suptitle('A sub-frame impact separates two non-contact field samples')
    save(fig,'velocity_jump_vs_gap')
    fig,axs=plt.subplots(3,1,figsize=(7.2,6),sharex=True,layout='constrained')
    for i,c in enumerate('xyz'):
        axs[i].plot(t,d['P_'+c]*1e6,color=BLUE,label='m ΔV')
        axs[i].plot(t,(d['Jmag_'+c]+d['Jhydro_'+c]+d['Jwall_'+c])*1e6,'--',color=ORANGE,label='Jmag + Jhydro + Jwall')
        axs[i].plot(t,(d['Jmag_'+c]+d['Jhydro_'+c])*1e6,':',color=GRAY,label='Jmag + Jhydro only')
        axs[i].set_ylabel(c.upper()+' impulse (µN s)');axs[i].axhline(0,color='black',lw=.4)
    axs[0].legend(loc='upper left');axs[-1].set_xlim(0,1.1);axs[-1].set_xlabel('Time (ms)');fig.suptitle('Independent wall-force output closes all three momentum components')
    save(fig,'momentum_balance_0_to_1p1ms')
    miss=np.linalg.norm(d[['P_'+c for c in 'xyz']].to_numpy()-d[['Jmag_'+c for c in 'xyz']].to_numpy()-d[['Jhydro_'+c for c in 'xyz']].to_numpy(),axis=1)
    fig,ax=plt.subplots(2,1,figsize=(7.2,4.8),sharex=True,layout='constrained')
    ax[0].plot(t,miss*1e6,color=BLUE,label='Missing impulse if wall is omitted');ax[0].set_ylabel('Impulse norm (µN s)');ax[0].legend()
    ax[1].plot(ev,g.gap_um,color=ORANGE);ax[1].axhline(0,color='black',lw=.6);ax[1].set_ylabel('Signed gap (µm)');ax[1].set_xlim(.998,1.01);ax[1].set_xlabel('Time (ms)');fig.suptitle('Momentum deficit grows during the near-wall rebound')
    visible=g.loc[(ev>=.998)&(ev<=1.01),'gap_um'];ax[1].set_ylim(visible.min()-.4,visible.max()+.4)
    save(fig,'missing_impulse_vs_nearwall')
    force=pd.read_csv(HERE/'first_impact_contact_force_history.csv');f=force[['Ftotal_robot_'+c+'_N' for c in 'xyz']].to_numpy()
    fig,ax=plt.subplots(2,1,figsize=(7.2,5.2),sharex=True,layout='constrained')
    ax[0].plot(force.time_s*1e3,np.linalg.norm(f,axis=1),color=ORANGE,label='Robot-side CFN + CFS (every increment)');ax[0].set_ylabel('Wall force (N)');ax[0].legend()
    ax[1].plot(t,np.linalg.norm(d[['Jwall_'+c for c in 'xyz']],axis=1)*1e6,color=ORANGE,label='Cumulative independent Jwall')
    ax[1].plot(t,np.linalg.norm(d[['P_'+c for c in 'xyz']],axis=1)*1e6,'--',color=BLUE,label='m ΔV norm')
    ax[1].set_ylabel('Impulse norm (µN s)');ax[1].set_xlabel('Time (ms)');ax[1].set_xlim(.999,1.01);ax[1].legend(loc='lower right');fig.suptitle('Microsecond wall impulse, not direct magnetic acceleration')
    save(fig,'first_impact_force_impulse')
    fig,ax=plt.subplots(3,1,figsize=(7.2,6),sharex=True,layout='constrained')
    ax[0].plot(ev,g.gap_um,color=ORANGE);ax[0].axhline(0,color='black',lw=.5);ax[0].set_ylabel('Signed gap (µm)')
    ax[1].plot(t,d.speed_mm_s,color=BLUE);ax[1].set_ylabel('COM speed (mm/s)')
    ax[2].plot(cp.time_s*1e3,cp.CPRESS_MPa,'o-',ms=3,color=GRAY);ax[2].set_ylabel('GC CPRESS (MPa)');ax[2].set_xlabel('Time (ms)');ax[2].set_xlim(.999,1.008)
    visible=g.loc[(ev>=.999)&(ev<=1.008),'gap_um'];ax[0].set_ylim(visible.min()-.4,visible.max()+.4)
    fig.suptitle('First impact: penetration, velocity jump and 0.5-us contact frames')
    save(fig,'first_impact_gap_velocity')
if __name__=='__main__':main()
