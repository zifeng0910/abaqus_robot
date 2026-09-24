"""Analyze the sole F80 run, including its intentionally preserved early-stop ODB."""
from __future__ import annotations

import csv, json, zipfile
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from analyze_f100_candidate import longest_negative
from analyze_f100_g2p20_final import geometry, kinematics, pose, setup_axis

ROOT=Path(__file__).resolve().parents[1]
JOB='TRUECEL_B0P11_G2P20_A14P5_F80_FAST'
F100_JOB='TRUECEL_B0P11_G2P20_A14P5_F100_FAST'
CASE=ROOT/'case'/JOB
F100=ROOT/'case'/F100_JOB
PERIOD=.0125

def summarize_cycle(t,s,v,k,period):
    lo=(k-1)*period; hi=k*period
    grid=np.r_[lo,t[(t>lo)&(t<hi)],hi]
    x=np.interp(grid,t,s); y=np.interp(grid,t,v)
    return {'cycle':k,'t_start_s':lo,'t_end_s':hi,
            'delta_s_mm':float(x[-1]-x[0]),
            'mean_v_s_mm_s':float((x[-1]-x[0])/period),
            'end_v_s_mm_s':float(y[-1]),
            'minimum_v_s_mm_s':float(y.min()),
            'maximum_backward_excursion_mm':float((np.maximum.accumulate(x)-x).max()),
            'negative_velocity_duration_s':float(longest_negative(grid,y))}

def draw_single(data,path):
    ident,t,u,ur,s,v=data; rp,c,n,rel,colors,pc=geometry(CASE,ident)
    frames=np.linspace(0,min(.025,t[-1]),81)
    fig,ax=plt.subplots(figsize=(11,4.3),layout='constrained')
    setup_axis(ax,ident,'TRUE-CEL F80 G=2.20 mT — early failure at Cycle 2')
    p=pose(0,t,u,ur,rp,rel)
    sc=ax.scatter(float(ident['s_start_mm'])+(p-rp)@c,(p-pc)@n,s=5,c=colors,linewidths=0)
    tx=ax.text(.985,.985,'',transform=ax.transAxes,va='top',ha='right',family='monospace',
               bbox={'boxstyle':'round','facecolor':'white','alpha':.9})
    def update(i):
        q=frames[i]; p=pose(q,t,u,ur,rp,rel)
        sc.set_offsets(np.c_[float(ident['s_start_mm'])+(p-rp)@c,(p-pc)@n])
        valid=t<=q; back=float((np.maximum.accumulate(s[valid])-s[valid]).max())
        tx.set_text(f't={q*1000:6.2f} ms  cycle={min(2,int(q/PERIOD)+1)}\n'
                    f'delta_s={np.interp(q,t,s):+.5f} mm  v_s={np.interp(q,t,v):+.2f} mm/s\n'
                    f'MAX_BACKTRACK={back:.5f} mm')
        return sc,tx
    FuncAnimation(fig,update,frames=len(frames),interval=60).save(path,PillowWriter(fps=16),dpi=96)
    plt.close(fig)

def draw_sync(f100,f80,path):
    d=[]
    fig,axes=plt.subplots(2,1,figsize=(11,7.2),layout='constrained')
    for ax,case,data,title in ((axes[0],F100,f100,'F100, G=2.20 mT'),
                               (axes[1],CASE,f80,'F80, G=2.20 mT')):
        ident,t,u,ur,s,v=data; rp,c,n,rel,colors,pc=geometry(case,ident)
        setup_axis(ax,ident,title)
        p=pose(0,t,u,ur,rp,rel)
        sc=ax.scatter(float(ident['s_start_mm'])+(p-rp)@c,(p-pc)@n,s=5,c=colors,linewidths=0)
        tx=ax.text(.985,.985,'',transform=ax.transAxes,va='top',ha='right',family='monospace',
                   bbox={'boxstyle':'round','facecolor':'white','alpha':.9})
        d.append((sc,tx,ident,t,u,ur,s,v,rp,c,n,rel,pc))
    # Both panels show exactly two drive cycles at the same local phase.
    phase=np.linspace(0,2,81)
    def update(i):
        q=phase[i]; out=[]
        for j,(sc,tx,ident,t,u,ur,s,v,rp,c,n,rel,pc) in enumerate(d):
            period=.01 if j==0 else .0125
            ti=min(float(q*period),t[-1]); p=pose(ti,t,u,ur,rp,rel)
            sc.set_offsets(np.c_[float(ident['s_start_mm'])+(p-rp)@c,(p-pc)@n])
            valid=t<=ti; back=float((np.maximum.accumulate(s[valid])-s[valid]).max())
            tx.set_text(f'phase={q:.2f} cycles  t={ti*1000:.2f} ms\n'
                        f'delta_s={np.interp(ti,t,s):+.5f} mm  v_s={np.interp(ti,t,v):+.2f} mm/s\n'
                        f'back={back:.5f} mm')
            out.extend((sc,tx))
        return out
    fig.suptitle('Equivalent cycle phase — same geometry and G')
    FuncAnimation(fig,update,frames=len(phase),interval=60).save(path,PillowWriter(fps=16),dpi=96)
    plt.close(fig)

def main():
    f80=kinematics(CASE,False); f100=kinematics(F100,True)
    _,t,u,ur,s,v=f80
    event=(CASE/'f80_event.txt').read_text(encoding='ascii').splitlines()[0].strip()
    if event!='F80_RECOIL_FAIL_CYCLE':raise RuntimeError(f'Unexpected stop reason {event}')
    if not .02499<=t[-1]<=.02501:raise RuntimeError(f'Unexpected ODB end {t[-1]}')
    cycles=[summarize_cycle(t,s,v,k,PERIOD) for k in (1,2)]
    back=float((np.maximum.accumulate(s)-s).max())
    metrics={'candidate':JOB,'classification':'F80_RECOIL_FAIL','event':event,
             'full_cycles':2,'total_time_s':float(t[-1]),'fresh_start':True,
             'only_physics_change':'frequency 100 -> 80 Hz',
             'cycle_metrics':cycles,'max_backtrack_mm':back,
             'mean_speed_mm_s':float((s[-1]-s[0])/(t[-1]-t[0])),
             'total_delta_s_mm':float(s[-1]-s[0]),'first_negative_cycle_index':2,
             'stage1_three_cycle_complete':False,'five_cycle_continuation':False,
             'solver_status':'successful controlled early halt'}
    csv_path=ROOT/f'{JOB}_CYCLE_SUMMARY.csv'
    with csv_path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(cycles[0]));w.writeheader();w.writerows(cycles)
    f100m=json.loads((ROOT/f'{F100_JOB}_METRICS.json').read_text())
    comparison=[]
    for label,m in [('F100 G2.20',f100m),('F80 G2.20',metrics)]:
        row={'candidate':label}
        for k in range(1,6):
            row[f'cycle{k}_delta_s_mm']=next((r['delta_s_mm'] for r in m['cycle_metrics'] if r['cycle']==k),None)
        row.update({'max_backtrack_mm':m['max_backtrack_mm'],'mean_speed_mm_s':m['mean_speed_mm_s'],
                    'first_negative_cycle_index':next((r['cycle'] for r in m['cycle_metrics'] if r['delta_s_mm']<0),None)})
        comparison.append(row)
    comp_csv=ROOT/'F100_G2P20_VS_F80_G2P20_COMPARISON.csv'
    with comp_csv.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(comparison[0]));w.writeheader();w.writerows(comparison)
    fig,ax=plt.subplots(figsize=(7.5,4.4),layout='constrained')
    for row in comparison:
        values=[row[f'cycle{k}_delta_s_mm'] for k in range(1,6)]
        x=[k for k,val in enumerate(values,1) if val is not None]
        y=[val for val in values if val is not None]
        ax.plot(x,y,'o-',label=row['candidate'])
    ax.axhline(0,color='gray',lw=.8);ax.set(xlim=(.7,5.3),xticks=range(1,6),
        xlabel='Complete cycle number',ylabel='Net axial displacement per cycle (mm)',
        title='F80 reversal occurs one cycle earlier than F100')
    ax.legend();ax.grid(alpha=.2)
    plot=ROOT/'F100_G2P20_VS_F80_G2P20_CYCLE_TREND.png';fig.savefig(plot,dpi=160);plt.close(fig)
    gif=ROOT/f'{JOB}.gif';sync=ROOT/'F100_G2P20_VS_F80_G2P20_PHASE_SYNC.gif'
    draw_single(f80,gif);draw_sync(f100,f80,sync)
    zip_path=ROOT/f'{JOB}_GIFS.zip'
    with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
        for p in (gif,sync,plot,csv_path,comp_csv):z.write(p,p.name)
    metrics.update({'gif':gif.name,'phase_sync_gif':sync.name,'gif_zip':zip_path.name,
                    'trend_plot':plot.name,'comparison_csv':comp_csv.name})
    (ROOT/f'{JOB}_METRICS.json').write_text(json.dumps(metrics,indent=2)+'\n')
    report=f'''# F80 / G2.20 coarse TRUE-CEL result\n\nClassification: **F80_RECOIL_FAIL**.\n\nThe fresh-start solve halted after the second complete 12.5-ms cycle. C1 advanced {cycles[0]['delta_s_mm']:+.6f} mm; C2 reversed {cycles[1]['delta_s_mm']:+.6f} mm. MAX_BACKTRACK reached {back:.6f} mm. The solver completed its controlled stop, so this is a physical-screen failure rather than numerical invalidity.\n\nF100 / G2.20 had C1/C2 positive and first negative cycle C3. Lowering to F80 moved reversal to C2; it did not eliminate late-cycle reversal. No further frequency tuning is authorized by this task.\n\nThe private contact history provides whole robot and pipe General Contact resultants, not a pair-isolated direct robot-wall resultant; these are not interchangeable.\n'''
    report+='''\nPre-run F100 rocking check: `ROCKING_REMAINS_PHASE_COHERENT_THROUGH_CYCLE3`. Fundamental rocking amplitudes C1/C2/C3 were 12.937°, 13.856°, 14.477°; C2−C1 and C3−C1 phase lags were +3.03° and +3.87°. RMS phase-matched alpha errors were 2.563° and 2.363°; p95 absolute rocking-rate differences were 259.06 and 225.73 rad/s. At F100 C3 negative-velocity onset, signed TAIL/HEAD wall gaps were +0.00218/+0.04441 mm; at maximum backtrack, −0.000082/+0.04476 mm. The saved history lacks a pair-isolated robot-wall force.\n\nScreen conditions: reduced-sound-speed coarse CEL (`c0=100000 mm/s`, 44×20×20, Explicit scale factor 0.4). The F80 run stopped at 25 ms by the complete-cycle failure gate. The planned third and five-cycle continuations were not run.\n'''
    (ROOT/f'{JOB}_REPORT.md').write_text(report)
    print(json.dumps(metrics,indent=2))

if __name__=='__main__':main()
