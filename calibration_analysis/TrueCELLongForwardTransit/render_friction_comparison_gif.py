from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
import numpy as np

root=Path(r'J:/abaqusfangzhen/abaqus_robot/calibration_analysis/TrueCELLongForwardTransit')
cases=[('INITFIX20','initialization_closure/robot_motion_frames.npz','#555555'),('BOTH FRIC0','initialization_closure/fluidfric0/robot_motion_frames.npz','#b34d2e'),('RF FRIC0','initialization_closure/rf_fric0/robot_motion_frames.npz','#176b87'),('WF FRIC0','initialization_closure/wf_fric0/robot_motion_frames.npz','#6a4c93')]
data=[]
for name,rel,color in cases:
 z=np.load(root/rel); t=z['time']; xyz=z['coords']; ax=z['axis']; n=z['n']; o=z['origin']; q=np.array([np.c_[(a-o)@ax,(a-o)@n] for a in xyz]); data.append((name,t,q,color))
T=min(len(x[1]) for x in data); fig,axs=plt.subplots(2,2,figsize=(9,6),constrained_layout=True); axs=axs.ravel()
for ax,(name,t,q,color) in zip(axs,data):
 ax.set_title(name); ax.set_aspect('equal'); ax.set_xlim(-9,2); ax.set_ylim(-.9,.9); ax.set_xlabel('axial s (mm)'); ax.set_ylabel('n (mm)'); ax.grid(alpha=.2); ax.axhline(.667345,color='#999',lw=.6,ls='--'); ax.axhline(-.667345,color='#999',lw=.6,ls='--')
lines=[]; dots=[]
for ax,(name,t,q,color) in zip(axs,data):
 ln,=ax.plot([],[],color=color,lw=1.0); dot=ax.scatter([],[],s=3,color=color); lines.append(ln); dots.append(dot)
label=fig.text(.5,.01,'',ha='center')
def update(i):
 for (name,t,q,color),ln,dot in zip(data,lines,dots):
  ln.set_data(q[:i+1,:,0].mean(1),q[:i+1,:,1].mean(1)); dot.set_offsets(q[i]);
 label.set_text(f't = {data[0][1][i]*1000:.3f} ms'); return lines+dots+[label]
ani=FuncAnimation(fig,update,frames=T,interval=70,blit=False); out=root/'CEL_FLUID_CONTACT_FRICTION_comparison.gif'; ani.save(out,PillowWriter(fps=12),dpi=90); plt.close(fig); print(out,out.stat().st_size)
