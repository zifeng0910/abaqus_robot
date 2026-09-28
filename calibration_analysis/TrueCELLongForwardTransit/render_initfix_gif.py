import sys, numpy as np
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
p=Path(sys.argv[1]); out=Path(sys.argv[2]); z=np.load(p); t=z['time']; xyz=z['coords']; ax=z['axis']; n=z['n']; o=z['origin']
def project(a):
 q=a-o; return np.c_[q@ax,q@n]
xy=np.array([project(a) for a in xyz]);
fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
for a,title in zip(axs,['Axial translation / transverse position','Robot geometry in global projection']): a.set_title(title); a.set_aspect('equal'); a.grid(alpha=.2)
qall=xy.reshape(-1,2); smin,smax=np.nanpercentile(qall[:,0],[1,99]); pad=.5
axs[0].set(xlabel='axial coordinate s (mm)',ylabel='transverse coordinate n (mm)',xlim=(smin-pad,smax+pad),ylim=(-.9,.9))
g=xyz.reshape(-1,3); xmin,xmax=np.nanpercentile(g[:,0],[1,99]); zmin,zmax=np.nanpercentile(g[:,2],[1,99])
axs[1].set(xlabel='global X (mm)',ylabel='global Z (mm)',xlim=(xmin-1,xmax+1),ylim=(zmin-1,zmax+1))
axs[0].axhline(.667345,color='#999999',lw=.8,ls='--'); axs[0].axhline(-.667345,color='#999999',lw=.8,ls='--')
sc=axs[0].scatter([],[],s=3,c='#176b87'); line,=axs[1].plot([],[],'.',ms=2,color='#b34d2e'); label=fig.text(.5,.01,'',ha='center')
def update(i):
 q=xy[i]; sc.set_offsets(q); line.set_data(xyz[i,:,0],xyz[i,:,2]); label.set_text(f't = {t[i]*1000:.3f} ms'); return sc,line,label
ani=FuncAnimation(fig,update,frames=len(t),interval=60,blit=False); ani.save(out,PillowWriter(fps=12),dpi=90); plt.close(fig); print(out, out.stat().st_size)
