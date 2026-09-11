"""Annotate a cone-frame-azimuth GIF without changing the fixed camera."""
from pathlib import Path
import sys
import pandas as pd
import numpy as np
from PIL import Image, ImageDraw
job=sys.argv[1] if len(sys.argv)>1 else 'Wobble_F30_ConeFrameAzimuth_Chi134_WallOn_Free_003'
base=Path('output')/(job+'_bend_validation')
frames=sorted(base.glob(job+'_bend_frame_*.png'))
dyn=pd.read_csv(base/'rp_centerline_contact_history.csv')
cp=pd.read_csv(base/'cpress_peaks.csv')
gap=pd.read_csv(job+'_exact_wall_penetration.csv')
cl=pd.read_csv('curvenew_CEL_xyrot56_exact_abaqus.csv')
def tan(s):
    a=cl.arclength_mm.to_numpy(float); x=cl[['x_mm','y_mm','z_mm']].to_numpy(float); i=int(np.clip(np.searchsorted(a,float(s))-1,0,len(a)-2)); v=x[i+1]-x[i]; return v/max(np.linalg.norm(v),1e-15)
vv=pd.read_csv(job+'_rp_velocity_history.csv'); vt_arr=np.einsum('ij,ij->i',vv[['V1','V2','V3']].to_numpy(float),np.asarray([tan(s) for s in dyn.s_mm]))
imgs=[]
for i,f in enumerate(frames):
    im=Image.open(f).convert('RGB'); dr=ImageDraw.Draw(im); j=min(i,len(dyn)-1)
    s=float(dyn.s_mm.iloc[j]); gp=float(gap.min_signed_gap_mm.iloc[min(i,len(gap)-1)])*1000; c=float(cp.cpress_max_mpa.iloc[min(i,len(cp)-1)])
    vtval=float(vt_arr[min(j,len(vt_arr)-1)]); vt=f"{vtval:+.1f}"
    txt=f'χ=134° | HEAD/TAIL axis | s={s:.3f} mm | Vt={vt} mm/s | exact gap={gp:+.3f} μm | CPRESS General Contact={c:.3f} MPa'
    dr.rectangle((10,im.height-35,im.width-10,im.height-8),fill=(255,255,255)); dr.text((16,im.height-31),txt,fill=(20,20,20)); imgs.append(im)
out=Path(job+'_CEL_FSI_bend_validation.gif'); imgs[0].save(out,save_all=True,append_images=imgs[1:],duration=70,loop=0,optimize=False); imgs[-1].save(job+'_CEL_FSI_bend_validation_final.png'); print('GIF',out.resolve(),'frames',len(imgs))
