"""Add pulse/phase/gap/Vt audit text to the fixed-camera bend frames."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

job=sys.argv[1] if len(sys.argv)>1 else 'Wobble_F30_PhaseDirectionalAntiImpact_WallOn_Free_003'
base=Path('output')/(job+'_bend_validation'); frames=sorted(base.glob(job+'_bend_frame_*.png'))
tele=pd.read_csv(job+'_telemetry.csv'); cp=pd.read_csv(base/'cpress_peaks.csv'); gap=pd.read_csv(job+'_exact_wall_penetration.csv'); vel=pd.read_csv(job+'_rp_velocity_history.csv')
defs=pd.read_csv('phase_pulse_definition.csv').iloc[0]; center=float(defs.phase_center_deg); half=float(defs.half_width_deg)
def wrap(d): return (d+180)%360-180
imgs=[]
for i,f in enumerate(frames):
    im=Image.open(f).convert('RGB'); dr=ImageDraw.Draw(im)
    j=min(i,len(tele)-1); phase=float(np.degrees(tele.iloc[j].instantaneous_phase_rad)%360); w=0.5*(1+np.cos(np.pi*wrap(phase-center)/half)) if abs(wrap(phase-center))<half else 0.0
    cpv=float(cp.iloc[min(i,len(cp)-1)].cpress_max_mpa); gp=float(gap.iloc[min(i,len(gap)-1)].min_signed_gap_mm); vt=float(tele.iloc[j].force_tangent_N*1e6)
    vv=vel.iloc[min(i,len(vel)-1)]; speed=float(np.linalg.norm([vv.V1,vv.V2,vv.V3])); state='ON' if w>0 else 'OFF'
    txt=f'phase={phase:6.2f} deg | tensor pulse {state} w={w:.3f} | |V|={speed:7.2f} mm/s | gap={gp*1000:+.3f} um | CPRESS={cpv:.3f} MPa'
    dr.rectangle((12,im.height-34,im.width-12,im.height-8),fill=(255,255,255)); dr.text((18,im.height-31),txt,fill=(20,20,20))
    imgs.append(im)
out=Path(job+'_CEL_FSI_bend_validation.gif'); imgs[0].save(out,save_all=True,append_images=imgs[1:],duration=70,loop=0,optimize=False); imgs[-1].save(job+'_CEL_FSI_bend_validation_final.png')
print('GIF',out.resolve(),'frames',len(imgs))
