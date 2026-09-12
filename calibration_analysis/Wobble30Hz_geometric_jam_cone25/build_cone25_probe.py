from pathlib import Path
import re
ROOT=Path(r'J:\\abaqusfangzhen')
src=ROOT/'WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003.inp'
job='Wobble_F30_G6L45_WallOn_Free_0083_Cone25_WobbleSurvival'
txt=src.read_text(encoding='utf-8',errors='ignore')
txt=re.sub(r'(?m)^\*Heading\s*$',f'*Heading\n** AUTHORIZED 8.333 ms CONE25 GEOMETRIC-JAM PROBE: {job}',txt,count=1)
txt=re.sub(r'(?m)^\*Dynamic, Explicit\s*\r?\n,\s*0\.003\s*$','*Dynamic, Explicit\n, 0.008333',txt,count=1)
txt=re.sub(r'time interval=1\.0e-4','time interval=5.0e-5',txt)
for a,b in [('duration=0.003 s','duration=0.008333 s'),('duration=0.03 s','duration=0.008333 s'),('duration=3 ms','duration=8.333 ms')]: txt=txt.replace(a,b)
(ROOT/f'{job}.inp').write_text(txt,encoding='utf-8')
print(f'WROTE {job}.inp')
