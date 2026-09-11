from pathlib import Path
from PIL import Image
import sys

job = sys.argv[1] if len(sys.argv) > 1 else 'WobbleCal_COMFixed_F40_Cone30_B10'
stride = max(int(sys.argv[2]) if len(sys.argv) > 2 else 4, 1)
paths = sorted((Path('output') / (job + '_bend_validation')).glob(job + '_bend_frame_*.png'))[::stride]
if not paths:
    raise SystemExit('no frame PNGs found')
imgs = [Image.open(p).convert('RGB') for p in paths]
out = Path(job + '_partial_CEL_FSI_bend_validation.gif')
imgs[0].save(out, save_all=True, append_images=imgs[1:], duration=70, loop=0, optimize=False)
final = Path(job + '_partial_CEL_FSI_bend_validation_final.png')
imgs[-1].save(final)
print('GIF', out.resolve(), 'frames', len(imgs), 'bytes', out.stat().st_size)
print('PNG', final.resolve())
