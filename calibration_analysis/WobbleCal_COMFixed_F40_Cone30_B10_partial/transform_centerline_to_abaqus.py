"""Map the Magpylib flat-XY centerline into Abaqus global coordinates."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parent
src = root / 'curvenew_CEL_xyrot56_exact.csv'
dst = root / 'curvenew_CEL_xyrot56_exact_abaqus.csv'
tf = json.loads((root / 'abaqus_magpylib_frame_transform.json').read_text())
R = np.asarray(tf['R_aba_to_mag'], float)
o = np.asarray(tf['origin_aba_mm'], float)
df = pd.read_csv(src)
p = df[['x_mm','y_mm','z_mm']].to_numpy(float)
pab = (R.T.dot(p.T)).T + o
out = df.copy()
out[['x_mm','y_mm','z_mm']] = pab
out.to_csv(dst, index=False)
print('WROTE', dst)
print('first_abaqus_mm', pab[0].tolist())
print('last_abaqus_mm', pab[-1].tolist())
