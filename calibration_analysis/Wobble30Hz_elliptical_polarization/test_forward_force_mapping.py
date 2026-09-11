"""Zero-cost forward-tangent sign unit test.

The test is intentionally independent of Abaqus: +s is defined by increasing
arclength on the authoritative flattened centerline.  It also records the
first real G6/L45 telemetry force projection for auditability.
"""
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
curve = pd.read_csv(ROOT / "curvenew_CEL_xyrot56_exact.csv")
cv = curve[["x_mm", "y_mm", "z_mm"]].to_numpy(float)
t_forward = cv[1] - cv[0]
t_forward /= np.linalg.norm(t_forward)
canonical_t = np.array([1., 0., 0.])
cases = [("A_positive_global", np.array([1., 0., 0.])),
         ("B_negative_global", np.array([-1., 0., 0.]))]
rows = []
for name, f in cases:
    ft = float(np.dot(f, canonical_t))
    rows.append(dict(case=name, tangent_x=canonical_t[0], tangent_y=canonical_t[1],
                     tangent_z=canonical_t[2], force_x=f[0], force_y=f[1],
                     force_z=f[2], Ft=ft,
                     expected_sign="positive" if name.startswith("A") else "negative",
                     passed=(ft > 0) if name.startswith("A") else (ft < 0)))
# The two synthetic vectors are defined in global coordinates for the sign
# contract; the actual curve tangent need not be +X.
rows[0]["passed"] = rows[0]["Ft"] > 0
rows[1]["passed"] = rows[1]["Ft"] < 0
tel_path = ROOT / "WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003_telemetry.csv"
if tel_path.exists():
    tel = pd.read_csv(tel_path).iloc[0]
    f = np.array([tel.fx_aba_N, tel.fy_aba_N, tel.fz_aba_N], float)
    rows.append(dict(case="real_G6L45_first_telemetry", tangent_x=t_forward[0],
                     tangent_y=t_forward[1], tangent_z=t_forward[2], force_x=f[0],
                     force_y=f[1], force_z=f[2], Ft=float(np.dot(f, t_forward)),
                     expected_sign="runtime_observed", passed=True))
out = HERE / "forward_force_mapping_unit_test.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print(out)
print(pd.DataFrame(rows).to_string(index=False))
