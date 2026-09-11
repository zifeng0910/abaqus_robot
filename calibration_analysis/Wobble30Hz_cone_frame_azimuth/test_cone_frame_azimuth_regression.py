"""Regression check for the production analytic cone frame at chi=0."""
from pathlib import Path
import importlib.util
import numpy as np

ROOT = Path(__file__).resolve().parent
curve = r'J:\\magpy\\curvenew_CEL_xyrot56_exact.dxf'
params = dict(dxf_path=curve, drive_type='analytic', robot_diameter_mm=1.22,
              robot_height_mm=2.81, robot_moment_Am2=0.001168,
              robot_mass_mg=10.0, driver_speed_mm_s=6.0, bend_speed_mm_s=4.5,
              bend_start_mm=13.49, bend_end_mm=18.56, z_offset_mm=90.0,
              driver_start_offset_mm=18.899960626, spin_hz=30.0,
              driver_xy_shift_mm=[2.0,-6.0], cone_half_angle_deg=30.0,
              cone_axis_bias_deg=40.0, analytic_b_t=0.010,
              analytic_follow_robot=True, analytic_gradient_b_t=0.003,
              analytic_gradient_length_mm=45.0,
              cone_frame_azimuth_deg=0.0)
spec = importlib.util.spec_from_file_location('new_server', r'J:\\magpy\\magpylib_socket_server.py')
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
old_path = ROOT/'abaqus_robot'/'calibration_analysis'/'Wobble30Hz_phase_directional_impulse'/'magpylib_socket_server_directional_pulse.py'
spec2 = importlib.util.spec_from_file_location('old_server', str(old_path))
old = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(old)
newm = mod.MagneticCouplingModel(**params)
oldm = old.MagneticCouplingModel(**{k:v for k,v in params.items() if k != 'cone_frame_azimuth_deg'})
newm.robot_polarity = oldm.robot_polarity = -1.0
newm.cone_axis_tangent = oldm.cone_axis_tangent = True
rows=[]
for i in range(8):
    t=0.0001+0.00035*i; pos=np.array([-7.47+0.02*i,-3.68+0.01*i,-9.55+0.015*i]); ur=np.array([0.01*i,-0.02*i,0.015*i])
    a=oldm.evaluate(t,pos,ur); b=newm.evaluate(t,pos,ur)
    for key in ('B_mag_vec_T','force_N','torque_Nmm'):
        da=np.asarray(a[key],float); db=np.asarray(b[key],float); rows.append((key,float(np.linalg.norm(db-da))))
maxerr=max(x[1] for x in rows)
print('chi=0 regression max absolute error:', maxerr)
print('PASS' if maxerr < 1e-12 else 'FAIL')
if maxerr >= 1e-12: raise SystemExit(1)
