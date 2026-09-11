"""Run the archived production evaluator in the configured MagPy environment.
This is a read-only B/F/T regression; no socket or Abaqus process is started."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0,'J:/magpy')
import magpylib_socket_server as s
ROOT=Path('J:/abaqusfangzhen'); OUT=Path(__file__).resolve().parent
tf=s.RigidFrameTransform(str(ROOT/'abaqus_magpylib_frame_transform.json'))
m=s.MagneticCouplingModel(str(Path('J:/magpy/curvenew_CEL_xyrot56_exact.dxf')),drive_type='analytic',robot_diameter_mm=1.22,robot_height_mm=2.81,robot_br_t=1.46,robot_moment_Am2=.001168,robot_mass_mg=10.,driver_speed_mm_s=6.,bend_speed_mm_s=4.5,bend_start_mm=13.49,bend_end_mm=18.56,z_offset_mm=90.,driver_start_offset_mm=18.899960626,spin_hz=30.,cone_half_angle_deg=30.,driver_xy_shift_mm=[2.,-6.],cone_axis_bias_deg=40.,ramp_time_s=.001,endpoint_taper_mm=1.,analytic_b_t=.010,analytic_follow_robot=True,analytic_gradient_b_t=.006,analytic_gradient_length_mm=45.,analytic_rotation_sense=1.)
m.robot_polarity=-1.
m.cone_axis_tangent=True; m.robot_axis_global=(m.curve_mm[1]-m.curve_mm[0])/np.linalg.norm(m.curve_mm[1]-m.curve_mm[0]); m.phase_offset_rad=np.deg2rad(248.)
mapped=s.FrameMappedMagneticModel(m,tf)
tel=pd.read_csv(ROOT/'WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003_telemetry.csv'); rows=[]
for r in tel.itertuples():
    pos=np.array([r.rp_x_aba_mm,r.rp_y_aba_mm,r.rp_z_aba_mm]); ur=np.array([r.ur1,r.ur2,r.ur3]); z=mapped.evaluate(float(r.t_s),pos,ur); b=np.array(z['B_aba_vec_T']); f=np.array(z['force_N']); tau=np.array(z['torque_Nmm']); bt=np.array([r.Bx_aba_T,r.By_aba_T,r.Bz_aba_T]); ft=np.array([r.fx_aba_N,r.fy_aba_N,r.fz_aba_N]); tt=np.array([r.tx_aba_Nmm,r.ty_aba_Nmm,r.tz_aba_Nmm]); rows.append(dict(time_s=r.t_s,B_model_x_T=b[0],B_model_y_T=b[1],B_model_z_T=b[2],B_tele_x_T=bt[0],B_tele_y_T=bt[1],B_tele_z_T=bt[2],B_abs_error_T=np.linalg.norm(b-bt),F_abs_error_N=np.linalg.norm(f-ft),T_abs_error_Nmm=np.linalg.norm(tau-tt),B_model_norm_T=np.linalg.norm(b),B_tele_norm_T=np.linalg.norm(bt),T_model_norm_Nmm=np.linalg.norm(tau),T_tele_norm_Nmm=np.linalg.norm(tt)))
df=pd.DataFrame(rows);df.to_csv(OUT/'runtime_exact_evaluator_replay.csv',index=False);print(df[['B_abs_error_T','F_abs_error_N','T_abs_error_Nmm']].agg(['max','mean','median']).to_string())
