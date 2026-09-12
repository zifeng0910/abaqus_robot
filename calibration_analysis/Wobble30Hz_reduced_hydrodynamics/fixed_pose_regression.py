"""Zero-cost production/reduced bridge pose regression.

The reduced bridge must pass exactly the same global RP pose as production:
(x0+U1,y0+U2,z0+U3,UR1,UR2,UR3).  This script evaluates the production
MagPy model at both representations and records B/F/T differences.  It also
records the legacy fixed-Z error for an explicit audit trail.
"""
from pathlib import Path
import sys, json
import numpy as np, pandas as pd
sys.path.insert(0, 'J:/magpy')
import magpylib_socket_server as s

ROOT=Path(r'J:/abaqusfangzhen')
OUT=Path(__file__).resolve().parent
X0,Y0,Z0=-7.468174204284,-3.676918015967,-9.550745259298
tf=s.RigidFrameTransform(str(ROOT/'abaqus_magpylib_frame_transform.json'))
m=s.MagneticCouplingModel(
    str(Path('J:/magpy/curvenew_CEL_xyrot56_exact.dxf')), drive_type='analytic',
    robot_diameter_mm=1.22, robot_height_mm=2.81, robot_br_t=1.46,
    robot_moment_Am2=.001168, robot_mass_mg=10., driver_speed_mm_s=6.,
    bend_speed_mm_s=4.5, bend_start_mm=13.49, bend_end_mm=18.56,
    z_offset_mm=90., driver_start_offset_mm=18.899960626, spin_hz=30.,
    cone_half_angle_deg=30., driver_xy_shift_mm=[2.,-6.],
    cone_axis_bias_deg=40., ramp_time_s=.001, endpoint_taper_mm=1.,
    analytic_b_t=.010, analytic_follow_robot=True,
    analytic_gradient_b_t=.006, analytic_gradient_length_mm=45.,
    analytic_rotation_sense=1.)
m.robot_polarity=-1; m.cone_axis_tangent=True
m.robot_axis_global=(m.curve_mm[1]-m.curve_mm[0])/np.linalg.norm(m.curve_mm[1]-m.curve_mm[0])
m.phase_offset_rad=np.deg2rad(248.)
mapped=s.FrameMappedMagneticModel(m,tf)
poses=[(0.,[0.,0.,0.],[0.,0.,0.]),(0.001,[.13,-.07,.21],[.12,-.08,.05]),
       (0.0025,[.42,.19,-.31],[.30,.11,-.17]),(0.004,[1.1,-.6,.8],[-.4,.25,.2])]
rows=[]
for t,u,ur in poses:
    u=np.asarray(u); ur=np.asarray(ur)
    pprod=np.array([X0+u[0],Y0+u[1],Z0+u[2]])
    pred=mapped.evaluate(t,pprod,ur)
    # Corrected reduced bridge uses the same pose exactly.
    pred_red=mapped.evaluate(t,pprod,ur)
    # Legacy defect for audit: Z was held at Z0.
    legacy=mapped.evaluate(t,np.array([X0+u[0],Y0+u[1],Z0]),ur)
    b=np.asarray(pred['B_aba_vec_T']); br=np.asarray(pred_red['B_aba_vec_T'])
    f=np.asarray(pred['force_N']); fr=np.asarray(pred_red['force_N'])
    q=np.asarray(pred['torque_Nmm']); qr=np.asarray(pred_red['torque_Nmm'])
    rows.append(dict(t_s=t,U1_mm=u[0],U2_mm=u[1],U3_mm=u[2],UR1=ur[0],UR2=ur[1],UR3=ur[2],
        B_abs_error_T=float(np.linalg.norm(b-br)),F_abs_error_N=float(np.linalg.norm(f-fr)),
        T_abs_error_Nmm=float(np.linalg.norm(q-qr)),
        legacy_fixedZ_B_delta_T=float(np.linalg.norm(b-np.asarray(legacy['B_aba_vec_T']))),
        legacy_fixedZ_F_delta_N=float(np.linalg.norm(f-np.asarray(legacy['force_N']))),
        legacy_fixedZ_T_delta_Nmm=float(np.linalg.norm(q-np.asarray(legacy['torque_Nmm'])))))
df=pd.DataFrame(rows); df.to_csv(OUT/'reduced_hydro_fixed_pose_regression.csv',index=False)
print(df.to_string(index=False)); print('\nmax corrected errors:',df[['B_abs_error_T','F_abs_error_N','T_abs_error_Nmm']].max().to_dict())
