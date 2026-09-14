"""Zero-solve old/new magnetic regression at identical pose and time."""
from pathlib import Path
import importlib.util
import json
import math
import sys

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SERVER=Path('J:/magpy/magpylib_socket_server.py')
spec=importlib.util.spec_from_file_location('magserver',SERVER)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)


def make(moment):
    model=mod.MagneticCouplingModel(
        'J:/magpy/curvenew_CEL_xyrot56_exact.dxf',drive_type='analytic',
        robot_diameter_mm=.815,robot_height_mm=1.8,robot_br_t=1.46,
        robot_moment_Am2=moment,robot_mass_mg=6.7536638391088974,
        analytic_b_t=.010,analytic_follow_robot=True,analytic_gradient_b_t=.006,
        analytic_gradient_length_mm=45.,analytic_gradient_profile='legacy',
        analytic_gradient_switch_start_s=0.,analytic_gradient_transition_s=.0005,
        driver_speed_mm_s=6.,bend_speed_mm_s=4.5,bend_start_mm=13.49,bend_end_mm=18.56,
        z_offset_mm=90.,driver_start_offset_mm=18.899960626,spin_hz=30.,
        cone_half_angle_deg=30.,driver_xy_shift_mm=[2.,-6.],cone_axis_bias_deg=40.,
        ramp_time_s=.001,endpoint_taper_mm=1.,adaptive_lead_mm=0.,analytic_rotation_sense=1.)
    model.robot_polarity=-1.;model.cone_axis_tangent=True
    axis=model.curve_mm[1]-model.curve_mm[0];model.robot_axis_global=axis/np.linalg.norm(axis)
    model.phase_offset_rad=math.radians(248.)
    return mod.FrameMappedMagneticModel(model,mod.RigidFrameTransform(str(ROOT/'abaqus_magpylib_frame_transform.json')))


def main():
    ident=json.loads((HERE/'freecad_preflight_identity.json').read_text())
    old_m=ident['old_moment_Am2'];new_m=ident['new_moment_Am2']
    pose=np.asarray(ident['RP_mm']);ur=np.zeros(3);t=.004321
    old=make(old_m).evaluate(t,pose,ur);new=make(new_m).evaluate(t,pose,ur)
    b0=np.asarray(old['B_aba_vec_T']);b1=np.asarray(new['B_aba_vec_T'])
    f0=np.asarray(old['force_N']);f1=np.asarray(new['force_N'])
    q0=np.asarray(old['torque_Nmm']);q1=np.asarray(new['torque_Nmm'])
    ratio=new_m/old_m
    row={'time_s':t,'B_difference_T':np.linalg.norm(b1-b0),'B_old_T':np.linalg.norm(b0),'B_new_T':np.linalg.norm(b1),
         'moment_ratio':ratio,'force_norm_ratio':np.linalg.norm(f1)/np.linalg.norm(f0),
         'torque_norm_ratio':np.linalg.norm(q1)/np.linalg.norm(q0)}
    for name,a in [('B_old',b0),('B_new',b1),('F_old',f0),('F_new',f1),('T_old',q0),('T_new',q1)]:
        for c,v in zip('xyz',a):row[name+'_'+c]=v
    row['pass']=bool(row['B_difference_T']<1e-15 and abs(row['force_norm_ratio']-ratio)<1e-10 and abs(row['torque_norm_ratio']-ratio)<1e-10)
    pd.DataFrame([row]).to_csv(HERE/'freecad_same_pose_magnetic_regression.csv',index=False)
    if not row['pass']:raise RuntimeError(row)
    print(json.dumps(row,indent=2))


if __name__=='__main__':main()
