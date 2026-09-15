"""Offline gate: audit wrapper must create the requested local cone exactly."""
import importlib.util
import math
from pathlib import Path

import numpy as np

from audit_common import DXF, TRANSFORM, RP0, unit


wrapper_path = Path(__file__).with_name("local_tangent_socket_server.py")
spec = importlib.util.spec_from_file_location("local_audit_server", str(wrapper_path))
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
model = module.LocalTangentConeModel(str(DXF), drive_type="analytic", robot_diameter_mm=.815, robot_height_mm=2.4,
    robot_moment_Am2=.0010876227522174417, analytic_b_t=.010, analytic_follow_robot=True,
    analytic_gradient_b_t=0, driver_speed_mm_s=6, bend_speed_mm_s=4.5, bend_start_mm=13.49, bend_end_mm=18.56,
    z_offset_mm=90, driver_start_offset_mm=18.899960626, spin_hz=30, cone_half_angle_deg=30,
    driver_xy_shift_mm=[2,-6], cone_axis_bias_deg=0, analytic_rotation_sense=1, ramp_time_s=.001)
model.robot_polarity=-1; model.cone_axis_tangent=True; model.robot_axis_global=unit(model.curve_mm[1]-model.curve_mm[0]); model.phase_offset_rad=math.radians(248)
transform=module.production.RigidFrameTransform(str(TRANSFORM)); wrapped=module.production.FrameMappedMagneticModel(model,transform)
position_mag=transform.position_to_mag(RP0); robot_arc=model._project_arc_mm(position_mag)
c_aba=unit(transform.vector_to_aba(model._tangent_at_distance(robot_arc)))
angles=[]; phases=[]
for time_s in np.linspace(0,1/30,721):
    result=wrapped.evaluate(time_s,RP0,np.zeros(3)); b=np.asarray(result['B_aba_vec_T']); b/=np.linalg.norm(b)
    angles.append(math.degrees(math.acos(np.clip(np.dot(b,c_aba),-1,1))))
    normal=unit(transform.vector_to_aba(model._z_reference_normal(model._tangent_at_distance(robot_arc))))
    side=unit(np.cross(c_aba,normal)); phases.append(math.atan2(np.dot(b,side),np.dot(b,normal)))
winding=(np.unwrap(phases)[-1]-np.unwrap(phases)[0])/(2*math.pi)
assert max(abs(np.asarray(angles)-30.0)) < 1e-8, max(abs(np.asarray(angles)-30.0))
assert abs(winding-1.0) < 1e-8, winding
assert abs(wrapped.evaluate(1/30,RP0,np.zeros(3))['drive_scale']-1.0) < 1e-12
model.ramp_time_s = .005
assert abs(wrapped.evaluate(.0025,RP0,np.zeros(3))['drive_scale']-.5) < 1e-12
assert abs(wrapped.evaluate(.005,RP0,np.zeros(3))['drive_scale']-1.0) < 1e-12
print("PASS LOCAL_TANGENT_CONE theta=30 deg winding=1.0 robot_arc=%.9g ramp5ms=0.5@2.5ms"%robot_arc)
