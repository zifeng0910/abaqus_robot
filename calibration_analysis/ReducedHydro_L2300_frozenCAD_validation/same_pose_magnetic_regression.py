"""Verify that L2300 changes magnetic loads only through the frozen moment."""

from pathlib import Path
import importlib.util
import json
import math

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SERVER = Path("J:/magpy/magpylib_socket_server.py")
spec = importlib.util.spec_from_file_location("magserver_l2300", SERVER)
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)


def model(height, moment, mass):
    inner = mod.MagneticCouplingModel(
        "J:/magpy/curvenew_CEL_xyrot56_exact.dxf", drive_type="analytic",
        robot_diameter_mm=.815, robot_height_mm=height, robot_br_t=1.46,
        robot_moment_Am2=moment, robot_mass_mg=mass,
        analytic_b_t=.010, analytic_follow_robot=True, analytic_gradient_b_t=.006,
        analytic_gradient_length_mm=45., analytic_gradient_profile="legacy",
        analytic_gradient_switch_start_s=0., analytic_gradient_transition_s=.0005,
        driver_speed_mm_s=6., bend_speed_mm_s=4.5, bend_start_mm=13.49, bend_end_mm=18.56,
        z_offset_mm=90., driver_start_offset_mm=18.899960626, spin_hz=30.,
        cone_half_angle_deg=30., driver_xy_shift_mm=[2., -6.], cone_axis_bias_deg=40.,
        ramp_time_s=.001, endpoint_taper_mm=1., adaptive_lead_mm=0., analytic_rotation_sense=1.)
    inner.robot_polarity = -1.; inner.cone_axis_tangent = True
    axis = inner.curve_mm[1] - inner.curve_mm[0]
    inner.robot_axis_global = axis / np.linalg.norm(axis)
    inner.phase_offset_rad = math.radians(248.)
    return mod.FrameMappedMagneticModel(inner, mod.RigidFrameTransform(str(ROOT / "abaqus_magpylib_frame_transform.json")))


def main():
    old = json.loads((HERE.parents[0] / "ReducedHydro_FreeCAD_L1800_validation/freecad_preflight_identity.json").read_text())
    new = json.loads((HERE / "L2300_deck_identity.json").read_text())
    m0 = old["new_moment_Am2"]; m1 = new["magnetic_moment_Am2"]
    pose = np.asarray(new["RP_mm"]); ur = np.zeros(3); time = .004321
    a = model(1.8, m0, old["mesh_mass_mg"]).evaluate(time, pose, ur)
    b = model(2.3, m1, new["mesh_mass_mg"]).evaluate(time, pose, ur)
    B0 = np.asarray(a["B_aba_vec_T"]); B1 = np.asarray(b["B_aba_vec_T"])
    F0 = np.asarray(a["force_N"]); F1 = np.asarray(b["force_N"])
    T0 = np.asarray(a["torque_Nmm"]); T1 = np.asarray(b["torque_Nmm"])
    ratio = m1 / m0
    row = {
        "time_s": time, "B_difference_T": np.linalg.norm(B1 - B0),
        "B_L1800_T": np.linalg.norm(B0), "B_L2300_T": np.linalg.norm(B1),
        "moment_ratio": ratio, "force_norm_ratio": np.linalg.norm(F1) / np.linalg.norm(F0),
        "torque_norm_ratio": np.linalg.norm(T1) / np.linalg.norm(T0),
    }
    for name, values in (("B_L1800", B0), ("B_L2300", B1), ("F_L1800", F0),
                         ("F_L2300", F1), ("T_L1800", T0), ("T_L2300", T1)):
        for component, value in zip("xyz", values): row[f"{name}_{component}"] = value
    row["pass"] = bool(row["B_difference_T"] < 1e-15
                       and abs(row["force_norm_ratio"] - ratio) < 1e-10
                       and abs(row["torque_norm_ratio"] - ratio) < 1e-10)
    pd.DataFrame([row]).to_csv(HERE / "L2300_same_pose_magnetic_regression.csv", index=False)
    if not row["pass"]: raise RuntimeError(row)
    print(json.dumps(row, indent=2))


if __name__ == "__main__":
    main()
