"""Audit the rigid flip, torque zero, rocking field, and forward gradient gate."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
SERVER = REPO / "calibration_analysis" / "FastStraightDynamicScreen" / "production" / "magpylib_socket_server_fast.py"
TRANSFORM = REPO.parent / "abaqus_magpylib_frame_transform.json"
SOURCE = REPO / "calibration_analysis" / "FastStraightDynamicScreen" / "cases" / "FAST_ELLIPTIC_20HZ_A14P343_X2P5_FORWARD"
CASES = ("S4_HEADFORWARD_G0", "S4_HEADFORWARD_G0P25", "S4_HEADFORWARD_G0P5", "S4_HEADFORWARD_G1P0")


def load_module(path):
    spec = importlib.util.spec_from_file_location("lowg_server", str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def unit(value):
    value = np.asarray(value, dtype=float)
    return value / np.linalg.norm(value)


def make_model(module, identity, gradient_mT, field_mode="ROBOT_LOCAL_ELLIPTIC_ROCKING",
               main_amp=14.343111711438091, cross_amp=2.5):
    c = unit(identity["canonical_plus_s_axis_aba"])
    moment = unit(identity["initial_magnetic_moment_axis_aba"])
    model = module.MagneticCouplingModel(
        str(SOURCE / "straight_control_centerline.dxf"), drive_type="analytic",
        robot_diameter_mm=0.815, robot_height_mm=2.4, robot_br_t=1.46,
        robot_moment_Am2=identity["robot_moment_Am2"], robot_mass_mg=identity["robot_mass_mg"],
        analytic_b_t=0.010, analytic_follow_robot=True,
        field_frame_mode=field_mode, rocking_amplitude_deg=main_amp,
        rocking_cross_amplitude_deg=cross_amp,
        rocking_frame_azimuth_deg=identity["routeA_gauge_for_flipped_body_deg"],
        analytic_gradient_b_t=gradient_mT / 1000.0, analytic_gradient_length_mm=45.0,
        analytic_gradient_profile="legacy", driver_speed_mm_s=6.0,
        bend_speed_mm_s=4.5, bend_start_mm=13.49, bend_end_mm=18.56,
        z_offset_mm=90.0, driver_start_offset_mm=18.899960626, spin_hz=20.0,
        cone_half_angle_deg=30.0, driver_xy_shift_mm=[2.0, -6.0],
        cone_axis_bias_deg=0.0, analytic_rotation_sense=1.0,
        ramp_time_s=0.001, endpoint_taper_mm=1.0, adaptive_lead_mm=0.0)
    model.set_driver_normal_offset(0.0)
    model.robot_polarity = 1.0
    model.cone_axis_tangent = True
    transform = module.RigidFrameTransform(str(TRANSFORM))
    moment_mag = unit(transform.vector_to_mag(moment))
    model.robot_axis_global = moment_mag / model.robot_polarity
    model.phase_offset_rad = 0.0
    return module.FrameMappedMagneticModel(model, transform), transform


def main():
    module = load_module(SERVER)
    identity = json.loads((OUT / "cases" / CASES[0] / "case_identity.json").read_text())
    c = unit(identity["canonical_plus_s_axis_aba"])
    n = unit(identity["n_routeA_aba"])
    b = unit(identity["b_routeA_aba"])
    rp0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    moment = unit(identity["initial_magnetic_moment_axis_aba"])

    zero_model, zero_transform = make_model(module, identity, 0.0, main_amp=0.0, cross_amp=0.0)
    zero = zero_transform and zero_model.evaluate(0.001, rp0, np.zeros(3))
    zero_torque = np.asarray(zero["torque_Nmm"], dtype=float)

    rocking_model, rocking_transform = make_model(module, identity, 0.0)
    rocking = []
    for t in (0.0, 0.0125, 0.025, 0.0375, 0.05):
        result = rocking_model.evaluate(t, rp0, np.zeros(3))
        B = unit(result["B_aba_vec_T"])
        rocking.append({
            "time_s": t,
            "alpha_main_deg": float(14.343111711438091 * np.sin(2.0 * np.pi * 20.0 * t)),
            "alpha_cross_deg": float(2.5 * np.cos(2.0 * np.pi * 20.0 * t)),
            "B_dot_c": float(np.dot(B, c)), "B_dot_n": float(np.dot(B, n)),
            "B_dot_b": float(np.dot(B, b)), "moment_dot_B": float(np.dot(moment, B)),
        })

    gates = []
    for gradient in (0.25, 0.50, 1.00):
        model, transform = make_model(module, identity, gradient)
        values = []
        for offset in (0.0, 1.0, 3.0, 5.0):
            result = transform and model.evaluate(0.001, rp0 + offset * c, np.zeros(3))
            force = np.asarray(result["force_N"], dtype=float)
            values.append({"offset_s_mm": offset,
                           "force_aba_N": force.tolist(),
                           "force_dot_canonical_plus_s_N": float(np.dot(force, c)),
                           "driver_arc_mm": result["driver_arc_mm"],
                           "robot_arc_mm": result["robot_arc_mm"]})
        gates.append({"gradient_mT": gradient, "samples": values,
                      "passed": all(v["force_dot_canonical_plus_s_N"] > 0.0 for v in values)})

    audit = {
        "canonical_plus_s_definition": "increasing centerline arclength; camera left-to-right",
        "mesh_flip": {
            "axis_aba": identity["mesh_flip_axis_aba"],
            "rotation_matrix": identity["mesh_flip_rotation_matrix"],
            "COM_aba_mm": identity["initial_center_aba_mm"],
            "flipped_head_tail_dot_plus_s_mm": identity["head_tail_dot_canonical_plus_s_mm"],
            "passed": identity["head_tail_dot_canonical_plus_s_mm"] > 0.0,
        },
        "ur0_magnetic_moment_axis_aba": moment.tolist(),
        "reduced_hydro_body_axis_aba": identity["reduced_hydro_body_axis_aba"],
        "zero_dynamics_torque": {
            "alpha_main_deg": 0.0, "alpha_cross_deg": 0.0,
            "torque_Nmm": zero_torque.tolist(),
            "norm_Nmm": float(np.linalg.norm(zero_torque)),
            "passed": bool(np.linalg.norm(zero_torque) < 1.0e-10),
        },
        "rocking_sign_check": {
            "mode": "ROBOT_LOCAL_ELLIPTIC_ROCKING; fixed RouteA gauge; no precession",
            "samples": rocking,
            "field_n_sign_matches_v_minus_tan_alpha_n": all(
                (abs(r["alpha_main_deg"]) < 1.0e-10 or np.sign(r["B_dot_n"]) == -np.sign(r["alpha_main_deg"]))
                for r in rocking),
        },
        "forward_gradient_gates": gates,
        "passed": bool(identity["head_tail_dot_canonical_plus_s_mm"] > 0.0 and
                       np.linalg.norm(zero_torque) < 1.0e-10 and all(g["passed"] for g in gates)),
    }
    (OUT / "lowg_setup_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")
    print(json.dumps(audit, indent=2))
    if not audit["passed"]:
        raise SystemExit("LOWG_SETUP_AUDIT_FAILED")


if __name__ == "__main__":
    main()
