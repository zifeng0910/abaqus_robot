"""Verify G6/L45 produces positive force along the current canonical +s axis."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
SERVER = OUT / "production" / "magpylib_socket_server_fast.py"
SOURCE = REPO / "calibration_analysis" / "LocalHeadTailRocking" / "case" / "PROD_LOCAL_ROCK10_5HZ_G0_STRAIGHT"
TRANSFORM = REPO.parent / "abaqus_magpylib_frame_transform.json"


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("fast_screen_server", str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    module = load_module(SERVER)
    identity = json.loads((SOURCE / "case_identity.json").read_text())
    c = np.asarray(identity["initial_axis_aba"], dtype=float)
    c /= np.linalg.norm(c)
    model = module.MagneticCouplingModel(
        str(SOURCE / "straight_control_centerline.dxf"), drive_type="analytic",
        robot_diameter_mm=0.815, robot_height_mm=2.4, robot_br_t=1.46,
        robot_moment_Am2=0.0010876227522174417, robot_mass_mg=9.207793514166587,
        analytic_b_t=0.010, analytic_follow_robot=True,
        field_frame_mode="ROBOT_LOCAL_ROCKING", rocking_amplitude_deg=12.0,
        rocking_frame_azimuth_deg=-61.37284757596327,
        analytic_gradient_b_t=0.006, analytic_gradient_length_mm=45.0,
        analytic_gradient_profile="legacy", driver_speed_mm_s=6.0,
        bend_speed_mm_s=4.5, bend_start_mm=13.49, bend_end_mm=18.56,
        z_offset_mm=90.0, driver_start_offset_mm=18.899960626, spin_hz=20.0,
        cone_half_angle_deg=30.0, driver_xy_shift_mm=[2.0, -6.0],
        cone_axis_bias_deg=0.0, analytic_rotation_sense=1.0,
        ramp_time_s=0.001, endpoint_taper_mm=1.0, adaptive_lead_mm=0.0)
    model.set_driver_normal_offset(0.0)
    model.robot_polarity = -1.0
    model.cone_axis_tangent = True
    transform = module.RigidFrameTransform(str(TRANSFORM))
    moment_axis_mag = transform.vector_to_mag(c)
    moment_axis_mag /= np.linalg.norm(moment_axis_mag)
    model.robot_axis_global = moment_axis_mag / model.robot_polarity
    model.phase_offset_rad = 0.0
    wrapped = module.FrameMappedMagneticModel(model, transform)
    rp0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    result = wrapped.evaluate(0.001, rp0, np.zeros(3))
    force = np.asarray(result["force_N"], dtype=float)
    ft = float(np.dot(force, c))
    model.field_frame_mode = "ROBOT_LOCAL_ELLIPTIC_ROCKING"
    model.rocking_amplitude_deg = 14.343111711438091
    model.rocking_cross_amplitude_deg = 2.5
    model.reset_robot_arc_continuity()
    n = np.asarray(identity["n_rock_aba"], dtype=float)
    b = np.asarray(identity["b_rock_aba"], dtype=float)
    field_errors = []
    magnitude_errors = []
    for sample_time in (0.0, 0.0125, 0.025, 0.0375, 0.05):
        sample = wrapped.evaluate(sample_time, rp0, np.zeros(3))
        phase = 2.0 * np.pi * 20.0 * sample_time
        alpha_main = np.radians(14.343111711438091) * np.sin(phase)
        alpha_cross = np.radians(2.5) * np.cos(phase)
        expected = c - np.tan(alpha_main) * n - np.tan(alpha_cross) * b
        expected = 0.010 * expected / np.linalg.norm(expected)
        actual = np.asarray(sample["B_aba_vec_T"], dtype=float)
        field_errors.append(float(np.max(np.abs(actual - expected))))
        magnitude_errors.append(abs(float(np.linalg.norm(actual)) - 0.010))
    elliptic_gate = {
        "formula": "normalize(c-tan(alpha_main)*n_routeA-tan(alpha_cross)*b_routeA)",
        "max_component_error_T": max(field_errors),
        "max_magnitude_error_T": max(magnitude_errors),
        "passed": max(field_errors) < 1e-12 and max(magnitude_errors) < 1e-12,
    }
    audit = {
        "definition": "canonical increasing centerline arclength +s",
        "gradient_mT": 6.0,
        "gradient_length_mm": 45.0,
        "time_s": 0.001,
        "c_hat_aba": c.tolist(),
        "force_aba_N": force.tolist(),
        "force_dot_canonical_plus_s_N": ft,
        "driver_arc_mm": float(result["driver_arc_mm"]),
        "robot_arc_mm": float(result["robot_arc_mm"]),
        "passed": ft > 0.0,
        "elliptic_field_gate": elliptic_gate,
    }
    (OUT / "forward_definition_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")
    print(json.dumps(audit, indent=2))
    if not audit["passed"] or not elliptic_gate["passed"]:
        raise SystemExit("FORWARD_DEFINITION_GATE_FAILED")


if __name__ == "__main__":
    main()
