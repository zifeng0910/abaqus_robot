"""Offline gates for the production robot-local magnetic field frame."""
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
VALIDATION = HERE.parent
REPO = VALIDATION.parents[1]
ROOT = REPO.parent
VENDORED = VALIDATION / "production" / "magpylib_socket_server.py"
LIVE = Path(r"J:\magpy\magpylib_socket_server.py")
SCREENING = REPO / "calibration_analysis" / "MagneticDriveFrameAudit" / "scripts" / "local_tangent_socket_server.py"
DXF = Path(r"J:\magpy\curvenew_CEL_xyrot56_exact.dxf")
TRANSFORM = ROOT / "abaqus_magpylib_frame_transform.json"
TRAJECTORY = REPO / "calibration_analysis" / "MagneticDriveFrameAudit" / "cases" / "FRAME_LOCAL30" / "private" / "rp_history_private.npz"
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
PHASE0 = math.radians(248.0)
MOMENT = 0.0010876227522174417


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def unit(vector):
    vector = np.asarray(vector, dtype=float)
    return vector / max(float(np.linalg.norm(vector)), 1.0e-30)


def model_kwargs(sense=1.0, bias=0.0, gradient=0.0):
    return dict(
        drive_type="analytic", robot_diameter_mm=0.815, robot_height_mm=2.4,
        robot_br_t=1.46, robot_moment_Am2=MOMENT,
        robot_mass_mg=9.207793514166587, analytic_b_t=0.010,
        analytic_follow_robot=True, analytic_gradient_b_t=gradient,
        analytic_gradient_length_mm=45.0, analytic_gradient_profile="legacy",
        driver_speed_mm_s=6.0, bend_speed_mm_s=4.5,
        bend_start_mm=13.49, bend_end_mm=18.56, z_offset_mm=90.0,
        driver_start_offset_mm=18.899960626, spin_hz=30.0,
        cone_half_angle_deg=30.0, driver_xy_shift_mm=[2.0, -6.0],
        cone_axis_bias_deg=bias, analytic_rotation_sense=sense,
        ramp_time_s=0.001, endpoint_taper_mm=1.0, adaptive_lead_mm=0.0)


def finish_model(module, model):
    model.robot_polarity = -1.0
    model.cone_axis_tangent = True
    model.robot_axis_global = unit(model.curve_mm[1] - model.curve_mm[0])
    model.phase_offset_rad = PHASE0
    transform = module.RigidFrameTransform(str(TRANSFORM))
    return module.FrameMappedMagneticModel(model, transform), transform


def build(module, field_mode=None, sense=1.0, bias=0.0, gradient=0.0):
    kwargs = model_kwargs(sense=sense, bias=bias, gradient=gradient)
    if field_mode is not None:
        kwargs["field_frame_mode"] = field_mode
    model = module.MagneticCouplingModel(str(DXF), **kwargs)
    wrapped, transform = finish_model(module, model)
    return model, wrapped, transform


def vector(result, key):
    return np.asarray(result[key], dtype=float)


def angle_deg(a, b):
    a, b = unit(a), unit(b)
    chord = min(2.0, float(np.linalg.norm(a - b)))
    return math.degrees(2.0 * math.asin(0.5 * chord))


def legacy_gate(live, corrected):
    maxima = {"B_T": 0.0, "F_N": 0.0, "T_Nmm": 0.0}
    positions = [RP0, RP0 + [0.17, -0.08, 0.11]]
    rotations = [np.zeros(3), np.array([0.17, -0.09, 0.21])]
    for sense in (-1.0, 1.0):
        old, old_w, _ = build(live, sense=sense, bias=40.0, gradient=0.006)
        new, new_w, _ = build(corrected, field_mode="LEGACY_DRIVER_FRAME", sense=sense, bias=40.0, gradient=0.006)
        for time_s in np.linspace(0.0, 1.0 / 30.0, 41):
            for position, rotation in zip(positions, rotations):
                a = old_w.evaluate(time_s, position, rotation)
                b = new_w.evaluate(time_s, position, rotation)
                maxima["B_T"] = max(maxima["B_T"], float(np.max(np.abs(vector(a, "B_aba_vec_T") - vector(b, "B_aba_vec_T")))))
                maxima["F_N"] = max(maxima["F_N"], float(np.max(np.abs(vector(a, "force_N") - vector(b, "force_N")))))
                maxima["T_Nmm"] = max(maxima["T_Nmm"], float(np.max(np.abs(vector(a, "torque_Nmm") - vector(b, "torque_Nmm")))))
        assert old.field_frame_mode if hasattr(old, "field_frame_mode") else True
    assert max(maxima.values()) < 1.0e-14, maxima
    return maxima


def projection_gate(model):
    edge = np.diff(model.curve_mm, axis=0)
    length = np.linalg.norm(edge, axis=1)
    tangent = edge / length[:, None]
    index = min(7, len(edge) // 2)
    midpoint = model.curve_mm[index] + 0.37 * edge[index]
    expected = model.arc_mm[index] + 0.37 * length[index]
    model.reset_robot_arc_continuity()
    interior = model._project_robot_arc_mm(midpoint)
    model.reset_robot_arc_continuity()
    before = model._project_robot_arc_mm(model.curve_mm[0] - 0.25 * tangent[0])
    model.reset_robot_arc_continuity()
    after = model._project_robot_arc_mm(model.curve_mm[-1] + 0.25 * tangent[-1])
    c0, _, _ = model._robot_local_frame(interior)
    c1, _, _ = model._robot_local_frame(min(interior + 0.2, model.arc_mm[-1]))
    assert abs(interior - expected) < 1.0e-10
    assert abs(before + 0.25) < 1.0e-10
    assert abs(after - (model.arc_mm[-1] + 0.25)) < 1.0e-10
    assert np.dot(c0, model.robot_polarity * model.robot_axis_global) > 0.0
    assert np.dot(c0, c1) > 0.0
    return {"interior_error_mm": abs(interior - expected), "before_s_mm": before,
            "after_excess_mm": after - model.arc_mm[-1], "tangent_continuity_dot": float(np.dot(c0, c1))}


def fixed_com_gate(corrected):
    model, wrapped, transform = build(corrected, field_mode="ROBOT_LOCAL_TANGENT")
    position_mag = transform.position_to_mag(RP0)
    times = np.linspace(0.0, 1.0 / 30.0, 721)
    fields, phases, theta = [], [], []
    reference_max = 0.0
    for time_s in times:
        result = wrapped.evaluate(time_s, RP0, np.zeros(3))
        b = vector(result, "B_aba_vec_T")
        c, e1, e2 = model._robot_local_frame(result["robot_arc_mm"])
        c = transform.vector_to_aba(c); e1 = transform.vector_to_aba(e1); e2 = transform.vector_to_aba(e2)
        fields.append(b)
        theta.append(angle_deg(b, c))
        phases.append(math.atan2(np.dot(b, e2), np.dot(b, e1)))
        phi = PHASE0 + 2.0 * math.pi * 30.0 * time_s
        reference = 0.010 * (math.cos(math.radians(30.0)) * unit(c) +
                            math.sin(math.radians(30.0)) * (math.cos(phi) * unit(e1) + math.sin(phi) * unit(e2)))
        reference_max = max(reference_max, float(np.max(np.abs(b - reference))))
    fields = np.asarray(fields)
    winding = float((np.unwrap(phases)[-1] - np.unwrap(phases)[0]) / (2.0 * math.pi))
    metrics = {
        "component_reference_max_abs_T": reference_max,
        "magnitude_max_abs_T": float(np.max(np.abs(np.linalg.norm(fields, axis=1) - 0.010))),
        "theta_max_abs_deg": float(np.max(np.abs(np.asarray(theta) - 30.0))),
        "winding_turn": winding,
        "projected_robot_arc_mm": float(model._project_robot_arc_mm(position_mag)),
    }
    assert metrics["magnitude_max_abs_T"] < 1.0e-14, metrics
    assert metrics["theta_max_abs_deg"] < 1.0e-8, metrics
    assert abs(winding - 1.0) < 1.0e-10, metrics
    return metrics, projection_gate(model)


def screening_gate(corrected):
    screening = load_module("screening_local_server", SCREENING)
    production_model, production_wrapped, _ = build(corrected, field_mode="ROBOT_LOCAL_TANGENT")
    compatibility_model, compatibility_wrapped, _ = build(corrected, field_mode="ROBOT_LOCAL_TANGENT")
    compatibility_model._project_robot_arc_mm = compatibility_model._project_arc_mm
    screening_model = screening.LocalTangentConeModel(str(DXF), **model_kwargs())
    screening_wrapped, _ = finish_model(screening.production, screening_model)
    history = np.load(TRAJECTORY)
    count = len(history["U1"])
    indices = np.unique(np.linspace(0, count - 1, 361).astype(int))
    maxima = {"B_component_T": 0.0, "B_angle_deg": 0.0, "B_magnitude_T": 0.0,
              "F_component_N": 0.0, "T_component_Nmm": 0.0}
    compatibility = {"B_component_T": 0.0, "B_angle_deg": 0.0,
                     "F_component_N": 0.0, "T_component_Nmm": 0.0}
    production_model.reset_robot_arc_continuity()
    for index in indices:
        time_s = float(history["U1"][index, 0])
        position = RP0 + np.array([history[key][index, 1] for key in ("U1", "U2", "U3")])
        rotation = np.array([history[key][index, 1] for key in ("UR1", "UR2", "UR3")])
        reference = screening_wrapped.evaluate(time_s, position, rotation)
        result = production_wrapped.evaluate(time_s, position, rotation)
        compat = compatibility_wrapped.evaluate(time_s, position, rotation)
        br, bp = vector(reference, "B_aba_vec_T"), vector(result, "B_aba_vec_T")
        maxima["B_component_T"] = max(maxima["B_component_T"], float(np.max(np.abs(br - bp))))
        maxima["B_angle_deg"] = max(maxima["B_angle_deg"], angle_deg(br, bp))
        maxima["B_magnitude_T"] = max(maxima["B_magnitude_T"], abs(float(np.linalg.norm(br) - np.linalg.norm(bp))))
        maxima["F_component_N"] = max(maxima["F_component_N"], float(np.max(np.abs(vector(reference, "force_N") - vector(result, "force_N")))))
        maxima["T_component_Nmm"] = max(maxima["T_component_Nmm"], float(np.max(np.abs(vector(reference, "torque_Nmm") - vector(result, "torque_Nmm")))))
        bc = vector(compat, "B_aba_vec_T")
        compatibility["B_component_T"] = max(compatibility["B_component_T"], float(np.max(np.abs(br - bc))))
        compatibility["B_angle_deg"] = max(compatibility["B_angle_deg"], angle_deg(br, bc))
        compatibility["F_component_N"] = max(compatibility["F_component_N"], float(np.max(np.abs(vector(reference, "force_N") - vector(compat, "force_N")))))
        compatibility["T_component_Nmm"] = max(compatibility["T_component_Nmm"], float(np.max(np.abs(vector(reference, "torque_Nmm") - vector(compat, "torque_Nmm")))))
    maxima["samples"] = int(len(indices))
    maxima["nearest_vertex_compatibility"] = compatibility
    maxima["pass"] = bool(
        maxima["B_component_T"] < 1.0e-10 and
        maxima["B_angle_deg"] < 1.0e-6 and
        maxima["B_magnitude_T"] < 1.0e-14 and
        maxima["F_component_N"] < 1.0e-14 and
        maxima["T_component_Nmm"] < 1.0e-10)
    assert max(compatibility.values()) < 1.0e-10, compatibility
    return maxima


def rejection_gate(corrected):
    for kwargs in ({"field_mode": "ROBOT_LOCAL_TANGENT", "bias": 1.0},):
        try:
            build(corrected, **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError("ROBOT_LOCAL_TANGENT accepted nonzero legacy bias")


def main():
    live = load_module("live_legacy_server", LIVE)
    corrected = load_module("corrected_production_server", VENDORED)
    print("RUN legacy preservation gate", flush=True)
    results = {"legacy": legacy_gate(live, corrected)}
    print("RUN fixed-COM and projection gates", flush=True)
    results["fixed_com"], results["projection"] = fixed_com_gate(corrected)
    rejection_gate(corrected)
    print("RUN moving FRAME_LOCAL30 B/F/T gate", flush=True)
    results["moving_FRAME_LOCAL30"] = screening_gate(corrected)
    output = VALIDATION / "field_regression_metrics.json"
    output.write_text(json.dumps(results, indent=2) + "\n", encoding="ascii")
    print(json.dumps(results, indent=2))
    if not results["moving_FRAME_LOCAL30"]["pass"]:
        raise SystemExit("PRODUCTION_LOCAL_FRAME_REGRESSION_FAILED")
    print("PASS PRODUCTION_LOCAL_FRAME_REGRESSION")


if __name__ == "__main__":
    main()
