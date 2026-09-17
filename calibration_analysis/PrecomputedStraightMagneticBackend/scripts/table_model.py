"""Shared straight-tube magnetic table and authoritative live-model helpers."""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
CONFIG_PATH = OUT / "straight_tube_case_config.json"
SERVER_PATH = REPO / "calibration_analysis" / "FastStraightDynamicScreen" / "production" / "magpylib_socket_server_fast.py"
DXF_PATH = REPO / "calibration_analysis" / "F60LowGBackgroundFlowScreen" / "case" / "S4_HEADFORWARD_F60_G0P10_FLOW" / "straight_control_centerline.dxf"
TRANSFORM_PATH = REPO.parent / "abaqus_magpylib_frame_transform.json"


def config():
    return json.loads(CONFIG_PATH.read_text(encoding="ascii"))


def unit(vector):
    vector = np.asarray(vector, dtype=float)
    return vector / np.linalg.norm(vector)


def rotation_from_ur(ur):
    ur = np.asarray(ur, dtype=float)
    angle = np.linalg.norm(ur)
    if angle < 1e-14:
        return np.eye(3)
    axis = ur / angle
    cross = np.array([[0.0, -axis[2], axis[1]],
                      [axis[2], 0.0, -axis[0]],
                      [-axis[1], axis[0], 0.0]])
    return np.eye(3) + math.sin(angle) * cross + (1.0 - math.cos(angle)) * cross.dot(cross)


def field_unit(phase_rad, cfg):
    c = unit(cfg["canonical_plus_s_axis_aba"])
    n = unit(cfg["n_routeA_aba"])
    b = unit(cfg["b_routeA_aba"])
    am = math.radians(cfg["rocking_main_amplitude_deg"]) * math.sin(phase_rad)
    ac = math.radians(cfg["rocking_cross_amplitude_deg"]) * math.cos(phase_rad)
    direction = c - math.tan(am) * n - math.tan(ac) * b
    return unit(direction)


def gradient_unit(s_eff_mm, cfg):
    c = unit(cfg["canonical_plus_s_axis_aba"])
    length = cfg["gradient_length_mm"]
    derivative_per_m = -(s_eff_mm / (length * length)) * math.exp(-0.5 * (s_eff_mm / length) ** 2) * 1000.0
    return np.outer(c, c) * derivative_per_m


def ramp(t_s, cfg):
    duration = cfg["ramp_time_s"]
    if duration <= 0.0 or t_s >= duration:
        return 1.0
    q = max(0.0, min(1.0, t_s / duration))
    return q * q * (3.0 - 2.0 * q)


def table_loads(t_s, position_mm, ur_rad, cfg, phase_table, s_table):
    c = unit(cfg["canonical_plus_s_axis_aba"])
    displacement = np.asarray(position_mm) - np.asarray(cfg["initial_center_aba_mm"])
    robot_arc = cfg["robot_initial_arc_mm"] + float(displacement.dot(c))
    driver_arc = cfg["driver_initial_arc_mm"] + cfg["driver_speed_mm_s"] * t_s
    s_eff = robot_arc - driver_arc
    phase_deg = (cfg["phase0_deg"] + cfg["drive_sense"] * 360.0 * cfg["frequency_Hz"] * t_s) % 360.0
    b_unit = np.array([np.interp(phase_deg, phase_table[:, 0], phase_table[:, i]) for i in range(1, 4)])
    grad_flat = np.array([np.interp(s_eff, s_table[:, 0], s_table[:, i]) for i in range(1, 10)])
    grad_unit_value = grad_flat.reshape(3, 3)
    moment0 = cfg["robot_moment_Am2"] * unit(cfg["initial_magnetic_moment_axis_aba"])
    moment = rotation_from_ur(ur_rad).dot(moment0)
    magnetic_field = cfg["B0_T"] * b_unit
    gradient = cfg["gradient_G_T"] * grad_unit_value
    scale = ramp(t_s, cfg)
    force = moment.dot(gradient) * scale
    torque = np.cross(moment, magnetic_field) * 1000.0 * scale
    return magnetic_field, gradient, force, torque, s_eff, phase_deg


def load_server_module():
    spec = importlib.util.spec_from_file_location("production_magnetic_server", str(SERVER_PATH))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_live_model(cfg):
    module = load_server_module()
    model = module.MagneticCouplingModel(
        str(DXF_PATH), drive_type="analytic", robot_diameter_mm=cfg["robot_diameter_mm"],
        robot_height_mm=cfg["robot_length_mm"], robot_br_t=1.46,
        robot_moment_Am2=cfg["robot_moment_Am2"], robot_mass_mg=cfg["robot_mass_mg"],
        analytic_b_t=cfg["B0_T"], analytic_follow_robot=True,
        field_frame_mode=cfg["field_frame_mode"],
        rocking_amplitude_deg=cfg["rocking_main_amplitude_deg"],
        rocking_cross_amplitude_deg=cfg["rocking_cross_amplitude_deg"],
        rocking_frame_azimuth_deg=61.37284757596327,
        analytic_gradient_b_t=cfg["gradient_G_T"],
        analytic_gradient_length_mm=cfg["gradient_length_mm"],
        analytic_gradient_profile="legacy", driver_speed_mm_s=cfg["driver_speed_mm_s"],
        bend_speed_mm_s=4.5, bend_start_mm=13.49, bend_end_mm=18.56,
        z_offset_mm=90.0, driver_start_offset_mm=cfg["driver_initial_arc_mm"],
        spin_hz=cfg["frequency_Hz"], cone_half_angle_deg=30.0,
        driver_xy_shift_mm=[2.0, -6.0], cone_axis_bias_deg=0.0,
        analytic_rotation_sense=cfg["drive_sense"], ramp_time_s=cfg["ramp_time_s"],
        endpoint_taper_mm=1.0, adaptive_lead_mm=0.0)
    model.set_driver_normal_offset(0.0)
    model.robot_polarity = 1.0
    model.cone_axis_tangent = True
    transform = module.RigidFrameTransform(str(TRANSFORM_PATH))
    moment_axis_mag = unit(transform.vector_to_mag(cfg["initial_magnetic_moment_axis_aba"]))
    model.robot_axis_global = moment_axis_mag / model.robot_polarity
    model.phase_offset_rad = math.radians(cfg["phase0_deg"])
    return module.FrameMappedMagneticModel(model, transform)


def read_table(path):
    with Path(path).open("r", encoding="ascii") as handle:
        header = [float(value) for value in handle.readline().split()]
        ns, nphase = int(header[0]), int(header[1])
        phase = np.array([[float(value) for value in handle.readline().split()] for _ in range(nphase)])
        spatial = np.array([[float(value) for value in handle.readline().split()] for _ in range(ns)])
    return phase, spatial
