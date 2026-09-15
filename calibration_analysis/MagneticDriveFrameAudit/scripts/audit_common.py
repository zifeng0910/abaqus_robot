"""Shared geometry and rotation utilities for the magnetic-frame audit."""
from pathlib import Path
import importlib.util
import json
import math

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent
REPO = AUDIT.parents[1]
ROOT = REPO.parent
SERVER = Path(r"J:\magpy\magpylib_socket_server.py")
DXF = Path(r"J:\magpy\curvenew_CEL_xyrot56_exact.dxf")
TRANSFORM = ROOT / "abaqus_magpylib_frame_transform.json"
CURVE_CSV = ROOT / "CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_TRUE_centerline_odb.csv"
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298], dtype=float)
A0 = np.array([0.9647382600216, -0.1188742372140, 0.2348382536499], dtype=float)
A0 /= np.linalg.norm(A0)


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / max(float(np.linalg.norm(v)), 1.0e-30)


def rodrigues(rotvec):
    vector = np.asarray(rotvec, dtype=float)
    angle = float(np.linalg.norm(vector))
    if angle < 1.0e-14:
        return np.eye(3)
    axis = vector / angle
    cross = np.array([[0.0, -axis[2], axis[1]],
                      [axis[2], 0.0, -axis[0]],
                      [-axis[1], axis[0], 0.0]])
    return np.eye(3) + math.sin(angle) * cross + (1.0 - math.cos(angle)) * cross.dot(cross)


def rotation_error_deg(a, b):
    delta = np.asarray(a).T.dot(np.asarray(b))
    return math.degrees(math.acos(np.clip((np.trace(delta) - 1.0) / 2.0, -1.0, 1.0)))


def load_production_module():
    spec = importlib.util.spec_from_file_location("production_magnetic_server", str(SERVER))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TransportedTubeFrame:
    """Continuous minimal-rotation frame on the authoritative Abaqus centerline."""

    def __init__(self):
        table = pd.read_csv(CURVE_CSV)
        self.points = table[["x_mm", "y_mm", "z_mm"]].to_numpy(float)
        edge = np.diff(self.points, axis=0)
        self.length = np.linalg.norm(edge, axis=1)
        raw_tangent = edge / self.length[:, None]
        sign = 1.0 if np.dot(raw_tangent[np.argmin(np.linalg.norm(self.points[:-1] - RP0, axis=1))], A0) > 0 else -1.0
        self.tangent = sign * raw_tangent
        self.arc = np.r_[0.0, np.cumsum(self.length)]
        first = self.tangent[0]
        e = np.array([0.0, 0.0, 1.0]) - first[2] * first
        if np.linalg.norm(e) < 1.0e-10:
            e = np.array([0.0, 1.0, 0.0]) - first[1] * first
        frames = [unit(e)]
        for previous, current in zip(self.tangent[:-1], self.tangent[1:]):
            cross = np.cross(previous, current)
            sine = float(np.linalg.norm(cross))
            cosine = float(np.clip(np.dot(previous, current), -1.0, 1.0))
            if sine > 1.0e-12:
                rotation = rodrigues(unit(cross) * math.atan2(sine, cosine))
                candidate = rotation.dot(frames[-1])
            else:
                candidate = frames[-1]
            candidate -= np.dot(candidate, current) * current
            candidate = unit(candidate)
            if np.dot(candidate, frames[-1]) < 0:
                candidate *= -1.0
            frames.append(candidate)
        self.e1 = np.asarray(frames)
        self.e2 = np.cross(self.tangent, self.e1)
        self.e2 /= np.linalg.norm(self.e2, axis=1)[:, None]

    def at(self, position):
        point = np.asarray(position, dtype=float)
        origin = self.points[:-1]
        edge = np.diff(self.points, axis=0)
        frac = np.clip(np.einsum("ij,ij->i", point - origin, edge) / np.einsum("ij,ij->i", edge, edge), 0.0, 1.0)
        closest = origin + frac[:, None] * edge
        index = int(np.argmin(np.linalg.norm(closest - point, axis=1)))
        c = self.tangent[index]
        e1 = self.e1[index] - np.dot(self.e1[index], c) * c
        e1 = unit(e1)
        e2 = unit(np.cross(c, e1))
        s = self.arc[index] + frac[index] * self.length[index]
        return float(s), closest[index], c, e1, e2


def build_production_model(gradient_mT=6.0, length_mm=2.4, moment_Am2=0.0010876227522174417):
    production = load_production_module()
    model = production.MagneticCouplingModel(
        str(DXF), drive_type="analytic", robot_diameter_mm=0.815,
        robot_height_mm=length_mm, robot_br_t=1.46,
        robot_moment_Am2=moment_Am2, robot_mass_mg=9.207793514166587,
        analytic_b_t=0.010, analytic_follow_robot=True,
        analytic_gradient_b_t=gradient_mT / 1000.0,
        analytic_gradient_length_mm=45.0, analytic_gradient_profile="legacy",
        analytic_gradient_switch_start_s=0.0, analytic_gradient_transition_s=0.0005,
        driver_speed_mm_s=6.0, bend_speed_mm_s=4.5,
        bend_start_mm=13.49, bend_end_mm=18.56, z_offset_mm=90.0,
        driver_start_offset_mm=18.899960626, spin_hz=30.0,
        cone_half_angle_deg=30.0, driver_xy_shift_mm=[2.0, -6.0],
        cone_axis_bias_deg=40.0, analytic_rotation_sense=1.0,
        ramp_time_s=0.001, endpoint_taper_mm=1.0, adaptive_lead_mm=0.0)
    model.robot_polarity = -1.0
    model.cone_axis_tangent = True
    model.robot_axis_global = unit(model.curve_mm[1] - model.curve_mm[0])
    model.phase_offset_rad = math.radians(248.0)
    transform = production.RigidFrameTransform(str(TRANSFORM))
    return production, model, production.FrameMappedMagneticModel(model, transform), transform


def local_tangent_field(position, time_s, alpha_deg, tube_frame, phase0_deg=248.0, sense=1.0, magnitude_T=0.010):
    _, _, c, e1, e2 = tube_frame.at(position)
    phase = math.radians(phase0_deg) + sense * 2.0 * math.pi * 30.0 * float(time_s)
    alpha = math.radians(alpha_deg)
    field = magnitude_T * (math.cos(alpha) * c + math.sin(alpha) * (math.cos(phase) * e1 + sense * math.sin(phase) * e2))
    return field, phase, c, e1, e2
