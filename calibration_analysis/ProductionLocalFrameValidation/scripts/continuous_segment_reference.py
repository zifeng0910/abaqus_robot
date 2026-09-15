"""Independent continuous-segment oracle for the robot-local magnetic frame."""
import json
import math

import ezdxf
import numpy as np


def unit(vector):
    vector = np.asarray(vector, dtype=float)
    return vector / max(float(np.linalg.norm(vector)), 1.0e-30)


def _point3(point):
    return np.array((point.x, point.y, point.z), dtype=float)


def read_curve(dxf_path, spacing_mm=0.1):
    pieces = []
    for entity in ezdxf.readfile(str(dxf_path)).modelspace():
        kind = entity.dxftype()
        if kind not in ("LINE", "ARC", "CIRCLE", "LWPOLYLINE", "POLYLINE", "SPLINE"):
            continue
        try:
            if kind == "POLYLINE":
                try:
                    children = list(entity.vertices())
                except TypeError:
                    children = list(entity.vertices)
                vertices = np.asarray([_point3(child.dxf.location) for child in children])
            else:
                vertices = np.asarray([_point3(point) for point in entity.flattening(0.01)])
        except (AttributeError, TypeError):
            if kind != "LINE":
                continue
            vertices = np.asarray([_point3(entity.dxf.start), _point3(entity.dxf.end)])
        if len(vertices) > 1 and np.linalg.norm(vertices[-1] - vertices[0]) > 1.0e-12:
            pieces.append(vertices)
    if not pieces:
        raise RuntimeError("Independent DXF reader found no curve")
    joined = pieces.pop(0)
    while pieces:
        choices = []
        for index, item in enumerate(pieces):
            choices.extend(((np.linalg.norm(joined[-1] - item[0]), index, 1, 0),
                            (np.linalg.norm(joined[-1] - item[-1]), index, 1, 1),
                            (np.linalg.norm(joined[0] - item[-1]), index, 0, 0),
                            (np.linalg.norm(joined[0] - item[0]), index, 0, 1)))
        gap, index, append, reverse = min(choices, key=lambda row: row[0])
        if gap > 0.02:
            break
        item = pieces.pop(index)
        if reverse:
            item = item[::-1]
        joined = (np.vstack((joined, item[1:])) if append
                  else np.vstack((item[:-1], joined)))
    raw_length = np.linalg.norm(np.diff(joined, axis=0), axis=1)
    raw_s = np.r_[0.0, np.cumsum(raw_length)]
    stations = np.r_[np.arange(0.0, raw_s[-1], spacing_mm), raw_s[-1]]
    return np.column_stack([np.interp(stations, raw_s, joined[:, axis]) for axis in range(3)])


def rodrigues(rotation_vector):
    vector = np.asarray(rotation_vector, dtype=float)
    angle = float(np.linalg.norm(vector))
    if angle < 1.0e-12:
        return np.eye(3)
    axis = vector / angle
    skew = np.array([[0.0, -axis[2], axis[1]],
                     [axis[2], 0.0, -axis[0]],
                     [-axis[1], axis[0], 0.0]])
    return np.eye(3) + math.sin(angle) * skew + (1.0 - math.cos(angle)) * skew.dot(skew)


class ContinuousSegmentReference:
    def __init__(self, dxf_path, transform_path, initial_axis_mag, polarity=-1.0,
                 moment_Am2=0.0010876227522174417, magnitude_T=0.010,
                 cone_deg=30.0, frequency_Hz=30.0, phase_deg=248.0, sense=1.0,
                 ramp_time_s=0.001):
        self.curve = read_curve(dxf_path)
        self.edge = np.diff(self.curve, axis=0)
        self.length = np.linalg.norm(self.edge, axis=1)
        self.raw_tangent = self.edge / self.length[:, None]
        self.arc = np.r_[0.0, np.cumsum(self.length)]
        self.initial_axis_mag = unit(initial_axis_mag)
        self.polarity = float(polarity)
        self.tangent_sign = (1.0 if np.dot(self.raw_tangent[0], self.polarity * self.initial_axis_mag) > 0.0 else -1.0)
        self.transport_e1 = self._transport(self.raw_tangent)
        transform = json.loads(open(transform_path, "r").read())
        self.origin_aba = np.asarray(transform["origin_aba_mm"], dtype=float)
        self.R = np.asarray(transform["R_aba_to_mag"], dtype=float)
        self.moment = float(moment_Am2)
        self.magnitude = float(magnitude_T)
        self.cone = math.radians(cone_deg)
        self.frequency = float(frequency_Hz)
        self.phase0 = math.radians(phase_deg)
        self.sense = float(sense)
        self.ramp_time = max(float(ramp_time_s), 0.0)
        self.previous_s = None

    @staticmethod
    def _transport(tangents):
        first = tangents[0]
        normal = np.array((0.0, 0.0, 1.0)) - first[2] * first
        if np.linalg.norm(normal) < 1.0e-12:
            normal = np.array((0.0, 1.0, 0.0)) - first[1] * first
        frames = [unit(normal)]
        for previous, current in zip(tangents[:-1], tangents[1:]):
            cross = np.cross(previous, current)
            sine = float(np.linalg.norm(cross))
            cosine = float(np.clip(np.dot(previous, current), -1.0, 1.0))
            if sine > 1.0e-12:
                axis = cross / sine
                angle = math.atan2(sine, cosine)
                candidate = (frames[-1] * math.cos(angle) +
                             np.cross(axis, frames[-1]) * math.sin(angle) +
                             axis * np.dot(axis, frames[-1]) * (1.0 - math.cos(angle)))
            else:
                candidate = frames[-1].copy()
            candidate = unit(candidate - np.dot(candidate, current) * current)
            if np.dot(candidate, frames[-1]) < 0.0:
                candidate *= -1.0
            frames.append(candidate)
        return np.asarray(frames)

    def reset(self):
        self.previous_s = None

    def position_to_mag(self, position_aba):
        return self.R.dot(np.asarray(position_aba, dtype=float) - self.origin_aba)

    def vector_to_aba(self, vector_mag):
        return self.R.T.dot(np.asarray(vector_mag, dtype=float))

    def project_mag(self, point_mag):
        point = np.asarray(point_mag, dtype=float)
        raw = np.einsum("ij,ij->i", point - self.curve[:-1], self.edge) / (self.length * self.length)
        fraction = np.clip(raw, 0.0, 1.0)
        fraction[0] = min(raw[0], 1.0)
        fraction[-1] = max(raw[-1], 0.0)
        centers = self.curve[:-1] + fraction[:, None] * self.edge
        distance2 = np.einsum("ij,ij->i", centers - point, centers - point)
        candidate_s = self.arc[:-1] + fraction * self.length
        candidates = np.flatnonzero(distance2 <= float(distance2.min()) + 1.0e-12)
        if self.previous_s is None or len(candidates) == 1:
            index = int(candidates[0])
        else:
            index = int(candidates[np.argmin(np.abs(candidate_s[candidates] - self.previous_s))])
        self.previous_s = float(candidate_s[index])
        return self.previous_s, centers[index]

    def tangent_at(self, s):
        distance = min(max(float(s), 0.0), float(self.arc[-1]))
        index = int(np.searchsorted(self.arc, distance, side="right") - 1)
        index = max(0, min(index, len(self.edge) - 1))
        first = self.raw_tangent[index]
        if index == len(self.edge) - 1:
            return first
        weight = np.clip((distance - self.arc[index]) / max(self.arc[index + 1] - self.arc[index], 1.0e-12), 0.0, 1.0)
        return unit((1.0 - weight) * first + weight * self.raw_tangent[index + 1])

    def frame_at(self, s):
        distance = min(max(float(s), 0.0), float(self.arc[-1]))
        tangent = self.tangent_sign * self.tangent_at(distance)
        index = int(np.searchsorted(self.arc, distance, side="right") - 1)
        index = max(0, min(index, len(self.transport_e1) - 1))
        if index < len(self.transport_e1) - 1:
            weight = np.clip((distance - self.arc[index]) / max(self.arc[index + 1] - self.arc[index], 1.0e-15), 0.0, 1.0)
            e1 = (1.0 - weight) * self.transport_e1[index] + weight * self.transport_e1[index + 1]
        else:
            e1 = self.transport_e1[index].copy()
        e1 = unit(e1 - np.dot(e1, tangent) * tangent)
        e2 = unit(np.cross(tangent, e1))
        return tangent, e1, e2

    def evaluate(self, time_s, position_aba, ur_aba):
        position_mag = self.position_to_mag(position_aba)
        s, center_mag = self.project_mag(position_mag)
        tangent, e1, e2 = self.frame_at(s)
        phase = self.phase0 + self.sense * 2.0 * math.pi * self.frequency * float(time_s)
        field_mag = self.magnitude * (math.cos(self.cone) * tangent +
            math.sin(self.cone) * (math.cos(phase) * e1 + self.sense * math.sin(phase) * e2))
        moment_mag = rodrigues(self.R.dot(np.asarray(ur_aba, dtype=float))).dot(
            self.polarity * self.moment * self.initial_axis_mag)
        torque_mag = np.cross(moment_mag, field_mag) * 1000.0
        if self.ramp_time > 0.0 and float(time_s) < self.ramp_time:
            fraction = np.clip(float(time_s) / self.ramp_time, 0.0, 1.0)
            torque_mag *= fraction * fraction * (3.0 - 2.0 * fraction)
        return {"s_mm": s, "center_mag_mm": center_mag, "tangent_mag": tangent,
                "e1_mag": e1, "e2_mag": e2, "B_aba_T": self.vector_to_aba(field_mag),
                "T_aba_Nmm": self.vector_to_aba(torque_mag), "F_aba_N": np.zeros(3)}
