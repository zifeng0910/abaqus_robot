"""SCREENING-ONLY server wrapper for a COM-local, polarity-correct tube-axis cone.

The production server remains unchanged. This wrapper subclasses its magnetic
model only for FRAME_LOCAL* cases and then calls the production CLI entrypoint.
"""
import importlib.util
import math
from pathlib import Path

import numpy as np


SERVER = Path(r"J:\magpy\magpylib_socket_server.py")
spec = importlib.util.spec_from_file_location("production_magnetic_server", str(SERVER))
production = importlib.util.module_from_spec(spec)
spec.loader.exec_module(production)
BaseModel = production.MagneticCouplingModel


def unit(vector):
    return vector / max(float(np.linalg.norm(vector)), 1.0e-15)


class LocalTangentConeModel(BaseModel):
    """Use robot COM arc, c_hat aligned with initial robot, and transported e1/e2."""

    def __init__(self, *args, **kwargs):
        super(LocalTangentConeModel, self).__init__(*args, **kwargs)
        raw = np.diff(self.curve_mm, axis=0)
        raw /= np.linalg.norm(raw, axis=1)[:, None]
        # Production a0/polarity maps to the negative DXF tangent.
        self._local_c = -raw
        e = np.array([0.0, 0.0, 1.0]) - self._local_c[0, 2] * self._local_c[0]
        if np.linalg.norm(e) < 1.0e-12:
            e = np.array([0.0, 1.0, 0.0]) - self._local_c[0, 1] * self._local_c[0]
        frames = [unit(e)]
        for previous, current in zip(self._local_c[:-1], self._local_c[1:]):
            cross = np.cross(previous, current); sine = float(np.linalg.norm(cross)); cosine = float(np.clip(np.dot(previous,current),-1,1))
            if sine > 1.0e-12:
                axis = cross / sine; angle = math.atan2(sine, cosine)
                candidate = frames[-1]*math.cos(angle) + np.cross(axis,frames[-1])*math.sin(angle) + axis*np.dot(axis,frames[-1])*(1-math.cos(angle))
            else:
                candidate = frames[-1]
            candidate = unit(candidate - np.dot(candidate, current) * current)
            if np.dot(candidate, frames[-1]) < 0: candidate *= -1
            frames.append(candidate)
        self._local_e1 = np.asarray(frames)

    def _tangent_at_distance(self, distance_mm):
        return -super(LocalTangentConeModel, self)._tangent_at_distance(distance_mm)

    def _z_reference_normal(self, tangent):
        tangent = unit(np.asarray(tangent, dtype=float))
        index = int(np.argmax(self._local_c.dot(tangent)))
        e1 = self._local_e1[index] - np.dot(self._local_e1[index], tangent) * tangent
        return unit(e1)

    def _driver_pose(self, t_s, robot_position_mm=None):
        position, direction, distance, robot_arc = super(LocalTangentConeModel, self)._driver_pose(t_s, robot_position_mm)
        if self.analytic_follow_robot and robot_position_mm is not None:
            robot_arc = self._project_arc_mm(robot_position_mm)
        return position, direction, distance, robot_arc

    def _drive_scale(self, t_s, distance_mm):
        # A COM-local analytic field has no physical driver endpoint. Retain
        # the established smooth startup, then hold the requested B0 exactly.
        if self.ramp_time_s <= 0.0:
            return 1.0
        x = min(1.0, max(0.0, float(t_s) / self.ramp_time_s))
        return x*x*(3.0-2.0*x)


production.MagneticCouplingModel = LocalTangentConeModel

if __name__ == "__main__":
    production.main()
