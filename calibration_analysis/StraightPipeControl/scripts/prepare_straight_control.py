"""Build the one authorized straight-pipe control from the solved production deck."""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

import ezdxf
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
BASE = REPO / "calibration_analysis" / "ProductionLocalFrameValidation"
sys.path.insert(0, str(BASE / "scripts"))
SOURCE_JOB = "PROD_LOCAL30_G0_3CYCLE"
JOB = "PROD_LOCAL30_G0_STRAIGHT_CTRL"
SOURCE_CASE = BASE / "case" / SOURCE_JOB
CASE = OUT / "case" / JOB
SOURCE_DXF = Path(r"J:\magpy\curvenew_CEL_xyrot56_exact.dxf")
TRANSFORM = Path(r"J:\abaqusfangzhen\abaqus_magpylib_frame_transform.json")
RADIUS_MM = 0.66734524
HALF_LENGTH_MM = 12.0
GAUGE_PREFIX_OFFSET_MM = 40.0
AXIAL_DIVISIONS = 120
CIRCUMFERENTIAL_DIVISIONS = 24


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unit(vector: np.ndarray) -> np.ndarray:
    return vector / np.linalg.norm(vector)


def replace_part(deck: str, part_name: str, replacement: str) -> str:
    pattern = re.compile(
        rf"(?ms)^\*Part, name={re.escape(part_name)}\s*$.*?^\*End Part\s*$"
    )
    result, count = pattern.subn(replacement.rstrip(), deck, count=1)
    if count != 1:
        raise RuntimeError(f"Expected exactly one {part_name} part; found {count}")
    return result


def cylinder_part(center: np.ndarray, tangent: np.ndarray, e1: np.ndarray, e2: np.ndarray):
    node_radius = RADIUS_MM / np.cos(np.pi / CIRCUMFERENTIAL_DIVISIONS)
    nodes = []
    for i, axial in enumerate(np.linspace(-HALF_LENGTH_MM, HALF_LENGTH_MM, AXIAL_DIVISIONS + 1)):
        for j, theta in enumerate(np.linspace(0.0, 2.0 * np.pi, CIRCUMFERENTIAL_DIVISIONS, endpoint=False)):
            point = center + axial * tangent + node_radius * (np.cos(theta) * e1 + np.sin(theta) * e2)
            nodes.append((i * CIRCUMFERENTIAL_DIVISIONS + j + 1, *point))
    elements = []
    for i in range(AXIAL_DIVISIONS):
        for j in range(CIRCUMFERENTIAL_DIVISIONS):
            jp = (j + 1) % CIRCUMFERENTIAL_DIVISIONS
            n1 = i * CIRCUMFERENTIAL_DIVISIONS + j + 1
            n2 = (i + 1) * CIRCUMFERENTIAL_DIVISIONS + j + 1
            n3 = (i + 1) * CIRCUMFERENTIAL_DIVISIONS + jp + 1
            n4 = i * CIRCUMFERENTIAL_DIVISIONS + jp + 1
            elements.append((len(elements) + 1, n1, n2, n3, n4))
    lines = ["*Part, name=Pipe_WALL_HELPER", "*Node"]
    lines.extend(f"{label}, {x:.12g}, {y:.12g}, {z:.12g}" for label, x, y, z in nodes)
    lines.append("*Element, type=R3D4")
    lines.extend(f"{label}, {n1}, {n2}, {n3}, {n4}" for label, n1, n2, n3, n4 in elements)
    lines.extend([
        "*Elset, elset=PIPE_WALL_HELPER_ALL, generate",
        f"1, {len(elements)}, 1",
        "*Surface, type=ELEMENT, name=PIPE_WALL_HELPER_SURF",
        "PIPE_WALL_HELPER_ALL, SPOS",
        "*End Part",
    ])
    return "\n".join(lines), len(nodes), len(elements)


def read_source_curve() -> np.ndarray:
    from continuous_segment_reference import read_curve
    return read_curve(SOURCE_DXF)


def write_control_dxf(center_aba: np.ndarray, tangent_aba: np.ndarray) -> dict:
    transform = json.loads(TRANSFORM.read_text())
    origin = np.asarray(transform["origin_aba_mm"], dtype=float)
    rotation = np.asarray(transform["R_aba_to_mag"], dtype=float)
    center_mag = rotation.dot(center_aba - origin)
    tangent_mag = unit(rotation.dot(tangent_aba))

    source = read_source_curve()
    source_edge = np.diff(source, axis=0)
    source_len = np.linalg.norm(source_edge, axis=1)
    source_s = np.r_[0.0, np.cumsum(source_len)]
    com_mag = rotation.dot(np.asarray([-7.468174204284, -3.676918015967, -9.550745259298]) - origin)
    nearest = int(np.argmin(np.linalg.norm(source - com_mag, axis=1)))
    prefix = source[: nearest + 1].copy()
    if len(prefix) < 3:
        raise RuntimeError("Source curve prefix is too short to preserve the transported basis")
    # The translated prefix only transports the original transverse gauge. It
    # ends 40 mm behind the robot, well outside the physical tube and reachable
    # region; a straight continuation then crosses the complete control tube.
    gauge_end = center_mag - GAUGE_PREFIX_OFFSET_MM * tangent_mag
    prefix += gauge_end - prefix[-1]
    line = np.asarray([gauge_end + distance * tangent_mag
                       for distance in np.arange(0.1, GAUGE_PREFIX_OFFSET_MM + HALF_LENGTH_MM + 0.05, 0.1)])
    curve = np.vstack((prefix, line))
    doc = ezdxf.new("R2010")
    doc.modelspace().add_polyline3d(curve.tolist())
    path = CASE / "straight_control_centerline.dxf"
    doc.saveas(path)
    return {
        "source_projection_s_mm": float(source_s[nearest]),
        "gauge_prefix_length_mm": float(np.sum(np.linalg.norm(np.diff(prefix, axis=0), axis=1))),
        "gauge_end_mag_mm": gauge_end.tolist(),
        "straight_wall_start_mag_mm": (center_mag - HALF_LENGTH_MM * tangent_mag).tolist(),
        "straight_center_mag_mm": center_mag.tolist(),
        "straight_tangent_mag": tangent_mag.tolist(),
        "dxf_sha256": sha256(path),
    }


def validate_deck(deck: str, source_deck: str, wall_elements: int) -> None:
    if deck.count("*Part, name=Robot_SOLID") != 1:
        raise RuntimeError("Robot part identity failed")
    source_robot = re.search(r"(?ms)^\*Part, name=Robot_SOLID\s*$.*?^\*End Part\s*$", source_deck).group(0)
    control_robot = re.search(r"(?ms)^\*Part, name=Robot_SOLID\s*$.*?^\*End Part\s*$", deck).group(0)
    if source_robot != control_robot:
        raise RuntimeError("Robot part or mesh changed")
    required = [
        "1.0e-7, 0.033333333", "0.03,", "0.5", "HYDRO_FX",
        "ROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF",
        f"1, {wall_elements}, 1",
    ]
    missing = [token for token in required if token not in deck]
    if missing:
        raise RuntimeError(f"Control deck validation failed; missing {missing}")


def main() -> None:
    CASE.mkdir(parents=True, exist_ok=True)
    pose = pd.read_csv(BASE / "three_cycle_pose.csv", nrows=1).iloc[0]
    # A straight control places the unchanged robot in the tube centerline
    # reference frame. The old curved-DXF projection is not the geometric
    # center of a compatible straight circular cross-section at both robot ends.
    center = np.asarray([-7.468174204284, -3.676918015967, -9.550745259298], dtype=float)
    tangent = unit(pose[["t_x", "t_y", "t_z"]].to_numpy(float))
    e1 = unit(pose[["e1_x", "e1_y", "e1_z"]].to_numpy(float))
    e2 = unit(np.cross(tangent, e1))

    source_path = SOURCE_CASE / f"{SOURCE_JOB}.inp"
    source_deck = source_path.read_text()
    wall_text, wall_nodes, wall_elements = cylinder_part(center, tangent, e1, e2)
    deck = replace_part(source_deck, "Pipe_WALL_HELPER", wall_text)
    deck = deck.replace(SOURCE_JOB, JOB)
    deck, duration_count = re.subn(
        r"(\*Dynamic, Explicit, DIRECT\s*\n)1\.0e-7,\s*0\.100000000",
        r"\g<1>1.0e-7, 0.033333333", deck, count=1,
    )
    deck, elset_count = re.subn(
        r"(\*Elset, elset=PIPE_WALL_HELPER_CEL_ALL, instance=Pipe_WALL_HELPER-1, generate\s*\n)1,\s*764,\s*1",
        rf"\g<1>1, {wall_elements}, 1", deck, count=1,
    )
    if duration_count != 1 or elset_count != 1:
        raise RuntimeError(f"Deck replacement counts duration={duration_count}, assembly_elset={elset_count}")
    header = (
        "** STRAIGHT-PIPE CONTROL; ONLY PIPE_WALL_HELPER AND FIELD CENTERLINE GEOMETRY CHANGED\n"
        "** SAME ROBOT/MESH/CONTACT/REDUCED-HYDRO/dt/B/f/cone/bias/G/ramp/phase AS 3CYCLE\n"
    )
    deck = header + deck
    validate_deck(deck, source_deck, wall_elements)
    deck_path = CASE / f"{JOB}.inp"
    deck_path.write_text(deck)
    shutil.copyfile(SOURCE_CASE / "vuforc_production_local.f", CASE / "vuforc_production_local.f")
    dxf_meta = write_control_dxf(center, tangent)
    identity = {
        "case_id": JOB, "status": "PREPARED", "control_type": "STRAIGHT_PIPE_ONLY",
        "source_job": SOURCE_JOB, "source_commit": "ba8b5071ce1707912b8682097f687e1ef33b83f3",
        "field_frame_mode": "ROBOT_LOCAL_TANGENT", "duration_s": 0.033333333,
        "direct_dt_s": 1.0e-7, "B0_mT": 10.0, "frequency_Hz": 30.0,
        "cone_deg": 30.0, "phase_deg": 248.0, "bias_deg": 0.0, "gradient_mT": 0.0,
        "ramp_time_s": 0.001, "mu": 0.03, "zeta": 0.5,
        "robot_node_count": 1588, "robot_element_count": 7302,
        "wall_radius_mm": RADIUS_MM, "wall_half_length_mm": HALF_LENGTH_MM,
        "wall_node_count": wall_nodes, "wall_element_count": wall_elements,
        "initial_center_aba_mm": center.tolist(), "straight_tangent_aba": tangent.tolist(),
        "rejected_curved_projection_center_aba_mm": pose[["center_x", "center_y", "center_z"]].to_numpy(float).tolist(),
        "initial_e1_aba": e1.tolist(), "initial_e2_aba": e2.tolist(),
        "source_input_sha256": sha256(source_path), "input_sha256": sha256(deck_path),
        "fortran_sha256": sha256(CASE / "vuforc_production_local.f"),
        "production_server_sha256": sha256(BASE / "production" / "magpylib_socket_server.py"),
        "straight_centerline": dxf_meta,
    }
    (CASE / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n")
    print(json.dumps(identity, indent=2))


if __name__ == "__main__":
    main()
