"""Prepare the single authorized 200 ms straight magnetic-rocking case."""
from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
SOURCE = REPO / "calibration_analysis" / "StraightPipeControl" / "case" / "PROD_LOCAL30_G0_STRAIGHT_CTRL"
SOURCE_JOB = "PROD_LOCAL30_G0_STRAIGHT_CTRL"
JOB = "PROD_LOCAL_ROCK10_5HZ_G0_STRAIGHT"
CASE = OUT / "case" / JOB
SERVER = REPO / "calibration_analysis" / "ProductionLocalFrameValidation" / "production" / "magpylib_socket_server.py"
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
A0 = np.array([0.9647382600215763, -0.11887423721399708, 0.2348382536498942])
ROUTEA_GAUGE_DEG = -61.37284757596327


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unit(v):
    return np.asarray(v, dtype=float) / np.linalg.norm(v)


def rotation_from_to(source, target):
    source, target = unit(source), unit(target)
    cross = np.cross(source, target)
    return Rotation.from_rotvec(unit(cross) * math.atan2(np.linalg.norm(cross), np.dot(source, target))).as_matrix()


def rotate_robot_nodes(deck, target_axis):
    pattern = re.compile(r"(?ms)(^\*Part, name=Robot_SOLID\s*$.*?^\*Node\s*$\n)(.*?)(?=^\*)")
    match = pattern.search(deck)
    if not match:
        raise RuntimeError("Robot node block not found")
    rows = []
    for line in match.group(2).splitlines():
        values = [item.strip() for item in line.split(",")]
        rows.append((int(values[0]), np.asarray([float(x) for x in values[1:4]])))
    rotation = rotation_from_to(A0, target_axis)
    rotated = [(label, RP0 + rotation.dot(point - RP0)) for label, point in rows]
    block = "\n".join(f"{label}, {p[0]:.12g}, {p[1]:.12g}, {p[2]:.12g}" for label, p in rotated) + "\n"
    result = deck[:match.start(2)] + block + deck[match.end(2):]
    return result, rotation, len(rows)


def write_pose_consistent_hydro(source, destination, target_axis):
    text = Path(source).read_text()
    pattern = re.compile(r"      DATA a0 /0\.9647382600216D0,-0\.1188742372140D0,\n     \*          0\.2348382536499D0/")
    replacement = ("      DATA a0 /{:.13f}D0,{:.13f}D0,\n"
                   "     *          {:.13f}D0/").format(*target_axis)
    text, count = pattern.subn(replacement, text, count=1)
    if count != 1:
        raise RuntimeError("Reduced-Hydro initial body-axis replacement failed")
    Path(destination).write_text(text)


def main():
    CASE.mkdir(parents=True, exist_ok=True)
    source_deck = SOURCE / f"{SOURCE_JOB}.inp"
    source = source_deck.read_text()
    straight_identity = json.loads((SOURCE / "case_identity.json").read_text())
    target = unit(straight_identity["straight_tangent_aba"])
    e1 = unit(straight_identity["initial_e1_aba"]); e2 = unit(straight_identity["initial_e2_aba"])
    chi = math.radians(ROUTEA_GAUGE_DEG)
    n_rock = unit(math.cos(chi)*e1 + math.sin(chi)*e2)
    b_rock = unit(np.cross(target, n_rock))
    deck, initial_rotation, node_count = rotate_robot_nodes(source, target)
    deck = deck.replace(SOURCE_JOB, JOB)
    deck, duration_count = re.subn(r"(\*Dynamic, Explicit, DIRECT\s*\n)1\.0e-7,\s*0\.033333333", r"\g<1>1.0e-7, 0.200000000", deck, count=1)
    deck, field_count = re.subn(r"\*Output, field, time interval=1\.0e-4, time marks=NO", "*Output, field, time interval=5.0e-4, time marks=NO", deck, count=1)
    if duration_count != 1 or field_count != 1:
        raise RuntimeError(f"Deck replacements duration={duration_count}, field={field_count}")
    header = (
        "** ROUTEA HEAD-TAIL MAGNETIC ROCKING EXISTENCE PROOF\n"
        "** ROBOT_LOCAL_ROCKING; B0=10mT; alpha_B0=10deg; f=5Hz; phase=0; G=0\n"
        "** CLEAN INITIAL POSE: robot HEAD-to-TAIL axis aligned with straight local tangent\n"
        "** ROBOT SHAPE/RIGIDITY, ONE-WALL GENERAL CONTACT, mu=0.03, zeta=0.50, dt=1e-7, AND REDUCED HYDRO UNCHANGED\n"
    )
    deck_path = CASE / f"{JOB}.inp"
    deck_path.write_text(header + deck)
    write_pose_consistent_hydro(SOURCE / "vuforc_production_local.f", CASE / "vuforc_production_local.f", target)
    shutil.copyfile(SOURCE / "straight_control_centerline.dxf", CASE / "straight_control_centerline.dxf")
    identity = {
        "case_id": JOB, "status": "PREPARED", "source_job": SOURCE_JOB,
        "field_frame_mode": "ROBOT_LOCAL_ROCKING", "duration_s": .2, "direct_dt_s": 1e-7,
        "B0_mT": 10.0, "frequency_Hz": 5.0, "rocking_field_amplitude_deg": 10.0,
        "phase_deg": 0.0, "gradient_mT": 0.0, "cone_parameter_used": False,
        "mu": .03, "zeta": .50, "robot_length_mm": 2.4, "robot_diameter_mm": .815,
        "robot_moment_Am2": .0010876227522174417, "robot_mass_mg": 9.207793514166587,
        "initial_pose": "HEAD_TO_TAIL_AXIS_ALIGNED_WITH_LOCAL_TANGENT",
        "initial_magnetic_moment_axis_aba": target.tolist(),
        "initial_center_aba_mm": RP0.tolist(), "initial_axis_aba": target.tolist(),
        "routeA_gauge_from_production_e1_deg": ROUTEA_GAUGE_DEG,
        "n_rock_aba": n_rock.tolist(), "b_rock_aba": b_rock.tolist(),
        "initial_mesh_rotation_matrix": initial_rotation.tolist(), "robot_node_count": node_count,
        "field_output_interval_s": 5e-4, "dense_RP_and_contact_history": True,
        "source_input_sha256": sha256(source_deck), "input_sha256": sha256(deck_path),
        "fortran_sha256": sha256(CASE / "vuforc_production_local.f"),
        "centerline_sha256": sha256(CASE / "straight_control_centerline.dxf"),
        "production_server_sha256": sha256(SERVER),
    }
    (CASE / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n")
    if "ROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF" not in deck:
        raise RuntimeError("One-wall General Contact topology missing")
    if "0.03," not in deck or "\n0.5\n" not in deck or "HYDRO_FX" not in deck:
        raise RuntimeError("Contact or Reduced-Hydro invariant missing")
    print(json.dumps(identity, indent=2))


if __name__ == "__main__":
    main()
