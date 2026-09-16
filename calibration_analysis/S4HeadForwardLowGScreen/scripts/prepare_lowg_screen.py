"""Prepare the four new head-forward low-gradient dynamics cases."""
from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
SOURCE_JOB = "FAST_ELLIPTIC_20HZ_A14P343_X2P5_FORWARD"
SOURCE = REPO / "calibration_analysis" / "FastStraightDynamicScreen" / "cases" / SOURCE_JOB
SERVER = REPO / "calibration_analysis" / "FastStraightDynamicScreen" / "production" / "magpylib_socket_server_fast.py"
CASES = (
    ("S4_HEADFORWARD_G0", 0.0),
    ("S4_HEADFORWARD_G0P25", 0.25),
    ("S4_HEADFORWARD_G0P5", 0.50),
    ("S4_HEADFORWARD_G1P0", 1.00),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unit(value):
    value = np.asarray(value, dtype=float)
    return value / np.linalg.norm(value)


def robot_nodes(deck: str):
    part = re.search(r"(?ms)^\*Part, name=Robot_SOLID\s*$.*?^\*End Part\s*$", deck).group(0)
    match = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", part)
    rows = []
    for line in match.group(1).splitlines():
        values = [item.strip() for item in line.split(",")]
        if values and values[0]:
            rows.append((int(values[0]), np.asarray([float(x) for x in values[1:4]])))
    return rows


def replace_robot_nodes(deck: str, rp0, rotation):
    pattern = re.compile(r"(?ms)(^\*Part, name=Robot_SOLID\s*$.*?^\*Node\s*$\n)(.*?)(?=^\*)")
    match = pattern.search(deck)
    if not match:
        raise RuntimeError("Robot_SOLID node block not found")
    rows = []
    for line in match.group(2).splitlines():
        values = [item.strip() for item in line.split(",")]
        if not values or not values[0]:
            continue
        label = int(values[0])
        point = np.asarray([float(x) for x in values[1:4]])
        point = rp0 + rotation.dot(point - rp0)
        rows.append(f"{label}, {point[0]:.12g}, {point[1]:.12g}, {point[2]:.12g}")
    block = "\n".join(rows) + "\n"
    return deck[:match.start(2)] + block + deck[match.end(2):], len(rows)


def replace_hydro_axis(text: str, axis):
    pattern = re.compile(
        r"      DATA a0 /0\.9762799602464D0,-0\.0618320247447D0,\n"
        r"     \*          0\.2074951564187D0/")
    replacement = ("      DATA a0 /{:.13f}D0,{:.13f}D0,\n"
                   "     *          {:.13f}D0/").format(*axis)
    text, count = pattern.subn(replacement, text, count=1)
    if count != 1:
        raise RuntimeError("ReducedHydro a0 replacement failed")
    return text


def main():
    source_identity = json.loads((SOURCE / "case_identity.json").read_text())
    source_deck_path = SOURCE / f"{SOURCE_JOB}.inp"
    source_deck = source_deck_path.read_text()
    if source_identity.get("status") != "SOLVED":
        raise RuntimeError("old S4 source is not solved")
    if "** Reduced-Hydro: CEL fluid part, material, initialization and CEL contact were removed." not in source_deck:
        raise RuntimeError("old S4 is not the verified non-CEL deck")

    rp0 = np.asarray(source_identity["initial_center_aba_mm"], dtype=float)
    c = unit(source_identity["initial_axis_aba"])
    n = unit(source_identity["n_routeA_aba"])
    b = unit(source_identity["b_routeA_aba"])
    # A proper 180 degree rigid rotation about the fixed RouteA side-view
    # vertical axis reverses the semantic head/tail order and leaves COM and
    # the circular tube relation unchanged.
    rotation = 2.0 * np.outer(n, n) - np.eye(3)
    if not np.allclose(rotation.dot(c), -c, atol=1e-12):
        raise RuntimeError("flip rotation does not reverse c_hat")
    target_moment = unit(rotation.dot(c))
    target_hydro_axis = target_moment.copy()

    source_rows = robot_nodes(source_deck)
    source_x = np.asarray([row[1] for row in source_rows])
    source_s = (source_x - rp0).dot(c)
    source_head = source_x[source_s <= source_s.min() + 0.25].mean(axis=0)
    source_tail = source_x[source_s >= source_s.max() - 0.25].mean(axis=0)
    flipped_deck, node_count = replace_robot_nodes(source_deck, rp0, rotation)
    flipped_x = rp0 + (rotation.dot((source_x - rp0).T)).T
    flipped_head = rp0 + rotation.dot(source_head - rp0)
    flipped_tail = rp0 + rotation.dot(source_tail - rp0)
    if float(np.dot(flipped_head - flipped_tail, c)) <= 0.0:
        raise RuntimeError("flipped semantic HEAD is not forward of TAIL")

    cases_root = OUT / "cases"
    cases_root.mkdir(parents=True, exist_ok=True)
    hydro_source = SOURCE / "vuforc_production_local.f"
    hydro_text = replace_hydro_axis(hydro_source.read_text(), target_hydro_axis)
    for job, gradient_mT in CASES:
        case = cases_root / job
        case.mkdir(parents=True, exist_ok=True)
        deck = flipped_deck.replace(SOURCE_JOB, job)
        header = (
            "** HEAD-FORWARD LOW-G DYNAMIC SCREEN\n"
            f"** {job}; f=20Hz; A_main=14.343111711438091deg; A_cross=2.5deg; "
            f"G={gradient_mT:.8g}mT; L=45mm; duration=0.075s\n"
            "** TRUE RIGID 180deg MESH FLIP ABOUT n_routeA; COM/TUBE/CONTACT/HYDRO COEFFICIENTS FROZEN\n"
            "** NON-CEL; ROBOT_LOCAL_ELLIPTIC_ROCKING; FIXED ROUTEA GAUGE; DT=1e-7s\n"
        )
        deck_path = case / f"{job}.inp"
        deck_path.write_text(header + deck, encoding="ascii")
        (case / "vuforc_production_local.f").write_text(hydro_text, encoding="ascii")
        shutil.copyfile(SOURCE / "straight_control_centerline.dxf", case / "straight_control_centerline.dxf")
        identity = {
            "case_id": job,
            "status": "PREPARED",
            "source_job": SOURCE_JOB,
            "field_frame_mode": "ROBOT_LOCAL_ELLIPTIC_ROCKING",
            "duration_s": 0.075,
            "direct_dt_s": 1e-7,
            "B0_mT": 10.0,
            "frequency_Hz": 20.0,
            "rocking_main_amplitude_deg": 14.343111711438091,
            "rocking_cross_amplitude_deg": 2.5,
            "gradient_mT": gradient_mT,
            "gradient_length_mm": 45.0,
            "gradient_profile": "legacy",
            "forward_definition": "canonical increasing centerline arclength +s",
            "mu": 0.03,
            "zeta": 0.50,
            "robot_length_mm": source_identity["robot_length_mm"],
            "robot_diameter_mm": source_identity["robot_diameter_mm"],
            "robot_moment_Am2": source_identity["robot_moment_Am2"],
            "robot_mass_mg": source_identity["robot_mass_mg"],
            "initial_center_aba_mm": rp0.tolist(),
            "canonical_plus_s_axis_aba": c.tolist(),
            "head_tail_axis_aba": c.tolist(),
            "initial_magnetic_moment_axis_aba": target_moment.tolist(),
            "reduced_hydro_body_axis_aba": target_hydro_axis.tolist(),
            "n_routeA_aba": n.tolist(),
            "b_routeA_aba": b.tolist(),
            "routeA_gauge_from_production_e1_deg": source_identity["routeA_gauge_from_production_e1_deg"],
            "routeA_gauge_for_flipped_body_deg": -source_identity["routeA_gauge_from_production_e1_deg"],
            "mesh_flip_axis_aba": n.tolist(),
            "mesh_flip_rotation_matrix": rotation.tolist(),
            "source_head_centroid_aba_mm": source_head.tolist(),
            "source_tail_centroid_aba_mm": source_tail.tolist(),
            "flipped_head_centroid_aba_mm": flipped_head.tolist(),
            "flipped_tail_centroid_aba_mm": flipped_tail.tolist(),
            "head_tail_dot_canonical_plus_s_mm": float(np.dot(flipped_head - flipped_tail, c)),
            "robot_node_count": node_count,
            "old_s4_reference_gradient_mT": 6.0,
            "source_input_sha256": sha256(source_deck_path),
            "input_sha256": sha256(deck_path),
            "fortran_sha256": sha256(case / "vuforc_production_local.f"),
            "centerline_sha256": sha256(case / "straight_control_centerline.dxf"),
            "fast_server_sha256": sha256(SERVER),
        }
        (case / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n", encoding="ascii")
        required = (
            "*Rigid Body, ref node=RP_ROBOT, elset=ROBOT_SOLID_CEL_ALL",
            "ROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF",
            "0.03,", "\n0.5\n", "HYDRO_FX")
        if not all(token in deck for token in required):
            raise RuntimeError(f"{job}: frozen invariant missing")
        if "*Eulerian" in deck or "*Eulerian Section" in deck:
            raise RuntimeError(f"{job}: unexpected true CEL definition")
        print(f"PREPARED {job} G={gradient_mT:g}mT head_tail_dot={identity['head_tail_dot_canonical_plus_s_mm']:.6f}mm")


if __name__ == "__main__":
    main()
