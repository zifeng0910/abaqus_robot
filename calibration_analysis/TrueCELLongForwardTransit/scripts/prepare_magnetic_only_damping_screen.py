"""Estimate rotational scales and prepare four single-parameter damping cases."""
from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PARENT_NAME = "TRUECEL_B0P11_G2P20_F100_NOFLUID_CONTROL"
PARENT = ROOT / "case" / PARENT_NAME
PREFIX = "TRUECEL_MAGNETIC_ONLY_DAMPING"
FACTORS = (0.1, 0.3, 1.0, 3.0)
ROBOT_DENSITY = 7.80906654321e-9  # N s^2 / mm^4 in the N-mm-s deck
WATER_VISCOSITY = 7.1e-10  # N s / mm^2, frozen full-CEL material


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def block(text: str, part: str, keyword: str) -> list[str]:
    section = re.search(r"(?ms)^\*Part, name=" + re.escape(part) +
                        r"\s*$.*?^\*End Part\s*$", text).group()
    body = re.search(r"(?ms)^\*" + re.escape(keyword) +
                     r"[^\n]*\n(.*?)(?=^\*)", section).group(1)
    return [line for line in body.splitlines() if line.strip()]


def inertia_from_tetrahedra(deck: str, rp: np.ndarray) -> tuple[float, np.ndarray]:
    nodes = {}
    for line in block(deck, "Robot_SOLID", "Node"):
        values = [x.strip() for x in line.split(",")]
        nodes[int(values[0])] = np.asarray([float(x) for x in values[1:4]]) - rp
    mass = 0.0
    inertia = np.zeros((3, 3))
    for line in block(deck, "Robot_SOLID", "Element"):
        ids = [int(x.strip()) for x in line.split(",")]
        p = np.stack([nodes[i] for i in ids[1:5]])
        volume = abs(np.linalg.det(np.stack((p[1]-p[0], p[2]-p[0], p[3]-p[0])))) / 6
        dm = ROBOT_DENSITY * volume
        second = dm / 20 * (p.T @ p + np.outer(p.sum(axis=0), p.sum(axis=0)))
        inertia += np.trace(second) * np.eye(3) - second
        mass += dm
    return mass, inertia


def main() -> None:
    parent_identity = json.loads((PARENT / "case_identity.json").read_text())
    deck_path = PARENT / f"{PARENT_NAME}.inp"
    source_path = PARENT / "vuamp_precomputed_truecel.f90"
    deck = deck_path.read_text(encoding="latin1")
    source = source_path.read_text(encoding="latin1")
    rp = np.asarray(parent_identity["initial_center_aba_mm"])
    mass, inertia = inertia_from_tetrahedra(deck, rp)
    rock = np.asarray(parent_identity["b_routeA_aba"])
    rock /= np.linalg.norm(rock)
    i_rock = float(rock @ inertia @ rock)
    omega_drive = 2 * math.pi * 100
    omega_command_peak = omega_drive * math.radians(14.5)
    # Low-Re slender-cylinder transverse rotational resistance in bulk water.
    length = float(parent_identity["robot_length_mm"])
    diameter = float(parent_identity["robot_diameter_mm"])
    c_bulk = math.pi * WATER_VISCOSITY * length**3 / (
        3 * (math.log(length / diameter) - 0.66))
    c_ref = i_rock * omega_drive
    mag = np.loadtxt(PARENT / "magnetic_increment_g2p20_f100.csv",
                     delimiter=",", skiprows=1)
    torque_rock = mag[:, 6:9] @ rock
    estimate = {
        "parent": PARENT_NAME, "robot_mass_N_s2_per_mm": mass,
        "robot_mass_mg": mass * 1e9,
        "inertia_tensor_N_mm_s2": inertia.tolist(),
        "rocking_inertia_N_mm_s2": i_rock,
        "magnetic_torque_rock_p95_abs_N_mm": float(np.percentile(abs(torque_rock), 95)),
        "magnetic_torque_rock_peak_abs_N_mm": float(abs(torque_rock).max()),
        "command_omega_peak_rad_s": omega_command_peak,
        "water_viscosity_N_s_per_mm2": WATER_VISCOSITY,
        "bulk_water_slender_rod_c_N_mm_s": c_bulk,
        "screen_reference_c_I_2pif_N_mm_s": c_ref,
        "bulk_water_to_reference_ratio": c_bulk / c_ref,
        "reference_definition": "I_rock * 2*pi*f; dynamic screening scale, not calibrated fluid drag",
        "damping_model": "isotropic rotational torque -c_rot * global angular velocity at RP; no translational force",
        "parent_input_sha256": sha(deck_path),
        "parent_fortran_sha256": sha(source_path),
        "magnetic_table_sha256": parent_identity["magnetic_table_sha256"],
        "cases": [],
    }
    for factor in FACTORS:
        tag = str(factor).replace(".", "P")
        name = f"{PREFIX}_{tag}X"
        dest = ROOT / "case" / name
        if dest.exists():
            raise RuntimeError(f"Existing candidate must not be overwritten: {dest}")
        dest.mkdir()
        candidate_deck = deck.replace(PARENT_NAME, name)
        candidate_source = source.replace(PARENT_NAME, name)
        if candidate_deck == deck or candidate_source == source:
            raise RuntimeError("Parent path replacement failed")
        c_rot = factor * c_ref
        declaration = "  real(8) :: loads(6),last_time,t,u(3),ur(3),v(3),vs"
        if candidate_source.count(declaration) != 1:
            raise RuntimeError("VUAMP declarations changed")
        candidate_source = candidate_source.replace(
            declaration,
            declaration + "\n  real(8) :: vr(3)\n" +
            f"  real(8), parameter :: c_rot={f'{c_rot:.16e}'.replace('e', 'd')}", 1)
        sensor_anchor = "    dotu=dot_product(u,c)"
        if candidate_source.count(sensor_anchor) != 1:
            raise RuntimeError("VUAMP sensor insertion point changed")
        vr_lines = "\n".join(
            f"    vr({i})=VGETSENSORVALUE('RP_VR{i}',jSensorLookUpTable,sensorValues)"
            for i in range(1, 4))
        candidate_source = candidate_source.replace(sensor_anchor, vr_lines + "\n" + sensor_anchor, 1)
        log_anchor = "    write(77,'(9(ES18.10,:,\",\"))') t,s_eff,phase_deg,loads"
        if candidate_source.count(log_anchor) != 1:
            raise RuntimeError("VUAMP magnetic-log insertion point changed")
        candidate_source = candidate_source.replace(
            log_anchor, log_anchor +
            "\n    loads(4:6)=loads(4:6)-c_rot*vr", 1)
        (dest / f"{name}.inp").write_text(candidate_deck, encoding="latin1")
        (dest / "vuamp_precomputed_truecel.f90").write_text(candidate_source, encoding="latin1")
        table = PARENT / "magnetic_field_gradient_table_B0P11_A14P5.dat"
        (dest / table.name).write_bytes(table.read_bytes())
        if sha(dest / table.name) != parent_identity["magnetic_table_sha256"]:
            raise RuntimeError("Candidate magnetic table changed")
        identity = dict(parent_identity)
        identity.update({"case_id": name, "source_case": PARENT_NAME,
                         "status": "PREPARED", "dynamics_run_count": 0,
                         "rotational_damping_c_N_mm_s": c_rot,
                         "rotational_damping_factor": factor,
                         "input_sha256": sha(dest / f"{name}.inp"),
                         "fortran_sha256": sha(dest / "vuamp_precomputed_truecel.f90")})
        (dest / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n")
        estimate["cases"].append({"name": name, "factor": factor,
                                  "c_rot_N_mm_s": c_rot,
                                  "ratio_to_bulk_water": c_rot / c_bulk})
    output = ROOT / f"{PREFIX}_OFFLINE_ESTIMATE.json"
    output.write_text(json.dumps(estimate, indent=2) + "\n")
    print(json.dumps(estimate, indent=2))


if __name__ == "__main__":
    main()
