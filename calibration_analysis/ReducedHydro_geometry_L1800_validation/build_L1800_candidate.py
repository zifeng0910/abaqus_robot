"""Build one physically scaled L=1.800 mm geometry candidate; never launches Abaqus."""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORKDIR = REPO.parent
BASE = REPO / "calibration_analysis" / "ReducedHydro_zeta050_8p333_validation"
SOURCE_JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_WallOn_Free_0083"
JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083"
TARGET_LENGTH_MM = 1.800
BASE_MOMENT_AM2 = 0.001168
AXIS = np.array([0.9647382600216, -0.1188742372140, 0.2348382536499], dtype=float)
AXIS /= np.linalg.norm(AXIS)

sys.path.insert(0, str(REPO / "calibration_analysis" / "ReducedHydro_hidden_impact_audit"))
from audit_stage_a import mesh_properties


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_robot_nodes(text, center, scale, expected_node_count):
    lines = text.splitlines()
    inside_part = False
    inside_nodes = False
    changed = 0
    output = []
    for line in lines:
        stripped = line.strip()
        lower = stripped.lower()
        if lower.startswith("*part,"):
            inside_part = "name=robot_solid" in lower.replace(" ", "")
            inside_nodes = False
        elif inside_part and lower.startswith("*end part"):
            inside_part = False
            inside_nodes = False
        elif inside_part and lower == "*node":
            inside_nodes = True
        elif inside_nodes and stripped.startswith("*"):
            inside_nodes = False
        if inside_nodes and stripped and not stripped.startswith("*"):
            fields = [item.strip() for item in line.split(",")]
            if len(fields) >= 4:
                label = int(fields[0])
                point = np.array([float(value) for value in fields[1:4]])
                offset = point - center
                parallel = np.dot(offset, AXIS) * AXIS
                scaled = center + (offset - parallel) + scale * parallel
                output.append(f"{label}, {scaled[0]:.12g}, {scaled[1]:.12g}, {scaled[2]:.12g}")
                changed += 1
                continue
        output.append(line)
    assert changed == expected_node_count, (changed, expected_node_count)
    return "\n".join(output) + "\n"


def dimensions(mesh):
    nodes = np.array(list(mesh["nodes"].values()), dtype=float) + mesh["shift"]
    axial = nodes @ AXIS
    center = mesh["com"]
    radial = np.linalg.norm((nodes - center) - ((nodes - center) @ AXIS)[:, None] * AXIS, axis=1)
    return float(np.ptp(axial)), float(2 * radial.max())


def main():
    source_path = BASE / f"{SOURCE_JOB}.inp"
    source = source_path.read_text()
    meshes, rp, _ = mesh_properties(source)
    base_mesh = meshes["Robot_SOLID"]
    base_length, base_diameter = dimensions(base_mesh)
    scale = TARGET_LENGTH_MM / base_length
    center_part = base_mesh["com"] - base_mesh["shift"]
    candidate = replace_robot_nodes(source, center_part, scale, len(base_mesh["nodes"]))
    candidate = candidate.replace(
        "** ROBOT TRANSVERSE SHRINK: PCA minor axes scaled to 0.90; long axis unchanged",
        "** ROBOT GEOMETRY CANDIDATE: transverse geometry unchanged; PCA long axis scaled to 1.800 mm",
        1,
    )
    local = HERE / f"{JOB}.inp"
    work = WORKDIR / f"{JOB}.inp"
    local.write_text(candidate)
    work.write_text(candidate)

    new_mesh = mesh_properties(candidate)[0]["Robot_SOLID"]
    new_length, new_diameter = dimensions(new_mesh)
    volume_ratio = new_mesh["mass"] / base_mesh["mass"]
    moment = BASE_MOMENT_AM2 * volume_ratio
    assert abs(new_length - TARGET_LENGTH_MM) < 1e-9
    assert abs(new_diameter - base_diameter) < 1e-9
    assert np.linalg.norm(new_mesh["com"] - base_mesh["com"]) < 1e-10
    identity = {
        "job": JOB,
        "source_job": SOURCE_JOB,
        "source_sha256": sha(source_path),
        "candidate_sha256": sha(local),
        "geometry_change_only": "PCA long-axis coordinates scaled about the volume COM",
        "dimension_definition": "all Robot_SOLID nodes projected onto the frozen PCA axis about the volume COM",
        "base_pca_axial_span_mm": base_length,
        "candidate_pca_axial_span_mm": new_length,
        "axial_scale": scale,
        "base_pca_transverse_diameter_mm": base_diameter,
        "candidate_pca_transverse_diameter_mm": new_diameter,
        "legacy_unshrunk_reported_length_mm": 2.9314146861447767,
        "legacy_unshrunk_reported_diameter_mm": 0.9030748963882974,
        "baseline_mesh_identity_correction": "the validated zeta050 INP already carries a 0.90 transverse PCA shrink; legacy 2.931 x 0.903 dimensions describe the pre-shrink CEL geometry",
        "density_tonne_per_mm3": 7.80906654321e-9,
        "base_mass_mg": base_mesh["mass"] * 1e9,
        "candidate_mass_mg": new_mesh["mass"] * 1e9,
        "mass_volume_ratio": volume_ratio,
        "base_magnetic_moment_Am2": BASE_MOMENT_AM2,
        "candidate_magnetic_moment_Am2": moment,
        "magnetization_rule": "constant Br; explicit moment scaled by mesh volume ratio",
        "base_inertia_tonne_mm2": base_mesh["I"].tolist(),
        "candidate_inertia_tonne_mm2": new_mesh["I"].tolist(),
        "COM_abaqus_mm": new_mesh["com"].tolist(),
        "RP_abaqus_mm": rp.tolist(),
        "frozen": {
            "duration_s": 0.008333,
            "direct_dt_s": 1e-7,
            "field_interval_s": 2.5e-5,
            "critical_damping_fraction": 0.50,
            "tangent_fraction": 0.0,
            "friction": 0.03,
            "B0_T": 0.010,
            "frequency_Hz": 30.0,
            "cone_deg": 30.0,
            "bias_deg": 40.0,
            "phase_deg": 248.0,
            "sense": 1,
            "gradient_T": 0.006,
            "gradient_length_mm": 45.0,
            "Cparallel_Ns_per_mm": 4e-9,
            "Cperp_Ns_per_mm": 1.2e-8,
            "Kspin_Nmm_s": 1e-9,
            "Kwobble_Nmm_s": 3e-9,
        },
    }
    (HERE / "L1800_input_identity.json").write_text(json.dumps(identity, indent=2))
    print(json.dumps(identity, indent=2))


if __name__ == "__main__":
    main()
