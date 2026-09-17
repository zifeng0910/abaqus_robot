"""Prepare the one accepted 10 ms refined-geometry strict-CEL gate."""
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
JOB = "REFINED_DUALEND_F60_A14P8_G0P15_TRUECEL_10MS"
CASE = OUT / "case" / JOB
SOURCE = REPO / "calibration_analysis" / "F60G015TrueCEL" / "case" / "S4_HEADFORWARD_F60_G0P15_TRUECEL"
SOURCE_INP = SOURCE / "S4_HEADFORWARD_F60_G0P15_TRUECEL.inp"
SOURCE_ID = SOURCE / "case_identity.json"
GEOMETRY = OUT / "selected_refined_geometry.json"

CELL = 0.15
S_MIN, S_MAX = -6.0, 6.0
Q_MIN, Q_MAX = -0.75, 0.75
LUMEN_RADIUS = 0.667345
LENGTH = 2.6
NOMINAL_DIAMETER = 0.815
OFFSET_N = 0.015
HEAD_CORRECTION = 0.0425
HEAD_CORRECTION_DIRECTION = -1.0
ROBOT_EXCLUSION_RADIUS = 0.500


def unit(values):
    values = np.asarray(values, dtype=float)
    return values / np.linalg.norm(values)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def chunks(values, size=16):
    values = list(values)
    return [", ".join(str(v) for v in values[i:i + size]) for i in range(0, len(values), size)]


def parse_nodes(part):
    match = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", part)
    labels, coords = [], []
    for line in match.group(1).splitlines():
        values = [value.strip() for value in line.split(",")]
        labels.append(int(values[0]))
        coords.append([float(value) for value in values[1:4]])
    return match, np.asarray(labels, dtype=int), np.asarray(coords, dtype=float)


def make_fluid_part(pipe_center, c, n, b, robot_center):
    ns = int(round((S_MAX - S_MIN) / CELL))
    nq = int(round((Q_MAX - Q_MIN) / CELL))

    def nid(i, j, k):
        return 1 + i + (ns + 1) * (j + (nq + 1) * k)

    lines = ["** REPAIRED TRUE-CEL DOMAIN: actual lumen coverage; h=0.15 mm",
             "*Part, name=FLUID_EULERIAN", "*Node"]
    for k in range(nq + 1):
        qb = Q_MIN + k * CELL
        for j in range(nq + 1):
            qn = Q_MIN + j * CELL
            for i in range(ns + 1):
                s = S_MIN + i * CELL
                xyz = pipe_center + s * c + qn * n + qb * b
                lines.append("{}, {:.12f}, {:.12f}, {:.12f}".format(nid(i, j, k), *xyz))

    lines.append("*Element, type=EC3D8R")
    filled, eid = [], 0
    center_n = float(np.dot(robot_center - pipe_center, n))
    center_b = float(np.dot(robot_center - pipe_center, b))
    half_length = 0.5 * LENGTH + 0.10
    half_cell = 0.5 * CELL
    for k in range(nq):
        qb = Q_MIN + (k + 0.5) * CELL
        for j in range(nq):
            qn = Q_MIN + (j + 0.5) * CELL
            for i in range(ns):
                s = S_MIN + (i + 0.5) * CELL
                eid += 1
                conn = (nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                        nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1))
                lines.append("{}, {}".format(eid, ", ".join(str(v) for v in conn)))
                in_lumen = qn * qn + qb * qb <= LUMEN_RADIUS ** 2
                dn = max(abs(qn - center_n) - half_cell, 0.0)
                db = max(abs(qb - center_b) - half_cell, 0.0)
                axial_overlap = abs(s) - half_cell < half_length
                intersects_robot = axial_overlap and dn * dn + db * db < (ROBOT_EXCLUSION_RADIUS - 1e-9) ** 2
                if in_lumen and not intersects_robot:
                    filled.append(eid)
    n_nodes = (ns + 1) * (nq + 1) * (nq + 1)
    lines += ["*Elset, elset=FLUID_CEL_BODY, generate", "1, {}, 1".format(eid),
              "*Nset, nset=FLUID_CEL_NODES, generate", "1, {}, 1".format(n_nodes),
              "*Eulerian Section, elset=FLUID_CEL_BODY", "MAT_FLUID_CEL, WATER", "*End Part", "**"]
    return "\n".join(lines) + "\n", filled, eid, n_nodes


def main():
    accepted = json.loads(GEOMETRY.read_text())
    near = accepted.get("best_near_candidate_not_selected", {})
    expected = np.asarray((near.get("length_mm"), near.get("diameter_mm"), near.get("offset_n_mm"),
                           near.get("head_radial_correction_mm")), dtype=float)
    if not np.allclose(expected, (LENGTH, NOMINAL_DIAMETER, OFFSET_N, HEAD_CORRECTION), atol=1e-12):
        raise RuntimeError("accepted geometry does not match audited near candidate")

    source_id = json.loads(SOURCE_ID.read_text())
    deck = SOURCE_INP.read_text(encoding="latin1")
    pipe_center = np.asarray(source_id["initial_center_aba_mm"], dtype=float)
    c = unit(source_id["canonical_plus_s_axis_aba"])
    n = unit(source_id["n_routeA_aba"])
    b = unit(source_id["b_routeA_aba"])
    robot_center = pipe_center + OFFSET_N * n

    robot_match = re.search(r"(?ms)^\*Part, name=Robot_SOLID\s*$.*?^\*End Part\s*$", deck)
    robot_part = robot_match.group(0)
    node_match, labels, nodes = parse_nodes(robot_part)
    relative = nodes - pipe_center
    axial = relative.dot(c)
    radial = relative - np.outer(axial, c)
    base_length = axial.max() - axial.min()
    x = (axial - axial.min()) / base_length
    blend = np.clip((x - 0.75) / 0.25, 0.0, 1.0)
    blend = blend * blend * (3.0 - 2.0 * blend)
    side_weight = np.clip(HEAD_CORRECTION_DIRECTION * radial.dot(n) / (0.5 * NOMINAL_DIAMETER), 0.0, 1.0)
    corrected = (pipe_center + OFFSET_N * n + np.outer(axial * LENGTH / base_length, c) + radial +
                 np.outer(blend * side_weight * HEAD_CORRECTION * HEAD_CORRECTION_DIRECTION, n))
    node_lines = ["{}, {:.12f}, {:.12f}, {:.12f}".format(label, *xyz) for label, xyz in zip(labels, corrected)]
    transformed_axial = (corrected - robot_center).dot(c)
    tail = labels[transformed_axial <= transformed_axial.min() + 0.25]
    head = labels[transformed_axial >= transformed_axial.max() - 0.25]
    body = labels[(transformed_axial > transformed_axial.min() + 0.25) &
                  (transformed_axial < transformed_axial.max() - 0.25)]
    region_text = "\n".join(["*Nset, nset=TAIL_REGION"] + chunks(tail) +
                            ["*Nset, nset=HEAD_REGION"] + chunks(head) +
                            ["*Nset, nset=BODY_REGION"] + chunks(body)) + "\n"
    robot_part = robot_part[:node_match.start(1)] + "\n".join(node_lines) + "\n" + robot_part[node_match.end(1):]
    robot_part = robot_part.replace("*End Part", region_text + "*End Part", 1)
    deck = deck[:robot_match.start()] + robot_part + deck[robot_match.end():]

    fluid_match = re.search(r"(?ms)^\*Part, name=FLUID_EULERIAN\s*$.*?^\*End Part\s*$\n\*\*\n", deck)
    fluid_part, filled, n_elem, n_nodes = make_fluid_part(pipe_center, c, n, b, robot_center)
    deck = deck[:fluid_match.start()] + fluid_part + deck[fluid_match.end():]
    deck = re.sub(r"(?ms)^\*Elset, elset=FLUID_CEL_ALL, instance=Fluid_EULERIAN-1, generate\n.*?^\*\*\n",
                  "*Elset, elset=FLUID_CEL_ALL, instance=Fluid_EULERIAN-1, generate\n1, {}, 1\n"
                  "*Elset, elset=FLUID_CEL_INIT_ALL, instance=Fluid_EULERIAN-1\n{}\n**\n".format(
                      n_elem, "\n".join(chunks(filled))), deck, count=1)

    old_rp = "      2, {:.12f}, {:.12f}, {:.12f}".format(*pipe_center)
    new_rp = "      2, {:.12f}, {:.12f}, {:.12f}".format(*robot_center)
    if old_rp not in deck:
        raise RuntimeError("RP_ROBOT coordinate record not found")
    deck = deck.replace(old_rp, new_rp, 1)
    deck = deck.replace(", 0.033333333333", ", 0.010000000000", 1)
    deck = deck.replace("S4_HEADFORWARD_F60_G0P15_TRUECEL", JOB)
    deck = deck.replace("A_main=14.343111711438091deg", "A_main=14.8deg")
    deck = deck.replace("duration=0.033333333333s", "duration=0.010000000000s")
    deck = deck.replace("*Contact Output, general contact\nCSTRESS, CFORCE",
                        "*Contact Output, general contact\nCSTRESS, CFORCE, CDISP", 1)

    CASE.mkdir(parents=True, exist_ok=True)
    inp = CASE / (JOB + ".inp")
    inp.write_text(deck, encoding="latin1")
    shutil.copy2(SOURCE / "straight_control_centerline.dxf", CASE / "straight_control_centerline.dxf")
    shutil.copy2(SOURCE / "vuforc_magnetic_only.f", CASE / "vuforc_magnetic_only.f")
    fortran_path = CASE / "vuforc_magnetic_only.f"
    fortran = fortran_path.read_text(encoding="latin1")
    old_center = ["-7.468174204284D0", "-3.676918015967D0", "-9.550745259298D0"]
    new_center = ["{:.12f}D0".format(value) for value in robot_center]
    for old, new in zip(old_center, new_center):
        if old not in fortran:
            raise RuntimeError("VUAMP RP center token missing: " + old)
        fortran = fortran.replace(old, new, 1)
    fortran_path.write_text(fortran, encoding="latin1")

    filled_volume = len(filled) * CELL ** 3
    identity = dict(source_id)
    identity.update({
        "case_id": JOB, "status": "PREPARED", "duration_s": 0.01, "dynamics_run_count": 0,
        "rocking_main_amplitude_deg": 14.8, "rocking_cross_amplitude_deg": 2.5,
        "robot_length_mm": LENGTH, "robot_diameter_mm": NOMINAL_DIAMETER,
        "initial_center_aba_mm": robot_center.tolist(), "radial_offset_n_mm": OFFSET_N,
        "head_side_correction_mm": HEAD_CORRECTION, "head_correction_direction_n": -1,
        "geometry_realization": "actual smooth one-sided HEAD mesh profile change over final 25% axial length",
        "lumen_radius_mm": LUMEN_RADIUS, "Eulerian_transverse_bounds_mm": [Q_MIN, Q_MAX],
        "Eulerian_element_count": n_elem, "fluid_element_count": len(filled),
        "fluid_filled_volume_mm3": filled_volume, "fluid_mass_mg": filled_volume * 1.0007,
        "contact_regions": {"HEAD_REGION_nodes": len(head), "TAIL_REGION_nodes": len(tail), "BODY_REGION_nodes": len(body)},
        "contact_output": "direct general-contact CSTRESS/CFORCE/CDISP fields partitioned by actual node regions; pair-isolated wall history unavailable",
        "input_sha256": sha256(inp), "fortran_sha256": sha256(CASE / "vuforc_magnetic_only.f"),
    })
    (CASE / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n", encoding="ascii")
    audit = {"strict_CEL": True, "EVF_initial": 1.0, "fluid_mass_mg": identity["fluid_mass_mg"],
             "filled_elements": len(filled), "total_elements": n_elem, "total_nodes": n_nodes,
             "lumen_radius_mm": LUMEN_RADIUS, "robot_exclusion_radius_mm": ROBOT_EXCLUSION_RADIUS,
             "ReducedHydro": "OFF", "Magpylib_socket": "ON", "U_flow_mm_s": 10.0}
    (OUT / "accepted_truecel_setup_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")
    print(json.dumps({"job": JOB, **audit, "case": str(CASE)}, indent=2))


if __name__ == "__main__":
    main()
