"""Prepare the sole A14.5 TAIL-gap locally refined TRUE-CEL comparison."""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
PARENT_JOB = "TRUECEL_A14P5_FORCE_FLUX_DIAG_2CYCLES"
JOB = "TRUECEL_A14P5_TAILGAP_REFINE_2CYCLES"
PARENT = OUT / "case" / PARENT_JOB
CASE = OUT / "case" / JOB
S_EDGES = np.linspace(-6.0, 6.0, 45)
BASE_Q = np.linspace(-0.75, 0.75, 21)
RADIUS = 0.667345
BUFFER_MM = 0.15

sys.path.insert(0, str(REPO / "calibration_analysis" / "TrueCELStability0p5ms" / "scripts"))
import prepare_stability_gate as stable


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def unit(values):
    values = np.asarray(values, dtype=float)
    return values / np.linalg.norm(values)


def part(text: str, name: str) -> str:
    match = re.search(rf"(?ms)^\*Part, name={re.escape(name)}\s*$.*?^\*End Part\s*$", text)
    if not match:
        raise RuntimeError(f"Part not found: {name}")
    return match.group(0)


def swept_tail_envelope(identity: dict):
    fields = np.load(PARENT / "private" / "robot_contact_fields_private.npz")
    motion = np.load(PARENT / "private" / "rp_history_private.npz")
    c = unit(identity["canonical_plus_s_axis_aba"])
    n = unit(identity["n_routeA_aba"])
    b = unit(identity["b_routeA_aba"])
    center0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    shift = float(identity["initial_axial_shift_mm"])
    pipe = center0 - shift * c - float(identity["radial_offset_n_mm"]) * n
    reference = fields["node_coordinates_mm"].astype(float)
    material_s = np.dot(reference - center0, c)
    tail_nodes = reference[material_s < -0.6]
    rows = []
    for time in fields["time"].astype(float):
        displacement = np.asarray([
            np.interp(time, motion[f"U{i}"][:, 0], motion[f"U{i}"][:, 1]) for i in (1, 2, 3)
        ])
        rotation = np.asarray([
            np.interp(time, motion[f"UR{i}"][:, 0], motion[f"UR{i}"][:, 1]) for i in (1, 2, 3)
        ])
        current = center0 + displacement + Rotation.from_rotvec(rotation).apply(tail_nodes - center0)
        relative = current - pipe
        coordinates = np.column_stack((relative @ c, relative @ n, relative @ b))
        radial = np.linalg.norm(coordinates[:, 1:], axis=1)
        near_wall = coordinates[radial >= RADIUS - 0.05]
        if len(near_wall):
            rows.append(near_wall)
    if not rows:
        raise RuntimeError("No TAIL nodes entered the 0.05-mm near-wall band")
    values = np.vstack(rows)
    raw_min = values.min(axis=0)
    raw_max = values.max(axis=0)
    buffered_min = np.maximum(raw_min - BUFFER_MM, [-6.0, -0.75, -0.75])
    buffered_max = np.minimum(raw_max + BUFFER_MM, [6.0, 0.75, 0.75])
    return {
        "raw_min_snb_mm": raw_min.tolist(),
        "raw_max_snb_mm": raw_max.tolist(),
        "buffer_mm": BUFFER_MM,
        "buffered_min_snb_mm": buffered_min.tolist(),
        "buffered_max_snb_mm": buffered_max.tolist(),
        "samples": int(len(values)),
    }


def smoothstep(value):
    value = np.clip(value, 0.0, 1.0)
    return value * value * (3.0 - 2.0 * value)


def transverse_grid_at_s(s, core_min, core_max):
    transition = 2.0 * (S_EDGES[1] - S_EDGES[0])
    if core_min <= s <= core_max:
        weight = 1.0
    elif core_min - transition < s < core_min:
        weight = smoothstep((s - (core_min - transition)) / transition)
    elif core_max < s < core_max + transition:
        weight = smoothstep(((core_max + transition) - s) / transition)
    else:
        weight = 0.0
    qn, qb = np.meshgrid(BASE_Q, BASE_Q, indexing="xy")
    radius = np.sqrt(qn*qn + qb*qb)
    old_r = np.asarray((0.0, 0.35, 0.475, 0.55, 0.625, 0.70, 0.85, np.sqrt(2.0)*0.75))
    new_r = np.asarray((0.0, 0.35, 0.50, 0.560, 0.610, 0.6120, 0.85, np.sqrt(2.0)*0.75))
    mapped_r = np.interp(radius, old_r, new_r)
    boundary_window = 1.0 - np.clip(np.maximum(np.abs(qn), np.abs(qb)) / 0.75, 0.0, 1.0)**8
    scale = np.ones_like(radius)
    nonzero = radius > 1.0e-12
    scale[nonzero] += 1.10 * weight * boundary_window[nonzero] * (
        mapped_r[nonzero] / radius[nonzero] - 1.0
    )
    return qn * scale, qb * scale, float(weight)


def build_fluid(pipe, c, n, b, robot_nodes, envelope):
    hull = ConvexHull(robot_nodes)
    plane_n, plane_d = hull.equations[:, :3], hull.equations[:, 3]
    core_min = float(envelope["buffered_min_snb_mm"][0])
    core_max = float(envelope["buffered_max_snb_mm"][0])
    transverse = [transverse_grid_at_s(s, core_min, core_max) for s in S_EDGES]
    weights = np.asarray([row[2] for row in transverse])
    nx, nn, nb = len(S_EDGES) - 1, 20, 20

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (nn + 1) * k)

    lines = [
        f"** TAIL-GAP LOCAL TRANSVERSE R-REFINE: {nx} x {nn} x {nb} EC3D8R",
        "** local near-wall spacing=0.0375 mm; two-layer smooth transition; bulk mesh unchanged",
        "*Part, name=FLUID_EULERIAN",
        "*Node",
    ]
    node_xyz = np.empty((nb + 1, nn + 1, nx + 1, 3), dtype=float)
    for k in range(nb + 1):
        for j in range(nn + 1):
            for i, s in enumerate(S_EDGES):
                qn = transverse[i][0][k, j]
                qb = transverse[i][1][k, j]
                xyz = pipe + s * c + qn * n + qb * b
                node_xyz[k, j, i] = xyz
                lines.append(f"{nid(i,j,k)}, {xyz[0]:.12f}, {xyz[1]:.12f}, {xyz[2]:.12f}")

    lines.append("*Element, type=EC3D8R")
    filled = []
    cell_volumes = []
    min_edge = np.inf
    max_transverse_edge = 0.0
    max_aspect = 0.0
    refined_target = 0
    eid = 0
    for k in range(nb):
        for j in range(nn):
            for i in range(nx):
                eid += 1
                conn = (nid(i,j,k), nid(i+1,j,k), nid(i+1,j+1,k), nid(i,j+1,k),
                        nid(i,j,k+1), nid(i+1,j,k+1), nid(i+1,j+1,k+1), nid(i,j+1,k+1))
                lines.append(f"{eid}, " + ", ".join(str(value) for value in conn))
                xyz = np.asarray((node_xyz[k,j,i], node_xyz[k,j,i+1],
                                  node_xyz[k,j+1,i+1], node_xyz[k,j+1,i],
                                  node_xyz[k+1,j,i], node_xyz[k+1,j,i+1],
                                  node_xyz[k+1,j+1,i+1], node_xyz[k+1,j+1,i]))
                cell_hull = ConvexHull(xyz)
                cell_volumes.append(float(cell_hull.volume))
                edge_pairs = ((0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),
                              (0,4),(1,5),(2,6),(3,7))
                edge_lengths = np.asarray([np.linalg.norm(xyz[a]-xyz[d]) for a,d in edge_pairs])
                transverse_lengths = edge_lengths[[1, 3, 5, 7, 8, 9, 10, 11]]
                min_edge = min(min_edge, float(edge_lengths.min()))
                max_transverse_edge = max(max_transverse_edge, float(transverse_lengths.max()))
                max_aspect = max(max_aspect, float(edge_lengths.max()/edge_lengths.min()))
                relative = xyz - pipe
                qn = relative @ n
                qb = relative @ b
                if np.max(qn*qn + qb*qb) > RADIUS**2:
                    continue
                intersects = linprog(np.zeros(3),
                                     A_ub=np.vstack((plane_n, cell_hull.equations[:, :3])),
                                     b_ub=np.r_[-plane_d, -cell_hull.equations[:, 3]], bounds=[(None,None)]*3,
                                     method="highs").success
                if not intersects:
                    filled.append(eid)
                s_mid = 0.5 * (S_EDGES[i] + S_EDGES[i+1])
                radial_max = float(np.sqrt(np.max(qn*qn + qb*qb)))
                if (core_min <= s_mid <= core_max and transverse_lengths.min() <= 0.050001
                        and radial_max >= RADIUS - 0.20):
                    refined_target += 1
    nnode = (nx + 1) * (nn + 1) * (nb + 1)
    lines += [
        "*Elset, elset=FLUID_CEL_BODY, generate", f"1, {eid}, 1",
        "*Nset, nset=FLUID_CEL_NODES, generate", f"1, {nnode}, 1",
        "*Eulerian Section, elset=FLUID_CEL_BODY", "MAT_FLUID_CEL, WATER", "*End Part", "**",
    ]
    quality = {
        "minimum_edge_length_mm": float(min_edge),
        "maximum_transverse_edge_length_mm": float(max_transverse_edge),
        "maximum_edge_aspect_ratio": float(max_aspect),
        "axial_node_refinement_weights": weights.tolist(),
        "core_s_range_mm": [core_min, core_max],
        "transition_width_each_side_mm": float(2.0 * (S_EDGES[1] - S_EDGES[0])),
    }
    return "\n".join(lines) + "\n", filled, np.asarray(cell_volumes), eid, nnode, refined_target, quality


def replace_mesh_and_initialization(deck, fluid, filled, nelem):
    deck, count = re.subn(
        r"(?ms)^\*Part, name=FLUID_EULERIAN\s*$.*?^\*End Part\s*$\n\*\*\n",
        fluid, deck, count=1,
    )
    if count != 1:
        raise RuntimeError("Fluid part replacement failed")
    assembly = (
        "*Instance, name=Fluid_EULERIAN-1, part=FLUID_EULERIAN\n*End Instance\n"
        f"*Elset, elset=FLUID_CEL_ALL, instance=Fluid_EULERIAN-1, generate\n1, {nelem}, 1\n"
        "*Elset, elset=FLUID_CEL_INIT_ALL, instance=Fluid_EULERIAN-1\n"
        + "\n".join(stable.old.chunks(filled)) + "\n**\n"
    )
    deck, count = re.subn(
        r"(?ms)^\*Instance, name=Fluid_EULERIAN-1.*?^\*\*\n(?=\*Elset, elset=PIPE_SOLID_CEL_ALL)",
        assembly, deck, count=1,
    )
    if count != 1:
        raise RuntimeError("Fluid assembly/initialization replacement failed")
    return deck


def main():
    if CASE.exists():
        raise RuntimeError(f"Refusing to overwrite existing case: {CASE}")
    identity = json.loads((PARENT / "case_identity.json").read_text())
    envelope = swept_tail_envelope(identity)
    parent_inp = PARENT / f"{PARENT_JOB}.inp"
    old_deck = parent_inp.read_text(encoding="latin1")
    deck = old_deck.replace(PARENT_JOB, JOB)
    c = unit(identity["canonical_plus_s_axis_aba"])
    n = unit(identity["n_routeA_aba"])
    b = unit(identity["b_routeA_aba"])
    center = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    pipe = center - float(identity["initial_axial_shift_mm"])*c - float(identity["radial_offset_n_mm"])*n
    robot_nodes, _ = stable.old.nodes_and_tets(part(deck, "Robot_SOLID"))
    fluid, filled, cell_volumes, nelem, nnode, refined_target, quality = build_fluid(
        pipe, c, n, b, robot_nodes, envelope
    )
    deck = replace_mesh_and_initialization(deck, fluid, filled, nelem)

    CASE.mkdir(parents=True)
    inp = CASE / f"{JOB}.inp"
    inp.write_text(deck, encoding="latin1")
    shutil.copy2(PARENT / "magnetic_field_gradient_table_A14P5.dat", CASE)
    shutil.copy2(PARENT / "magnetic_field_gradient_table_A14P5.json", CASE)
    source_sub = (PARENT / "vuamp_precomputed_truecel.f90").read_text(encoding="ascii")
    source_sub = source_sub.replace(
        str(OUT / "case" / "TRUECEL_A14P5_DIAG_2CYCLES" / "magnetic_increment.csv"),
        str(CASE / "magnetic_increment.csv"),
    )
    sub = CASE / "vuamp_precomputed_truecel.f90"
    sub.write_text(source_sub, encoding="ascii")

    water_volume = float(np.sum(cell_volumes[np.asarray(filled, dtype=int) - 1]))
    ds = np.diff(S_EDGES)
    min_cell = quality["minimum_edge_length_mm"]
    max_aspect = quality["maximum_edge_aspect_ratio"]
    for key in ("wallclock_s", "cpus", "socket_calls", "failure"):
        identity.pop(key, None)
    identity.update({
        "case_id": JOB,
        "status": "PREPARED",
        "physical_parent": PARENT_JOB,
        "dynamics_run_count": 0,
        "single_change": "TAIL-gap transverse CEL mesh resolution only",
        "physics_changed_from_parent": False,
        "solver_controls_changed_from_parent": False,
        "duration_s": 2.0 / 120.0,
        "Eulerian_dimensions": [44, 20, 20],
        "Eulerian_element_count": nelem,
        "Eulerian_node_count": nnode,
        "Eulerian_spacing_mm": [float(ds[0]), min_cell, min_cell],
        "Eulerian_spacing_ranges_mm": {
            "s": [float(ds.min()), float(ds.max())],
            "n": [min_cell, quality["maximum_transverse_edge_length_mm"]],
            "b": [min_cell, quality["maximum_transverse_edge_length_mm"]],
        },
        "tail_swept_envelope": envelope,
        "mesh_refinement_type": "LOCAL_RADIAL_R_REFINEMENT_WITH_FIXED_ELEMENT_COUNT",
        "mesh_quality": quality,
        "conforming_mesh_extension": (
            "No hanging nodes or global h-refinement: existing 20x20 transverse nodes are smoothly "
            "redistributed toward the wall only over the TAIL swept s corridor and two transition layers."
        ),
        "target_corridor_refined_elements": refined_target,
        "old_element_count": int(json.loads((PARENT / "case_identity.json").read_text())["Eulerian_element_count"]),
        "new_element_count": nelem,
        "minimum_cell_size_mm": min_cell,
        "maximum_local_aspect_ratio": max_aspect,
        "fluid_containing_element_count": len(filled),
        "fluid_initialized_volume_mm3": water_volume,
        "fluid_mass_mg": water_volume * float(identity["fluid_density_tonne_mm3"]) * 1.0e9,
        "initialization_audit_status": "PENDING",
        "robot_fluid_initial_overlap_mm3": None,
        "initial_EVF_max_geometrically_inside_robot": None,
        "input_sha256": sha(inp),
        "fortran_sha256": sha(sub),
        "parent_magnetic_table_sha256_verified": sha(PARENT / "magnetic_field_gradient_table_A14P5.dat"),
        "magnetic_table_sha256": sha(CASE / "magnetic_field_gradient_table_A14P5.dat"),
        "frozen_for_tailgap_refine": [
            "B0=12mT", "G=2.0mT", "f=120Hz", "A_main=14.5deg", "A_cross=2.5deg",
            "PRECOMPUTED_TABLE magnetic values", "robot geometry and pose", "robot-wall contact law",
            "mu and contact damping", "EOS c0=100000mm/s", "water density and viscosity",
            "axial CEL boundaries", "startup ramp", "Explicit scale factor=0.4", "two cycles",
        ],
    })
    (CASE / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n", encoding="ascii")
    setup = {
        "case": JOB,
        "single_change": identity["single_change"],
        "tail_swept_envelope": envelope,
        "mesh_dimensions": identity["Eulerian_dimensions"],
        "old_element_count": identity["old_element_count"],
        "new_element_count": nelem,
        "target_corridor_refined_elements": refined_target,
        "minimum_cell_size_mm": min_cell,
        "maximum_local_aspect_ratio": max_aspect,
        "water_volume_mm3": water_volume,
        "water_mass_mg": identity["fluid_mass_mg"],
        "conforming_mesh_extension": identity["conforming_mesh_extension"],
        "physics_frozen": True,
        "dynamics_run_count": 0,
        "input_sha256": sha(inp),
        "fortran_sha256": sha(sub),
    }
    (OUT / f"{JOB}_Setup_Audit.json").write_text(json.dumps(setup, indent=2) + "\n", encoding="ascii")
    print(json.dumps(setup, indent=2))


if __name__ == "__main__":
    main()
