"""Fail-closed T0 audit for the TAIL-gap refined Eulerian mesh."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from scipy.optimize import linprog
from scipy.spatial import ConvexHull


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
JOB = "TRUECEL_A14P5_TAILGAP_REFINE_2CYCLES"
CASE = OUT / "case" / JOB
sys.path.insert(0, str(REPO / "calibration_analysis" / "TrueCELStability0p5ms" / "scripts"))
import prepare_stability_gate as stable


def unit(values):
    values = np.asarray(values, dtype=float)
    return values / np.linalg.norm(values)


def parse_nodes(block):
    body = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", block).group(1)
    rows = [[float(x.strip()) for x in line.split(",")[:4]] for line in body.splitlines() if line.strip()]
    return np.asarray([int(row[0]) for row in rows]), np.asarray([row[1:] for row in rows])


def parse_hexes(block):
    body = re.search(r"(?ms)^\*Element, type=EC3D8R\s*$\n(.*?)(?=^\*)", block).group(1)
    result = {}
    for line in body.splitlines():
        fields = [int(x.strip()) for x in line.split(",") if x.strip()]
        result[fields[0]] = fields[1:9]
    return result


def parse_filled(deck):
    match = re.search(
        r"(?ms)^\*Elset, elset=FLUID_CEL_INIT_ALL, instance=Fluid_EULERIAN-1\s*$\n(.*?)(?=^\*)",
        deck,
    )
    return [int(x) for x in re.findall(r"\d+", match.group(1))]


def main():
    inp = CASE / f"{JOB}.inp"
    identity_path = CASE / "case_identity.json"
    deck = inp.read_text(encoding="latin1")
    identity = json.loads(identity_path.read_text())
    c = unit(identity["canonical_plus_s_axis_aba"])
    n = unit(identity["n_routeA_aba"])
    b = unit(identity["b_routeA_aba"])
    center = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    pipe = center - float(identity["initial_axial_shift_mm"])*c - float(identity["radial_offset_n_mm"])*n

    robot, _ = stable.old.nodes_and_tets(stable.old.part(deck, "Robot_SOLID"))
    fluid_block = stable.old.part(deck, "FLUID_EULERIAN")
    labels, nodes = parse_nodes(fluid_block)
    index = {int(label): i for i, label in enumerate(labels)}
    hexes = parse_hexes(fluid_block)
    filled = parse_filled(deck)
    hull = ConvexHull(robot)
    plane_n, plane_d = hull.equations[:, :3], hull.equations[:, 3]
    intersections = []
    volume = 0.0
    slice_polygons = []
    for eid in filled:
        xyz = nodes[[index[node] for node in hexes[eid]]]
        cell_hull = ConvexHull(xyz)
        center_cell = xyz.mean(axis=0)
        ds = float(np.ptp(xyz @ c))
        dn = float(np.ptp(xyz @ n))
        db = float(np.ptp(xyz @ b))
        volume += float(cell_hull.volume)
        if linprog(np.zeros(3), A_ub=np.vstack((plane_n, cell_hull.equations[:, :3])),
                   b_ub=np.r_[-plane_d, -cell_hull.equations[:, 3]], bounds=[(None,None)]*3,
                   method="highs").success:
            intersections.append(eid)
        rel = center_cell - pipe
        if abs(rel.dot(b)) <= db * 0.51:
            s = rel.dot(c)
            qn = rel.dot(n)
            slice_polygons.append([(s-ds/2, qn-dn/2), (s+ds/2, qn-dn/2),
                                   (s+ds/2, qn+dn/2), (s-ds/2, qn+dn/2)])

    relative_robot = robot - pipe
    radial = np.sqrt((relative_robot @ n)**2 + (relative_robot @ b)**2)
    wall_gap = float(identity["lumen_radius_mm"] - radial.max())
    wall_penetration = max(0.0, -wall_gap)
    overlap_upper = 0.0
    for eid in intersections:
        xyz = nodes[[index[node] for node in hexes[eid]]]
        overlap_upper += float(ConvexHull(xyz).volume)
    max_evf_inside = 1.0 if intersections else 0.0
    mass = volume * float(identity["fluid_density_tonne_mm3"]) * 1.0e9
    passed = max_evf_inside <= 1e-6 and overlap_upper <= 1e-12 and wall_penetration <= 1e-9

    audit = {
        "case": JOB,
        "initialization_gate_passed": bool(passed),
        "audit_basis": "generated nonuniform input mesh; exact convex-hull/cell linear feasibility",
        "fluid_containing_element_count": len(filled),
        "initial_water_volume_mm3": volume,
        "initial_water_mass_mg": mass,
        "intersecting_filled_elements": len(intersections),
        "intersecting_element_ids": intersections[:100],
        "max_EVF_geometrically_inside_robot": max_evf_inside,
        "estimated_water_overlap_volume_upper_bound_mm3": overlap_upper,
        "robot_wall_initial_penetration_mm": wall_penetration,
        "minimum_robot_wall_gap_mm": wall_gap,
        "canonical_plus_s": "left_to_right",
    }
    path = OUT / f"{JOB}_PreSolve_Initialization_Audit.json"
    path.write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")
    identity.update({
        "initialization_audit_status": "PASSED" if passed else "FAILED",
        "initialization_audit_file": path.name,
        "robot_fluid_initial_overlap_mm3": overlap_upper,
        "initial_EVF_max_geometrically_inside_robot": max_evf_inside,
        "fluid_initialized_volume_mm3": volume,
        "fluid_mass_mg": mass,
    })
    identity_path.write_text(json.dumps(identity, indent=2) + "\n", encoding="ascii")

    rs = relative_robot @ c
    rn = relative_robot @ n
    outline = ConvexHull(np.column_stack((rs, rn)))
    order = np.r_[outline.vertices, outline.vertices[0]]
    fig, ax = plt.subplots(figsize=(12, 4.2))
    ax.add_collection(PolyCollection(slice_polygons, facecolor="#4d9fd6", edgecolor="white",
                                     linewidth=0.15, alpha=0.7))
    ax.plot(rs[order], rn[order], color="#25282c", lw=2.0)
    ax.axhline(identity["lumen_radius_mm"], color="#5996a5", lw=1.2)
    ax.axhline(-identity["lumen_radius_mm"], color="#5996a5", lw=1.2)
    tail = int(np.argmin(rs)); head = int(np.argmax(rs))
    ax.scatter([rs[tail], rs[head]], [rn[tail], rn[head]], c=["#2864a8", "#d1493f"], s=36)
    ax.text(rs[tail], rn[tail]+0.08, "TAIL", color="#2864a8", ha="center", weight="bold")
    ax.text(rs[head], rn[head]+0.08, "HEAD", color="#d1493f", ha="center", weight="bold")
    ax.set(xlim=(-6, 6), ylim=(-0.8, 0.8), aspect="equal", xlabel="canonical +s (mm)",
           ylabel="n_routeA (mm)", title=f"{JOB} T0 EVF audit | overlap={overlap_upper:.3g} mm3 | {'PASS' if passed else 'FAIL'}")
    fig.tight_layout()
    fig.savefig(OUT / "TRUECEL_TAILREFINE_T0_EVF_AUDIT.png", dpi=190)
    plt.close(fig)
    print(json.dumps(audit, indent=2))
    if not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
