"""Render the user-stopped partial F120 strict-CEL result."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
JOB = "TAIL_STICK_F120_ZETA100_8P333MS_TRUECEL"
CASE = OUT / "case" / JOB
PRIVATE = CASE / "private"
sys.path.insert(0, str(REPO / "calibration_analysis" / "S4HeadForwardLowGScreen" / "scripts"))
import analyze_and_render_lowg as base

BLUE, RED, BODY = "#2864a8", "#d1493f", "#383c42"


def interp3(archive, prefix, time):
    return np.column_stack([np.interp(time, archive[prefix + str(i)][:, 0],
                                      archive[prefix + str(i)][:, 1]) for i in (1, 2, 3)])


def main():
    identity = json.loads((CASE / "case_identity.json").read_text())
    deck = (CASE / (JOB + ".inp")).read_text(encoding="latin1")
    c = base.unit(identity["canonical_plus_s_axis_aba"])
    n = base.unit(identity["n_routeA_aba"])
    b = base.unit(identity["b_routeA_aba"])
    rp0 = np.asarray(identity["initial_center_aba_mm"], float)
    pipe_center = rp0 - identity["radial_offset_n_mm"] * n
    nodes = base.part_nodes(deck, "Robot_SOLID")
    wall = base.part_nodes(deck, "Pipe_WALL_HELPER")
    normals, wall_radius = base.wall_planes(wall, pipe_center, c, n, b)
    s_ref = (nodes - rp0).dot(c)
    region = np.full(len(nodes), "BODY", dtype=object)
    region[s_ref <= s_ref.min() + .25] = "TAIL"
    region[s_ref >= s_ref.max() - .25] = "HEAD"

    rp = np.load(PRIVATE / "rp_history_private.npz")
    time = rp["U1"][:, 0]
    U, UR = interp3(rp, "U", time), interp3(rp, "UR", time)
    rp.close()
    cfield = np.load(PRIVATE / "robot_contact_fields_private.npz")
    ctime = cfield["time"].astype(float)
    node_force = cfield["CNORMF"].astype(float) + cfield["CSHEARF"].astype(float)
    cfield.close()

    positions, tail_gap = [], []
    for frame_time in ctime:
        j = int(np.clip(np.searchsorted(time, frame_time), 0, len(time) - 1))
        pts = rp0 + U[j] + Rotation.from_rotvec(UR[j]).apply(nodes - rp0)
        positions.append(pts)
        rel = pts - pipe_center
        transverse = rel - np.outer(rel.dot(c), c)
        gap = wall_radius - transverse.dot(normals.T).max(axis=1)
        tail_gap.append(float(gap[region == "TAIL"].min()))
    positions, tail_gap = np.asarray(positions), np.asarray(tail_gap)
    tail_force = np.linalg.norm(node_force[:, region == "TAIL", :].sum(axis=1), axis=1)
    contact = (tail_gap <= .010) & (tail_force > 1e-8)
    axial = U.dot(c)
    body_axis = Rotation.from_rotvec(UR).apply(np.broadcast_to(c, UR.shape))
    rocking = np.degrees(np.arctan2(body_axis.dot(n), body_axis.dot(c)))
    energy = np.load(PRIVATE / "energy_history_private.npz")["ETOTAL"]
    field = np.load(PRIVATE / "truecel_field_private.npz")
    fcoords, fvel = field["fluid_node_coordinates_mm"], field["fluid_velocity_mm_s"].astype(float)
    fs = (fcoords - pipe_center).dot(c)
    fr = fcoords - pipe_center - np.outer(fs, c)
    near = (np.abs(fs) <= 2.0) & (np.linalg.norm(fr, axis=1) <= wall_radius + .10)
    max_fluid = float(np.nanmax(np.linalg.norm(fvel[:, near], axis=2)))
    field.close()

    metrics = {
        "case": JOB, "status": "USER_STOPPED_PARTIAL_RESULT",
        "requested_duration_ms": identity["duration_s"] * 1e3,
        "available_field_time_ms": float(ctime[-1] * 1e3),
        "available_commanded_cycles": float(ctime[-1] * identity["frequency_Hz"]),
        "frequency_Hz": identity["frequency_Hz"], "zeta": identity["zeta"],
        "delta_s_mm": float(axial[-1] - axial[0]),
        "rocking_min_deg": float(rocking.min()), "rocking_max_deg": float(rocking.max()),
        "TAIL_min_gap_mm": float(tail_gap.min()),
        "TAIL_contact_frame_count": int(contact.sum()),
        "TAIL_total_general_contact_resultant_max_N": float(tail_force.max()),
        "max_abs_energy_drift_Nmm": float(np.max(np.abs(energy[:, 1] - energy[0, 1]))),
        "near_robot_fluid_velocity_max_mm_s": max_fluid,
        "interpretation_limit": "Only 0.54 commanded cycle is available; sustained periodic rocking cannot be established.",
    }
    (OUT / "TAIL_STICK_F120_PARTIAL_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="ascii")

    colors = np.where(region == "HEAD", RED, np.where(region == "TAIL", BLUE, BODY))
    fig, ax = plt.subplots(figsize=(10, 3))
    def draw(k):
        ax.clear(); pts = positions[k]
        x, y = (pts - pipe_center).dot(c), (pts - pipe_center).dot(n)
        ax.axhline(wall_radius, color="#5996a5"); ax.axhline(-wall_radius, color="#5996a5")
        ax.scatter(x, y, c=colors, s=6, linewidths=0)
        tail, head = pts[region == "TAIL"].mean(0), pts[region == "HEAD"].mean(0)
        ax.text((tail-pipe_center).dot(c)-.08, (tail-pipe_center).dot(n), "TAIL", color=BLUE, weight="bold", ha="right")
        ax.text((head-pipe_center).dot(c)+.08, (head-pipe_center).dot(n), "HEAD", color=RED, weight="bold", ha="left")
        state = "CONTACT / SUPPORT" if contact[k] else "FREE"
        ax.text(.02, .92, state, transform=ax.transAxes, color=RED if contact[k] else BODY, weight="bold")
        ax.text(.98, .92, "USER-STOPPED PARTIAL", transform=ax.transAxes, ha="right", color="#8a5a00", weight="bold")
        ax.set(xlim=(-3, 3), ylim=(-.9, .9), aspect="equal",
               xlabel="canonical +s (mm), left to right", ylabel="n_routeA (mm)")
        ax.set_title("{} | t={:.2f} ms | f=120 Hz | zeta=1.00 | {:.2f} cycles available".format(
            JOB, ctime[k]*1e3, ctime[k]*120.0), fontsize=8)
        ax.grid(axis="x", alpha=.15); fig.tight_layout()
    FuncAnimation(fig, draw, frames=len(ctime), interval=100).save(
        OUT / "TAIL_STICK_F120_ZETA100_PARTIAL.gif", writer=PillowWriter(fps=10), dpi=100)
    plt.close(fig)

    identity.update({"status": "USER_STOPPED_PARTIAL_RESULT",
                     "simulated_time_s": float(ctime[-1]),
                     "available_commanded_cycles": float(ctime[-1] * identity["frequency_Hz"]),
                     "stop_reason": "USER_REQUESTED_STOP_AND_EXPORT",
                     "solver_end_state": "INTERRUPTED_WITH_READABLE_PARTIAL_ODB"})
    (CASE / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n", encoding="ascii")
    report = """# F120 Critically Damped TRUE-CEL Partial Screen

This run was stopped on user request and exported from the readable partial ODB.
It contains {time:.3f} ms ({cycles:.3f} commanded cycles), so it cannot establish
sustained periodic rocking. The GIF and metrics must not be treated as a completed
8.333 ms gate.

- Frequency: 120 Hz
- Contact damping: zeta=1.00
- TAIL total General Contact resultant max: {force:.9g} N
- Energy drift max: {energy:.9g} N mm
- Near-robot fluid speed max: {fluid:.9g} mm/s
""".format(time=metrics["available_field_time_ms"], cycles=metrics["available_commanded_cycles"],
           force=metrics["TAIL_total_general_contact_resultant_max_N"],
           energy=metrics["max_abs_energy_drift_Nmm"], fluid=max_fluid)
    (OUT / "TAIL_STICK_F120_PARTIAL_Report.md").write_text(report, encoding="ascii")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
