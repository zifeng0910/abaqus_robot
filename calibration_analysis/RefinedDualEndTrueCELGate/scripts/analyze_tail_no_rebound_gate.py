"""Analyze the zeta=0.80, 8 ms strict-CEL TAIL no-rebound gate."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
JOB = "TAIL_NO_REBOUND_ZETA080_8MS_TRUECEL"
CASE = OUT / "case" / JOB
PRIVATE = CASE / "private"
sys.path.insert(0, str(REPO / "calibration_analysis" / "S4HeadForwardLowGScreen" / "scripts"))
import analyze_and_render_lowg as base

BLUE, RED, BODY = "#2864a8", "#d1493f", "#383c42"


def interp3(archive, prefix, time):
    return np.column_stack([np.interp(time, archive[prefix + str(i)][:, 0],
                                      archive[prefix + str(i)][:, 1]) for i in (1, 2, 3)])


def runs(active):
    result, start = [], None
    for i, value in enumerate(active):
        if value and start is None:
            start = i
        if start is not None and (not value or i == len(active) - 1):
            end = i if value and i == len(active) - 1 else i - 1
            result.append((start, end))
            start = None
    return result


def longest(time, active):
    return max((time[j] - time[i] for i, j in runs(active)), default=0.0)


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
    region[s_ref <= s_ref.min() + 0.25] = "TAIL"
    region[s_ref >= s_ref.max() - 0.25] = "HEAD"

    rp = np.load(PRIVATE / "rp_history_private.npz")
    time = rp["U1"][:, 0]
    U, UR = interp3(rp, "U", time), interp3(rp, "UR", time)
    rp.close()
    axial = U.dot(c)
    body_axis = Rotation.from_rotvec(UR).apply(np.broadcast_to(c, UR.shape))
    rocking = np.degrees(np.arctan2(body_axis.dot(n), body_axis.dot(c)))

    contact = np.load(PRIVATE / "robot_contact_fields_private.npz")
    ctime = contact["time"].astype(float)
    if not np.array_equal(contact["node_labels"], np.arange(1, len(nodes) + 1)):
        raise RuntimeError("robot node ordering changed")
    node_force = contact["CNORMF"].astype(float) + contact["CSHEARF"].astype(float)
    contact.close()
    positions, gaps = [], {key: [] for key in ("HEAD", "TAIL", "BODY")}
    for frame_time in ctime:
        j = int(np.clip(np.searchsorted(time, frame_time), 0, len(time) - 1))
        pts = rp0 + U[j] + Rotation.from_rotvec(UR[j]).apply(nodes - rp0)
        positions.append(pts)
        rel = pts - pipe_center
        transverse = rel - np.outer(rel.dot(c), c)
        node_gap = wall_radius - transverse.dot(normals.T).max(axis=1)
        for key in gaps:
            gaps[key].append(float(node_gap[region == key].min()))
    positions = np.asarray(positions)
    gaps = {key: np.asarray(value) for key, value in gaps.items()}
    force = {key: node_force[:, region == key, :].sum(axis=1) for key in gaps}
    force["TOTAL"] = node_force.sum(axis=1)

    proximity = {key: gaps[key] <= 0.005 for key in gaps}
    # The analytic wall-plane gap differs from the faceted General Contact surface
    # by about 0.002 mm here. Require both a 0.010 mm geometric neighborhood and
    # a directly reported contact resultant so the first real support event is not
    # discarded solely by the plane approximation.
    tail_contact = ((gaps["TAIL"] <= 0.010) &
                    (np.linalg.norm(force["TAIL"], axis=1) > 1.0e-8))
    tail_runs = runs(tail_contact)
    tail_vgap = np.gradient(gaps["TAIL"], ctime)
    first = tail_runs[0] if tail_runs else None
    impact_i = first[0] if first else None
    release_i = first[1] + 1 if first and first[1] + 1 < len(ctime) else None
    if first:
        peak_i = first[0] + int(np.argmin(gaps["TAIL"][first[0]:first[1] + 1]))
        approach_start = max(0, first[0] - 3)
        rebound_end = tail_runs[1][0] if len(tail_runs) > 1 else len(ctime) - 1
        impact_v = float(np.min(tail_vgap[approach_start:peak_i + 1]))
        rebound_v = float(np.max(tail_vgap[peak_i:rebound_end + 1]))
    else:
        impact_v = rebound_v = None
    rv = abs(rebound_v) / abs(impact_v) if impact_v not in (None, 0.0) and rebound_v is not None else None
    rock_c = np.interp(ctime, time, rocking)
    center_i = len(ctime) - 1
    if impact_i is not None:
        sign0 = np.sign(rock_c[impact_i])
        candidates = np.flatnonzero(np.sign(rock_c[impact_i + 1:]) != sign0)
        if len(candidates):
            center_i = impact_i + 1 + int(candidates[0])
    same_halfcycle = [episode for episode in tail_runs if episode[0] <= center_i]
    second_touch = len(same_halfcycle) >= 2

    energy = np.load(PRIVATE / "energy_history_private.npz")["ETOTAL"]
    drift = energy[:, 1] - energy[0, 1]
    field = np.load(PRIVATE / "truecel_field_private.npz")
    fcoords = field["fluid_node_coordinates_mm"]
    fvel = field["fluid_velocity_mm_s"].astype(float)
    fevf = field["fluid_evf"].astype(float)
    fs = (fcoords - pipe_center).dot(c)
    fr = fcoords - pipe_center - np.outer(fs, c)
    near = (np.abs(fs) <= 2.0) & (np.linalg.norm(fr, axis=1) <= wall_radius + 0.10)
    max_fluid = float(np.nanmax(np.linalg.norm(fvel[:, near], axis=2)))
    evf_v_present = bool(np.any(fevf[0] > 0) and np.any(fevf[-1] > 0) and np.isfinite(fvel).any())
    field.close()

    telemetry = pd.read_csv(CASE / (JOB + "_telemetry.csv")).drop_duplicates("t_s").sort_values("t_s")
    sta = (CASE / (JOB + ".sta")).read_text(encoding="latin1")
    rows = re.findall(r"(?m)^\s*(\d+)\s+([0-9.E+-]+)\s+[0-9.E+-]+\s+\S+\s+([0-9.E+-]+)\s+\d+", sta)
    dts = np.asarray([float(row[2]) for row in rows])
    max_contact = float(np.linalg.norm(force["TOTAL"], axis=1).max())
    min_gap = min(float(value.min()) for value in gaps.values())
    both = proximity["HEAD"] & proximity["TAIL"]
    tumble = bool(np.any(body_axis.dot(c) < 0.0))
    gates = {
        "energy_drift_below_0p01_Nmm": float(np.max(np.abs(drift))) < 0.01,
        "near_robot_fluid_velocity_below_10000_mm_s": max_fluid < 10000.0,
        "robot_total_general_contact_resultant_below_0p1_N": max_contact < 0.1,
        "no_tumble": not tumble,
        "no_gross_penetration": min_gap > -0.05,
        "no_persistent_both_bridge": longest(ctime, both) <= 0.001,
        "stable_dt_no_collapse": float(dts.min()) > 1e-8,
        "EVF_and_V_present": evf_v_present,
        "socket_continuous": len(telemetry) > 10 and telemetry.t_s.max() >= 0.0079,
    }
    gates = {key: bool(value) for key, value in gates.items()}
    if second_touch:
        classification = "TRUECEL_CONTACT_GATE_FAILED_REBOUND"
    elif not all(gates.values()):
        classification = "CONTACT_REBOUND_FIXED_BUT_CEL_NUMERICS_FAIL"
    else:
        classification = "TRUECEL_CONTACT_GATE_PASSED_NO_REBOUND"

    metrics = {
        "case": JOB, "classification": classification, "simulated_time_ms": float(time[-1] * 1e3),
        "tail_gap_threshold_mm": 0.005,
        "tail_contact_event_definition": "direct TAIL General Contact resultant >1e-8 N and analytic wall-plane gap <=0.010 mm; 0.010 mm accounts for the measured faceted-wall/plane offset",
        "first_TAIL_touch_ms": None if impact_i is None else float(ctime[impact_i] * 1e3),
        "first_TAIL_release_ms": None if release_i is None else float(ctime[release_i] * 1e3),
        "second_TAIL_touch_ms": None if len(same_halfcycle) < 2 else float(ctime[same_halfcycle[1][0]] * 1e3),
        "N_tail_contacts_same_halfcycle": len(same_halfcycle),
        "TAIL_impact_gap_velocity_mm_s": impact_v,
        "TAIL_rebound_gap_velocity_mm_s": rebound_v,
        "R_v": rv,
        "rocking_min_deg": float(rocking.min()), "rocking_max_deg": float(rocking.max()),
        "delta_s_mm": float(axial[-1] - axial[0]),
        "HEAD_min_gap_mm": float(gaps["HEAD"].min()), "TAIL_min_gap_mm": float(gaps["TAIL"].min()),
        "BODY_min_gap_mm": float(gaps["BODY"].min()),
        "both_bridge_longest_ms": float(longest(ctime, both) * 1e3),
        "max_abs_energy_drift_Nmm": float(np.max(np.abs(drift))),
        "near_robot_fluid_velocity_max_mm_s": max_fluid,
        "robot_total_general_contact_resultant_max_N": max_contact,
        "wall_only_force_directly_available": False,
        "contact_force_limitation": "ODB exposes robot-total General Contact; robot-wall and robot-fluid are not pair-isolated.",
        "minimum_stable_dt_s": float(dts.min()), "total_increments": max(int(row[0]) for row in rows),
        "socket_last_time_ms": float(telemetry.t_s.max() * 1e3), "gates": gates,
    }
    (OUT / "TAIL_NO_REBOUND_8MS_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="ascii")

    pd.DataFrame({"time_s": ctime, "TAIL_gap_mm": gaps["TAIL"], "TAIL_gap_velocity_mm_s": tail_vgap,
                  "TAIL_wall_proximity": proximity["TAIL"].astype(int),
                  "TAIL_contact_event": tail_contact.astype(int),
                  "TAIL_total_general_contact_resultant_N": np.linalg.norm(force["TAIL"], axis=1),
                  "rocking_deg": rock_c}).to_csv(OUT / "TAIL_NO_REBOUND_8MS_metrics.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    ax.plot(energy[:, 0] * 1e3, drift, color="#2f5d62")
    ax.axhline(0.01, color=RED, ls="--"); ax.axhline(-0.01, color=RED, ls="--")
    ax.set(xlabel="time (ms)", ylabel="ETOTAL drift (N mm)", title="8 ms TRUE-CEL energy gate")
    ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(OUT / "TAIL_NO_REBOUND_ENERGY.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(ctime * 1e3, gaps["TAIL"], color=BLUE, marker="o", ms=3, label="TAIL geometric gap")
    ax.axhline(.005, color=RED, ls="--", label="0.005 mm nominal proximity threshold")
    ax.axhline(.010, color="#9a6b00", ls=":", label="0.010 mm faceted-wall contact neighborhood")
    ax.fill_between(ctime * 1e3, ax.get_ylim()[0], ax.get_ylim()[1], where=tail_contact, color=RED, alpha=.14)
    ax.set(xlabel="time (ms)", ylabel="TAIL minimum wall gap (mm)", title="TAIL impact / release / recontact diagnosis")
    ax.grid(alpha=.2); ax.legend(); fig.tight_layout(); fig.savefig(OUT / "TAIL_NO_REBOUND_CONTACT_DIAGNOSIS.png", dpi=180); plt.close(fig)

    frame_ids = np.arange(len(ctime))
    colors = np.where(region == "HEAD", RED, np.where(region == "TAIL", BLUE, BODY))
    fig, ax = plt.subplots(figsize=(10, 3))
    def draw(k):
        ax.clear(); pts = positions[int(frame_ids[k])]
        x, y = (pts - pipe_center).dot(c), (pts - pipe_center).dot(n)
        ax.axhline(wall_radius, color="#5996a5"); ax.axhline(-wall_radius, color="#5996a5")
        ax.scatter(x, y, c=colors, s=6, linewidths=0)
        tail, head = pts[region == "TAIL"].mean(0), pts[region == "HEAD"].mean(0)
        ax.text((tail-pipe_center).dot(c), (tail-pipe_center).dot(n), "TAIL", color=BLUE, weight="bold", ha="right")
        ax.text((head-pipe_center).dot(c), (head-pipe_center).dot(n), "HEAD", color=RED, weight="bold", ha="left")
        ax.set(xlim=(-3, 3), ylim=(-.9, .9), aspect="equal", xlabel="canonical +s (mm), left to right", ylabel="n_routeA (mm)")
        ax.set_title("{} | t={:.2f} ms | zeta=0.80 | strict CEL".format(JOB, ctime[k]*1e3), fontsize=8)
        ax.grid(axis="x", alpha=.15); fig.tight_layout()
    FuncAnimation(fig, draw, frames=len(frame_ids), interval=100).save(OUT / "TAIL_NO_REBOUND_8MS.gif", writer=PillowWriter(fps=10), dpi=100)
    plt.close(fig)

    failed = [key for key, value in gates.items() if not value]
    report = """# TAIL No-Rebound zeta=0.80, 8 ms TRUE-CEL Gate

Final classification: `{classification}`

The sole dynamics run completed successfully. TAIL contact was diagnosed from
the actual transformed TAIL mesh geometry, wall-normal gap velocity, and direct
robot-total General Contact fields. The nominal 0.005 mm plane-gap threshold
misses the first directly loaded event because the wall is faceted; the event
gate therefore requires direct contact force and a documented 0.010 mm
faceted-wall neighborhood. The ODB does not pair-isolate robot-wall from
robot-fluid contact.

- First touch: {touch} ms
- First release: {release} ms
- Second touch in the same half-cycle: {second} ms
- Same-half-cycle TAIL contacts: {count}
- Impact / rebound gap velocity: {vi} / {vr} mm/s
- R_v: {rv}
- Failed CEL numerical gates: {failed}

No FWD/REV interpretation is made.
""".format(classification=classification, touch=metrics["first_TAIL_touch_ms"],
           release=metrics["first_TAIL_release_ms"], second=metrics["second_TAIL_touch_ms"],
           count=metrics["N_tail_contacts_same_halfcycle"], vi=impact_v, vr=rebound_v,
           rv=rv, failed=", ".join(failed) if failed else "none")
    (OUT / "TAIL_NO_REBOUND_8MS_Report.md").write_text(report, encoding="ascii")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
