"""Analyze and render the accepted refined-geometry 10 ms strict-CEL gate."""
from __future__ import annotations

import json
import math
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
JOB = "REFINED_DUALEND_F60_A14P8_G0P15_TRUECEL_10MS"
CASE = OUT / "case" / JOB
PRIVATE = CASE / "private"
sys.path.insert(0, str(REPO / "calibration_analysis" / "S4HeadForwardLowGScreen" / "scripts"))
import analyze_and_render_lowg as base

BLUE, RED, BODY = "#2864a8", "#d1493f", "#383c42"


def interp_history(archive, prefix, time):
    return np.column_stack([np.interp(time, archive[prefix + str(i)][:, 0], archive[prefix + str(i)][:, 1])
                            for i in (1, 2, 3)])


def longest(time, active):
    best = start = None
    for i, value in enumerate(active):
        if value and start is None:
            start = i
        if start is not None and (not value or i == len(active) - 1):
            end = i if value and i == len(active) - 1 else i - 1
            duration = time[end] - time[start]
            best = duration if best is None else max(best, duration)
            start = None
    return 0.0 if best is None else float(best)


def episodes(time, active, force, gap, region):
    rows, start = [], None
    for i, value in enumerate(active):
        if value and start is None:
            start = i
        if start is not None and (not value or i == len(active) - 1):
            end = i if value and i == len(active) - 1 else i - 1
            rows.append({"region": region, "start_s": time[start], "end_s": time[end],
                         "duration_ms": (time[end] - time[start]) * 1e3,
                         "max_total_general_contact_resultant_N": float(np.linalg.norm(force[start:end + 1], axis=1).max()),
                         "min_geometric_gap_mm": float(np.min(gap[start:end + 1])),
                         "reopened": bool(end < len(active) - 1 and not active[end + 1])})
            start = None
    return rows


def main():
    identity = json.loads((CASE / "case_identity.json").read_text())
    deck = (CASE / (JOB + ".inp")).read_text(encoding="latin1")
    c = base.unit(identity["canonical_plus_s_axis_aba"])
    n = base.unit(identity["n_routeA_aba"])
    b = base.unit(identity["b_routeA_aba"])
    rp0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    pipe_center = rp0 - identity["radial_offset_n_mm"] * n
    nodes = base.part_nodes(deck, "Robot_SOLID")
    wall = base.part_nodes(deck, "Pipe_WALL_HELPER")
    normals, wall_radius = base.wall_planes(wall, pipe_center, c, n, b)
    s_ref = (nodes - rp0).dot(c)
    regions = np.full(len(nodes), "BODY", dtype=object)
    regions[s_ref <= s_ref.min() + 0.25] = "TAIL"
    regions[s_ref >= s_ref.max() - 0.25] = "HEAD"

    rp_archive = np.load(PRIVATE / "rp_history_private.npz")
    time = rp_archive["U1"][:, 0]
    U = interp_history(rp_archive, "U", time)
    UR = interp_history(rp_archive, "UR", time)
    V = interp_history(rp_archive, "V", time)
    rp_archive.close()
    axial = U.dot(c)
    body_axis = Rotation.from_rotvec(UR).apply(np.broadcast_to(c, UR.shape))
    rocking = np.degrees(np.arctan2(body_axis.dot(n), body_axis.dot(c)))

    contact = np.load(PRIVATE / "robot_contact_fields_private.npz")
    ctime = contact["time"]
    labels = contact["node_labels"]
    if not np.array_equal(labels, np.arange(1, len(nodes) + 1)):
        raise RuntimeError("robot node ordering changed")
    cnorm = contact["CNORMF"].astype(float)
    cshear = contact["CSHEARF"].astype(float)
    total_node_force = cnorm + cshear
    force = {}
    gaps = {name: [] for name in ("HEAD", "TAIL", "BODY")}
    positions = []
    for i, frame_time in enumerate(ctime):
        j = int(np.clip(np.searchsorted(time, frame_time), 0, len(time) - 1))
        points = rp0 + U[j] + Rotation.from_rotvec(UR[j]).apply(nodes - rp0)
        positions.append(points)
        rel = points - pipe_center
        transverse = rel - np.outer(rel.dot(c), c)
        node_gap = wall_radius - transverse.dot(normals.T).max(axis=1)
        for name in gaps:
            gaps[name].append(float(node_gap[regions == name].min()))
    positions = np.asarray(positions)
    gaps = {name: np.asarray(values) for name, values in gaps.items()}
    for name in ("HEAD", "TAIL", "BODY"):
        force[name] = total_node_force[:, regions == name, :].sum(axis=1)
    force["TOTAL"] = total_node_force.sum(axis=1)
    contact.close()
    threshold = 1e-8
    active = {name: np.linalg.norm(force[name], axis=1) > threshold for name in ("HEAD", "TAIL", "BODY")}
    wall_proximity = {name: gaps[name] <= 0.005 for name in ("HEAD", "TAIL", "BODY")}
    both = wall_proximity["HEAD"] & wall_proximity["TAIL"]
    episode_rows = []
    for name in ("HEAD", "TAIL", "BODY"):
        episode_rows += episodes(ctime, active[name], force[name], gaps[name], name)
    pd.DataFrame(episode_rows).to_csv(CASE / (JOB + "_contact_episodes.csv"), index=False)
    pd.DataFrame({"time_s": ctime, "TAIL_gap_mm": gaps["TAIL"],
                  "TAIL_wall_proximity": wall_proximity["TAIL"].astype(int),
                  "TAIL_total_contact_force_N": np.linalg.norm(force["TAIL"], axis=1),
                  "TAIL_total_contact_axial_N": force["TAIL"].dot(c),
                  "rocking_deg": np.interp(ctime, time, rocking)}).to_csv(
                      CASE / (JOB + "_TAIL_contact_diagnosis.csv"), index=False)

    energy = np.load(PRIVATE / "energy_history_private.npz")
    etotal = energy["ETOTAL"]
    energy.close()
    energy_drift = etotal[:, 1] - etotal[0, 1]

    field = np.load(PRIVATE / "truecel_field_private.npz")
    ftime = field["time"]
    fcoords = field["fluid_node_coordinates_mm"]
    fvelocity = field["fluid_velocity_mm_s"].astype(float)
    fevf = field["fluid_evf"].astype(float)
    fs = (fcoords - pipe_center).dot(c)
    frvec = fcoords - pipe_center - np.outer(fs, c)
    near = (np.abs(fs) <= 2.0) & (np.linalg.norm(frvec, axis=1) <= wall_radius + 0.10)
    near_speed = np.linalg.norm(fvelocity[:, near], axis=2)
    max_fluid_speed = float(np.nanmax(near_speed))
    max_frame = int(np.nanargmax(np.nanmax(near_speed, axis=1)))
    evf_present = bool(np.any(fevf[0] > 0.0) and np.any(fevf[-1] > 0.0))

    telemetry = pd.read_csv(CASE / (JOB + "_telemetry.csv")).drop_duplicates("t_s").sort_values("t_s")
    fmag = telemetry[["fx_aba_N", "fy_aba_N", "fz_aba_N"]].to_numpy(float).dot(c)
    jmag = float(np.trapz(fmag, telemetry.t_s.to_numpy(float)))
    total_contact_s = force["TOTAL"].dot(c)
    jcontact = float(np.trapz(total_contact_s, ctime))

    sta = (CASE / (JOB + ".sta")).read_text(encoding="latin1")
    rows = re.findall(r"(?m)^\s*(\d+)\s+([0-9.E+-]+)\s+[0-9.E+-]+\s+\S+\s+([0-9.E+-]+)\s+\d+", sta)
    dts = np.asarray([float(row[2]) for row in rows])
    increments = max(int(row[0]) for row in rows)
    socket_continuous = len(telemetry) > 10 and telemetry.t_s.max() >= 0.0099
    max_contact = float(np.linalg.norm(force["TOTAL"], axis=1).max())
    min_gap = min(float(values.min()) for values in gaps.values())
    tumble = bool(np.any(body_axis.dot(c) < 0.0))
    persistent_bridge = longest(ctime, both) > 0.001
    gates = {
        "no_tumble": not tumble,
        "no_gross_penetration": min_gap > -0.05,
        "no_persistent_both_bridge": not persistent_bridge,
        "energy_drift_below_0p01_Nmm": float(np.max(np.abs(energy_drift))) < 0.01,
        "near_robot_fluid_velocity_below_10000_mm_s": max_fluid_speed < 10000.0,
        "robot_total_general_contact_resultant_below_0p1_N": max_contact < 0.1,
        "stable_dt_no_collapse": float(dts.min()) > 1e-8,
        "socket_continuous": socket_continuous,
        "EVF_and_V_present": evf_present and np.isfinite(fvelocity).any(),
    }
    gates = {name: bool(value) for name, value in gates.items()}
    passed = all(gates.values())

    metrics = {
        "case": JOB, "classification": "STRICT_CEL_FLUID_PRESENT", "simulated_time_ms": float(time[-1] * 1e3),
        "geometry": {"L_mm": 2.6, "D_mm": 0.815, "L_over_D": 2.6 / 0.815,
                     "radial_offset_n_mm": 0.015, "HEAD_side_correction_mm": 0.0425,
                     "HEAD_touch_deg": 14.387837515678257, "TAIL_touch_deg": 13.423312820959836},
        "rocking_command_deg": 14.8, "rocking_min_deg": float(rocking.min()), "rocking_max_deg": float(rocking.max()),
        "delta_s_mm": float(axial[-1] - axial[0]), "max_radial_COM_displacement_mm": float(np.linalg.norm(U - np.outer(axial, c), axis=1).max()),
        "energy_initial_Nmm": float(etotal[0, 1]), "energy_final_Nmm": float(etotal[-1, 1]),
        "max_abs_energy_drift_Nmm": float(np.max(np.abs(energy_drift))),
        "near_robot_fluid_velocity_max_mm_s": max_fluid_speed,
        "robot_total_general_contact_resultant_max_N": max_contact,
        "wall_only_force_directly_available": False,
        "contact_force_limitation": "Abaqus ODB exposes robot total General Contact; robot-wall and robot-fluid are not pair-isolated.",
        "HEAD_min_gap_mm": float(gaps["HEAD"].min()), "TAIL_min_gap_mm": float(gaps["TAIL"].min()),
        "BODY_min_gap_mm": float(gaps["BODY"].min()),
        "HEAD_total_general_contact_episodes": int(sum(row["region"] == "HEAD" for row in episode_rows)),
        "TAIL_total_general_contact_episodes": int(sum(row["region"] == "TAIL" for row in episode_rows)),
        "BODY_total_general_contact_episodes": int(sum(row["region"] == "BODY" for row in episode_rows)),
        "HEAD_reopened": bool(np.any(active["HEAD"][:-1] & ~active["HEAD"][1:])),
        "TAIL_reopened": bool(np.any(active["TAIL"][:-1] & ~active["TAIL"][1:])),
        "BODY_reopened": bool(np.any(active["BODY"][:-1] & ~active["BODY"][1:])),
        "both_bridge_longest_ms": longest(ctime, both) * 1e3,
        "J_robot_total_contact_s_Ns": jcontact, "Jmag_s_Ns": jmag,
        "minimum_stable_dt_s": float(dts.min()), "total_increments": increments,
        "socket_last_time_ms": float(telemetry.t_s.max() * 1e3), "gates": gates,
        "gate_passed": passed, "final_status": "TRUECEL_GATE_PASSED" if passed else "TRUECEL_GATE_FAILED",
    }
    (OUT / (JOB + "_metrics.json")).write_text(json.dumps(metrics, indent=2) + "\n", encoding="ascii")

    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    ax.plot(etotal[:, 0] * 1e3, energy_drift, color="#2f5d62", lw=1.5)
    ax.axhline(0.01, color=RED, ls="--"); ax.axhline(-0.01, color=RED, ls="--")
    ax.set(xlabel="time (ms)", ylabel="ETOTAL drift (N mm)", title="10 ms TRUE-CEL energy gate")
    ax.grid(alpha=0.2); fig.tight_layout(); fig.savefig(OUT / (JOB + "_energy.png"), dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    for name, color in (("HEAD", RED), ("TAIL", BLUE), ("BODY", BODY)):
        ax.plot(ctime * 1e3, np.linalg.norm(force[name], axis=1), marker="o", ms=2.5, label=name, color=color)
    ax.axhline(0.1, color=RED, ls="--", label="0.1 N gate")
    ax.set(xlabel="time (ms)", ylabel="direct total General Contact resultant (N)", title="Region contact-force gate")
    ax.grid(alpha=0.2); ax.legend(); fig.tight_layout(); fig.savefig(OUT / (JOB + "_contact_force.png"), dpi=180); plt.close(fig)

    speed = np.linalg.norm(fvelocity[max_frame], axis=1)
    local_s = fs; local_n = (fcoords - pipe_center).dot(n)
    fig, ax = plt.subplots(figsize=(10.0, 3.2))
    sc = ax.scatter(local_s[near], local_n[near], c=speed[near], s=9, cmap="viridis")
    j = int(np.clip(np.searchsorted(ctime, ftime[max_frame]), 0, len(ctime) - 1))
    pts = positions[j]
    ax.scatter((pts - pipe_center).dot(c), (pts - pipe_center).dot(n),
               c=np.where(regions == "HEAD", RED, np.where(regions == "TAIL", BLUE, BODY)), s=5)
    ax.axhline(wall_radius, color="#5996a5"); ax.axhline(-wall_radius, color="#5996a5")
    ax.set(xlim=(-3, 3), ylim=(-0.9, 0.9), xlabel="canonical +s (mm)", ylabel="n_routeA (mm)",
           title="TRUE-CEL velocity snapshot at t={:.3f} ms".format(ftime[max_frame] * 1e3))
    fig.colorbar(sc, ax=ax, label="fluid nodal speed (mm/s)"); fig.tight_layout()
    fig.savefig(OUT / (JOB + "_CEL_velocity_snapshot.png"), dpi=180); plt.close(fig)

    frame_ids = np.arange(0, len(ctime), 2)
    fig, ax = plt.subplots(figsize=(10.0, 3.0))
    colors = np.where(regions == "HEAD", RED, np.where(regions == "TAIL", BLUE, BODY))
    def draw(k):
        ax.clear(); i = int(frame_ids[k]); pts = positions[i]
        x = (pts - pipe_center).dot(c); y = (pts - pipe_center).dot(n)
        ax.axhline(wall_radius, color="#5996a5"); ax.axhline(-wall_radius, color="#5996a5")
        ax.scatter(x, y, c=colors, s=6, linewidths=0)
        tail = pts[regions == "TAIL"].mean(axis=0); head = pts[regions == "HEAD"].mean(axis=0)
        ax.text((tail-pipe_center).dot(c), (tail-pipe_center).dot(n), "TAIL", color=BLUE, weight="bold", ha="right")
        ax.text((head-pipe_center).dot(c), (head-pipe_center).dot(n), "HEAD", color=RED, weight="bold", ha="left")
        ax.set(xlim=(-3, 3), ylim=(-0.9, 0.9), aspect="equal", xlabel="canonical +s (mm), left to right", ylabel="n_routeA (mm)")
        ax.set_title("{} | t={:.2f} ms | f=60 Hz | A=14.8 deg | G=0.15 mT | U0=+10 mm/s".format(JOB, ctime[i]*1e3), fontsize=8)
        ax.grid(axis="x", alpha=0.15); fig.tight_layout()
    anim = FuncAnimation(fig, draw, frames=len(frame_ids), interval=100)
    anim.save(OUT / (JOB + "_2D.gif"), writer=PillowWriter(fps=10), dpi=100); plt.close(fig)

    failed = [name for name, value in gates.items() if not value]
    report = """# Refined Dual-End 10 ms TRUE-CEL Gate

This was one strict-CEL dynamics run. It completed 10.0 ms with real water EVF,
initial `+10 mm/s` flow, Magpylib socket loading, General Contact, and no
ReducedHydro load.

The accepted physical mesh change is a smooth one-sided 0.0425 mm HEAD profile
increase over the final 25% of the body, plus a +0.015 mm radial pose offset.
It is an actual node-coordinate geometry change, not a label or measurement hack.

## Gate result

Final status: `{status}`

Failed hard gates: `{failed}`.

| quantity | value |
|---|---:|
| ETOTAL max absolute drift | {energy:.9g} N mm |
| near-robot fluid velocity max | {fluid:.9g} mm/s |
| robot total General Contact resultant max | {contact:.9g} N |
| minimum stable dt | {dt:.9g} s |
| rocking min / max | {rmin:.6f} / {rmax:.6f} deg |
| HEAD / TAIL / BODY minimum gap | {hg:.6f} / {tg:.6f} / {bg:.6f} mm |
| HEAD / TAIL / BODY total-General-Contact episodes | {he} / {te} / {be} |
| longest simultaneous HEAD+TAIL wall proximity | {bridge:.6f} ms |
| delta_s | {ds:+.9g} mm |
| Jmag_s | {jmag:+.9g} N s |
| J robot-total-contact_s | {jc:+.9g} N s |

Abaqus directly provides total robot General Contact fields (`CNORMF`,
`CSHEARF`, `COPEN`). It does not pair-isolate robot-wall from robot-fluid in
this ODB, so no momentum residual or geometry-filtered proxy is labeled as a
direct wall-only force. The energy hard gate failed at the endpoint; therefore
no propulsion interpretation and no FWD/REV pair are permitted.
""".format(status=metrics["final_status"], failed=", ".join(failed), energy=metrics["max_abs_energy_drift_Nmm"],
           fluid=max_fluid_speed, contact=max_contact, dt=metrics["minimum_stable_dt_s"],
           rmin=metrics["rocking_min_deg"], rmax=metrics["rocking_max_deg"],
           hg=metrics["HEAD_min_gap_mm"], tg=metrics["TAIL_min_gap_mm"], bg=metrics["BODY_min_gap_mm"],
           he=metrics["HEAD_total_general_contact_episodes"], te=metrics["TAIL_total_general_contact_episodes"], be=metrics["BODY_total_general_contact_episodes"],
           bridge=metrics["both_bridge_longest_ms"], ds=metrics["delta_s_mm"], jmag=jmag, jc=jcontact)
    (OUT / (JOB + "_Gate_Report.md")).write_text(report, encoding="ascii")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
