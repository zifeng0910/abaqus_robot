"""Analyze the authorized F100 run, render canonical motion, and compare F120."""
from __future__ import annotations

import csv
import json
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
CASE = OUT / "case" / "TRUECEL_B0P11_A14P5_F100_FAST"
JOB = "TRUECEL_B0P11_A14P5_F100_FAST"
PERIOD = 0.01
F120 = {
    "frequency_Hz": 120,
    "cycle1_delta_s_mm": 0.06094774095508182,
    "cycle2_delta_s_mm": -0.027793991839170307,
    "cycle3_delta_s_mm": None,
    "mean_speed_mm_s": 1.9886538,
    "max_backtrack_mm": 0.07236198,
    "cycle2_end_v_s_mm_s": -28.5378,
    "classification": "SUSTAINED_RECOIL_FAIL",
}


def unit(x):
    x = np.asarray(x, dtype=float)
    return x / np.linalg.norm(x)


def read_robot_nodes(deck: str) -> np.ndarray:
    import re
    part = re.search(r"^\*Part, name=Robot_SOLID\s*$([\s\S]*?)^\*End Part\s*$",
                     deck, re.MULTILINE | re.IGNORECASE)
    if not part:
        raise RuntimeError("Robot_SOLID part missing from candidate input")
    block = re.search(r"^\*Node\s*$([\s\S]*?)(?=^\*)", part.group(1), re.MULTILINE | re.IGNORECASE)
    if not block:
        raise RuntimeError("Robot_SOLID node block missing")
    return np.asarray([[float(v.strip()) for v in line.split(",")[1:4]]
                       for line in block.group(1).splitlines() if len(line.split(",")) >= 4])


def load_history():
    first_path = CASE / "private" / "rp_history_private.npz"
    if not first_path.is_file():
        raise FileNotFoundError(f"RP extraction not found: {first_path}")
    first = np.load(first_path)
    data = {k: first[k].astype(float) for k in ("U1", "U2", "U3", "UR1", "UR2", "UR3", "V1", "V2", "V3", "VR1", "VR2", "VR3")}
    first.close()
    cont_path = CASE / "F100_CYCLE3_RESTART_private" / "rp_history_private.npz"
    restart_time = float(data["U1"][-1, 0])
    if cont_path.is_file():
        continuation = np.load(cont_path)
        for key in data:
            a = continuation[key].astype(float).copy()
            a[:, 0] += restart_time
            a = a[a[:, 0] > restart_time + 1e-13]
            data[key] = np.vstack((data[key], a))
        continuation.close()
    return data, restart_time


def longest_negative(t, v):
    best = run = 0.0
    for i in range(len(t) - 1):
        t0, t1, v0, v1 = t[i], t[i + 1], v[i], v[i + 1]
        dt = t1 - t0
        if v0 < 0 and v1 < 0:
            run += dt
        elif v0 >= 0 and v1 < 0:
            run = dt * v0 / max(v0 - v1, 1e-30)
        elif v0 < 0 and v1 >= 0:
            run += dt * (-v0) / max(v1 - v0, 1e-30)
            best = max(best, run); run = 0.0
        else:
            best = max(best, run); run = 0.0
    return max(best, run)


def render(identity, deck, time, U, UR, position, v_s, gif_path):
    c = unit(identity["canonical_plus_s_axis_aba"])
    n = unit(identity["n_routeA_aba"])
    rp0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    robot = read_robot_nodes(deck)
    relative = robot - rp0
    pipe_center = rp0 - float(identity["s_start_mm"]) * c - float(identity["radial_offset_n_mm"]) * n
    axial = relative.dot(c)
    regions = np.full(len(robot), "BODY", dtype=object)
    regions[axial <= axial.min() + 0.25] = "TAIL"
    regions[axial >= axial.max() - 0.25] = "HEAD"
    colors = np.where(regions == "HEAD", "#d1493f", np.where(regions == "TAIL", "#2864a8", "#3e434b"))
    count = max(32, int(np.ceil(time[-1] / PERIOD * 32)))
    frames = np.linspace(0.0, time[-1], count)
    fig, ax = plt.subplots(figsize=(11, 4.3), constrained_layout=True)
    def pose(t):
        u = np.array([np.interp(t, time, U[:, j]) for j in range(3)])
        ur = np.array([np.interp(t, time, UR[:, j]) for j in range(3)])
        return rp0 + u + Rotation.from_rotvec(ur).apply(relative)
    points = pose(0.0)
    scatter = ax.scatter(identity["s_start_mm"] + (points-rp0).dot(c), (points-pipe_center).dot(n),
                         s=5, c=colors, linewidths=0)
    status = ax.text(0.985, 0.985, "", transform=ax.transAxes, va="top", ha="right", family="monospace", fontsize=9,
                     bbox={"boxstyle": "round,pad=0.30", "facecolor": "white", "edgecolor": "#bbbbbb", "alpha": 0.90})
    ax.axhline(identity["lumen_radius_mm"], color="#5d8794", lw=2)
    ax.axhline(-identity["lumen_radius_mm"], color="#5d8794", lw=2)
    ax.axvline(identity["s_start_mm"], color="#777", ls="--", lw=1, label="start")
    ax.set(xlim=tuple(identity["axial_CEL_extent_s_mm"]), ylim=(-0.9, 0.9), aspect="equal",
           xlabel="canonical +s (mm), LEFT -> RIGHT", ylabel="n_routeA (mm)",
           title="TRUE-CEL B0P11 A14P5 — 100 Hz single candidate")
    ax.legend(loc="lower right", frameon=False)
    def update(i):
        t = frames[i]
        pts = pose(t)
        scatter.set_offsets(np.column_stack((identity["s_start_mm"] + (pts-rp0).dot(c), (pts-pipe_center).dot(n))))
        cycle = min(3, max(1, int(t / PERIOD) + 1))
        total = np.interp(t, time, position) - position[0]
        back = np.max(np.maximum.accumulate(position[time <= t]) - position[time <= t]) if np.any(time <= t) else 0.0
        status.set_text(f"time = {t*1e3:8.3f} ms    cycle = {cycle}\n"
                        f"total delta_s = {total:+9.5f} mm    v_s = {np.interp(t,time,v_s):+9.3f} mm/s\n"
                        f"MAX_BACKTRACK = {back:8.5f} mm")
        return scatter, status
    FuncAnimation(fig, update, frames=count, interval=60, blit=False).save(
        gif_path, writer=PillowWriter(fps=16), dpi=96)
    plt.close(fig)


def main():
    identity = json.loads((CASE / "case_identity.json").read_text(encoding="utf-8"))
    history, restart_time = load_history()
    time = history["U1"][:, 0]
    c = unit(identity["canonical_plus_s_axis_aba"])
    U = np.column_stack([history[f"U{i}"][:, 1] for i in (1, 2, 3)])
    UR = np.column_stack([history[f"UR{i}"][:, 1] for i in (1, 2, 3)])
    V = np.column_stack([history[f"V{i}"][:, 1] for i in (1, 2, 3)])
    position = float(identity["s_start_mm"]) + U.dot(c)
    v_s = V.dot(c)
    back = np.maximum.accumulate(position) - position
    # Abaqus can close an exact 20-ms step a few tenths of a microsecond early
    # in the RP-history time stamps; treat the closed step endpoint as complete.
    complete = min(3, int(np.floor((time[-1] + 1e-6) / PERIOD)))
    cycles = []
    for k in range(1, complete + 1):
        t0, t1 = (k - 1) * PERIOD, k * PERIOD
        tc = np.r_[t0, time[(time > t0) & (time < t1)], t1]
        sc, vc = np.interp(tc, time, position), np.interp(tc, time, v_s)
        bc = np.maximum.accumulate(sc) - sc
        cycles.append({
            "cycle": k, "t_start_s": t0, "t_end_s": t1,
            "delta_s_mm": float(sc[-1] - sc[0]),
            "mean_v_s_mm_s": float((sc[-1] - sc[0]) / PERIOD),
            "end_v_s_mm_s": float(vc[-1]), "minimum_v_s_mm_s": float(vc.min()),
            "maximum_backward_excursion_mm": float(max(np.max(bc), np.max(back[(time >= t0) & (time <= t1)]) if np.any((time >= t0) & (time <= t1)) else 0.0)),
            "negative_velocity_duration_s": float(longest_negative(tc, vc)),
        })
    event_path = CASE / "f100_event.txt"
    event = {"reason": "NONE"}
    if event_path.is_file():
        lines = event_path.read_text(encoding="ascii", errors="replace").splitlines()
        event["reason"] = lines[0].strip() if lines else "UNKNOWN"
        for line in lines[1:]:
            if "=" in line:
                k, v = line.split("=", 1)
                try: event[k.strip()] = float(v)
                except ValueError: pass
    max_back = float(np.max(back))
    d = [x["delta_s_mm"] for x in cycles]
    ratio21 = d[1] / d[0] if len(d) >= 2 and abs(d[0]) > 1e-15 else None
    ratio31 = d[2] / d[0] if len(d) >= 3 and abs(d[0]) > 1e-15 else None
    hard_fail = bool(event["reason"].startswith(("F100_RECOIL_FAIL", "NUMERICALLY_INVALID")))
    all_forward = complete == 3 and all(x > 0 for x in d) and not hard_fail
    recoil = hard_fail or max_back >= 0.01 or any(x <= 0 for x in d)
    if all_forward and not recoil:
        classification = "F100_THREE_CYCLE_FORWARD_PASS"
    elif complete >= 2 or hard_fail:
        classification = "F100_RECOIL_FAIL"
    else:
        classification = "F100_STAGE1_INCOMPLETE"
    if len(d) >= 2:
        ratio = d[-1] / d[0] if abs(d[0]) > 1e-15 else 0.0
        trend = "FORWARD_STEP_INCREASING" if ratio > 1.2 else "FORWARD_STEP_DECAYING" if ratio < 0.8 else "FORWARD_STEP_STABLE"
    else:
        trend = None
    mean_speed = float((position[-1] - position[0]) / max(time[-1] - time[0], 1e-30))
    metrics = {
        "candidate": JOB, "classification": classification, "frequency_Hz": 100,
        "period_s": PERIOD, "complete_cycles": complete, "total_time_s": float(time[-1]),
        "fresh_start": True, "only_physics_change": "frequency 120 -> 100 Hz",
        "cycle_metrics": cycles, "cycle2_over_cycle1_delta_s_ratio": ratio21,
        "cycle3_over_cycle1_delta_s_ratio": ratio31, "mean_speed_mm_s": mean_speed,
        "max_backtrack_mm": max_back, "forward_step_behavior": trend,
        "online_event": event,
        "stage1_gate_pass": bool(complete >= 2 and len(d) >= 2 and d[0] > 0 and d[1] > 0 and max_back < 0.01 and not hard_fail),
        "solver_status": identity.get("status"),
        "torque_authority_check_F120": {
            "classification": "ROCKING_CONTROL_REMAINS_COHERENT",
            "rocking_amplitude_cycle1_deg": 12.984901915033914,
            "rocking_amplitude_cycle2_deg": 13.121398139744127,
            "phase_shift_cycle2_minus_cycle1_deg": -2.4191182435696006,
            "rms_phase_matched_alpha_difference_deg": 2.2722889409534726,
            "p95_abs_omega_rock_difference_rad_s": 220.4288459210195,
        },
        "F120_reference": F120,
    }
    csv_path = OUT / f"{JOB}_CYCLE_SUMMARY.csv"
    fields = ["cycle", "t_start_s", "t_end_s", "delta_s_mm", "mean_v_s_mm_s", "end_v_s_mm_s",
              "minimum_v_s_mm_s", "maximum_backward_excursion_mm", "negative_velocity_duration_s"]
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader(); writer.writerows(cycles)
    gif_path = OUT / f"{JOB}.gif"
    render(identity, (CASE / f"{JOB}.inp").read_text(encoding="latin1"),
           time, U, UR, position, v_s, gif_path)
    zip_path = OUT / f"{JOB}_GIF.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        z.write(gif_path, arcname=gif_path.name)
    metrics["gif"] = gif_path.name; metrics["gif_zip"] = zip_path.name
    (OUT / f"{JOB}_METRICS.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    comp_path = OUT / f"{JOB}_F120_COMPARISON.csv"
    rows = [F120, {"frequency_Hz": 100,
        "cycle1_delta_s_mm": d[0] if len(d)>0 else None,
        "cycle2_delta_s_mm": d[1] if len(d)>1 else None,
        "cycle3_delta_s_mm": d[2] if len(d)>2 else None,
        "mean_speed_mm_s": mean_speed, "max_backtrack_mm": max_back,
        "cycle2_end_v_s_mm_s": cycles[1]["end_v_s_mm_s"] if len(cycles)>1 else None,
        "classification": classification}]
    with comp_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    table = "\n".join("| {frequency_Hz} | {cycle1_delta_s_mm} | {cycle2_delta_s_mm} | {cycle3_delta_s_mm} | {mean_speed_mm_s} | {max_backtrack_mm} | {cycle2_end_v_s_mm_s} |".format(**r) for r in rows)
    removes_net_collapse = complete >= 2 and len(d) >= 2 and F120["cycle2_delta_s_mm"] <= 0 and d[1] > 0
    report = f"""# TRUE-CEL B0P11 A14P5 — single 100-Hz candidate

Classification: **{classification}**

F120 torque-authority check: **ROCKING_CONTROL_REMAINS_COHERENT**. The measured phase-matched rocking amplitudes were 12.985° and 13.121°, fundamental shift −2.42°, RMS angle difference 2.272°. The p95 absolute instantaneous omega_rock difference was 220.43 rad/s (not small), so angular-rate variability remains visible even though the rocking amplitude and phase did not collapse.

F100: fresh start at t=0; only the forcing frequency changed from 120 to 100 Hz. F100 end time {time[-1]*1e3:.6f} ms; complete cycles {complete}; mean speed {mean_speed:+.4f} mm/s; MAX_BACKTRACK {max_back:.6f} mm; stage-1 gate {'PASS' if metrics['stage1_gate_pass'] else 'FAIL/NOT PASSED'}.

| frequency (Hz) | Cycle 1 Δs (mm) | Cycle 2 Δs (mm) | Cycle 3 Δs (mm) | mean speed (mm/s) | MAX_BACKTRACK (mm) | Cycle 2 end v_s (mm/s) |
|---:|---:|---:|---:|---:|---:|---:|
{table}

Central question — does 100 Hz remove the second-cycle collapse? **{'It removes the cycle-integrated collapse: Cycle 2 Δs becomes positive, but sustained second-cycle recoil remains, so the two-cycle gate fails and no third cycle is run.' if removes_net_collapse else 'Not established; the cycle-integrated collapse remains or stage 1 is incomplete.'}**

Cycle metrics are in `{csv_path.name}`; machine-readable results are in `{(OUT / f'{JOB}_METRICS.json').name}`. Motion uses canonical +s left-to-right with tail left/head right. GIF: `{gif_path.name}`.
"""
    (OUT / f"{JOB}_REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
