"""Summarize and render the single coarse TRUE-CEL fast full-transit screen."""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
JOB = "TRUECEL_B0P11_A14P5_FAST_FULL_TRANSIT"
CASE = OUT / "case" / JOB
PRIVATE = CASE / "private"
FREQUENCY_HZ = 120.0


def unit(values):
    vector = np.asarray(values, dtype=float)
    return vector / np.linalg.norm(vector)


def vec(data, prefix, time):
    return np.column_stack(
        [np.interp(time, data[f"{prefix}{axis}"][:, 0], data[f"{prefix}{axis}"][:, 1]) for axis in (1, 2, 3)]
    )


def runs(mask):
    edges = np.diff(np.r_[False, mask, False].astype(np.int8))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1) - 1
    return list(zip(starts, ends))


def longest_duration(time, mask):
    return max((float(time[j] - time[i]) for i, j in runs(mask)), default=0.0)


def maximum_backward_excursion(position):
    return float(np.max(np.maximum.accumulate(position) - position))


def part_nodes(deck, name):
    part = re.search(
        rf"^\*Part, name={re.escape(name)}\s*$([\s\S]*?)^\*End Part\s*$",
        deck,
        re.MULTILINE | re.IGNORECASE,
    )
    block = re.search(r"^\*Node\s*$([\s\S]*?)(?=^\*)", part.group(1), re.MULTILINE | re.IGNORECASE)
    rows = []
    for line in block.group(1).splitlines():
        fields = [item.strip() for item in line.split(",")]
        if len(fields) >= 4:
            rows.append([float(x) for x in fields[1:4]])
    return np.asarray(rows, dtype=float)


def read_event():
    path = CASE / "fast_screen_event.txt"
    if not path.exists():
        return {"reason": "NONE"}
    lines = path.read_text().splitlines()
    event = {"reason": lines[0].strip()}
    for line in lines[1:]:
        if "=" in line:
            key, value = line.split("=", 1)
            event[key.strip()] = float(value)
    return event


def classify(event, reached_finish, max_back, longest_negative, solver_ok):
    reason = event["reason"]
    if not solver_ok and reason == "NONE":
        return "NUMERICALLY_INVALID"
    if reason.startswith("NUMERICALLY_INVALID"):
        return "NUMERICALLY_INVALID"
    if reason.startswith("SUSTAINED_RECOIL_FAIL"):
        return "SUSTAINED_RECOIL_FAIL"
    if reason == "FORWARD_PROGRESS_STALL":
        return "FORWARD_PROGRESS_STALL"
    if reached_finish:
        if max_back <= 0.02 and longest_negative <= 0.00025:
            return "FULL_TRANSIT_FORWARD_PASS"
        return "FULL_TRANSIT_WITH_MINOR_RECOIL"
    return "FORWARD_PROGRESS_STALL"


def main():
    identity = json.loads((CASE / "case_identity.json").read_text())
    deck = (CASE / f"{JOB}.inp").read_text(encoding="latin1")
    c, n, b = (unit(identity[key]) for key in ("canonical_plus_s_axis_aba", "n_routeA_aba", "b_routeA_aba"))
    rp0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    s_start = float(identity["s_start_mm"])
    s_finish = float(identity["s_finish_mm"])

    rp = np.load(PRIVATE / "rp_history_private.npz")
    time = rp["U1"][:, 0].astype(float)
    U, UR, V = (vec(rp, prefix, time) for prefix in ("U", "UR", "V"))
    rp.close()
    displacement = U.dot(c)
    position = s_start + displacement
    v_s = V.dot(c)
    max_back = maximum_backward_excursion(position)
    longest_negative = longest_duration(time, v_s < 0.0)
    reached_finish = bool(np.max(position) >= s_finish - 1e-6)
    net = float(position[-1] - position[0])
    duration = float(time[-1] - time[0])
    mean_speed = net / duration if duration > 0 else 0.0

    period = 1.0 / FREQUENCY_HZ
    completed_cycles = int(np.floor((time[-1] + 1e-12) / period))
    cycles = []
    for cycle in range(1, completed_cycles + 1):
        t0, t1 = (cycle - 1) * period, cycle * period
        mask = (time >= t0) & (time <= t1)
        tc = np.r_[t0, time[mask], t1]
        sc = np.interp(tc, time, position)
        vc = np.interp(tc, time, v_s)
        delta = float(sc[-1] - sc[0])
        cycles.append({
            "cycle": cycle,
            "start_time_ms": t0 * 1e3,
            "end_time_ms": t1 * 1e3,
            "delta_s_mm": delta,
            "mean_v_s_mm_s": delta / period,
            "end_v_s_mm_s": float(vc[-1]),
            "minimum_v_s_mm_s": float(np.min(vc)),
            "maximum_backward_excursion_mm": maximum_backward_excursion(sc),
        })

    event = read_event()
    event_time = float(event.get("time_s", time[-1]))
    event_position = float(np.interp(event_time, time, position))
    event_position_2ms = float(np.interp(max(time[0], event_time - 0.002), time, position))
    rolling_2ms_gain = event_position - event_position_2ms
    sta_path = CASE / f"{JOB}.sta"
    msg_path = CASE / f"{JOB}.msg"
    sta = sta_path.read_text(encoding="latin1", errors="replace") if sta_path.exists() else ""
    msg = msg_path.read_text(encoding="latin1", errors="replace") if msg_path.exists() else ""
    solver_ok = "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in sta and "***ERROR" not in msg.upper()
    classification = classify(event, reached_finish, max_back, longest_negative, solver_ok)

    metrics = {
        "classification": classification,
        "reached_finish": "YES" if reached_finish else "NO",
        "travel_distance_mm": float(identity["travel_distance_mm"]),
        "achieved_net_displacement_mm": net,
        "maximum_s_mm": float(np.max(position)),
        "s_start_mm": s_start,
        "s_finish_mm": s_finish,
        "transit_time_ms": float(time[-1] * 1e3),
        "termination_time_ms": event_time * 1e3,
        "rolling_2ms_gain_at_termination_mm": rolling_2ms_gain,
        "mean_speed_mm_s": mean_speed,
        "minimum_v_s_mm_s": float(np.min(v_s)),
        "final_v_s_mm_s": float(v_s[-1]),
        "max_backward_excursion_mm": max_back,
        "longest_negative_velocity_duration_ms": longest_negative * 1e3,
        "full_rocking_cycles_completed": completed_cycles,
        "online_event": event,
        "solver_completed_successfully": solver_ok,
        "cycle_metrics": cycles,
        "screen_label": "REDUCED_SOUND_SPEED_COARSE_CEL_FAST_SCREEN",
    }
    (OUT / f"{JOB}_Metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="ascii")

    csv_path = OUT / f"{JOB}_Cycle_Summary.csv"
    with csv_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cycles[0]) if cycles else [
            "cycle", "start_time_ms", "end_time_ms", "delta_s_mm", "mean_v_s_mm_s",
            "end_v_s_mm_s", "minimum_v_s_mm_s", "maximum_backward_excursion_mm",
        ])
        writer.writeheader()
        writer.writerows(cycles)

    robot = part_nodes(deck, "Robot_SOLID")
    relative = robot - rp0
    pipe_center = rp0 - s_start * c - float(identity["radial_offset_n_mm"]) * n
    axial_ref = relative.dot(c)
    region = np.full(len(robot), "BODY", dtype=object)
    region[axial_ref <= axial_ref.min() + 0.25] = "TAIL"
    region[axial_ref >= axial_ref.max() - 0.25] = "HEAD"
    colors = np.where(region == "HEAD", "#d1493f", np.where(region == "TAIL", "#2864a8", "#3e434b"))
    frame_count = min(180, max(48, int(np.ceil(time[-1] * 120.0 * 12.0))))
    frame_times = np.linspace(time[0], time[-1], frame_count)
    fig, ax = plt.subplots(figsize=(11, 3.8), constrained_layout=True)
    initial_points = rp0 + U[0] + Rotation.from_rotvec(UR[0]).apply(relative)
    initial_x = s_start + (initial_points - rp0).dot(c)
    initial_y = (initial_points - pipe_center).dot(n)
    scatter = ax.scatter(initial_x, initial_y, s=5, c=colors, linewidths=0)
    status = ax.text(0.015, 0.965, "", transform=ax.transAxes, va="top", ha="left", family="monospace", fontsize=10)
    ax.axhline(identity["lumen_radius_mm"], color="#5d8794", lw=2)
    ax.axhline(-identity["lumen_radius_mm"], color="#5d8794", lw=2)
    ax.axvline(s_start, color="#777777", ls="--", lw=1, label="start")
    ax.axvline(s_finish, color="#2a9d55", ls="--", lw=1.5, label="finish")
    ax.set(
        xlim=tuple(identity["axial_CEL_extent_s_mm"]),
        ylim=(-0.9, 0.9),
        aspect="equal",
        xlabel="canonical +s (mm), LEFT -> RIGHT",
        ylabel="n_routeA (mm)",
        title="TRUE-CEL B0=11 mT, A_main=14.5 deg — FAST FULL TRANSIT",
    )
    ax.legend(loc="lower right", frameon=False)

    def update(frame):
        t = frame_times[frame]
        i = int(np.clip(np.searchsorted(time, t), 0, len(time) - 1))
        points = rp0 + U[i] + Rotation.from_rotvec(UR[i]).apply(relative)
        x = s_start + (points - rp0).dot(c)
        y = (points - pipe_center).dot(n)
        scatter.set_offsets(np.column_stack((x, y)))
        status.set_text(
            f"t = {time[i]*1e3:8.3f} ms    cycle = {time[i]*FREQUENCY_HZ:6.2f}\n"
            f"delta_s = {displacement[i]-displacement[0]:+8.4f} mm    v_s = {v_s[i]:+8.3f} mm/s"
        )
        return scatter, status

    animation = FuncAnimation(fig, update, frames=frame_count, interval=60, blit=False)
    gif_path = OUT / f"{JOB}.gif"
    animation.save(gif_path, writer=PillowWriter(fps=16), dpi=100)
    plt.close(fig)

    cycle_lines = [
        "| cycle | delta_s (mm) | mean v_s (mm/s) | end v_s | min v_s | max back (mm) |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    cycle_lines += [
        f"| {row['cycle']} | {row['delta_s_mm']:+.6f} | {row['mean_v_s_mm_s']:+.4f} | "
        f"{row['end_v_s_mm_s']:+.4f} | {row['minimum_v_s_mm_s']:+.4f} | "
        f"{row['maximum_backward_excursion_mm']:.6f} |"
        for row in cycles
    ]
    report = f"""# TRUE-CEL B0P11 A14P5 FAST full-transit screen

Classification: **{classification}**

- reached_finish: {'YES' if reached_finish else 'NO'}
- travel_distance_mm: {identity['travel_distance_mm']:.9f}
- achieved_net_displacement_mm: {net:.9f}
- transit_time_ms: {time[-1]*1e3:.6f}
- termination_time_ms: {event_time*1e3:.6f}
- rolling_2ms_gain_at_termination_mm: {rolling_2ms_gain:.9f}
- mean_speed_mm_s: {mean_speed:.6f}
- minimum_v_s_mm_s: {np.min(v_s):.6f}
- max_backward_excursion_mm: {max_back:.9f}
- longest_negative_velocity_duration_ms: {longest_negative*1e3:.6f}
- full_rocking_cycles_completed: {completed_cycles}
- online_event: {event['reason']}

This is a `REDUCED_SOUND_SPEED_COARSE_CEL_FAST_SCREEN`, not final quantitative FSI validation.

## Frozen parent configuration

- parent: `TRUECEL_NO_RECOIL_B0P11_2CYCLES`
- B0: 11 mT
- G: 2 mT
- frequency: 120 Hz
- A_main: 14.8 -> 14.5 deg (only physics change)
- A_cross: 2.5 deg
- c0: 100000 mm/s
- CEL mesh: 44 x 20 x 20
- contact: parent one-wall General Contact unchanged; mu=0.03; zeta=1.0
- Explicit: automatic stable increment, scale factor 0.4, no mass scaling
- s_start / s_finish: {s_start:.9f} / {s_finish:.9f} mm

## Cycle-by-cycle

{chr(10).join(cycle_lines)}
"""
    (OUT / f"{JOB}_Report.md").write_text(report, encoding="utf-8")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
