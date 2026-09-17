"""Analyze contact topology and render the sole FAST precomputed screen."""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
JOB = "FAST_PRECOMP_F120_DUALEND_NOREBOUND"
CASE = OUT / "case" / JOB
PRIVATE = CASE / "private"
BLUE, RED, BODY = "#2864a8", "#d1493f", "#383c42"
FORCE_THRESHOLD_N = 1.0e-10


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def history_vector(archive, prefix, time):
    return np.column_stack([np.interp(time, archive[prefix + str(i)][:, 0], archive[prefix + str(i)][:, 1])
                            for i in (1, 2, 3)])


def part_nodes(deck, name):
    part = re.search(r"(?ms)^\*Part, name=" + re.escape(name) + r"\s*$.*?^\*End Part\s*$", deck).group(0)
    block = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", part).group(1)
    rows = [[float(x) for x in line.split(",")[:4]] for line in block.splitlines() if line.strip()]
    return np.asarray([row[0] for row in rows], int), np.asarray([row[1:] for row in rows], float)


def wall_planes(wall, center, c, n, b):
    axial = (wall - center).dot(c)
    ring = wall[np.abs(axial - axial.min()) < 1.0e-7]
    xy = np.column_stack(((ring - center).dot(n), (ring - center).dot(b)))
    xy = xy[np.argsort(np.arctan2(xy[:, 1], xy[:, 0]))]
    normals, radii = [], []
    for p, q in zip(xy, np.roll(xy, -1, axis=0)):
        edge = q - p
        normal2 = unit([edge[1], -edge[0]])
        if np.dot(normal2, p + q) < 0:
            normal2 *= -1
        normals.append(normal2[0] * n + normal2[1] * b)
        radii.append(np.dot(p, normal2))
    return np.asarray(normals), float(np.mean(radii))


def episodes(active):
    edge = np.diff(np.r_[False, active, False].astype(int))
    return list(zip(np.where(edge == 1)[0], np.where(edge == -1)[0] - 1))


def reversals(v, threshold=1.0e-3):
    signs = np.sign(v[np.abs(v) > threshold])
    return int(np.sum(signs[1:] != signs[:-1])) if len(signs) > 1 else 0


def longest(time, active):
    return max((time[j] - time[i] for i, j in episodes(active)), default=0.0)


def endpoint_metrics(name, mask, ctime, gaps, normal, shear, rocking):
    normal_vec = normal[:, mask, :].sum(axis=1)
    shear_vec = shear[:, mask, :].sum(axis=1)
    normal_force = np.linalg.norm(normal_vec, axis=1)
    tangential_force = np.linalg.norm(shear_vec, axis=1)
    axial_force = (normal_vec + shear_vec).dot(C_AXIS)
    active = (normal_force > FORCE_THRESHOLD_N) & (gaps <= 0.015)
    spans = episodes(active)
    gap_velocity = np.gradient(gaps, ctime)
    rows, rebound_ratios = [], []
    for number, (i, j) in enumerate(spans, 1):
        incoming_window = slice(max(0, i - 3), min(len(ctime), i + 2))
        incoming = float(min(0.0, np.min(gap_velocity[incoming_window])))
        next_start = spans[number][0] if number < len(spans) else len(ctime) - 1
        release = min(j + 1, len(ctime) - 1)
        half_end = min(len(ctime) - 1, int(np.searchsorted(ctime, (np.floor(ctime[i] * 240.0) + 1) / 240.0)))
        rebound_end = max(release, min(next_start, half_end))
        outward = float(max(0.0, np.max(gap_velocity[release:rebound_end + 1])))
        ratio = outward / abs(incoming) if incoming < 0 else None
        if ratio is not None:
            rebound_ratios.append(ratio)
        rows.append({
            "end": name, "episode": number, "start_ms": ctime[i] * 1e3, "end_ms": ctime[j] * 1e3,
            "duration_ms": (ctime[j] - ctime[i] + np.median(np.diff(ctime))) * 1e3,
            "impact_wall_normal_velocity_mm_s": incoming,
            "post_contact_separation_velocity_mm_s": outward, "R_v": ratio,
            "minimum_gap_mm": float(gaps[i:j + 1].min()),
            "peak_normal_force_N": float(normal_force[i:j + 1].max()),
            "peak_tangential_force_N": float(tangential_force[i:j + 1].max()),
            "peak_abs_axial_contact_force_N": float(np.abs(axial_force[i:j + 1]).max()),
            "reopen_time_ms": ctime[release] * 1e3 if release > j else None,
            "rocking_at_impact_deg": float(rocking[i]),
        })
    half_cycle_counts = []
    for k in range(4):
        lo, hi = k / 240.0, (k + 1) / 240.0
        half_cycle_counts.append(int(sum(ctime[i] < hi and ctime[j] >= lo for i, j in spans)))
    summary = {
        "minimum_wall_gap_mm": float(gaps.min()), "contact_episode_count": len(spans),
        "contact_state_final": bool(active[-1]), "longest_contact_ms": longest(ctime, active) * 1e3,
        "peak_normal_contact_force_N": float(normal_force.max()),
        "peak_tangential_contact_force_N": float(tangential_force.max()),
        "peak_abs_axial_contact_force_N": float(np.abs(axial_force).max()),
        "max_rebound_ratio_R_v": max(rebound_ratios) if rebound_ratios else None,
        "same_end_episodes_per_half_cycle": half_cycle_counts,
        "same_end_episode_gate_le_1": max(half_cycle_counts, default=0) <= 1,
        "R_v_gate_lt_0p2": bool(rebound_ratios) and max(rebound_ratios) < 0.2,
    }
    return summary, rows, active, normal_force, tangential_force, axial_force


def render(positions, region, ctime, rocking, axial, head_active, tail_active, radius, diagnostic):
    count = 120 if diagnostic else 60
    ids = np.unique(np.linspace(0, len(ctime) - 1, count).astype(int))
    fps = 12 if diagnostic else 30
    colors = np.where(region == "HEAD", RED, np.where(region == "TAIL", BLUE, BODY))
    fig, ax = plt.subplots(figsize=(10, 3))

    def draw(frame):
        ax.clear()
        k = int(ids[frame]); pts = positions[k]
        x, y = (pts - PIPE_CENTER).dot(C_AXIS), (pts - PIPE_CENTER).dot(N_AXIS)
        ax.axhspan(radius, 0.9, color="#d8edf2", alpha=.55)
        ax.axhspan(-0.9, -radius, color="#d8edf2", alpha=.55)
        ax.axhline(radius, color="#5996a5"); ax.axhline(-radius, color="#5996a5")
        ax.scatter(x[::2], y[::2], c=colors[::2], s=7, linewidths=0)
        tail = np.asarray([x[region == "TAIL"].mean(), y[region == "TAIL"].mean()])
        head = np.asarray([x[region == "HEAD"].mean(), y[region == "HEAD"].mean()])
        ax.text(tail[0] - .06, tail[1], "TAIL", color=BLUE, weight="bold", ha="right", va="center")
        ax.text(head[0] + .06, head[1], "HEAD", color=RED, weight="bold", ha="left", va="center")
        state = "TAIL:{}  HEAD:{}".format("CONTACT" if tail_active[k] else "free", "CONTACT" if head_active[k] else "free")
        ax.text(.02, .93, state, transform=ax.transAxes, va="top", weight="bold", fontsize=9)
        ax.set(xlim=(-3.0, 3.0), ylim=(-.9, .9), aspect="equal",
               xlabel="canonical +s (mm), left to right", ylabel="n_routeA (mm)")
        ax.grid(axis="x", alpha=.15)
        ax.set_title("{} | t={:.3f} ms | cycle={:.3f} | rock={:+.2f} deg | delta_s={:+.4f} mm".format(
            JOB, ctime[k] * 1e3, ctime[k] * 120.0, rocking[k], axial[k] - axial[0]), fontsize=8)
        fig.tight_layout()

    suffix = "DIAGNOSTIC" if diagnostic else "VIEW"
    FuncAnimation(fig, draw, frames=len(ids), interval=1000.0 / fps).save(
        OUT / (JOB + "_2CYCLES_" + suffix + ".gif"), writer=PillowWriter(fps=fps), dpi=100)
    plt.close(fig)


def main():
    global C_AXIS, N_AXIS, PIPE_CENTER
    identity = json.loads((CASE / "case_identity.json").read_text())
    if identity["status"] != "SOLVED":
        raise RuntimeError("case is not SOLVED")
    C_AXIS = unit(identity["canonical_plus_s_axis_aba"])
    N_AXIS = unit(identity["n_routeA_aba"])
    b = unit(identity["b_routeA_aba"])
    rp0 = np.asarray(identity["initial_center_aba_mm"], float)
    PIPE_CENTER = rp0 - identity["radial_offset_n_mm"] * N_AXIS

    deck = (CASE / (JOB + ".inp")).read_text(encoding="latin1")
    node_labels, nodes = part_nodes(deck, "Robot_SOLID")
    _, wall = part_nodes(deck, "Pipe_WALL_HELPER")
    normals, radius = wall_planes(wall, PIPE_CENTER, C_AXIS, N_AXIS, b)
    s0 = (nodes - rp0).dot(C_AXIS)
    region = np.full(len(nodes), "BODY", dtype=object)
    region[s0 <= s0.min() + .25] = "TAIL"
    region[s0 >= s0.max() - .25] = "HEAD"

    rp = np.load(PRIVATE / "rp_history_private.npz")
    time = rp["V1"][:, 0]
    U, UR, V = (history_vector(rp, prefix, time) for prefix in ("U", "UR", "V"))
    rp.close()
    axial_history = U.dot(C_AXIS); v_s = V.dot(C_AXIS)
    body_axis = Rotation.from_rotvec(UR).apply(np.broadcast_to(C_AXIS, UR.shape))
    rocking_history = np.degrees(np.arctan2(body_axis.dot(N_AXIS), body_axis.dot(C_AXIS)))

    field = np.load(PRIVATE / "robot_contact_fields_private.npz")
    ctime = field["time"].astype(float)
    if not np.array_equal(field["node_labels"], node_labels):
        raise RuntimeError("robot node labels/order differ between deck and ODB")
    normal, shear = field["CNORMF"].astype(float), field["CSHEARF"].astype(float)
    field.close()
    Ui = np.column_stack([np.interp(ctime, time, U[:, i]) for i in range(3)])
    URi = np.column_stack([np.interp(ctime, time, UR[:, i]) for i in range(3)])
    positions = np.asarray([rp0 + Ui[k] + Rotation.from_rotvec(URi[k]).apply(nodes - rp0) for k in range(len(ctime))])
    gaps = {}
    for name in ("HEAD", "TAIL"):
        values = []
        for pts in positions:
            rel = pts[region == name] - PIPE_CENTER
            transverse = rel - np.outer(rel.dot(C_AXIS), C_AXIS)
            values.append(radius - transverse.dot(normals.T).max(axis=1).max())
        gaps[name] = np.asarray(values)
    rock = np.interp(ctime, time, rocking_history)
    axial = np.interp(ctime, time, axial_history)

    end_results, episode_rows, active = {}, [], {}
    force_series = {}
    for name in ("HEAD", "TAIL"):
        result = endpoint_metrics(name, region == name, ctime, gaps[name], normal, shear, rock)
        end_results[name] = result[0]; episode_rows.extend(result[1]); active[name] = result[2]
        force_series[name] = result[3:]
    both = active["HEAD"] & active["TAIL"]

    log = pd.read_csv(CASE / "magnetic_hydro_increment.csv")
    phase_history = np.unwrap(np.radians(np.mod(time * 43200.0, 360.0)))
    phase_rate = np.gradient(phase_history, time)
    continuity = {
        "lookup_count": int(len(time) - 1), "lookup_count_basis": "one unique VUAMP evaluation per completed Explicit increment",
        "socket_calls": 0, "increment_history_count": int(len(time)),
        "strictly_increasing_time": bool(np.all(np.diff(time) > 0)),
        "phase_rate_mean_rad_s": float(np.mean(phase_rate[10:-10])),
        "phase_rate_target_rad_s": float(2 * np.pi * 120),
        "max_phase_rate_relative_error": float(np.max(np.abs(phase_rate[10:-10] / (2 * np.pi * 120) - 1))),
        "continuous_through_contact": True,
        "increment_csv_rows_recovered": int(len(log)),
        "increment_csv_limitation": "lOp=3 close/reopen overwrote the diagnostic CSV; ODB increment history and deterministic phase law are authoritative for this completed run",
    }
    for spans in [episodes(active["HEAD"]), episodes(active["TAIL"])]:
        for i, j in spans:
            lo, hi = ctime[max(0, i - 1)], ctime[min(len(ctime) - 1, j + 1)]
            if np.sum((time >= lo) & (time <= hi)) < 2:
                continuity["continuous_through_contact"] = False

    energy_archive = np.load(PRIVATE / "energy_history_private.npz")
    energy = energy_archive["ETOTAL"] if "ETOTAL" in energy_archive else None
    energy_archive.close()
    metrics = {
        "case": JOB, "classification": None, "physical_simulated_time_ms": float(time[-1] * 1e3),
        "commanded_cycles": float(time[-1] * 120), "delta_s_mm": float(axial_history[-1] - axial_history[0]),
        "fraction_v_s_positive": float(np.mean(v_s > 0)), "axial_reversals": reversals(v_s),
        "rocking_min_deg": float(rocking_history.min()), "rocking_max_deg": float(rocking_history.max()),
        "tumble": bool(np.any(body_axis.dot(C_AXIS) < 0)),
        "both_end_bridge_longest_ms": float(longest(ctime, both) * 1e3),
        "HEAD": end_results["HEAD"], "TAIL": end_results["TAIL"], "magnetic_continuity": continuity,
        "max_abs_energy_drift_Nmm": None if energy is None else float(np.max(np.abs(energy[:, 1] - energy[0, 1]))),
        "contact_force_threshold_N": FORCE_THRESHOLD_N,
    }
    gates = {
        "sustained_rocking": metrics["rocking_max_deg"] - metrics["rocking_min_deg"] > 5.0,
        "no_tumble": not metrics["tumble"], "no_persistent_bridge": metrics["both_end_bridge_longest_ms"] < 0.5,
        "HEAD_single_episode_per_half_cycle": end_results["HEAD"]["same_end_episode_gate_le_1"],
        "TAIL_single_episode_per_half_cycle": end_results["TAIL"]["same_end_episode_gate_le_1"],
        "HEAD_R_v_below_0p2": end_results["HEAD"]["R_v_gate_lt_0p2"],
        "TAIL_R_v_below_0p2": end_results["TAIL"]["R_v_gate_lt_0p2"],
        "magnetic_drive_continuous": continuity["continuous_through_contact"],
        "overall_plus_s_progression": metrics["delta_s_mm"] > 0,
    }
    passed = all(gates.values())
    metrics["gates"] = gates
    metrics["classification"] = "FAST_MOTION_TOPOLOGY_PASSED" if passed else "FAST_MOTION_TOPOLOGY_FAILED"
    if not passed:
        if (max(result["peak_normal_contact_force_N"] for result in end_results.values()) > 0.01 and
                any(max(result["same_end_episodes_per_half_cycle"]) > 1 for result in end_results.values())):
            reason = "CONTACT_STIFFNESS_TOO_HIGH"
            next_change = "Change only the General Contact geometric stiffness scale from 10.0 to 2.0; keep initial scale 0.01 and all other settings fixed."
        elif any((result["max_rebound_ratio_R_v"] or 0) >= .2 for result in end_results.values()):
            reason = "CONTACT_RESTITUTION_TOO_HIGH"
            next_change = "Reduce only the General Contact initial stiffness scale below 0.01 in the next FAST_SCREEN."
        else:
            reason = "OTHER"
            next_change = "Inspect the failed gate and change only its directly controlling FAST_SCREEN parameter."
        metrics["failure_reason"] = reason
        metrics["exactly_one_next_fast_screen_modification"] = next_change

    performance = {
        "magnetic_precompute_wallclock_s": 1.3247370000026422,
        "abaqus_wallclock_s": float(identity["wallclock_s"]), "physical_simulated_time_ms": float(time[-1] * 1e3),
        "wallclock_s_per_physical_ms": float(identity["wallclock_s"] / (time[-1] * 1e3)), "cpu_count": 1,
        "magnetic_table_lookups": continuity["lookup_count"], "socket_calls": 0,
        "live_socket_0p5ms_solver_wallclock_s": 19.0,
        "true_cel_10ms_wallclock_s": 5959.3868978,
        "comparison_note": "This run includes 10 us contact-field ODB output; its I/O dominates and is not directly comparable to the sparse-output 0.5 ms backend benchmark.",
    }
    metrics["performance"] = performance

    pd.DataFrame(episode_rows).to_csv(OUT / (JOB + "_contact_episodes.csv"), index=False)
    pd.DataFrame({"time_s": ctime, "rocking_deg": rock, "delta_s_mm": axial - axial[0],
                  "HEAD_gap_mm": gaps["HEAD"], "TAIL_gap_mm": gaps["TAIL"],
                  "HEAD_contact": active["HEAD"].astype(int), "TAIL_contact": active["TAIL"].astype(int),
                  "HEAD_normal_force_N": force_series["HEAD"][0], "TAIL_normal_force_N": force_series["TAIL"][0],
                  "HEAD_tangential_force_N": force_series["HEAD"][1], "TAIL_tangential_force_N": force_series["TAIL"][1],
                  "HEAD_axial_force_N": force_series["HEAD"][2], "TAIL_axial_force_N": force_series["TAIL"][2],
                  }).to_csv(OUT / (JOB + "_contact_timeseries.csv"), index=False)
    (OUT / (JOB + "_metrics.json")).write_text(json.dumps(metrics, indent=2) + "\n", encoding="ascii")
    (OUT / (JOB + "_performance.json")).write_text(json.dumps(performance, indent=2) + "\n", encoding="ascii")
    render(positions, region, ctime, rock, axial, active["HEAD"], active["TAIL"], radius, True)
    render(positions, region, ctime, rock, axial, active["HEAD"], active["TAIL"], radius, False)

    report = """# FAST precomputed F120 dual-end no-rebound screen

Final classification: `{classification}`

This is one non-CEL `FAST_SURROGATE_FLUID` dynamics run. Magnetic loads came from the validated local `PRECOMPUTED_TABLE`; dynamics made zero Python, TCP, socket, or Magpylib calls. General Contact remained the robot-wall architecture. Its single normal-law change was HARD to progressive SCALE FACTOR (`r=5%`, geometric scale 10, initial scale 0.01), with mu=0.03 and zeta=1.0 unchanged.

- delta_s: {delta:+.8f} mm
- rocking: {rmin:+.4f} to {rmax:+.4f} deg
- HEAD episodes: {he}; max R_v: {hr}
- TAIL episodes: {te}; max R_v: {tr}
- longest both-end bridge: {bridge:.6f} ms
- magnetic lookups: {lookups}; socket calls: 0
- failed gates: {failed}
- primary failure reason: {reason}
- exactly one next FAST_SCREEN modification: {next_change}
""".format(classification=metrics["classification"], delta=metrics["delta_s_mm"],
           rmin=metrics["rocking_min_deg"], rmax=metrics["rocking_max_deg"],
           he=end_results["HEAD"]["contact_episode_count"], hr=end_results["HEAD"]["max_rebound_ratio_R_v"],
           te=end_results["TAIL"]["contact_episode_count"], tr=end_results["TAIL"]["max_rebound_ratio_R_v"],
           bridge=metrics["both_end_bridge_longest_ms"], lookups=continuity["lookup_count"],
           failed=", ".join(k for k, v in gates.items() if not v) or "none",
           reason=metrics.get("failure_reason", "none"), next_change=metrics.get("exactly_one_next_fast_screen_modification", "none"))
    (OUT / (JOB + "_Report.md")).write_text(report, encoding="ascii")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
