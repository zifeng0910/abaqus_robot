"""Audit the completed uninterrupted 3x rotational-damping FULL CEL run."""
from __future__ import annotations

import csv
import json
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

from analyze_f100_g2p20_clean30 import motion
from analyze_f100_g2p20_clean50 import load_rp
from analyze_f100_g2p20_final import geometry, pose, setup_axis
from analyze_magnetic_only_damping_screen import wall_metrics


ROOT = Path(__file__).resolve().parents[1]
JOB = "TRUECEL_B0P11_G2P20_A14P5_F100_CROT3_FULLCEL50"
CASE = ROOT / "case" / JOB
BASE = ROOT / "case" / "TRUECEL_B0P11_G2P20_A14P5_F100_CLEAN50"
PREFIX = "F100_G2P20_CROT3_FULLCEL50"


def cycle_rows(m: dict) -> list[dict]:
    rows = []
    for cycle in range(1, 6):
        start, end = (cycle - 1) * .01, cycle * .01
        t = np.r_[start, m["t"][(m["t"] > start) & (m["t"] < end)], end]
        s = np.interp(t, m["t"], m["s"])
        v = np.interp(t, m["t"], m["v_s"])
        angle = np.interp(t, m["t"], m["rocking_angle"])
        omega = np.interp(t, m["t"], m["omega_rock"])
        basis = np.column_stack((np.sin(2 * np.pi * 100 * t),
                                 np.cos(2 * np.pi * 100 * t), np.ones(len(t))))
        a, b, offset = np.linalg.lstsq(basis, angle, rcond=None)[0]
        fitted = basis @ np.array([a, b, offset])
        amp = float(np.hypot(a, b))
        rows.append({
            "cycle": cycle, "delta_s_mm": float(s[-1] - s[0]),
            "mean_v_s_mm_s": float((s[-1] - s[0]) / .01),
            "end_v_s_mm_s": float(v[-1]),
            "minimum_v_s_mm_s": float(v.min()),
            "max_backtrack_mm": float((np.maximum.accumulate(s) - s).max()),
            "rocking_half_range_deg": float(np.ptp(angle) / 2),
            "rocking_fundamental_amplitude_deg": amp,
            "phase_lag_deg": float((-np.degrees(np.arctan2(b, a)) + 180) % 360 - 180),
            "theta_offset_deg": float(offset),
            "sinusoid_residual_rms_over_amplitude": float(
                np.sqrt(np.mean((angle - fitted) ** 2)) / max(amp, 1e-12)),
            "omega_peak_abs_rad_s": float(np.max(np.abs(omega))),
        })
    return rows


def plot(m: dict, baseline: dict, rows: list[dict], output: Path) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(12, 9), constrained_layout=True)
    t = m["t"] * 1000
    tb = baseline["t"] * 1000
    for ax, key, label in ((axes[0, 0], "s", "s (mm)"),
                           (axes[0, 1], "v_s", "v_s (mm/s)"),
                           (axes[1, 0], "rocking_angle", "theta (deg)"),
                           (axes[1, 1], "omega_rock", "omega (rad/s)")):
        ax.plot(tb, baseline[key], color="#999999", lw=.8, label="clean FULL CEL, 0x")
        ax.plot(t, m[key], color="#176b60", lw=.9, label="FULL CEL, 3x")
        ax.set(xlabel="time (ms)", ylabel=label, xlim=(0, 50))
        for boundary in (10, 20, 30, 40):
            ax.axvline(boundary, color="#cccccc", lw=.6, ls="--")
        ax.grid(alpha=.15)
    axes[1, 0].plot(t, 14.5 * np.sin(2 * np.pi * 100 * m["t"]),
                    color="#b4663d", lw=.6, ls="--", label="14.5 deg command")
    axes[0, 0].legend(fontsize=8)
    axes[1, 0].legend(fontsize=8)
    cycles = [row["cycle"] for row in rows]
    axes[2, 0].plot(cycles, [r["rocking_fundamental_amplitude_deg"] for r in rows],
                    "o-", color="#176b60")
    axes[2, 1].plot(cycles, [r["phase_lag_deg"] for r in rows],
                    "o-", color="#a34c3f")
    axes[2, 0].set(xlabel="cycle", ylabel="first-harmonic amplitude (deg)", xticks=cycles)
    axes[2, 1].set(xlabel="cycle", ylabel="phase lag (deg)", xticks=cycles)
    axes[2, 0].grid(alpha=.15)
    axes[2, 1].grid(alpha=.15)
    fig.savefig(output, dpi=160)
    plt.close(fig)


def gif(m: dict, ident: dict, rows: list[dict], output: Path) -> None:
    rp0, c, n, rel, colors, pipe = geometry(CASE, ident)
    fig, (ax, trace) = plt.subplots(2, 1, figsize=(10, 6.5),
                                    gridspec_kw={"height_ratios": [2, 1]},
                                    constrained_layout=True)
    setup_axis(ax, ident, "FULL CEL | 3x rotational damping")
    p0 = pose(0, m["t"], m["u"], m["ur"], rp0, rel)
    robot = ax.scatter(ident["s_start_mm"] + (p0-rp0) @ c,
                       (p0-pipe) @ n, s=5, c=colors, linewidths=0)
    trace.plot(m["t"] * 1000, m["rocking_angle"], color="#b9cbc8", lw=.8)
    trace.plot(m["t"] * 1000, 14.5*np.sin(2*np.pi*100*m["t"]),
               color="#bd784d", lw=.7, ls="--")
    point, = trace.plot([], [], "o", color="#176b60", ms=5)
    trace.set(xlim=(0, 50), xlabel="time (ms)", ylabel="theta / command (deg)")
    trace.grid(alpha=.15)
    label = ax.text(.98, .98, "", transform=ax.transAxes, ha="right", va="top",
                    family="monospace", fontsize=8,
                    bbox={"facecolor": "white", "alpha": .9})
    times = np.linspace(0, .05, 121)

    def update(i):
        ti = times[i]
        p = pose(ti, m["t"], m["u"], m["ur"], rp0, rel)
        robot.set_offsets(np.c_[ident["s_start_mm"]+(p-rp0)@c, (p-pipe)@n])
        angle = float(np.interp(ti, m["t"], m["rocking_angle"]))
        displacement = float(np.interp(ti, m["t"], m["s"]))
        cycle = min(5, int(ti/.01)+1)
        point.set_data([ti*1000], [angle])
        label.set_text(f"t={ti*1000:5.2f} ms  C{cycle}\n"
                       f"theta={angle:+.1f} deg  lag={rows[cycle-1]['phase_lag_deg']:+.1f} deg\n"
                       f"s={displacement:+.4f} mm")
        return robot, point, label

    FuncAnimation(fig, update, frames=len(times), interval=67).save(
        output, PillowWriter(fps=15), dpi=90)
    plt.close(fig)


def main() -> None:
    ident = json.loads((CASE / "case_identity.json").read_text(encoding="utf-8-sig"))
    base_ident = json.loads((BASE / "case_identity.json").read_text(encoding="utf-8-sig"))
    sta = (CASE / f"{JOB}.sta").read_text(encoding="latin1")
    if "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" not in sta or "***ERROR" in sta:
        raise RuntimeError("FULL CEL solve did not pass completion gate")
    m = motion(load_rp(CASE), ident)
    baseline = motion(load_rp(BASE), base_ident)
    if m["t"][0] > 1e-10 or m["t"][-1] < .05 - 1e-8:
        raise RuntimeError("Incomplete uninterrupted RP history")
    if not all(np.all(np.isfinite(m[k])) for k in
               ("s", "v_s", "rocking_angle", "omega_rock")):
        raise RuntimeError("Nonfinite RP history")
    rows = cycle_rows(m)
    wall = wall_metrics(CASE, ident, m)
    with np.load(CASE / "private" / "contact_history_private.npz") as z:
        key = next(k for k in z.files if "ROBOT_SOLID" in k and "|CFNM " in k)
        contact = z[key].astype(float)
    contact_summary = {"whole_robot_general_contact_CFNM_peak_N": float(contact[:, 1].max()),
                       "whole_robot_general_contact_nonzero_fraction": float(
                           np.mean(contact[:, 1] > 1e-9)),
                       "interpretation": "whole robot General Contact; not pair-isolated robot-wall force"}
    writes = [float(line.rsplit(" at ", 1)[1]) * 1000 for line in sta.splitlines()
              if "Restart Number" in line and " at " in line]
    if not np.allclose(writes, np.arange(5, 51, 5), atol=.001, rtol=0):
        raise RuntimeError(f"Unexpected restart write cadence: {writes}")
    nofluid = json.loads((ROOT / "TRUECEL_B0P11_G2P20_F100_NOFLUID_CONTROL_METRICS.json").read_text())
    screen = json.loads((ROOT / "TRUECEL_MAGNETIC_ONLY_DAMPING_METRICS.json").read_text())
    nofluid_3x = next(x for x in screen["results"] if x["factor"] == 3.0)
    baseline_report = json.loads((ROOT / "F100_G2P20_CLEAN50_EARLY_STOP_METRICS.json").read_text())
    result = {"classification": "CROT3_FULLCEL50_RECOIL_CONFIRMED",
              "solver_completed": True, "t_start_s": float(m["t"][0]),
              "t_end_s": float(m["t"][-1]), "single_step": True,
              "restart_read": False, "restart_writes_ms": writes,
              "c_rot_N_mm_s": ident["rotational_damping_c_N_mm_s"],
              "cycles": rows, "net_displacement_mm": float(m["s"][-1]-m["s"][0]),
              "whole_run_max_backtrack_mm": float((np.maximum.accumulate(m["s"])-m["s"]).max()),
              "wall_geometry": wall, "general_contact": contact_summary,
              "clean_full_CEL_baseline_cycles_1_4_delta_s_mm": [
                  r["delta_s_mm"] for r in baseline_report["cycles"][:4]],
              "no_fluid_baseline_cycles_delta_s_mm": [
                  r["delta_s_mm"] for r in nofluid["cycle_metrics"]],
              "no_fluid_3x_rocking_cycles": nofluid_3x["cycles"],
              "contact_limit": "geometry overlap is sampled at 0.05 ms; whole-robot CFNM includes fluid contact"}
    metrics = ROOT / f"{PREFIX}_METRICS.json"
    metrics.write_text(json.dumps(result, indent=2) + "\n")
    cycles = ROOT / f"{PREFIX}_CYCLES.csv"
    with cycles.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    chart = ROOT / f"{PREFIX}_ROCKING_MOTION.png"
    animation = ROOT / f"{PREFIX}.gif"
    plot(m, baseline, rows, chart)
    gif(m, ident, rows, animation)
    report = ROOT / f"{PREFIX}_REPORT.md"
    lines = ["# F100/G2P20 3x rotational damping, full CEL", "",
             "**CROT3_FULLCEL50_RECOIL_CONFIRMED.** Completed one Explicit step from t=0 to 50 ms; restart snapshots were written at 5-ms intervals and never read.", "",
             "| Cycle | delta s (mm) | mean v (mm/s) | end v (mm/s) | minimum v (mm/s) | max backtrack (mm) | rocking amplitude (deg) | phase lag (deg) |",
             "|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r['cycle']} | {r['delta_s_mm']:+.6f} | {r['mean_v_s_mm_s']:+.3f} | {r['end_v_s_mm_s']:+.3f} | {r['minimum_v_s_mm_s']:+.3f} | {r['max_backtrack_mm']:.6f} | {r['rocking_fundamental_amplitude_deg']:.3f} | {r['phase_lag_deg']:+.2f} |")
    lines += ["", f"Final net displacement: {result['net_displacement_mm']:+.6f} mm; whole-run max backtrack: {result['whole_run_max_backtrack_mm']:.6f} mm.",
              f"Geometric wall overlap: {wall['virtual_wall_crossing_events']} sampled episodes, minimum signed gap {wall['signed_min_wall_gap_mm']:.6f} mm. This is a sampled geometry indicator, not an isolated contact-force count.",
              f"Whole-robot General Contact CFNM peak: {contact_summary['whole_robot_general_contact_CFNM_peak_N']:.6g} N; this signal includes fluid contact.",
              "", "The clean FULL CEL baseline had C3 +0.013300 mm and C4 -0.564697 mm before user termination at 40.75 ms. Magnetic-only 3x had nearly constant 14.728-deg amplitude and 34.17-deg lag in C2-C5, but its passive-wall audit also found wall overlap. The coupled 3x case therefore does not establish a physical wall-free rocking attractor.",
              "", "Next bounded test: reduce only rotational damping to 1x (3.5006005000754733e-6 N mm s), whose magnetic-only screen showed stable C2-C5 amplitude and phase; hold all other model terms fixed."]
    report.write_text("\n".join(lines) + "\n")
    archive = ROOT / f"{PREFIX}_GIFS.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
        for path in (animation, chart, report, metrics, cycles):
            z.write(path, path.name)
    ident.update({"classification": result["classification"],
                  "last_ODB_RP_time_ms": float(m["t"][-1] * 1000)})
    (CASE / "case_identity.json").write_text(json.dumps(ident, indent=2) + "\n")
    print(json.dumps({"classification": result["classification"],
                      "net_displacement_mm": result["net_displacement_mm"],
                      "wall": wall, "archive": str(archive)}, indent=2))


if __name__ == "__main__":
    main()
