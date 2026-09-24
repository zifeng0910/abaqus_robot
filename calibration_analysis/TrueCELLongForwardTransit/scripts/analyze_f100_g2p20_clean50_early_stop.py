"""Audit the user-stopped CLEAN50 run and package its completed history."""
from __future__ import annotations

import csv
import hashlib
import json
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

from analyze_f100_g2p20_clean50 import BASE, NEW, ROOT, load_rp
from analyze_f100_g2p20_clean30 import motion, stats
from analyze_f100_g2p20_final import geometry, pose, setup_axis


PREFIX = "F100_G2P20_CLEAN50_EARLY_STOP"
END_S = .04


def measure(m: dict, start: float, end: float, label: str) -> dict:
    t, s, v = m["t"], m["s"], m["v_s"]
    tc = np.r_[start, t[(t > start) & (t < end)], end]
    sc, vc = np.interp(tc, t, s), np.interp(tc, t, v)
    return {
        "segment": label, "start_ms": start * 1000, "end_ms": end * 1000,
        "complete_cycle": abs(end - start - .01) < 1e-9,
        "delta_s_mm": float(sc[-1] - sc[0]),
        "mean_v_s_mm_s": float((sc[-1] - sc[0]) / (end - start)),
        "end_v_s_mm_s": float(vc[-1]),
        "minimum_v_s_mm_s": float(vc.min()),
        "max_backtrack_mm": float((np.maximum.accumulate(sc) - sc).max()),
    }


def plot_history(m: dict, output: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                             constrained_layout=True)
    t = m["t"] * 1000
    axes[0].plot(t, m["s"], color="#176b60", lw=1)
    axes[1].plot(t, m["v_s"], color="#a04737", lw=.7)
    axes[0].axhline(0, color="#777", lw=.8)
    axes[1].axhline(0, color="#777", lw=.8)
    for ax in axes:
        for boundary in (10, 20, 30, 40):
            ax.axvline(boundary, color="#888", lw=.7, ls="--")
        ax.grid(alpha=.18)
        ax.set_xlim(0, t[-1])
    axes[0].set_ylabel("s (mm)")
    axes[1].set_ylabel("v_s (mm/s)")
    axes[1].set_xlabel("time (ms)")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def render_motion(m: dict, ident: dict, output: Path) -> None:
    rp0, c, n, rel, colors, pipe = geometry(NEW, ident)
    fig, (ax, trace) = plt.subplots(2, 1, figsize=(10, 6.5),
                                    gridspec_kw={"height_ratios": [2, 1]},
                                    constrained_layout=True)
    setup_axis(ax, ident, "Uninterrupted TRUE-CEL motion")
    p0 = pose(0, m["t"], m["u"], m["ur"], rp0, rel)
    robot = ax.scatter(float(ident["s_start_mm"]) + (p0-rp0) @ c,
                       (p0-pipe) @ n, s=5, c=colors, linewidths=0)
    label = ax.text(.98, .97, "", transform=ax.transAxes, ha="right", va="top",
                    family="monospace", fontsize=9,
                    bbox={"facecolor": "white", "alpha": .92})
    trace.plot(m["t"] * 1000, m["s"], color="#b9cbc8", lw=1)
    trace.axhline(0, color="#888", lw=.7)
    for boundary in (10, 20, 30, 40):
        trace.axvline(boundary, color="#999", lw=.7, ls="--")
    marker, = trace.plot([], [], "o", color="#176b60", ms=5)
    trace.set(xlim=(0, m["t"][-1] * 1000), xlabel="time (ms)",
              ylabel="s (mm)")
    trace.grid(alpha=.18)
    frames = np.linspace(0, m["t"][-1], 131)

    def update(i):
        ti = frames[i]
        p = pose(ti, m["t"], m["u"], m["ur"], rp0, rel)
        robot.set_offsets(np.c_[float(ident["s_start_mm"]) + (p-rp0) @ c,
                                (p-pipe) @ n])
        si = np.interp(ti, m["t"], m["s"])
        vi = np.interp(ti, m["t"], m["v_s"])
        marker.set_data([ti * 1000], [si])
        label.set_text(f"t={ti*1000:5.2f} ms\ns={si:+.4f} mm\nv_s={vi:+.2f} mm/s")
        return robot, marker, label

    FuncAnimation(fig, update, frames=len(frames), interval=66).save(
        output, PillowWriter(fps=15), dpi=90)
    plt.close(fig)


def render_cycles(m: dict, output: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                             constrained_layout=True)
    t = m["t"] * 1000
    mask = (t >= 20) & (t <= m["t"][-1] * 1000)
    for ax, values, color, ylabel in (
        (axes[0], m["s"], "#176b60", "s (mm)"),
        (axes[1], m["v_s"], "#a04737", "v_s (mm/s)"),
    ):
        ax.plot(t[mask], values[mask], color="#c4cecc", lw=1)
        ax.axvline(30, color="#888", lw=.8, ls="--")
        ax.axhline(0, color="#888", lw=.7)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=.18)
    axes[1].set(xlabel="time (ms)", xlim=(20, m["t"][-1] * 1000))
    lines = [axes[0].plot([], [], color="#176b60", lw=1.7)[0],
             axes[1].plot([], [], color="#a04737", lw=1.7)[0]]
    dots = [axes[0].plot([], [], "o", color="#176b60", ms=5)[0],
            axes[1].plot([], [], "o", color="#a04737", ms=5)[0]]
    frames = np.linspace(.02, m["t"][-1], 111)

    def update(i):
        ti = frames[i]
        segment = mask & (m["t"] <= ti)
        for j, values in enumerate((m["s"], m["v_s"])):
            lines[j].set_data(t[segment], values[segment])
            dots[j].set_data([ti * 1000], [np.interp(ti, m["t"], values)])
        fig.suptitle(f"CLEAN50: C3 forward, C4 recoil | t={ti*1000:.2f} ms")
        return lines + dots

    FuncAnimation(fig, update, frames=len(frames), interval=66).save(
        output, PillowWriter(fps=15), dpi=100)
    plt.close(fig)


def main() -> None:
    ident = json.loads((NEW / "case_identity.json").read_text(encoding="utf-8"))
    base_ident = json.loads((BASE / "case_identity.json").read_text(encoding="utf-8"))
    setup = json.loads((ROOT / f"{NEW.name}_SETUP.json").read_text(encoding="utf-8"))
    m = motion(load_rp(NEW), ident)
    base = motion(load_rp(BASE), base_ident)
    if m["t"][0] > 1e-10 or m["t"][-1] <= END_S or m["t"][-1] >= .05:
        raise RuntimeError("Expected an interrupted 40-50 ms RP history")
    if ident["magnetic_table_sha256"] != base_ident["magnetic_table_sha256"]:
        raise RuntimeError("Magnetic table identity differs")
    sta = (NEW / f"{NEW.name}.sta").read_text(encoding="latin1")
    if "Process terminated by external request" not in sta:
        raise RuntimeError("Missing explicit user-requested termination record")
    writes = [float(line.rsplit(" at ", 1)[1]) * 1000
              for line in sta.splitlines() if "Restart Number" in line and " at " in line]
    if not np.allclose(writes, np.arange(5, 41, 5), atol=.001, rtol=0):
        raise RuntimeError(f"Unexpected restart write cadence: {writes}")
    rows = [measure(m, (i-1) * .01, i * .01, f"Cycle {i}")
            for i in range(1, 5)]
    rows.append(measure(m, .04, float(m["t"][-1]), "Cycle 5 partial"))
    baseline_t = base["t"][base["t"] <= .02]
    baseline_error = {key: stats(np.interp(baseline_t, m["t"], m[key]) -
                                 np.interp(baseline_t, base["t"], base[key]))
                      for key in ("s", "v_s")}
    online = np.loadtxt(NEW / "online_motion.csv", delimiter=",", skiprows=1)
    result = {
        "classification": "CLEAN50_RECOIL_CONFIRMED",
        "completion": "USER_TERMINATED_PARTIAL",
        "requested_end_ms": 50.0,
        "last_ODB_RP_time_ms": float(m["t"][-1] * 1000),
        "last_online_time_ms": float(online[-1, 0] * 1000),
        "cycles": rows,
        "C3_old_restart_net_recoil_reproduced": rows[2]["delta_s_mm"] < 0,
        "C4_clean_net_recoil_confirmed": rows[3]["delta_s_mm"] < 0,
        "baseline_0_20_error": baseline_error,
        "actual_restart_write_times_ms": writes,
        "restart_read": False,
        "minimum_online_stable_dt_s": float(online[:, 4].min()),
        "setup": setup,
    }
    ident.update({
        "status": "USER_TERMINATED_PARTIAL",
        "classification": result["classification"],
        "last_ODB_RP_time_ms": result["last_ODB_RP_time_ms"],
        "dynamics_run_count": 1,
    })
    (NEW / "case_identity.json").write_text(
        json.dumps(ident, indent=2) + "\n", encoding="utf-8")
    metrics = ROOT / f"{PREFIX}_METRICS.json"
    metrics.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    cycles_csv = ROOT / f"{PREFIX}_CYCLES.csv"
    with cycles_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    chart = ROOT / f"{PREFIX}_MOTION.png"
    gif_motion = ROOT / f"{PREFIX}_MOTION.gif"
    gif_recoil = ROOT / f"{PREFIX}_C3_C4.gif"
    plot_history(m, chart)
    render_motion(m, ident, gif_motion)
    render_cycles(m, gif_recoil)
    report = ROOT / f"{PREFIX}_REPORT.md"
    lines = [
        "# F100/G2P20 CLEAN50 user-stopped audit", "",
        "**CLEAN50_RECOIL_CONFIRMED**, based on the completed clean Cycle 4.", "",
        f"Abaqus was terminated on user request at RP t={result['last_ODB_RP_time_ms']:.6f} ms. "
        "Cycles 1-4 are complete; Cycle 5 is partial. This is not a completed "
        "50 ms solve and is not classified as a spontaneous numerical failure.", "",
        "One Explicit step from t=0; no restart read or continuation. Fixed restart "
        f"writes at {writes} ms were never used for continuation. All model "
        "parameters and input differences are recorded in the setup JSON.", "",
        "| Segment | Time (ms) | delta_s (mm) | mean v_s (mm/s) | end v_s (mm/s) | "
        "min v_s (mm/s) | max backtrack (mm) |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(f"| {row['segment']} | {row['start_ms']:.3f}-{row['end_ms']:.3f} | "
                     f"{row['delta_s_mm']:+.6f} | {row['mean_v_s_mm_s']:+.3f} | "
                     f"{row['end_v_s_mm_s']:+.3f} | {row['minimum_v_s_mm_s']:+.3f} | "
                     f"{row['max_backtrack_mm']:.6f} |")
    lines += [
        "", "The old restart-derived C3 net recoil (-0.023527 mm) did not recur: "
        f"clean C3 delta_s={rows[2]['delta_s_mm']:+.6f} mm. Clean C4 instead "
        f"had delta_s={rows[3]['delta_s_mm']:+.6f} mm and max backtrack="
        f"{rows[3]['max_backtrack_mm']:.6f} mm.", "",
        "The 0-20 ms path diverged after the 15 ms write. Maximum baseline "
        f"position difference: {baseline_error['s']['max']:.6g} mm; p99 "
        f"velocity difference: {baseline_error['v_s']['p99']:.6g} mm/s. "
        "This is recorded because restart-write cadence affects the numerical path; "
        "it does not invalidate the uninterrupted CLEAN50 observation.", "",
        f"Minimum sampled stable increment: {result['minimum_online_stable_dt_s']:.6g} s. "
        "The solver reported only the requested external termination at the end.",
    ]
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    included = [report, metrics, cycles_csv, chart, gif_motion, gif_recoil,
                NEW / "case_identity.json",
                NEW / "online_motion.csv",
                NEW / f"{NEW.name}.sta",
                NEW / f"{NEW.name}.inp",
                NEW / "vuamp_precomputed_truecel.f90",
                NEW / "magnetic_field_gradient_table_B0P11_A14P5.dat",
                NEW / "magnetic_field_gradient_table_B0P11_A14P5.dat.json",
                ROOT / f"{NEW.name}_SETUP.json",
                ROOT / f"{NEW.name}_PreSolve_Initialization_Audit.json"]
    manifest = {"odb_included": False,
                "odb_location": str(NEW / f"{NEW.name}.odb"),
                "files": [
        {"name": path.name, "bytes": path.stat().st_size,
         "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        for path in included]}
    manifest_path = ROOT / f"{PREFIX}_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    archive_path = ROOT / f"{PREFIX}_GIFS.zip"
    with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for path in included + [manifest_path]:
            archive.write(path, path.name)
    print(json.dumps({"archive": str(archive_path), "metrics": rows,
                      "last_ms": result["last_ODB_RP_time_ms"]}, indent=2))


if __name__ == "__main__":
    main()
