"""Compare the four clean FULL CEL damping trajectories after 0.3x completes."""
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
from scipy.spatial import ConvexHull

from analyze_f100_g2p20_clean30 import motion
from analyze_f100_g2p20_clean50 import load_rp
from analyze_f100_g2p20_final import geometry, pose, setup_axis
from analyze_magnetic_only_damping_screen import wall_facets
from analyze_f100_crot3_fullcel import cycle_rows


ROOT = Path(__file__).resolve().parents[1]
CASES = (
    ("0x", "TRUECEL_B0P11_G2P20_A14P5_F100_CLEAN50", 0.0),
    ("0.3x", "TRUECEL_B0P11_G2P20_A14P5_F100_CROT0P3_FULLCEL50", .3),
    ("1x", "TRUECEL_B0P11_G2P20_A14P5_F100_CROT1_FULLCEL50", 1.0),
    ("3x", "TRUECEL_B0P11_G2P20_A14P5_F100_CROT3_FULLCEL50", 3.0),
)
PREFIX = "F100_G2P20_CROT0_CROT03_CROT1_CROT3_COMPARISON"


def sampled_wall(case: Path, ident: dict, m: dict) -> dict:
    rp0, _, _, rel, _, _ = geometry(case, ident)
    rel = rel[ConvexHull(rel).vertices]
    pipe, basis, planes = wall_facets(case, ident)
    end = min(.05, float(m["t"][-1]))
    times = np.arange(0, end + 1e-10, .00005)
    gaps = np.empty(len(times))
    for index, ti in enumerate(times):
        p = pose(ti, m["t"], m["u"], m["ur"], rp0, rel)
        gaps[index] = -np.max((p-pipe) @ basis @ planes[:, :2].T + planes[:, 2])
    outside = gaps < -1e-6
    return {"minimum_signed_gap_mm": float(gaps.min()),
            "negative_gap_episodes": int(outside[0]) + int(np.count_nonzero(outside[1:] & ~outside[:-1])),
            "negative_gap_duration_ms": float(outside.sum() * .05),
            "sample_end_ms": float(times[-1] * 1000)}


def load_case(label: str, name: str, factor: float) -> dict:
    case = ROOT / "case" / name
    ident = json.loads((case / "case_identity.json").read_text(encoding="utf-8-sig"))
    m = motion(load_rp(case), ident)
    complete = label in ("1x", "3x")
    sta = (case / f"{name}.sta").read_text(encoding="latin1")
    if complete and ("THE ANALYSIS HAS COMPLETED SUCCESSFULLY" not in sta or
                     m["t"][-1] < .05 - 1e-8):
        raise RuntimeError(f"Incomplete case: {name}")
    if label == "0x" and (m["t"][-1] <= .04 or m["t"][-1] >= .05):
        raise RuntimeError("Expected user-stopped 0x baseline ending in C5")
    if label == "0.3x" and (m["t"][-1] <= .03 or m["t"][-1] >= .04 or
                            "Process terminated by external request" not in sta):
        raise RuntimeError("Expected user-stopped 0.3x case in C4")
    # The generic helper interpolates to 50 ms; retain only complete cycles.
    rows = cycle_rows(m)[:min(5, int((float(m["t"][-1]) + 1e-8) / .01))]
    with np.load(case / "private" / "energy_history_private.npz") as archive:
        energy = {key: archive[key].astype(float) for key in
                  ("ETOTAL", "ALLPW", "ALLFD", "ALLIE", "ALLKE")}
    online = np.loadtxt(case / "online_motion.csv", delimiter=",", skiprows=1)
    wall = sampled_wall(case, ident, m)
    amps = np.array([r["rocking_fundamental_amplitude_deg"] for r in rows[1:]])
    phase = np.unwrap(np.deg2rad([r["phase_lag_deg"] for r in rows[1:]]))
    cycle_energy = []
    for cycle in range(1, 6):
        end = cycle * .01
        if end > m["t"][-1] + 1e-8:
            break
        sample = {"case": label, "cycle": str(cycle), "end_ms": end*1000,
                  "complete_cycle": True}
        for key, array in energy.items():
            sample[f"{key}_N_mm"] = float(np.interp(end, array[:, 0], array[:, 1]))
        dt = online[(online[:, 0] > end-.01) & (online[:, 0] <= end), 4]
        sample["minimum_stable_dt_s"] = float(dt.min()) if len(dt) else None
        cycle_energy.append(sample)
    if not complete:
        sample = {"case": label, "cycle": f"C{len(rows)+1} partial",
                  "end_ms": float(energy["ETOTAL"][-1, 0]*1000),
                  "complete_cycle": False}
        for key, array in energy.items():
            sample[f"{key}_N_mm"] = float(array[-1, 1])
        dt = online[online[:, 0] > len(rows)*.01, 4]
        sample["minimum_stable_dt_s"] = float(dt.min()) if len(dt) else None
        cycle_energy.append(sample)
    return {"label": label, "name": name, "factor": factor, "case": case,
            "identity": ident, "motion": m, "cycles": rows, "energy": energy,
            "cycle_energy": cycle_energy, "wall": wall,
            "max_backtrack_mm": float((np.maximum.accumulate(m["s"])-m["s"]).max()),
            "amplitude_spread_C2_end_fraction": float(np.ptp(amps)/np.mean(amps)),
            "phase_spread_C2_end_deg": float(np.ptp(np.rad2deg(phase))),
            "end_ETOTAL_N_mm": float(energy["ETOTAL"][-1, 1]),
            "end_ALLPW_N_mm": float(energy["ALLPW"][-1, 1]),
            "end_time_ms": float(m["t"][-1]*1000),
            "minimum_stable_dt_s": float(online[:, 4].min())}


def target_gif(target: dict, output: Path) -> None:
    m, ident = target["motion"], target["identity"]
    rp0, c, n, rel, colors, pipe = geometry(target["case"], ident)
    fig, (ax, trace) = plt.subplots(2, 1, figsize=(10, 6.5),
                                    gridspec_kw={"height_ratios": [2, 1]},
                                    constrained_layout=True)
    setup_axis(ax, ident, "FULL CEL | 0.3x damping | user stopped after recoil")
    p0 = pose(0, m["t"], m["u"], m["ur"], rp0, rel)
    robot = ax.scatter(ident["s_start_mm"]+(p0-rp0)@c,
                       (p0-pipe)@n, s=5, c=colors, linewidths=0)
    trace.plot(m["t"]*1000, m["rocking_angle"], color="#176b60", lw=.8,
               label="measured theta")
    trace.plot(m["t"]*1000, 14.5*np.sin(2*np.pi*100*m["t"]),
               color="#bd784d", lw=.7, ls="--", label="command")
    trace.set(xlim=(0, 50), xlabel="time (ms)", ylabel="theta / command (deg)")
    trace.grid(alpha=.15)
    trace.legend(fontsize=8)
    marker, = trace.plot([], [], "o", color="#176b60", ms=5)
    label = ax.text(.98, .98, "", transform=ax.transAxes, ha="right", va="top",
                    family="monospace", fontsize=8,
                    bbox={"facecolor": "white", "alpha": .92})
    times = np.linspace(0, m["t"][-1], 121)

    def update(index):
        ti = times[index]
        p = pose(ti, m["t"], m["u"], m["ur"], rp0, rel)
        robot.set_offsets(np.c_[ident["s_start_mm"]+(p-rp0)@c, (p-pipe)@n])
        angle = float(np.interp(ti, m["t"], m["rocking_angle"]))
        command = 14.5*np.sin(2*np.pi*100*ti)
        cycle = min(5, int(ti/.01)+1)
        marker.set_data([ti*1000], [angle])
        label.set_text(f"t={ti*1000:5.2f} ms  C{cycle}\n"
                       f"s={np.interp(ti,m['t'],m['s']):+.4f} mm\n"
                       f"v_s={np.interp(ti,m['t'],m['v_s']):+.2f} mm/s\n"
                       f"theta={angle:+.1f}  cmd={command:+.1f} deg")
        return robot, marker, label

    FuncAnimation(fig, update, frames=len(times), interval=67).save(
        output, PillowWriter(fps=15), dpi=90)
    plt.close(fig)


def comparison_plot(cases: list[dict], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    colors = ("#777777", "#1c7a70", "#ae7541", "#ac4d48")
    for result, color in zip(cases, colors):
        x = np.array([r["cycle"] for r in result["cycles"]])
        axes[0, 0].plot(x, [r["delta_s_mm"] for r in result["cycles"]],
                        "o-", label=result["label"], color=color)
        axes[0, 1].plot([result["factor"]], [result["max_backtrack_mm"]],
                        "o", label=result["label"], color=color, ms=8)
        axes[1, 0].plot(x, [r["rocking_fundamental_amplitude_deg"] for r in result["cycles"]],
                        "o-", label=result["label"], color=color)
        axes[1, 1].plot(x, [r["phase_lag_deg"] for r in result["cycles"]],
                        "o-", label=result["label"], color=color)
    axes[0, 0].axhline(0, color="#999999", lw=.7)
    axes[0, 0].set(xlabel="cycle", ylabel="delta s (mm)", xticks=range(1, 6),
                   title="Cycle displacement")
    axes[0, 1].set(xlabel="damping factor (x I 2pi f)",
                   ylabel="whole-run max backtrack (mm)", title="Maximum positional backtrack")
    axes[1, 0].set(xlabel="cycle", ylabel="first-harmonic amplitude (deg)",
                   xticks=range(1, 6), title="Rocking amplitude")
    axes[1, 1].set(xlabel="cycle", ylabel="phase lag (deg)",
                   xticks=range(1, 6), title="Rocking phase lag")
    for ax in axes.flat:
        ax.grid(alpha=.15)
    axes[0, 0].legend(title="FULL CEL")
    fig.savefig(output, dpi=170)
    plt.close(fig)


def comparison_gif(cases: list[dict], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 7), constrained_layout=True)
    artists = []
    for ax, result in zip(axes.flat, cases):
        ident, m = result["identity"], result["motion"]
        rp0, c, n, rel, colors, pipe = geometry(result["case"], ident)
        setup_axis(ax, ident, f"{result['label']} rotational damping")
        p = pose(0, m["t"], m["u"], m["ur"], rp0, rel)
        robot = ax.scatter(ident["s_start_mm"]+(p-rp0)@c,
                           (p-pipe)@n, s=3, c=colors, linewidths=0)
        label = ax.text(.99, .98, "", transform=ax.transAxes, ha="right", va="top",
                        family="monospace", fontsize=7,
                        bbox={"facecolor": "white", "alpha": .92})
        artists.append((result, rp0, c, n, rel, pipe, robot, label))
    times = np.linspace(0, .05, 121)

    def update(index):
        ti = times[index]
        changed = []
        for result, rp0, c, n, rel, pipe, robot, label in artists:
            m, ident = result["motion"], result["identity"]
            if ti > m["t"][-1] + 1e-8:
                robot.set_offsets(np.empty((0, 2)))
                label.set_text(f"t={ti*1000:.2f} ms\nno 0x data after {m['t'][-1]*1000:.2f} ms")
            else:
                p = pose(ti, m["t"], m["u"], m["ur"], rp0, rel)
                robot.set_offsets(np.c_[ident["s_start_mm"]+(p-rp0)@c, (p-pipe)@n])
                cycle = min(5, int(ti/.01)+1)
                angle = np.interp(ti, m["t"], m["rocking_angle"])
                command = 14.5*np.sin(2*np.pi*100*ti)
                label.set_text(f"t={ti*1000:5.2f} ms  C{cycle}\n"
                               f"s={np.interp(ti,m['t'],m['s']):+.4f} mm\n"
                               f"v_s={np.interp(ti,m['t'],m['v_s']):+.2f} mm/s\n"
                               f"theta={angle:+.1f}  cmd={command:+.1f} deg")
            changed.extend((robot, label))
        return changed

    FuncAnimation(fig, update, frames=len(times), interval=67).save(
        output, PillowWriter(fps=15), dpi=85)
    plt.close(fig)


def main() -> None:
    cases = [load_case(*entry) for entry in CASES]
    target, one = cases[1], cases[2]
    target_rows = target["cycles"]
    short_window_rocking_coherent = (
        target["amplitude_spread_C2_end_fraction"] < .10 and
        target["phase_spread_C2_end_deg"] < 5 and
        all(r["sinusoid_residual_rms_over_amplitude"] < .20
            for r in target_rows[1:]))
    energy_suspect = (abs(target["end_ETOTAL_N_mm"]) >= .1*abs(one["end_ETOTAL_N_mm"]) or
                      target["end_ALLPW_N_mm"] >= .1*one["end_ALLPW_N_mm"])
    energy_class = ("ENERGY_BEHAVIOR_NUMERICALLY_SUSPECT" if energy_suspect else
                    "ENERGY_BEHAVIOR_ACCEPTABLE_FOR_SCREEN")
    if not np.isfinite(target["end_ETOTAL_N_mm"]) or target["minimum_stable_dt_s"] <= 0:
        raise RuntimeError("Partial trajectory contains nonfinite energy or invalid timestep")
    primary = "CROT03_EARLY_RECOIL_FAIL"
    completion = "USER_STOPPED_AFTER_DECISIVE_RECOIL"
    branch = "ROTATIONAL_DAMPING_DOES_NOT_REMOVE_COUPLED_RECOIL"
    mechanism = ("SHORT_WINDOW_ROCKING_COHERENT_AXIAL_RECOIL_PRESENT"
                 if short_window_rocking_coherent else "MECHANISM_INCONCLUSIVE")
    m = target["motion"]
    start, end = .03, float(m["t"][-1])
    t = np.r_[start, m["t"][(m["t"] > start) & (m["t"] < end)], end]
    s = np.interp(t, m["t"], m["s"])
    v = np.interp(t, m["t"], m["v_s"])
    angle = np.interp(t, m["t"], m["rocking_angle"])
    omega = np.interp(t, m["t"], m["omega_rock"])
    partial_c4 = {"cycle": "C4 partial", "start_ms": 30.0, "end_ms": end*1000,
                  "delta_s_mm": float(s[-1]-s[0]),
                  "mean_v_s_mm_s": float((s[-1]-s[0])/(end-start)),
                  "end_v_s_mm_s": float(v[-1]), "minimum_v_s_mm_s": float(v.min()),
                  "max_backtrack_mm": float((np.maximum.accumulate(s)-s).max()),
                  "rocking_half_range_deg": float(np.ptp(angle)/2),
                  "rocking_fundamental_amplitude_deg": None,
                  "phase_lag_deg": None, "theta_offset_deg": None,
                  "sinusoid_residual_rms_over_amplitude": None,
                  "omega_peak_abs_rad_s": float(np.max(np.abs(omega)))}
    online = np.loadtxt(target["case"] / "online_motion.csv", delimiter=",", skiprows=1)
    if target["max_backtrack_mm"] <= .010 or target_rows[2]["delta_s_mm"] >= 0:
        raise RuntimeError("User-stop classification requires decisive C3 recoil")
    result = {"classification": primary, "completion": completion,
              "branch_conclusion": branch, "mechanism": mechanism,
              "energy_classification": energy_class,
              "energy_suspect_rule": "0.3x partial |ETOTAL| or ALLPW within one order of magnitude of completed 1x",
              "rocking_C2_C3_coherent_provisional": short_window_rocking_coherent,
              "rocking_C2_C5_stable": None, "axial_translation_stable": False,
              "partial_C4_ODB": partial_c4,
              "last_online_time_ms": float(online[-1, 0]*1000),
              "last_online_max_backtrack_mm": float(online[-1, 3]),
              "comparison": [{k: v for k, v in case.items() if k in (
                  "label", "name", "factor", "cycles", "cycle_energy", "wall",
                  "max_backtrack_mm", "amplitude_spread_C2_end_fraction",
                  "phase_spread_C2_end_deg", "end_ETOTAL_N_mm", "end_ALLPW_N_mm",
                  "end_time_ms", "minimum_stable_dt_s")}
                  for case in cases]}
    ident = target["identity"]
    ident.update({"status": completion, "classification": primary,
                  "branch_conclusion": branch, "restart_read": False,
                  "last_ODB_RP_time_ms": target["end_time_ms"],
                  "last_online_time_ms": result["last_online_time_ms"]})
    (target["case"] / "case_identity.json").write_text(json.dumps(ident, indent=2)+"\n")
    metrics = ROOT / f"{PREFIX}_METRICS.json"
    metrics.write_text(json.dumps(result, indent=2) + "\n")
    energy_csv = ROOT / f"{PREFIX}_CYCLE_ENERGY.csv"
    energy_rows = [row for case in cases for row in case["cycle_energy"]]
    with energy_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(energy_rows[0]))
        writer.writeheader()
        writer.writerows(energy_rows)
    chart = ROOT / f"{PREFIX}.png"
    animation = ROOT / f"{PREFIX}.gif"
    target_animation = ROOT / "F100_G2P20_CROT03_FULLCEL50.gif"
    comparison_plot(cases, chart)
    comparison_gif(cases, animation)
    target_gif(target, target_animation)
    report = ROOT / f"{PREFIX}_REPORT.md"
    lines = ["# F100/G2P20 rotational damping, FULL CEL", "",
             f"**{primary} / {completion}**. Branch: **{branch}**. Energy: **{energy_class}**.", "",
             f"Abaqus terminated on user request after completed C3 net recoil {target_rows[2]['delta_s_mm']:+.6f} mm and whole-run online backtrack {result['last_online_max_backtrack_mm']:.6f} mm, above the predeclared 0.010 mm gate. This is not a completed 50-ms solve or a spontaneous solver failure. ODB RP history ends at {target['end_time_ms']:.6f} ms; the online motion file ends at {result['last_online_time_ms']:.6f} ms. C4 is partial and C5 unavailable.", "",
             "0x was separately stopped by the user at 40.75 ms. The 1x and 3x cases completed 50 ms. All cases began at t=0 with one Explicit step and no restart continuation.", "",
             "## 0.3x cycle evidence", "",
             "| Cycle | delta s (mm) | mean v (mm/s) | end v (mm/s) | min v (mm/s) | max backtrack (mm) | fundamental amp (deg) | half range (deg) | phase lag (deg) | theta offset (deg) | residual / amp | omega peak (rad/s) |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in target_rows:
        lines.append(f"| C{row['cycle']} | {row['delta_s_mm']:+.6f} | {row['mean_v_s_mm_s']:+.3f} | {row['end_v_s_mm_s']:+.3f} | {row['minimum_v_s_mm_s']:+.3f} | {row['max_backtrack_mm']:.6f} | {row['rocking_fundamental_amplitude_deg']:.3f} | {row['rocking_half_range_deg']:.3f} | {row['phase_lag_deg']:+.2f} | {row['theta_offset_deg']:+.2f} | {row['sinusoid_residual_rms_over_amplitude']:.3f} | {row['omega_peak_abs_rad_s']:.2f} |")
    lines += [f"| C4 partial ({partial_c4['end_ms']:.3f} ms) | {partial_c4['delta_s_mm']:+.6f} | {partial_c4['mean_v_s_mm_s']:+.3f} | {partial_c4['end_v_s_mm_s']:+.3f} | {partial_c4['minimum_v_s_mm_s']:+.3f} | {partial_c4['max_backtrack_mm']:.6f} | NA | {partial_c4['rocking_half_range_deg']:.3f} | NA | NA | NA | {partial_c4['omega_peak_abs_rad_s']:.2f} |",
              "| C5 unavailable | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |",
              "", "## Four-level comparison", "",
             "| Metric | 0x | 0.3x | 1x | 3x |", "|---|---:|---:|---:|---:|"]
    def cells(fn, fmt=".4f"):
        return " | ".join("NA" if (v := fn(c)) is None else format(v, fmt) for c in cases)
    for cycle in range(1, 6):
        lines.append(f"| C{cycle} delta s (mm) | {cells(lambda c: next((r['delta_s_mm'] for r in c['cycles'] if r['cycle']==cycle), None), '+.6f')} |")
    lines += [f"| Max backtrack (mm) | {cells(lambda c: c['max_backtrack_mm'], '.6f')} |",
              f"| Rocking amplitude spread C2-last complete (%) | {cells(lambda c: 100*c['amplitude_spread_C2_end_fraction'], '.2f')} |",
              f"| Phase spread C2-last complete (deg) | {cells(lambda c: c['phase_spread_C2_end_deg'], '.2f')} |",
              f"| Minimum wall gap (mm) | {cells(lambda c: c['wall']['minimum_signed_gap_mm'], '.6f')} |",
              f"| Negative-gap episodes | {cells(lambda c: c['wall']['negative_gap_episodes'], 'd')} |",
              f"| Negative-gap duration (ms) | {cells(lambda c: c['wall']['negative_gap_duration_ms'], '.2f')} |",
              f"| ETOTAL at last available time (N mm) | {cells(lambda c: c['end_ETOTAL_N_mm'], '+.5g')} |",
              f"| ALLPW at last available time (N mm) | {cells(lambda c: c['end_ALLPW_N_mm'], '.5g')} |",
              f"| Last available time (ms) | {cells(lambda c: c['end_time_ms'], '.2f')} |",
              "", "0.3x rocking stability over C2-C5 cannot be established from the stopped run. C2-C3 coherence is provisional and does not explain away axial recoil. Negative wall gap is a 0.05-ms sampled geometric overlap; whole General Contact includes robot-fluid interaction and is not a pair-isolated robot-wall force. Cycle-wise ETOTAL, ALLPW, ALLFD, ALLIE, ALLKE and minimum stable timestep are in the accompanying CSV. Partial case endpoints are not equal-duration comparisons.",
              "", "Energy is marked suspect because the partial 0.3x |ETOTAL| or ALLPW is within one order of magnitude of the completed 1x case. This flags large progressive drift, without identifying its cause. The 0x energy is available only through 40.75 ms.",
              "", "Rotational-damping tuning stops here. The next investigation is a separate no-fluid control with real robot-wall contact isolated; it was not started in this branch."]
    report.write_text("\n".join(lines) + "\n")
    archive = ROOT / f"{PREFIX}_GIFS.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
        for path in (target_animation, animation, chart, report, metrics, energy_csv):
            z.write(path, path.name)
    print(json.dumps({"classification": primary, "completion": completion,
                      "branch_conclusion": branch, "mechanism": mechanism,
                      "energy": energy_class, "archive": str(archive)}, indent=2))


if __name__ == "__main__":
    main()
