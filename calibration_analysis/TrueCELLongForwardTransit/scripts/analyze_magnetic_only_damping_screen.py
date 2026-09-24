"""Compare 50-ms rotational-damping runs against a passive wall geometry."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
from scipy.spatial import ConvexHull

from analyze_f100_g2p20_clean30 import KEYS, motion
from analyze_f100_g2p20_final import geometry, pose, setup_axis


ROOT = Path(__file__).resolve().parents[1]
PREFIX = "TRUECEL_MAGNETIC_ONLY_DAMPING"
PERIOD = .01
FREQUENCY = 100.0


def nodes(deck: str, part: str) -> np.ndarray:
    match = re.search(r"(?ms)^\*Part, name=" + re.escape(part) +
                      r"\s*$.*?^\*End Part\s*$", deck)
    body = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", match.group()).group(1)
    return np.asarray([[float(x) for x in line.split(",")[1:4]]
                       for line in body.splitlines() if line.strip()])


def wall_facets(case: Path, ident: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    deck = (case / f"{case.name}.inp").read_text(encoding="latin1")
    wall = nodes(deck, "Pipe_WALL_HELPER")
    c = np.asarray(ident["canonical_plus_s_axis_aba"])
    n = np.asarray(ident["n_routeA_aba"])
    b = np.asarray(ident["b_routeA_aba"])
    rp0 = np.asarray(ident["initial_center_aba_mm"])
    pipe = rp0 - ident["s_start_mm"] * c - ident["radial_offset_n_mm"] * n
    s = (wall - pipe) @ c
    ring = wall[abs(s-s.min()) < 1e-6]
    if len(ring) < 6:
        raise RuntimeError("Passive wall ring could not be identified")
    cross = np.column_stack(((ring-pipe) @ n, (ring-pipe) @ b))
    planes = ConvexHull(cross).equations
    return pipe, np.column_stack((n, b)), planes


def wall_metrics(case: Path, ident: dict, m: dict) -> dict:
    rp0, _, _, rel, _, _ = geometry(case, ident)
    rel = rel[ConvexHull(rel).vertices]
    pipe, basis, planes = wall_facets(case, ident)
    ts = np.linspace(0, .05, 1001)
    gaps = np.empty(len(ts))
    for i, ti in enumerate(ts):
        p = pose(ti, m["t"], m["u"], m["ur"], rp0, rel)
        cross = (p-pipe) @ basis
        gaps[i] = -np.max(cross @ planes[:, :2].T + planes[:, 2])
    outside = gaps < -1e-6
    events = int(outside[0]) + int(np.count_nonzero(outside[1:] & ~outside[:-1]))
    return {"signed_min_wall_gap_mm": float(gaps.min()),
            "virtual_wall_crossing_events": events,
            "virtual_wall_penetration_duration_ms": float(outside.sum() * .05),
            "wall_gap_sample_interval_ms": .05,
            "wall_gap_definition": "minimum signed gap of rigid robot convex hull to passive wall polygon; negative is geometric overlap, not contact force"}


def longest_true_duration(t: np.ndarray, mask: np.ndarray) -> float:
    starts = np.flatnonzero(mask & ~np.r_[False, mask[:-1]])
    ends = np.flatnonzero(mask & ~np.r_[mask[1:], False])
    return float(max((t[j] - t[i] for i, j in zip(starts, ends)), default=0.0))


def analyze_case(entry: dict) -> tuple[dict, dict, dict]:
    case = ROOT / "case" / entry["name"]
    ident = json.loads((case / "case_identity.json").read_text())
    sta = (case / f"{case.name}.sta").read_text(encoding="latin1")
    with np.load(case / "private" / "rp_history_private.npz") as archive:
        history = {key: archive[key].astype(float) for key in KEYS}
    m = motion(history, ident)
    valid = ("THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in sta and
             m["t"][0] <= 1e-10 and m["t"][-1] >= .05-1e-8 and
             all(np.all(np.isfinite(m[key])) for key in
                 ("u", "ur", "s", "v_s", "rocking_angle", "omega_rock")))
    if not valid:
        raise RuntimeError(f"Incomplete or nonfinite history: {case.name}")
    rows = []
    for cycle in range(1, 6):
        t0, t1 = (cycle-1)*PERIOD, cycle*PERIOD
        t = np.r_[t0, m["t"][(m["t"] > t0) & (m["t"] < t1)], t1]
        theta = np.interp(t, m["t"], m["rocking_angle"])
        omega = np.interp(t, m["t"], m["omega_rock"])
        coeff = np.column_stack((np.sin(2*np.pi*FREQUENCY*t),
                                 np.cos(2*np.pi*FREQUENCY*t), np.ones(len(t))))
        a, b, offset = np.linalg.lstsq(coeff, theta, rcond=None)[0]
        fitted = coeff @ np.asarray([a, b, offset])
        amplitude = float(np.hypot(a, b))
        phase = float((-np.degrees(np.arctan2(b, a))+180) % 360-180)
        rows.append({"cycle": cycle,
                     "rocking_half_range_deg": float(np.ptp(theta)/2),
                     "rocking_fundamental_amplitude_deg": amplitude,
                     "phase_lag_deg": phase,
                     "theta_offset_deg": float(offset),
                     "sinusoid_residual_rms_over_amplitude": float(
                         np.sqrt(np.mean((theta-fitted)**2))/max(amplitude, 1e-12)),
                     "omega_peak_abs_rad_s": float(abs(omega).max()),
                     "theta_peak_abs_deg": float(abs(theta).max())})
    amps = np.asarray([r["rocking_fundamental_amplitude_deg"] for r in rows])
    phases = np.asarray([r["phase_lag_deg"] for r in rows])
    steady = slice(1, None)  # exclude the unchanged 1-ms magnetic start-up ramp
    amp_spread = float(np.ptp(amps[steady])/max(np.mean(amps[steady]), 1e-12))
    phase_spread = float(np.ptp(np.rad2deg(np.unwrap(np.deg2rad(phases[steady])))))
    transverse = abs(m["rocking_angle"]) >= 70
    hover = transverse & (abs(m["omega_rock"]) <= 20)
    wall = wall_metrics(case, ident, m)
    gates = {
        "amplitude_spread_C2_C5_below_10pct": amp_spread < .10,
        "phase_spread_C2_C5_below_10deg": phase_spread < 10,
        "sinusoid_residual_C2_C5_below_20pct": all(
            r["sinusoid_residual_rms_over_amplitude"] < .20 for r in rows[1:]),
        "no_virtual_wall_crossings": wall["virtual_wall_crossing_events"] == 0,
        "no_transverse_hover_over_0p5ms": longest_true_duration(m["t"], hover) < .0005,
    }
    result = {"case": case.name, "c_rot_N_mm_s": entry["c_rot_N_mm_s"],
              "factor": entry["factor"], "ratio_to_bulk_water": entry["ratio_to_bulk_water"],
              "cycles": rows, "amplitude_spread_C2_C5": amp_spread,
              "phase_spread_C2_C5_deg": phase_spread,
              "transverse_time_ms": float(np.trapz(transverse.astype(float), m["t"])*1000),
              "transverse_hover_time_ms": float(np.trapz(hover.astype(float), m["t"])*1000),
              "longest_transverse_hover_ms": longest_true_duration(m["t"], hover)*1000,
              "max_abs_omega_rad_s": float(abs(m["omega_rock"]).max()),
              "wall": wall, "gates": gates, "all_dynamics_gates_pass": all(gates.values())}
    return result, ident, m


def make_plots(results: list[dict], motions: list[dict]) -> None:
    colors = ("#16806d", "#d19628", "#ac5047", "#425b95")
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
    for result, m, color in zip(results, motions, colors):
        name = f"{result['factor']:g}x"
        cycles = [r["cycle"] for r in result["cycles"]]
        axes[0, 0].plot(cycles, [r["rocking_fundamental_amplitude_deg"] for r in result["cycles"]], "o-", label=name, color=color)
        axes[0, 1].plot(cycles, [r["phase_lag_deg"] for r in result["cycles"]], "o-", label=name, color=color)
        axes[1, 0].plot(m["t"]*1000, m["rocking_angle"], lw=.9, label=name, color=color)
        axes[1, 1].plot(m["t"]*1000, m["omega_rock"], lw=.9, label=name, color=color)
    axes[0, 0].set(xlabel="cycle", ylabel="first-harmonic rocking amplitude (deg)")
    axes[0, 1].set(xlabel="cycle", ylabel="phase lag (deg)")
    axes[1, 0].set(xlabel="time (ms)", ylabel="theta (deg)")
    axes[1, 1].set(xlabel="time (ms)", ylabel="omega_rock (rad/s)")
    for ax in axes.flat:
        ax.grid(alpha=.2)
    axes[0, 0].legend(title="c / (I 2pi f)")
    fig.savefig(ROOT / f"{PREFIX}_ROCKING_SCREEN.png", dpi=160)
    plt.close(fig)


def make_gif(result: dict, ident: dict, m: dict) -> None:
    case = ROOT / "case" / result["case"]
    rp0, c, n, rel, colors, pipe = geometry(case, ident)
    fig, (ax, trace) = plt.subplots(2, 1, figsize=(10, 6.5),
                                    gridspec_kw={"height_ratios": [2, 1]},
                                    constrained_layout=True)
    setup_axis(ax, ident, f"Magnetic-only rotational damping | {result['factor']:g}x")
    p0 = pose(0, m["t"], m["u"], m["ur"], rp0, rel)
    robot = ax.scatter(ident["s_start_mm"]+(p0-rp0)@c, (p0-pipe)@n,
                       s=5, c=colors, linewidths=0)
    label = ax.text(.98, .97, "", transform=ax.transAxes, ha="right", va="top",
                    family="monospace", fontsize=9,
                    bbox={"facecolor": "white", "alpha": .92})
    trace.plot(m["t"]*1000, m["rocking_angle"], color="#b5c8c4", lw=1)
    command = 14.5*np.sin(2*np.pi*FREQUENCY*m["t"])
    trace.plot(m["t"]*1000, command, color="#c77b4c", lw=.8, ls="--")
    point, = trace.plot([], [], "o", color="#16806d", ms=5)
    trace.set(xlim=(0, 50), xlabel="time (ms)", ylabel="theta / command (deg)")
    trace.grid(alpha=.2)
    times = np.linspace(0, .05, 121)

    def update(i):
        ti = times[i]
        p = pose(ti, m["t"], m["u"], m["ur"], rp0, rel)
        robot.set_offsets(np.c_[ident["s_start_mm"]+(p-rp0)@c, (p-pipe)@n])
        theta = float(np.interp(ti, m["t"], m["rocking_angle"]))
        omega = float(np.interp(ti, m["t"], m["omega_rock"]))
        cycle = min(5, int(ti/PERIOD)+1)
        phase = result["cycles"][cycle-1]["phase_lag_deg"]
        label.set_text(f"t={ti*1000:5.2f} ms  C{cycle}\n"
                       f"theta={theta:+.1f} deg  omega={omega:+.0f} rad/s\n"
                       f"C{cycle} lag={phase:+.1f} deg")
        point.set_data([ti*1000], [theta])
        return robot, point, label

    FuncAnimation(fig, update, frames=len(times), interval=67).save(
        ROOT / f"{PREFIX}.gif", PillowWriter(fps=15), dpi=90)
    plt.close(fig)


def main() -> None:
    estimate = json.loads((ROOT / f"{PREFIX}_OFFLINE_ESTIMATE.json").read_text())
    baseline, _, baseline_motion = analyze_case({
        "name": estimate["parent"], "c_rot_N_mm_s": 0.0,
        "factor": 0.0, "ratio_to_bulk_water": 0.0})
    analyzed = [analyze_case(entry) for entry in estimate["cases"]]
    results = [item[0] for item in analyzed]
    passing = [i for i, r in enumerate(results) if r["all_dynamics_gates_pass"]]
    selected = passing[0] if passing else min(
        range(len(results)), key=lambda i: (
            results[i]["wall"]["virtual_wall_crossing_events"],
            results[i]["amplitude_spread_C2_C5"],
            results[i]["phase_spread_C2_C5_deg"]))
    classification = ("MAGNETIC_ROCKING_LIMIT_CYCLE_ESTABLISHED" if passing else
                      "MAGNETIC_ROCKING_LIMIT_CYCLE_NOT_FOUND")
    output = {"classification": classification,
              "selected_visualization_case": results[selected]["case"],
              "gate_basis": "C2-C5 after 1-ms start-up ramp: amplitude spread <10%, phase spread <10 deg, harmonic residual <20%; no geometric wall crossing; no >=0.5-ms transverse hover (|theta|>=70 deg and |omega|<=20 rad/s)",
              "physical_plausibility_note": "Candidate c values are dynamic screening scales, 13-402 times the bulk-water slender-rod estimate. A passing dynamics gate alone does not calibrate the drag or prove hydrodynamic plausibility.",
              "offline_estimate": estimate, "baseline": baseline, "results": results}
    (ROOT / f"{PREFIX}_METRICS.json").write_text(json.dumps(output, indent=2)+"\n")
    with (ROOT / f"{PREFIX}_CYCLES.csv").open("w", newline="") as handle:
        fields = ["case", "factor", "c_rot_N_mm_s"] + list(results[0]["cycles"][0])
        writer = csv.DictWriter(handle, fields)
        writer.writeheader()
        for result in results:
            for row in result["cycles"]:
                writer.writerow({"case": result["case"], "factor": result["factor"],
                                 "c_rot_N_mm_s": result["c_rot_N_mm_s"], **row})
    with (ROOT / f"{PREFIX}_THETA_OMEGA_10US.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("case", "factor", "time_s", "theta_deg", "omega_rock_rad_s",
                         "command_theta_deg"))
        for result, m in [(baseline, baseline_motion)] + [
                (item[0], item[2]) for item in analyzed]:
            sample_t = np.linspace(0, .05, 5001)
            sample_theta = np.interp(sample_t, m["t"], m["rocking_angle"])
            sample_omega = np.interp(sample_t, m["t"], m["omega_rock"])
            writer.writerows((result["case"], result["factor"], float(ti),
                              float(theta), float(omega),
                              float(14.5*np.sin(2*np.pi*FREQUENCY*ti)))
                             for ti, theta, omega in zip(
                                 sample_t, sample_theta, sample_omega))
    make_plots(results, [item[2] for item in analyzed])
    make_gif(*analyzed[selected])
    lines = ["# Magnetic-only rotational damping screen", "", f"**{classification}**", "",
             "The passive wall has no contact interaction. Its crossing count is a geometric overlap diagnostic, not a simulated collision count.", "",
             "| c/(I 2pi f) | c (N mm s) | bulk-water multiple | amplitude spread C2-C5 | phase spread (deg) | virtual wall events | min gap (mm) | longest transverse hover (ms) | gates |",
             "|---:|---:|---:|---:|---:|---:|---:|---:|---|" ]
    for r in results:
        lines.append(f"| {r['factor']:g} | {r['c_rot_N_mm_s']:.3e} | {r['ratio_to_bulk_water']:.1f} | "
                     f"{r['amplitude_spread_C2_C5']*100:.1f}% | {r['phase_spread_C2_C5_deg']:.1f} | "
                     f"{r['wall']['virtual_wall_crossing_events']} | {r['wall']['signed_min_wall_gap_mm']:+.4f} | "
                     f"{r['longest_transverse_hover_ms']:.3f} | {'PASS' if r['all_dynamics_gates_pass'] else 'FAIL'} |")
    lines += ["", f"Undamped baseline: {baseline['wall']['virtual_wall_crossing_events']} virtual wall crossings; "
              f"minimum signed gap {baseline['wall']['signed_min_wall_gap_mm']:+.4f} mm.", "",
              output["gate_basis"], "", output["physical_plausibility_note"],
              f"The input-mesh mass integral is {estimate['robot_mass_mg']:.4f} mg, whereas the inherited case metadata records 9.2078 mg. The inertia and screening scale use the actual input mesh and density; the metadata discrepancy needs resolution before any calibrated hydrodynamic claim.", "",
              f"The GIF shows {results[selected]['case']}; it is the first passing candidate, or the smallest-wall-crossing diagnostic candidate if none passes.", "",
              "No translational damping, propulsion, position constraint, magnetic parameter change, fluid domain, or active wall contact was added. All cases start at t=0 and use one uninterrupted Explicit step. The magnetic-only torque log records magnetic torque before rotational damping is subtracted in VUAMP."]
    (ROOT / f"{PREFIX}_REPORT.md").write_text("\n".join(lines)+"\n")
    print(json.dumps({"classification": classification, "selected": selected,
                      "summary": [{"factor": r["factor"], "gates": r["gates"]}
                                  for r in results]}, indent=2))


if __name__ == "__main__":
    main()
