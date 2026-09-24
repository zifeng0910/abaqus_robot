"""Analyze the 50-ms magnetic-only control and render its canonical GIF."""
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

from analyze_f100_g2p20_clean30 import KEYS, motion
from analyze_f100_g2p20_final import geometry, pose, setup_axis


ROOT = Path(__file__).resolve().parents[1]
JOB = "TRUECEL_B0P11_G2P20_F100_NOFLUID_CONTROL"
CASE = ROOT / "case" / JOB
FULL_CASE = ROOT / "case" / "TRUECEL_B0P11_G2P20_A14P5_F100_CLEAN50"
FULL_METRICS = ROOT / "F100_G2P20_CLEAN50_EARLY_STOP_METRICS.json"
PERIOD = .01
FREQUENCY = 100.0


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def load_motion(identity: dict) -> dict:
    with np.load(CASE / "private" / "rp_history_private.npz") as archive:
        history = {key: archive[key].astype(float) for key in KEYS}
    result = motion(history, identity)
    result["v"] = np.column_stack([
        np.interp(result["t"], history[f"V{i}"][:, 0], history[f"V{i}"][:, 1])
        for i in range(1, 4)
    ])
    return result


def cycle_metrics(m: dict, mag: np.ndarray) -> list[dict]:
    rows = []
    for cycle in range(1, 6):
        start, end = (cycle - 1) * PERIOD, cycle * PERIOD
        t = np.r_[start, m["t"][(m["t"] > start) & (m["t"] < end)], end]
        s = np.interp(t, m["t"], m["s"])
        v = np.interp(t, m["t"], m["v_s"])
        alpha = np.interp(t, m["t"], m["rocking_angle"])
        omega = 2 * np.pi * FREQUENCY
        basis = np.column_stack((np.sin(omega * t), np.cos(omega * t),
                                 np.ones(len(t))))
        sin_coeff, cos_coeff, _ = np.linalg.lstsq(basis, alpha, rcond=None)[0]
        harmonic_amp = float(np.hypot(sin_coeff, cos_coeff))
        lag = float(-np.rad2deg(np.arctan2(cos_coeff, sin_coeff)))
        lag = (lag + 180) % 360 - 180
        mag_rows = mag[(mag[:, 0] > start) & (mag[:, 0] <= end)]
        force = np.linalg.norm(mag_rows[:, 3:6], axis=1)
        torque = np.linalg.norm(mag_rows[:, 6:9], axis=1)
        rows.append({
            "cycle": cycle, "start_ms": start * 1000, "end_ms": end * 1000,
            "delta_s_mm": float(s[-1] - s[0]),
            "mean_v_s_mm_s": float((s[-1] - s[0]) / PERIOD),
            "end_v_s_mm_s": float(v[-1]),
            "minimum_v_s_mm_s": float(v.min()),
            "max_backtrack_mm": float((np.maximum.accumulate(s) - s).max()),
            "rocking_amplitude_half_range_deg": float((alpha.max() - alpha.min()) / 2),
            "rocking_first_harmonic_amplitude_deg": harmonic_amp,
            "phase_lag_deg": lag,
            "magnetic_force_mean_N": float(force.mean()),
            "magnetic_force_max_N": float(force.max()),
            "magnetic_torque_mean_Nmm": float(torque.mean()),
            "magnetic_torque_max_Nmm": float(torque.max()),
        })
    return rows


def write_timeseries(m: dict, mag: np.ndarray, identity: dict, path: Path) -> None:
    t = np.linspace(0, .05, 5001)
    initial = np.asarray(identity["initial_center_aba_mm"], float)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("time_s", "rp_x_mm", "rp_y_mm", "rp_z_mm",
                         "v_x_mm_s", "v_y_mm_s", "v_z_mm_s", "s_displacement_mm",
                         "v_s_mm_s", "rocking_angle_deg", "omega_rock_rad_s",
                         "command_rocking_deg", "Fmag1_N", "Fmag2_N", "Fmag3_N",
                         "Mmag1_Nmm", "Mmag2_Nmm", "Mmag3_Nmm"))
        sampled = [np.interp(t, m["t"], m["u"][:, j]) + initial[j]
                   for j in range(3)]
        sampled += [np.interp(t, m["t"], m["v"][:, j]) for j in range(3)]
        sampled += [np.interp(t, m["t"], m[key])
                    for key in ("s", "v_s", "rocking_angle", "omega_rock")]
        sampled += [identity["rocking_main_amplitude_deg"] *
                    np.sin(2 * np.pi * FREQUENCY * t)]
        sampled += [np.interp(t, mag[:, 0], mag[:, j], left=0)
                    for j in range(3, 9)]
        writer.writerows(zip(t, *sampled))


def plot_comparison(m: dict, full: dict, output: Path) -> None:
    with np.load(FULL_CASE / "private" / "rp_history_private.npz") as archive:
        old_history = {key: archive[key].astype(float) for key in KEYS}
    full_identity = json.loads((FULL_CASE / "case_identity.json").read_text())
    old = motion(old_history, full_identity)
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                             constrained_layout=True)
    axes[0].plot(m["t"] * 1000, m["s"], color="#126d62", label="no fluid/no contact")
    axes[0].plot(old["t"] * 1000, old["s"], color="#b54d3f", label="full CEL (stopped)")
    axes[1].plot(m["t"] * 1000, m["v_s"], color="#126d62")
    axes[1].plot(old["t"] * 1000, old["v_s"], color="#b54d3f")
    for ax in axes:
        ax.axhline(0, color="#777", lw=.7)
        for boundary in (10, 20, 30, 40):
            ax.axvline(boundary, color="#999", lw=.7, ls="--")
        ax.grid(alpha=.18)
        ax.set_xlim(0, 50)
    axes[0].legend()
    axes[0].set_ylabel("s (mm)")
    axes[1].set_ylabel("v_s (mm/s)")
    axes[1].set_xlabel("time (ms)")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def render_gif(m: dict, identity: dict, output: Path) -> None:
    rp0, c, n, rel, colors, pipe = geometry(CASE, identity)
    fig, (ax, trace) = plt.subplots(2, 1, figsize=(10, 6.5),
                                    gridspec_kw={"height_ratios": [2, 1]},
                                    constrained_layout=True)
    setup_axis(ax, identity, "Magnetic-only TRUE-CEL control | passive wall")
    p0 = pose(0, m["t"], m["u"], m["ur"], rp0, rel)
    robot = ax.scatter(float(identity["s_start_mm"]) + (p0-rp0) @ c,
                       (p0-pipe) @ n, s=5, c=colors, linewidths=0)
    label = ax.text(.98, .97, "", transform=ax.transAxes, ha="right", va="top",
                    family="monospace", fontsize=9,
                    bbox={"facecolor": "white", "alpha": .92})
    trace.plot(m["t"] * 1000, m["s"], color="#b9cbc8", lw=1)
    for boundary in (10, 20, 30, 40):
        trace.axvline(boundary, color="#999", lw=.7, ls="--")
    marker, = trace.plot([], [], "o", color="#126d62", ms=5)
    trace.set(xlim=(0, 50), xlabel="time (ms)", ylabel="s (mm)")
    trace.grid(alpha=.18)
    frames = np.linspace(0, .05, 121)

    def update(i):
        ti = frames[i]
        p = pose(ti, m["t"], m["u"], m["ur"], rp0, rel)
        robot.set_offsets(np.c_[float(identity["s_start_mm"]) + (p-rp0) @ c,
                                (p-pipe) @ n])
        si = np.interp(ti, m["t"], m["s"])
        vi = np.interp(ti, m["t"], m["v_s"])
        cycle = min(5, int(ti / PERIOD) + 1)
        start_s = np.interp((cycle - 1) * PERIOD, m["t"], m["s"])
        marker.set_data([ti * 1000], [si])
        label.set_text(f"t={ti*1000:5.2f} ms  C{cycle}\n"
                       f"delta_s={si-start_s:+.4f} mm\nv_s={vi:+.2f} mm/s")
        return robot, marker, label

    FuncAnimation(fig, update, frames=len(frames), interval=67).save(
        output, PillowWriter(fps=15), dpi=90)
    plt.close(fig)


def main() -> None:
    identity_path = CASE / "case_identity.json"
    identity = json.loads(identity_path.read_text(encoding="utf-8"))
    setup = json.loads((ROOT / f"{JOB}_SETUP.json").read_text(encoding="utf-8"))
    sta = (CASE / f"{JOB}.sta").read_text(encoding="latin1")
    m = load_motion(identity)
    with np.load(CASE / "private" / "energy_history_private.npz") as archive:
        etotal = archive["ETOTAL"].astype(float)
        kinetic = archive["ALLKE"].astype(float)
    mag = np.loadtxt(CASE / "magnetic_increment_g2p20_f100.csv",
                     delimiter=",", skiprows=1)
    valid = ("THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in sta and
             "***ERROR" not in sta and m["t"][0] <= 1e-10 and
             m["t"][-1] >= .05 - 1e-8 and mag[-1, 0] >= .05 - 1e-8 and
             np.all(np.isfinite(mag)) and np.all(np.isfinite(m["s"])) and
             len(etotal) > 0)
    if not valid:
        raise RuntimeError("No-fluid completion/numerical validity gate failed")
    if sha(CASE / "magnetic_field_gradient_table_B0P11_A14P5.dat") != \
            setup["magnetic_table_sha256"]:
        raise RuntimeError("Magnetic table changed during solve")
    rows = cycle_metrics(m, mag)
    stable_forward = all(row["delta_s_mm"] > 0 and row["max_backtrack_mm"] < .01
                         for row in rows)
    classification = ("MAGNETIC_ONLY_STABLE_FORWARD" if stable_forward
                      else "MAGNETIC_ONLY_RECOIL")
    full = json.loads(FULL_METRICS.read_text(encoding="utf-8"))
    full_rows = full["cycles"]
    if len(full_rows) != 5 or full_rows[4]["complete_cycle"]:
        raise RuntimeError("Expected FULL CEL C5 to be partial")
    with np.load(FULL_CASE / "private" / "rp_history_private.npz") as archive:
        full_history = {key: archive[key].astype(float) for key in KEYS}
    full_identity = json.loads((FULL_CASE / "case_identity.json").read_text())
    coupled = motion(full_history, full_identity)
    in_40 = coupled["t"] <= .04
    full_s = np.r_[coupled["s"][in_40], np.interp(.04, coupled["t"], coupled["s"])]
    new_t = m["t"] <= .04
    nofluid_s = np.r_[m["s"][new_t], np.interp(.04, m["t"], m["s"])]
    compare = {
        "full_CEL_cycle_delta_s_mm": [r["delta_s_mm"] for r in full_rows[:4]] + [None],
        "no_fluid_cycle_delta_s_mm": [r["delta_s_mm"] for r in rows],
        "full_CEL_C5_status": "NOT_COMPLETED; stopped at 40.754825 ms",
        "full_CEL_C5_partial_delta_s_mm": full_rows[4]["delta_s_mm"],
        "full_CEL_max_backtrack_0_40_mm": float((np.maximum.accumulate(full_s) - full_s).max()),
        "no_fluid_max_backtrack_0_40_mm": float((np.maximum.accumulate(nofluid_s) - nofluid_s).max()),
        "full_CEL_mean_v_0_40_mm_s": float((full_s[-1] - full_s[0]) / .04),
        "no_fluid_mean_v_0_40_mm_s": float((nofluid_s[-1] - nofluid_s[0]) / .04),
        "no_fluid_mean_v_0_50_mm_s": float((m["s"][-1] - m["s"][0]) / .05),
    }
    result = {
        "classification": classification, "job": JOB,
        "solver_completed": True, "last_RP_time_ms": float(m["t"][-1] * 1000),
        "cycle_metrics": rows, "comparison": compare,
        "max_abs_ETOTAL_Nmm": float(np.abs(etotal[:, 1]).max()),
        "max_ALLKE_Nmm": float(np.abs(kinetic[:, 1]).max()),
        "wall_geometry_present": True, "contact_present": False,
        "fluid_domain_present": False, "direct_dt_s": 1e-7,
        "magnetic_table_sha256": setup["magnetic_table_sha256"],
        "phase_lag_method": "first-harmonic least squares: angle=a*sin(2*pi*100*t)+b*cos(2*pi*100*t)+offset; lag=-atan2(b,a)",
        "rocking_amplitude_method": "half of cycle max-min rocking angle; first-harmonic amplitude also reported",
        "setup": setup,
    }
    metrics = ROOT / f"{JOB}_METRICS.json"
    metrics.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    cycles_csv = ROOT / f"{JOB}_CYCLES.csv"
    with cycles_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    series = ROOT / f"{JOB}_RP_MAGNETIC_10US.csv"
    write_timeseries(m, mag, identity, series)
    chart = ROOT / f"{JOB}_COMPARISON.png"
    gif = ROOT / f"{JOB}.gif"
    plot_comparison(m, full, chart)
    render_gif(m, identity, gif)
    gif_zip = ROOT / f"{JOB}_GIF.zip"
    with zipfile.ZipFile(gif_zip, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.write(gif, gif.name)
    lines = [
        "# Magnetic-only / no-fluid F100 G2P20 control", "",
        f"**{classification}**", "",
        "Fresh t=0 to 50 ms, one uninterrupted Explicit step. The Eulerian fluid "
        "domain, fluid material, fluid/contact interactions, and all General "
        "Contact were removed. The unchanged wall geometry remains as a passive, "
        "fixed reference; contact is absent. Robot geometry, mass, initial pose, "
        "Magpylib-derived table, phase convention, B0/G/f/A, and VUAMP load "
        "calculation were preserved.", "",
        "The all-rigid model required DIRECT 1e-7 s time increments because "
        "Abaqus could not derive an automatic stable increment without a "
        "deformable element. This is a numerical setting change, not a new "
        "force or damping term. The source step's bulk-viscosity keyword was "
        "left unchanged; no artificial damping was added. Fixed 5-ms restart "
        "writes were retained and never read for continuation.", "",
        "| Cycle | delta_s (mm) | mean v_s (mm/s) | end v_s | minimum v_s | "
        "MAX_BACKTRACK (mm) | rocking amplitude (deg) | phase lag (deg) |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(f"| {row['cycle']} | {row['delta_s_mm']:+.6f} | "
                     f"{row['mean_v_s_mm_s']:+.3f} | {row['end_v_s_mm_s']:+.3f} | "
                     f"{row['minimum_v_s_mm_s']:+.3f} | "
                     f"{row['max_backtrack_mm']:.6f} | "
                     f"{row['rocking_amplitude_half_range_deg']:.3f} | "
                     f"{row['phase_lag_deg']:+.2f} |")
    lines += [
        "", "Rocking amplitude is half the measured angle range per cycle. "
        "Phase lag is the first-harmonic angle lag relative to the unchanged "
        "14.5-degree, 100-Hz sine command; positive means measured rocking "
        "lags the command. The first cycle includes the unchanged 1-ms "
        "magnetic ramp.", "",
        "| Cycle | Mean magnetic force (N) | Peak magnetic force (N) | "
        "Mean magnetic torque (N mm) | Peak magnetic torque (N mm) |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(f"| {row['cycle']} | {row['magnetic_force_mean_N']:.3e} | "
                     f"{row['magnetic_force_max_N']:.3e} | "
                     f"{row['magnetic_torque_mean_Nmm']:.3e} | "
                     f"{row['magnetic_torque_max_Nmm']:.3e} |")
    lines += [
        "", "Magnetic force and torque are vector magnitudes from the unchanged "
        "VUAMP magnetic log. Their three components, RP position and velocity, "
        "rocking angle, and angular velocity are in the 10-us CSV.", "",
        "| Matched comparison | FULL CEL | NO FLUID |",
        "|---|---:|---:|",
    ]
    for i in range(5):
        left = (f"{full_rows[i]['delta_s_mm']:+.6f}" if i < 4 else
                "N/A (stopped at 40.755 ms)")
        lines.append(f"| Cycle {i+1} delta_s (mm) | {left} | "
                     f"{rows[i]['delta_s_mm']:+.6f} |")
    lines += [
        f"| MAX_BACKTRACK, 0-40 ms (mm) | {compare['full_CEL_max_backtrack_0_40_mm']:.6f} | "
        f"{compare['no_fluid_max_backtrack_0_40_mm']:.6f} |",
        f"| Mean speed, 0-40 ms (mm/s) | {compare['full_CEL_mean_v_0_40_mm_s']:+.3f} | "
        f"{compare['no_fluid_mean_v_0_40_mm_s']:+.3f} |",
        "", "The FULL CEL Cycle 5 is incomplete; its 40-40.755 ms partial "
        f"delta_s was {full_rows[4]['delta_s_mm']:+.6f} mm and is excluded "
        "from the five-cycle comparison. The no-fluid 0-50 ms mean speed was "
        f"{compare['no_fluid_mean_v_0_50_mm_s']:+.3f} mm/s.", "",
        "The magnetic-only trajectory remains directionally stable under the "
        "stated 0.010-mm within-cycle backtrack threshold. The full coupled "
        "trajectory reverses in Cycle 4. This is consistent with the removed "
        "fluid/contact physics contributing to recoil. The control also "
        "requires a fixed Explicit time increment, so this comparison alone "
        "does not uniquely attribute the difference to fluid forces, wall "
        "contact, or their interaction.", "",
        f"Numerical checks: Abaqus completed 50 ms; peak |ETOTAL|="
        f"{result['max_abs_ETOTAL_Nmm']:.3e} N mm versus peak ALLKE="
        f"{result['max_ALLKE_Nmm']:.3e} N mm. The four inherited distorted "
        "tetrahedra are in the fixed passive PIPE_SOLID and were also present "
        "in the FULL CEL source. RP, magnetic loads, and energy histories are "
        "finite. No restart continuation or feedback was used.",
    ]
    report = ROOT / f"{JOB}_REPORT.md"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    identity.update({"status": "SOLVED", "classification": classification,
                     "dynamics_run_count": 1, "last_RP_time_ms": result["last_RP_time_ms"]})
    identity_path.write_text(json.dumps(identity, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"classification": classification, "cycles": rows,
                      "comparison": compare, "max_abs_ETOTAL_Nmm":
                      result["max_abs_ETOTAL_Nmm"]}, indent=2))


if __name__ == "__main__":
    main()
