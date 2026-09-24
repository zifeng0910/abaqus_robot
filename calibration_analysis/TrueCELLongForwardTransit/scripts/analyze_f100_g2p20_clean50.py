"""Validate the uninterrupted 0-50 ms run and summarize five 100-Hz cycles."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_f100_g2p20_clean30 import KEYS, motion, stats


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "case" / "TRUECEL_B0P11_G2P20_A14P5_F100_FAST"
NEW = ROOT / "case" / "TRUECEL_B0P11_G2P20_A14P5_F100_CLEAN50"
PERIOD = .01


def load_rp(folder: Path) -> dict[str, np.ndarray]:
    with np.load(folder / "private" / "rp_history_private.npz") as z:
        return {key: z[key].astype(float) for key in KEYS}


def compare_contact(t: np.ndarray) -> dict:
    with np.load(BASE / "private" / "contact_history_private.npz") as a, \
            np.load(NEW / "private" / "contact_history_private.npz") as b:
        keys = sorted(k for k in set(a.files) & set(b.files)
                      if "ROBOT_SOLID" in k and any(f"|CFN{i} " in k for i in (1, 2, 3))
                      or ("ROBOT_SOLID" in k and "|CFNM " in k))
        if len(keys) != 4:
            raise RuntimeError(f"Expected four whole robot General Contact outputs: {keys}")
        return {k.split("|")[1]: stats(np.interp(t, b[k][:, 0], b[k][:, 1]) -
                                        np.interp(t, a[k][:, 0], a[k][:, 1]))
                for k in keys}


def compare_magnetic(t: np.ndarray) -> dict:
    a = np.loadtxt(BASE / "magnetic_increment_g2p20_f100.csv", delimiter=",", skiprows=1)
    b = np.loadtxt(NEW / "magnetic_increment_g2p20_f100.csv", delimiter=",", skiprows=1)
    out = {}
    for col, key in enumerate(("Fmag1_N", "Fmag2_N", "Fmag3_N",
                               "Mmag1_Nmm", "Mmag2_Nmm", "Mmag3_Nmm"), 3):
        out[key] = stats(np.interp(t, b[:, 0], b[:, col]) -
                         np.interp(t, a[:, 0], a[:, col]))
    return out


def cycles(m: dict) -> list[dict]:
    rows = []
    t, s, vs = m["t"], m["s"], m["v_s"]
    for cycle in range(1, 6):
        start, end = (cycle - 1)*PERIOD, cycle*PERIOD
        tc = np.r_[start, t[(t > start) & (t < end)], end]
        sc = np.interp(tc, t, s)
        vc = np.interp(tc, t, vs)
        rows.append({
            "cycle": cycle, "start_ms": start*1000, "end_ms": end*1000,
            "delta_s_mm": float(sc[-1] - sc[0]),
            "mean_v_s_mm_s": float((sc[-1] - sc[0])/PERIOD),
            "end_v_s_mm_s": float(vc[-1]),
            "minimum_v_s_mm_s": float(vc.min()),
            "max_backtrack_mm": float((np.maximum.accumulate(sc) - sc).max()),
        })
    return rows


def plot(m: dict, output: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11, 6), sharex=True, constrained_layout=True)
    axes[0].plot(m["t"]*1000, m["s"], color="#126b66", lw=1)
    axes[1].plot(m["t"]*1000, m["v_s"], color="#9a4b37", lw=.8)
    axes[0].set_ylabel("s (mm)")
    axes[1].set_ylabel("v_s (mm/s)")
    axes[1].set_xlabel("time (ms)")
    for ax in axes:
        for boundary in (10, 20, 30, 40):
            ax.axvline(boundary, color="#888", lw=.7, ls="--")
        ax.grid(alpha=.2)
        ax.set_xlim(0, 50)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    identity = json.loads((NEW / "case_identity.json").read_text(encoding="utf-8"))
    base_ident = json.loads((BASE / "case_identity.json").read_text(encoding="utf-8"))
    setup = json.loads((ROOT / f"{NEW.name}_SETUP.json").read_text(encoding="utf-8"))
    if identity["magnetic_table_sha256"] != base_ident["magnetic_table_sha256"]:
        raise RuntimeError("Magnetic table identity differs")
    sta = (NEW / f"{NEW.name}.sta").read_text(encoding="latin1")
    solver_completed = "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in sta
    actual_writes = []
    for line in sta.splitlines():
        if "Restart Number" in line and " at " in line:
            actual_writes.append(float(line.rsplit(" at ", 1)[1])*1000)
    expected_writes = list(range(5, 51, 5))
    baseline_h = load_rp(BASE)
    clean_h = load_rp(NEW)
    clean = motion(clean_h, identity)
    baseline = motion(baseline_h, base_ident)
    completed_to_50 = clean["t"][0] <= 1e-10 and clean["t"][-1] >= .05 - 1e-8
    valid_solver = (solver_completed and completed_to_50 and
                    len(actual_writes) == len(expected_writes) and
                    np.allclose(actual_writes, expected_writes, atol=.001, rtol=0))

    t = baseline["t"][baseline["t"] <= .02]
    t = t[t <= clean["t"][-1]]
    baseline_error = {key: stats(np.interp(t, clean["t"], clean[key]) -
                                 np.interp(t, baseline["t"], baseline[key]))
                      for key in ("s", "v_s", "rocking_angle", "omega_rock")}
    baseline_error["magnetic"] = compare_magnetic(t)
    baseline_error["whole_general_contact"] = compare_contact(t)
    gates = {
        "position_max_mm": baseline_error["s"]["max"] < .001,
        "v_s_p99_mm_s": baseline_error["v_s"]["p99"] < .1,
        "rocking_angle_max_deg": baseline_error["rocking_angle"]["max"] < .02,
        "omega_rock_p99_rad_s": baseline_error["omega_rock"]["p99"] < .5,
        "magnetic_force_max_N": max(baseline_error["magnetic"][k]["max"]
                                    for k in ("Fmag1_N", "Fmag2_N", "Fmag3_N")) < 1e-8,
        "magnetic_torque_max_Nmm": max(baseline_error["magnetic"][k]["max"]
                                       for k in ("Mmag1_Nmm", "Mmag2_Nmm", "Mmag3_Nmm")) < 1e-5,
        "contact_p99_N": max(v["p99"] for v in
                             baseline_error["whole_general_contact"].values()) < 1e-5,
    }
    baseline_reproduced = bool(t[-1] >= .02-1e-8 and all(gates.values()))
    rows = cycles(clean) if valid_solver else []
    if not valid_solver:
        classification = "CLEAN50_NUMERICALLY_INVALID"
    elif all(r["delta_s_mm"] > 0 and r["max_backtrack_mm"] < .01 for r in rows):
        classification = "CLEAN50_FORWARD_SUSTAINED"
    else:
        classification = "CLEAN50_RECOIL_CONFIRMED"

    result = {
        "classification": classification, "job": NEW.name,
        "solver_completed": solver_completed,
        "last_RP_time_ms": float(clean["t"][-1]*1000),
        "actual_restart_write_times_ms": actual_writes,
        "expected_restart_write_times_ms": expected_writes,
        "restart_read": False,
        "baseline_0_20_reproduced": baseline_reproduced,
        "baseline_0_20_error": baseline_error,
        "baseline_0_20_gates": gates,
        "cycle_metrics": rows,
        "old_restart_C3_status": "RESTART_CONTAMINATED_C3_RESULT",
        "old_restart_C3_delta_s_mm": -.023527,
        "clean_C3_net_recoil": rows[2]["delta_s_mm"] < 0 if rows else None,
        "setup": setup,
    }
    (ROOT / "F100_G2P20_CLEAN50_METRICS.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if rows:
        with (ROOT / "F100_G2P20_CLEAN50_CYCLES.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    plot(clean, ROOT / "F100_G2P20_CLEAN50_MOTION.png")
    write_report(result)
    identity.update({"status": "SOLVED" if valid_solver else "FAILED",
                     "dynamics_run_count": 1, "classification": classification,
                     "cpus": 1, "precision": "BOTH"})
    (NEW / "case_identity.json").write_text(
        json.dumps(identity, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"classification": classification,
                      "baseline_reproduced": baseline_reproduced,
                      "last_RP_time_ms": result["last_RP_time_ms"],
                      "cycles": rows}, indent=2))


def write_report(result: dict) -> None:
    lines = ["# F100/G2.20 uninterrupted CLEAN50", "",
             f"**{result['classification']}**", "",
             "One fresh-start 0-50 ms Explicit step, no restart read or "
             "cycle-boundary steps. B0=11 mT, G=2.20, f=100 Hz, "
             "A_main=14.5 deg, A_cross=2.5 deg, 44x20x20 CEL, "
             "c0=100000 mm/s, original contact/fluid/EOS, one CPU, "
             "double=both and no mass scaling. Physical recoil early stop "
             "was disabled; numerical timestep-collapse stop was retained.", "",
             "The fixed restart WRITE keyword is `number interval=10, time "
             "marks=YES`, giving writes at 5,10,...,50 ms. These states were "
             "never read for continuation. Actual write times (ms): "
             f"{result['actual_restart_write_times_ms']}.", "",
             f"Solver completed: {result['solver_completed']}; final RP time "
             f"{result['last_RP_time_ms']:.6f} ms. The original 0-20 ms history "
             "was compared to document any numerical-path divergence.", "",
             "| Baseline difference | RMS | p99 absolute | max absolute |",
             "|---|---:|---:|---:|"]
    for key in ("s", "v_s", "rocking_angle", "omega_rock"):
        x = result["baseline_0_20_error"][key]
        lines.append(f"| {key} | {x['rms']:.6g} | {x['p99']:.6g} | {x['max']:.6g} |")
    lines.extend(["", "The metrics JSON also records magnetic-load and whole "
                  "General Contact RMS/p99/max errors and every reproduction "
                  "gate. Baseline reproduction passed: "
                  f"{result['baseline_0_20_reproduced']}. This comparison does "
                  "not gate interpretation of the uninterrupted CLEAN50 run; "
                  "restart-write cadence is known to affect the numerical path.", ""])
    if result["cycle_metrics"]:
        lines.extend(["| Cycle | Time (ms) | delta_s (mm) | mean v_s (mm/s) | "
                      "end v_s (mm/s) | min v_s (mm/s) | max backtrack (mm) |",
                      "|---:|---:|---:|---:|---:|---:|---:|"])
        for r in result["cycle_metrics"]:
            lines.append(f"| {r['cycle']} | {r['start_ms']:.0f}-{r['end_ms']:.0f} | "
                         f"{r['delta_s_mm']:+.6f} | {r['mean_v_s_mm_s']:+.4f} | "
                         f"{r['end_v_s_mm_s']:+.4f} | {r['minimum_v_s_mm_s']:+.4f} | "
                         f"{r['max_backtrack_mm']:.6f} |")
        lines.extend(["", "The historical C3 delta_s=-0.023527 mm came from "
                      "restart continuation and remains "
                      "`RESTART_CONTAMINATED_C3_RESULT`. CLEAN50 Cycle 3 "
                      f"net recoil reproduced: {result['clean_C3_net_recoil']}. "
                      "Classification treats any non-forward cycle or >=0.010 "
                      "mm within-cycle backtrack as recoil; a sustained-forward "
                      "result requires all five cycles to pass."])
    else:
        lines.extend(["No physical Cycle 1-5 interpretation was issued because "
                      "the solver validity gate failed."])
    (ROOT / "F100_G2P20_CLEAN50_REPORT.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
