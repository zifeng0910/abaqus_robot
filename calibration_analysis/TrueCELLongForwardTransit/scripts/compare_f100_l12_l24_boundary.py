"""Compare matched F100 axial-domain replays after both pass kinematic gates."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
L12 = "F100_G2P20_REALWALL_FIXEDVVR_FULLCEL40"
L24 = L12 + "_L24"
PREFIX = "F100_G2P20_L12_L24"


def rows(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def force(job: str) -> dict[str, np.ndarray]:
    records = rows(ROOT / f"{job}_AXIAL_FORCE_10US.csv")
    return {key: np.asarray([float(row[key]) for row in records]) for key in records[0]}


def correlation(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.corrcoef(a, b)[0, 1]) if np.std(a) > 0 and np.std(b) > 0 else float("nan")


def main() -> None:
    gate = {}
    metrics = {}
    for job in (L12, L24):
        gate[job] = json.loads((ROOT / f"{job}_KINEMATIC_GATE.json").read_text())
        if not gate[job]["passed"]:
            raise RuntimeError(f"{job}: KINEMATIC_REPLAY_INVALID; force comparison forbidden")
        metrics[job] = json.loads((ROOT / f"{job}_FORCE_METRICS.json").read_text())
    field = {job: rows(ROOT / f"{job}_FLUID_BOUNDARY_DIAGNOSTICS.csv") for job in (L12, L24)}
    source = {job: force(job) for job in (L12, L24)}
    t = source[L12]["time_s"]
    # Adaptive stable increments produce slightly different history timestamps.
    # Compare both runs on the L12 grid while preserving the native L24 data.
    t24 = source[L24]["time_s"]
    if t24[0] > t[0] + 1e-12 or t24[-1] < t[-1] - 1e-9:
        raise RuntimeError("Force histories do not cover the common 0-40 ms interval")
    for key, values in list(source[L24].items()):
        if key != "time_s":
            source[L24][key] = np.interp(t, t24, values)
    source[L24]["time_s"] = t
    if any(row["pressure_metric_source"] != field[L12][0]["pressure_metric_source"]
           for records in field.values() for row in records):
        raise RuntimeError("L12/L24 pressure metrics have different sources")
    cycles = []
    for index in range(4):
        a, b = metrics[L12]["cycles"][index], metrics[L24]["cycles"][index]
        fa, fb = a["fluid_J_s_Ns"], b["fluid_J_s_Ns"]
        wa, wb = a["wall_J_s_Ns"], b["wall_J_s_Ns"]
        use = (t >= index*.01) & (t <= (index+1)*.01)
        va = source[L12]["fluid_N"][use]
        vb = source[L24]["fluid_N"][use]
        cycles.append({"cycle": index+1,
                       "J_CEL_L12_Ns": fa, "J_CEL_L24_Ns": fb,
                       "D_CEL": abs(fb-fa)/max(abs(fa), abs(fb), 1e-30),
                       "J_wall_L12_Ns": wa, "J_wall_L24_Ns": wb,
                       "D_wall": abs(wb-wa)/max(abs(wa), abs(wb), 1e-30),
                       "fluid_force_correlation": correlation(va, vb),
                       "fluid_force_normalized_rms_difference": float(
                           np.sqrt(np.mean((va-vb)**2)) /
                           max(np.sqrt(np.mean(va**2)), np.sqrt(np.mean(vb**2)), 1e-30)),
                       "fluid_force_p95_abs_L12_N": float(np.percentile(np.abs(va), 95)),
                       "fluid_force_p95_abs_L24_N": float(np.percentile(np.abs(vb), 95)),
                       "fluid_force_peak_abs_L12_N": float(np.max(np.abs(va))),
                       "fluid_force_peak_abs_L24_N": float(np.max(np.abs(vb)))})
    field_comparison = []
    for region in ("near_robot", "low_s_end", "high_s_end"):
        a = [row for row in field[L12] if row["region"] == region]
        b = [row for row in field[L24] if row["region"] == region]
        if len(a) != len(b) or any(abs(float(x["time_s"])-float(y["time_s"])) > 1e-6
                                   for x, y in zip(a, b)):
            raise RuntimeError(f"Sparse field times differ in {region}")
        for x, y in zip(a, b):
            field_comparison.append({"region": region, "time_ms": float(x["time_s"])*1000,
                                     "pressure_proxy_p95_abs_L12_N_mm2": float(x["pressure_metric_p95_abs_N_mm2"]),
                                     "pressure_proxy_p95_abs_L24_N_mm2": float(y["pressure_metric_p95_abs_N_mm2"]),
                                     "pressure_proxy_mean_L12_N_mm2": float(x["pressure_metric_mean_N_mm2"]),
                                     "pressure_proxy_mean_L24_N_mm2": float(y["pressure_metric_mean_N_mm2"]),
                                     "axial_velocity_mean_L12_mm_s": float(x["axial_velocity_mean_mm_s"]),
                                     "axial_velocity_mean_L24_mm_s": float(y["axial_velocity_mean_mm_s"]),
                                     "speed_p95_L12_mm_s": float(x["speed_p95_mm_s"]),
                                     "speed_p95_L24_mm_s": float(y["speed_p95_mm_s"]),
                                     "EVF_volume_proxy_L12_mm3": float(x["represented_fluid_volume_proxy_mm3"]),
                                     "EVF_volume_proxy_L24_mm3": float(y["represented_fluid_volume_proxy_mm3"])})
    contact = {}
    for job in (L12, L24):
        items = [row for row in rows(ROOT / f"{job}_CONTACT_EVENTS.csv")
                 if row["case"] == "FULL_CEL_FIXED_REPLAY"]
        contact[job] = {end: {"events": sum(row["end"] == end for row in items),
                              "same_end_recontacts": sum(row["end"] == end and
                                                         row["same_end_recontact"] == "True" for row in items),
                              "minimum_signed_gap_mm": min(float(row["minimum_signed_gap_mm"])
                                                           for row in items if row["end"] == end)}
                        for end in ("HEAD", "TAIL")}
    comparison = {"cases": [L12, L24], "cycles": cycles, "contact": contact,
                  "pressure_metric_source": field[L12][0]["pressure_metric_source"],
                  "fluid_field_sampling": "six native field times over 0-40 ms; not a resolved pressure waveform",
                  "flux_status": "NO_NATIVE_FACE_MASS_FLUX; EVF volume is a represented-volume proxy",
                  "force_label": "INFERRED ROBOT-CEL FORCE = whole robot General Contact minus direct robot-wall pair",
                  "decision": "PENDING_EVIDENCE_REVIEW"}
    (ROOT / f"{PREFIX}_COMPARISON.json").write_text(json.dumps(comparison, indent=2) + "\n")
    with (ROOT / f"{PREFIX}_FIELD_COMPARISON.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(field_comparison[0]))
        writer.writeheader()
        writer.writerows(field_comparison)
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True, constrained_layout=True)
    for job, label, color in ((L12, "L12", "#426ea8"), (L24, "L24", "#bd563f")):
        axes[0].plot(t*1000, source[job]["fluid_N"], color=color, label=label, lw=.7)
        axes[1].plot(t*1000, source[job]["fluid_J_Ns"], color=color, label=label, lw=1.2)
    for ax in axes:
        for boundary in (10, 20, 30):
            ax.axvline(boundary, color="0.65", ls="--", lw=.6)
        ax.grid(alpha=.2)
        ax.legend()
    axes[0].set(ylabel="inferred robot-CEL axial force (N)")
    axes[1].set(xlabel="time (ms)", ylabel="cumulative axial impulse (N s)")
    fig.savefig(ROOT / f"{PREFIX}_FORCE_IMPULSE.png", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.5), sharex=True, constrained_layout=True)
    specs = (("near_robot", "pressure_proxy_p95_abs", axes[0, 0], "near-robot p95 |mean stress| (N/mm2)"),
             ("near_robot", "speed_p95", axes[0, 1], "near-robot p95 speed (mm/s)"),
             ("low_s_end", "pressure_proxy_mean", axes[1, 0], "low-s end mean pressure proxy (N/mm2)"),
             ("high_s_end", "axial_velocity_mean", axes[1, 1], "high-s end mean axial velocity (mm/s)"))
    for region, metric, ax, ylabel in specs:
        selected = [row for row in field_comparison if row["region"] == region]
        for label, color in (("L12", "#426ea8"), ("L24", "#bd563f")):
            ax.plot([row["time_ms"] for row in selected],
                    [row[f"{metric}_{label}_N_mm2"] if metric.startswith("pressure") else
                     row[f"{metric}_{label}_mm_s"] for row in selected],
                    marker="o", color=color, label=label)
        ax.set(xlabel="time (ms)", ylabel=ylabel)
        ax.grid(alpha=.2)
        ax.legend()
    fig.suptitle("Sparse native CEL field samples; pressure is -trace(S_water)/3 proxy")
    fig.savefig(ROOT / f"{PREFIX}_PRESSURE_VELOCITY.png", dpi=180)
    plt.close(fig)
    print(json.dumps({"D_C3": cycles[2]["D_CEL"], "D_C4": cycles[3]["D_CEL"],
                      "C2_C4_signs": [(np.sign(row["J_CEL_L12_Ns"]), np.sign(row["J_CEL_L24_Ns"]))
                                      for row in cycles[1:]]}, indent=2))


if __name__ == "__main__":
    main()
