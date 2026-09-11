"""Analyze the available COM-fixed wobble calibration output.

This is read-only post-processing.  The run may have stopped before the
requested 75 ms final frame; the report therefore records the actual ODB and
telemetry end times and never labels a partial run as successful.
"""
from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


JOB = "WobbleCal_COMFixed_F40_Cone30_B10"
TARGET_S = 0.075
root = Path(__file__).resolve().parent
out = root / (JOB + "_analysis")
out.mkdir(exist_ok=True)

ur_path = root / (JOB + "_UR_highres.csv")
ur = pd.read_csv(ur_path)
ur = ur.sort_values("time_s").drop_duplicates("time_s")

telemetry_paths = [root / (JOB + "_telemetry.csv"),
                   root / (JOB + "_telemetry_recovery.csv")]
parts = []
for p in telemetry_paths:
    if p.exists() and p.stat().st_size:
        parts.append(pd.read_csv(p))
if parts:
    tel = pd.concat(parts, ignore_index=True)
    tel = tel.sort_values("t_s").drop_duplicates("t_s", keep="last")
else:
    tel = pd.DataFrame()

# Merge the independent data streams on nearest output time.
if not tel.empty:
    ur2 = ur.rename(columns={"time_s": "t_s"})
    merged = pd.merge_asof(tel, ur2, on="t_s", direction="nearest",
                           tolerance=7.0e-5, suffixes=("", "_odb"))
else:
    merged = ur.rename(columns={"time_s": "t_s"})

for c in ["tx_aba_Nmm", "ty_aba_Nmm", "tz_aba_Nmm",
          "fx_aba_N", "fy_aba_N", "fz_aba_N",
          "B_aba_T", "commanded_frequency_Hz", "instantaneous_phase_rad"]:
    if c not in merged:
        merged[c] = np.nan
merged["torque_norm_Nmm"] = np.sqrt(sum(merged[c].astype(float)**2
                                          for c in ["tx_aba_Nmm", "ty_aba_Nmm", "tz_aba_Nmm"]))
merged["force_norm_N"] = np.sqrt(sum(merged[c].astype(float)**2
                                       for c in ["fx_aba_N", "fy_aba_N", "fz_aba_N"]))
merged.to_csv(out / "runtime_merged.csv", index=False)

def finite(x):
    return x[np.isfinite(x)]

time_end_odb = float(ur.time_s.max()) if len(ur) else float("nan")
time_end_tel = float(tel.t_s.max()) if not tel.empty else float("nan")
time_end = max(time_end_odb, time_end_tel)
complete = bool(time_end >= TARGET_S - 5.0e-5)

phase = finite(merged["instantaneous_phase_rad"].to_numpy(float))
phase_t = merged.loc[np.isfinite(merged["instantaneous_phase_rad"]), "t_s"].to_numpy(float)
phase_freq = float("nan")
if len(phase) >= 3:
    phase_u = np.unwrap(phase)
    phase_freq = float(np.polyfit(phase_t, phase_u, 1)[0] / (2.0 * math.pi))
else:
    phase_u = phase

rows = []
for i, (a, b) in enumerate([(0.0, 0.025), (0.025, 0.050), (0.050, 0.075)], 1):
    u = ur[(ur.time_s >= a) & (ur.time_s < min(b, time_end + 1e-12))]
    q = merged[(merged.t_s >= a) & (merged.t_s < min(b, time_end + 1e-12))]
    if len(u):
        ur_values = u[["UR1", "UR2", "UR3"]].to_numpy(float)
        ur_norm = np.linalg.norm(ur_values, axis=1)
        row = {
            "cycle": i, "window_start_s": a, "window_end_s": min(b, time_end),
            "frames_odb": int(len(u)),
            "ur1_min_rad": float(u.UR1.min()), "ur1_max_rad": float(u.UR1.max()),
            "ur1_peak_abs_rad": float(np.abs(u.UR1).max()),
            "ur_norm_peak_rad": float(ur_norm.max()),
            "torque_peak_Nmm": float(q.torque_norm_Nmm.max()) if len(q) else float("nan"),
            "force_peak_N": float(q.force_norm_N.max()) if len(q) else float("nan"),
        }
        rows.append(row)
cycle = pd.DataFrame(rows)
cycle.to_csv(out / "cycle_metrics.csv", index=False)

summary = {
    "job": JOB,
    "requested_duration_s": TARGET_S,
    "odb_end_s": time_end_odb,
    "telemetry_end_s": time_end_tel,
    "actual_end_s": time_end,
    "completion_status": "COMPLETE" if complete else "PARTIAL_STOP_OR_STALE_LOCK",
    "completion_fraction": time_end / TARGET_S if TARGET_S else float("nan"),
    "odb_frames": int(len(ur)),
    "telemetry_rows_merged": int(len(merged)),
    "phase_fit_frequency_hz": phase_freq,
    "phase_fit_expected_hz": 40.0,
    "ur1_min_rad": float(ur.UR1.min()), "ur1_max_rad": float(ur.UR1.max()),
    "ur2_min_rad": float(ur.UR2.min()), "ur2_max_rad": float(ur.UR2.max()),
    "ur3_min_rad": float(ur.UR3.min()), "ur3_max_rad": float(ur.UR3.max()),
    "ur_vector_peak_rad": float(np.linalg.norm(ur[["UR1", "UR2", "UR3"]].to_numpy(float), axis=1).max()),
    "torque_peak_Nmm": float(merged.torque_norm_Nmm.max()),
    "torque_median_Nmm": float(merged.torque_norm_Nmm.median()),
    "B_median_T": float(merged.B_aba_T.median()),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
ax[0].plot(ur.time_s * 1000, ur.UR1, label="UR1")
ax[0].plot(ur.time_s * 1000, ur.UR2, label="UR2")
ax[0].plot(ur.time_s * 1000, ur.UR3, label="UR3")
ax[0].set_ylabel("rotation (rad)")
ax[0].legend(loc="upper right")
ax[0].set_title(JOB + " — actual available output")
ax[1].plot(merged.t_s * 1000, merged.torque_norm_Nmm, label="|T| Abaqus")
ax[1].set_ylabel("|T| (N mm)")
ax[1].legend(loc="upper right")
ax[2].plot(merged.t_s * 1000, merged.instantaneous_phase_rad, label="field phase")
ax[2].set_ylabel("phase (rad)")
ax[2].set_xlabel("time (ms)")
ax[2].legend(loc="upper left")
fig.tight_layout()
fig.savefig(out / "runtime_diagnostic.png", dpi=180)
plt.close(fig)

report = f"""# {JOB} — partial calibration analysis

## Status

The requested duration was **75.0 ms**.  The ODB contains frames through
**{time_end_odb*1000:.3f} ms** and the merged Socket telemetry through
**{time_end_tel*1000:.3f} ms** ({time_end/TARGET_S*100:.2f}% of the request).
The lock file remains present but no `explicit_dp.exe` process is visible and
no final `ANALYSIS HAS COMPLETED SUCCESSFULLY` record was written.  This is
therefore classified as **`{summary['completion_status']}`**, not as a
successful complete run.

## Drive continuity

- Linear fit to the transmitted analytic phase: **{phase_freq:.3f} Hz**
  (configured 40 Hz).
- Median transmitted field magnitude: **{summary['B_median_T']:.6g} T**.
- Telemetry was recovered after the original server stopped; the main and
  recovery CSV streams were merged by time and deduplicated.

## Rotation response in the available window

- UR1 range: **{summary['ur1_min_rad']:.4f} to {summary['ur1_max_rad']:.4f} rad**
- UR2 range: **{summary['ur2_min_rad']:.4f} to {summary['ur2_max_rad']:.4f} rad**
- UR3 range: **{summary['ur3_min_rad']:.4f} to {summary['ur3_max_rad']:.4f} rad**
- Peak three-axis rotation-vector norm: **{summary['ur_vector_peak_rad']:.4f} rad**
- Applied torque norm: median **{summary['torque_median_Nmm']:.6g} N mm**, peak
  **{summary['torque_peak_Nmm']:.6g} N mm**.

The COM-fixed diagnostic intentionally suppresses translation, so this run
can only calibrate the rotational response and field continuity.  It cannot
validate forward transport, wall interaction, or a 75 ms cycle-to-cycle
steady state because the final 1.4 ms is absent and the job did not close
normally.

## Files

- `runtime_merged.csv` — merged main/recovery telemetry and ODB UR history
- `cycle_metrics.csv` — available cycle-window metrics
- `summary.json` — machine-readable summary
- `runtime_diagnostic.png` — UR, torque, and analytic phase plots

## Decision

**Do not use this partial run as a final calibration baseline.**  Preserve the
ODB and lock file for forensic inspection; after confirming the stale lock and
server shutdown cause, rerun only if a complete 75 ms record is required.
"""
(out / "report.md").write_text(report, encoding="utf-8")
print(json.dumps(summary, indent=2))
