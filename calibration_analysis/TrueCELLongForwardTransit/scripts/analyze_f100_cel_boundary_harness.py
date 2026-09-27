from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PROBES = ["PROBE_CENTER", "PROBE_LOW_QUARTER", "PROBE_HIGH_QUARTER", "PROBE_LOW_END_NEAR", "PROBE_HIGH_END_NEAR"]


def read_csv(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    data = {key: np.array([float(row[key]) for row in rows], dtype=float) for key in rows[0]}
    data["time_us"] = data.pop("time_s") * 1.0e6
    return data


def metric(t, x, incident, return_window):
    iwin = (t >= incident[0]) & (t <= incident[1])
    rwin = (t >= return_window[0]) & (t <= return_window[1])
    iidx = np.flatnonzero(iwin)[np.argmax(np.abs(x[iwin]))]
    ridx = np.flatnonzero(rwin)[np.argmax(np.abs(x[rwin]))]
    inc = float(abs(x[iidx]))
    ret = float(abs(x[ridx]))
    return {"incident_time_us": float(t[iidx]), "return_time_us": float(t[ridx]), "incident_peak": float(x[iidx]), "return_peak": float(x[ridx]), "R": ret / inc if inc else None}


def main():
    if len(sys.argv) not in (2, 3):
        raise SystemExit("usage: python analyze_f100_cel_boundary_harness.py output_dir [length_mm]")
    out = Path(sys.argv[1])
    length_mm = float(sys.argv[2]) if len(sys.argv) == 3 else 12.0
    quarter_us = (length_mm / 4.0) / 100000.0 * 1.0e6
    end_us = (length_mm / 2.0) / 100000.0 * 1.0e6
    return_quarter_us = (3.0 * length_mm / 4.0) / 100000.0 * 1.0e6
    incident_window = (max(5.0, 0.5 * quarter_us), quarter_us + 25.0)
    return_window = (return_quarter_us - 20.0, return_quarter_us + 35.0)
    h0 = read_csv(out / "HARNESS_DEFAULT_FREE.csv")
    h1 = read_csv(out / "HARNESS_NONREFLECTING.csv")
    metrics = {"official_boundary_reference": "https://help-3dexperience.aesvietnam.com/English/SIMA3DXKEYRefMap/simakey-r-eulerianboundary.htm", "length_mm": length_mm, "c0_mm_s": 100000.0, "expected_center_to_end_us": end_us, "expected_round_trip_us": 2.0 * end_us, "incident_window_us": incident_window, "return_window_us": return_window, "pressure_metric": "MEAN_NORMAL_STRESS_PRESSURE_PROXY", "cases": {}}
    for name, d in [("HARNESS_DEFAULT_FREE", h0), ("HARNESS_NONREFLECTING", h1)]:
        metrics["cases"][name] = {"low_quarter_pressure": metric(d["time_us"], d["PROBE_LOW_QUARTER_p_proxy"], incident_window, return_window), "high_quarter_pressure": metric(d["time_us"], d["PROBE_HIGH_QUARTER_p_proxy"], incident_window, return_window), "low_quarter_velocity": metric(d["time_us"], d["PROBE_LOW_QUARTER_v_axial"], incident_window, return_window), "high_quarter_velocity": metric(d["time_us"], d["PROBE_HIGH_QUARTER_v_axial"], incident_window, return_window)}
    with (out / "HARNESS_REFLECTION_METRICS.json").open("w") as stream:
        json.dump(metrics, stream, indent=2)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for ax, probe, field, title in [(axes[0,0], "PROBE_LOW_QUARTER", "p_proxy", "Low quarter pressure proxy"), (axes[0,1], "PROBE_HIGH_QUARTER", "p_proxy", "High quarter pressure proxy"), (axes[1,0], "PROBE_LOW_QUARTER", "v_axial", "Low quarter axial velocity"), (axes[1,1], "PROBE_HIGH_QUARTER", "v_axial", "High quarter axial velocity")]:
        ax.plot(h0["time_us"], h0[probe + "_" + field], label="H0 free")
        ax.plot(h1["time_us"], h1[probe + "_" + field], label="H1 nonreflecting", linestyle="--")
        ax.axvline(end_us, color="0.6", linewidth=0.8)
        ax.axvline(2.0 * end_us, color="0.6", linewidth=0.8)
        ax.set_title(title)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    axes[1,0].set_xlabel("time (us)")
    axes[1,1].set_xlabel("time (us)")
    fig.tight_layout()
    fig.savefig(out / "HARNESS_REFLECTION_COMPARE.png", dpi=180)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for name, d in [("H0 free", h0), ("H1 nonreflecting", h1)]:
        z = np.vstack([d[p + "_p_proxy"] for p in PROBES])
        axes[0].plot(d["time_us"], z[0], label=name + " center")
        axes[0].plot(d["time_us"], z[1], label=name + " low quarter", linestyle="--")
        axes[0].plot(d["time_us"], z[2], label=name + " high quarter", linestyle=":")
        axes[1].plot(d["time_us"], d["PROBE_LOW_END_NEAR_p_proxy"], label=name + " low end")
        axes[1].plot(d["time_us"], d["PROBE_HIGH_END_NEAR_p_proxy"], label=name + " high end", linestyle="--")
    for ax in axes:
        ax.axvline(end_us, color="0.6", linewidth=0.8)
        ax.axvline(2.0 * end_us, color="0.6", linewidth=0.8)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8, ncol=2)
    axes[0].set_ylabel("p_proxy")
    axes[1].set_ylabel("p_proxy")
    axes[1].set_xlabel("time (us)")
    fig.tight_layout()
    fig.savefig(out / "HARNESS_TIME_DISTANCE_PROBES.png", dpi=180)


if __name__ == "__main__":
    main()
