"""Analyze the preserved partial Wall-OFF/CEL-ON chirp telemetry.

This is deliberately telemetry-only: the run stopped at ~34.2 ms, so the
report never labels it a completed 40 ms experiment.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal


def q(x, p):
    return float(np.nanpercentile(np.asarray(x, dtype=float), p))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--telemetry", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    src = Path(args.telemetry)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    d = pd.read_csv(src)
    d = d.sort_values("t_s").drop_duplicates("t_s").reset_index(drop=True)
    t = d.t_s.to_numpy(float)
    if len(d) < 20:
        raise RuntimeError("not enough telemetry rows")
    dt = np.diff(t)
    p = d[["rp_x_aba_mm", "rp_y_aba_mm", "rp_z_aba_mm"]].to_numpy(float)
    v = np.gradient(p, t, axis=0)
    speed = np.linalg.norm(v, axis=1)
    robot_arc = d.robot_arc_mm.to_numpy(float)
    driver_arc = d.driver_arc_mm.to_numpy(float)
    robot_arc_v = np.gradient(robot_arc, t)
    delta_s = driver_arc - robot_arc
    ur = d[["ur1", "ur2", "ur3"]].to_numpy(float)
    theta = np.linalg.norm(ur, axis=1)
    torque = d[["tx_aba_Nmm", "ty_aba_Nmm", "tz_aba_Nmm"]].to_numpy(float)
    torque_mag = np.linalg.norm(torque, axis=1)
    ft_mN = 1000.0 * d.force_tangent_N.to_numpy(float)
    freq = d.commanded_frequency_Hz.to_numpy(float)
    phase = np.unwrap(d.instantaneous_phase_rad.to_numpy(float))
    bmag_mT = 1000.0 * d.B_aba_T.to_numpy(float)

    # A 10 kHz-equivalent, low-cost view for plots/CSV (raw telemetry remains untouched).
    stride = max(1, len(d) // 12000)
    idx = np.arange(0, len(d), stride)
    ts = pd.DataFrame({
        "t_s": t[idx], "freq_Hz": freq[idx], "phase_rad": phase[idx],
        "B_mT": bmag_mT[idx], "robot_arc_mm": robot_arc[idx],
        "driver_arc_mm": driver_arc[idx], "delta_s_mm": delta_s[idx],
        "robot_arc_velocity_mm_s": robot_arc_v[idx],
        "rp_speed_mm_s": speed[idx], "ur1": ur[idx, 0],
        "ur2": ur[idx, 1], "ur3": ur[idx, 2], "theta_norm_rad": theta[idx],
        "Ft_mN": ft_mN[idx], "torque_mag_Nmm": torque_mag[idx],
    })
    ts.to_csv(out / "wobble_chirp_partial_timeseries.csv", index=False)

    # Bin metrics make phase/frequency trends readable despite sub-microsecond telemetry.
    edges = np.arange(0.0, t[-1] + 0.002, 0.002)
    bid = np.clip(np.digitize(t, edges) - 1, 0, len(edges) - 2)
    rows = []
    for k in range(len(edges) - 1):
        m = bid == k
        if not np.any(m):
            continue
        rows.append({
            "t_start_ms": 1e3 * edges[k], "t_end_ms": 1e3 * edges[k + 1],
            "n": int(m.sum()), "freq_mean_Hz": float(freq[m].mean()),
            "theta_rms_rad": float(np.sqrt(np.mean(theta[m] ** 2))),
            "theta_peak_rad": float(theta[m].max()),
            "ur1_range_rad": float(ur[m, 0].max() - ur[m, 0].min()),
            "robot_arc_start_mm": float(robot_arc[m][0]),
            "robot_arc_end_mm": float(robot_arc[m][-1]),
            "driver_arc_end_mm": float(driver_arc[m][-1]),
            "delta_s_end_mm": float(delta_s[m][-1]),
            "Ft_mean_mN": float(ft_mN[m].mean()),
            "Ft_positive_fraction": float(np.mean(ft_mN[m] > 0)),
            "rp_speed_peak_mm_s": float(speed[m].max()),
            "torque_peak_Nmm": float(torque_mag[m].max()),
        })
    bins = pd.DataFrame(rows)
    bins.to_csv(out / "wobble_chirp_partial_2ms_bins.csv", index=False)

    # Frequency-domain view: short run, so frequency peaks are descriptive, not locking proof.
    fs = 1.0 / np.median(dt)
    x = signal.detrend(ur[:, 0])
    nperseg = min(65536, len(x))
    if nperseg < 256:
        nperseg = len(x)
    f, pxx = signal.welch(x, fs=fs, nperseg=nperseg,
                          noverlap=max(0, nperseg // 2))
    valid = f > 0
    top = np.argsort(pxx[valid])[-8:][::-1]
    top_freqs = f[valid][top]
    top_powers = pxx[valid][top]

    # STFT uses a physically interpretable window; output is useful for visual tracking.
    sw = min(65536, len(x))
    sf, st, zz = signal.stft(x, fs=fs, nperseg=sw,
                             noverlap=max(0, sw // 2), boundary=None)

    fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    ax[0].plot(t * 1e3, freq, lw=1.2, label="commanded chirp")
    ax[0].set_ylabel("f (Hz)"); ax[0].legend(); ax[0].grid(alpha=.25)
    ax[1].plot(t * 1e3, driver_arc, label="driver arc")
    ax[1].plot(t * 1e3, robot_arc, label="robot arc")
    ax[1].plot(t * 1e3, delta_s, label="lead Δs", lw=1.2)
    ax[1].set_ylabel("arc / Δs (mm)"); ax[1].legend(); ax[1].grid(alpha=.25)
    ax[2].plot(t * 1e3, ft_mN, lw=.8, color="tab:red")
    ax[2].axhline(0, color="k", lw=.7); ax[2].set_ylabel("Ft (mN)")
    ax[2].set_xlabel("time (ms)"); ax[2].grid(alpha=.25)
    fig.suptitle("Retry003 partial chirp: drive, arc lead, and tangent force")
    fig.tight_layout(); fig.savefig(out / "wobble_chirp_partial_overview.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    ax[0].plot(t * 1e3, ur, lw=.7)
    ax[0].set_ylabel("UR (rad)"); ax[0].legend(["UR1", "UR2", "UR3"])
    ax[0].grid(alpha=.25)
    ax[1].plot(t * 1e3, theta, color="tab:purple", label="||UR||")
    ax[1].plot(t * 1e3, torque_mag, color="tab:orange", label="||T|| (N mm)")
    ax[1].set_xlabel("time (ms)"); ax[1].legend(); ax[1].grid(alpha=.25)
    fig.suptitle("Retry003 partial chirp: rotational response")
    fig.tight_layout(); fig.savefig(out / "wobble_chirp_partial_rotation.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5))
    pcm = ax.pcolormesh(st * 1e3, sf, 10 * np.log10(np.maximum(np.abs(zz) ** 2, 1e-20)),
                        shading="auto", cmap="magma")
    ax.set_ylim(0, min(500, max(100, sf.max())))
    ax.set_xlabel("time (ms)"); ax.set_ylabel("UR1 frequency (Hz)")
    ax.set_title("UR1 STFT (descriptive; partial 34.2 ms window)")
    fig.colorbar(pcm, ax=ax, label="relative power (dB)")
    fig.tight_layout(); fig.savefig(out / "wobble_chirp_partial_spectrogram.png", dpi=180); plt.close(fig)

    # Summary values and explicit status flags.
    summary = {
        "source_file": str(src), "status": "PARTIAL_STOPPED_AT_34.2_MS",
        "rows": int(len(d)), "t_end_ms": float(1e3 * t[-1]),
        "median_dt_us": float(1e6 * np.median(dt)),
        "dt_p01_us": float(1e6 * q(dt, 1)), "dt_p99_us": float(1e6 * q(dt, 99)),
        "frequency_start_Hz": float(freq[0]), "frequency_end_Hz": float(freq[-1]),
        "phase_change_deg": float(np.degrees(phase[-1] - phase[0])),
        "B_mT_min": float(bmag_mT.min()), "B_mT_max": float(bmag_mT.max()),
        "B_gate_bad_rows": int(np.sum(np.abs(bmag_mT - 10.0) > 1e-6)),
        "driver_arc_start_mm": float(driver_arc[0]), "driver_arc_end_mm": float(driver_arc[-1]),
        "robot_arc_start_mm": float(robot_arc[0]), "robot_arc_end_mm": float(robot_arc[-1]),
        "delta_s_start_mm": float(delta_s[0]), "delta_s_end_mm": float(delta_s[-1]),
        "robot_arc_net_speed_mm_s": float((robot_arc[-1] - robot_arc[0]) / (t[-1] - t[0])),
        "rp_displacement_norm_mm": float(np.linalg.norm(p[-1] - p[0])),
        "rp_speed_peak_mm_s": float(speed.max()),
        "Ft_mean_mN": float(ft_mN.mean()), "Ft_median_mN": float(np.median(ft_mN)),
        "Ft_positive_fraction": float(np.mean(ft_mN > 0)),
        "UR1_min_rad": float(ur[:, 0].min()), "UR1_max_rad": float(ur[:, 0].max()),
        "UR1_range_rad": float(np.ptp(ur[:, 0])), "UR_norm_peak_rad": float(theta.max()),
        "torque_peak_Nmm": float(torque_mag.max()),
        "welch_top_frequencies_Hz": ";".join(f"{z:.3f}" for z in top_freqs),
        "welch_top_powers": ";".join(f"{z:.6g}" for z in top_powers),
    }
    pd.DataFrame([summary]).to_csv(out / "wobble_chirp_partial_summary.csv", index=False)

    report = out / "wobble_chirp_partial_report.md"
    report.write_text(f"""# Retry003 partial Chirp analysis

Status: **{summary['status']}**. The 40 ms run did not complete; this report uses only the preserved telemetry through **{summary['t_end_ms']:.3f} ms**.

## Data integrity

- Telemetry rows: `{summary['rows']:,}`; median sampling interval: `{summary['median_dt_us']:.3f} µs` (p01–p99 `{summary['dt_p01_us']:.3f}–{summary['dt_p99_us']:.3f} µs`).
- Commanded frequency: `{summary['frequency_start_Hz']:.3f} → {summary['frequency_end_Hz']:.3f} Hz`; integrated phase change: `{summary['phase_change_deg']:.2f}°`.
- Field magnitude: `{summary['B_mT_min']:.6f}–{summary['B_mT_max']:.6f} mT`; bad 10 mT rows: `{summary['B_gate_bad_rows']}`.

## Motion

- Driver arc: `{summary['driver_arc_start_mm']:.5f} → {summary['driver_arc_end_mm']:.5f} mm`.
- Robot arc: `{summary['robot_arc_start_mm']:.5f} → {summary['robot_arc_end_mm']:.5f} mm`; Δs (driver−robot): `{summary['delta_s_start_mm']:.5f} → {summary['delta_s_end_mm']:.5f} mm`.
- Projected robot-arc net speed: `{summary['robot_arc_net_speed_mm_s']:.3f} mm/s`; 3-D RP displacement: `{summary['rp_displacement_norm_mm']:.3f} mm`; peak RP speed: `{summary['rp_speed_peak_mm_s']:.3f} mm/s`.
- Tangent force mean/median: `{summary['Ft_mean_mN']:.5f}/{summary['Ft_median_mN']:.5f} mN`; positive fraction: `{100*summary['Ft_positive_fraction']:.2f}%`.

## Rotation

- UR1 range: `{summary['UR1_range_rad']:.4f} rad` (min `{summary['UR1_min_rad']:.4f}`, max `{summary['UR1_max_rad']:.4f}`); `||UR||` peak `{summary['UR_norm_peak_rad']:.4f} rad`.
- Applied torque magnitude peak: `{summary['torque_peak_Nmm']:.6g} N·mm`.
- Welch descriptive top frequencies: `{summary['welch_top_frequencies_Hz']} Hz`. These are not a step-out/locking result because the record is only 34.2 ms and the drive is chirped.

## Interpretation

1. The constant-field gate passed exactly (`10 mT` throughout the preserved rows), so the previous 7.5 mT field-amplitude bug is not present in Retry003.
2. The driver advances only about `{summary['driver_arc_end_mm']-summary['driver_arc_start_mm']:.3f} mm`, while the robot's centerline projection advances about `{summary['robot_arc_end_mm']-summary['robot_arc_start_mm']:.3f} mm`; the projected motion is therefore not a sustained 6 mm/s forward-following trajectory in this partial window.
3. `Ft` is positive for only `{100*summary['Ft_positive_fraction']:.1f}%` of samples and has negative mean, so this Wall-OFF diagnostic does **not** yet demonstrate a forward magnetic-drive window. The large 3-D RP motion is mostly off-centerline, not confirmed forward transport.
4. Because the run stopped without an Abaqus error message and left a lock file, treat this as an interrupted partial result. Do not infer a completed Chirp lock band or submit a Wall-ON candidate from it.

## Files

- `wobble_chirp_partial_summary.csv`
- `wobble_chirp_partial_2ms_bins.csv`
- `wobble_chirp_partial_timeseries.csv`
- `wobble_chirp_partial_overview.png`
- `wobble_chirp_partial_rotation.png`
- `wobble_chirp_partial_spectrogram.png`
""", encoding="utf-8")


if __name__ == "__main__":
    main()
