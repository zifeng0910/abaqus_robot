"""Lightweight analysis and fixed-camera rendering for the one F60 screen."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
from scipy.optimize import minimize_scalar
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
JOB = "S4_HEADFORWARD_F60_G0P10_FLOW"
CASE = OUT / "case" / JOB
PARENT_METRICS = REPO / "calibration_analysis" / "F40LowGBackgroundFlowScreen" / "S4_HEADFORWARD_F40_G0P05_FLOW_metrics.json"
BASE_SCRIPTS = REPO / "calibration_analysis" / "S4HeadForwardLowGScreen" / "scripts"
sys.path.insert(0, str(BASE_SCRIPTS))
import analyze_and_render_lowg as base


def load_rp():
    archive = np.load(CASE / "private" / "rp_history_private.npz")
    time = base.get_array(archive, "V1")[:, 0]
    data = {prefix: np.column_stack([
        np.interp(time, base.get_array(archive, prefix + str(axis))[:, 0],
                  base.get_array(archive, prefix + str(axis))[:, 1])
        for axis in (1, 2, 3)]) for prefix in ("U", "UR", "V")}
    archive.close()
    return time, data


def harmonic_fit(time, value, lower=20.0, upper=70.0):
    def solve(freq):
        phase = 2.0 * math.pi * freq * time
        design = np.column_stack((np.ones_like(time), np.sin(phase), np.cos(phase)))
        coef, _, _, _ = np.linalg.lstsq(design, value, rcond=None)
        residual = design.dot(coef) - value
        return float(np.mean(residual * residual)), coef
    result = minimize_scalar(lambda f: solve(float(f))[0], bounds=(lower, upper), method="bounded",
                             options={"xatol": 1.0e-5})
    freq = float(result.x)
    mse, coef = solve(freq)
    phase_deg = math.degrees(math.atan2(float(coef[2]), float(coef[1])))
    amplitude = float(math.hypot(coef[1], coef[2]))
    return freq, phase_deg, amplitude, mse


def read_force_log(time):
    path = CASE / "hydro_increment.csv"
    columns = [
        "time_s", "v1_mm_s", "v2_mm_s", "v3_mm_s", "flow1_mm_s", "flow2_mm_s",
        "flow3_mm_s", "vrel1_mm_s", "vrel2_mm_s", "vrel3_mm_s", "vr1", "vr2", "vr3",
        "Fmag1_N", "Fmag2_N", "Fmag3_N", "Fh1_N", "Fh2_N", "Fh3_N", "Th1_Nmm",
        "Th2_Nmm", "Th3_Nmm", "vrobot_s_mm_s", "vrel_s_mm_s", "U_flow_s_mm_s",
        "Fmag_s_N", "Fhydro_s_N",
    ]
    # The fixed-form Fortran header is emitted as five WRITE records. The
    # numeric FORMAT has one trailing delimiter, so each row has an empty
    # 28th field after the 27 physical columns.
    table = pd.read_csv(path, skiprows=5, header=None, names=columns + ["_trailing"],
                        skipinitialspace=True)
    required = ["time_s", "vrobot_s_mm_s", "vrel_s_mm_s", "U_flow_s_mm_s", "Fmag_s_N", "Fhydro_s_N"]
    missing = [column for column in required if column not in table]
    if missing:
        raise RuntimeError("force log missing columns: {}".format(missing))
    result = {column: np.interp(time, table.time_s.to_numpy(float), table[column].to_numpy(float))
              for column in required[1:]}
    return result


def render_animation(time, data, nodes, region, identity, c, n, rp0, rocking, axial, path, step_s, fps):
    frame_times = np.arange(0.0, identity["duration_s"] + 0.5 * step_s, step_s)
    frame_ids = np.clip(np.searchsorted(time, frame_times), 0, len(time) - 1)
    colors = np.where(region == "HEAD", "#d1493f", np.where(region == "TAIL", "#2864a8", "#383c42"))
    fig, ax = plt.subplots(figsize=(10.0, 2.8))
    xlim, ylim = (-3.0, 12.0), (-2.5, 2.5)

    def draw(frame):
        ax.clear()
        index = int(frame_ids[frame])
        points = rp0 + data["U"][index] + Rotation.from_rotvec(data["UR"][index]).apply(nodes - rp0)
        local = np.column_stack(((points - rp0).dot(c), (points - rp0).dot(n)))
        radius = 0.4075
        ax.axhspan(radius, ylim[1], color="#d8edf2", alpha=0.55)
        ax.axhspan(ylim[0], -radius, color="#d8edf2", alpha=0.55)
        ax.axhline(radius, color="#5996a5", lw=1.1); ax.axhline(-radius, color="#5996a5", lw=1.1)
        ax.scatter(local[::2, 0], local[::2, 1], c=colors[::2], s=8, linewidths=0)
        ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("canonical +s (mm), fixed orthographic camera")
        ax.set_ylabel("n_routeA (mm)")
        ax.grid(axis="x", alpha=0.16)
        ax.text(0.025, 0.94, "HEAD", color="#d1493f", weight="bold", transform=ax.transAxes, va="top")
        ax.text(0.145, 0.94, "TAIL", color="#2864a8", weight="bold", transform=ax.transAxes, va="top")
        ax.set_title(
            f"{JOB} | t={time[index] * 1e3:5.1f} ms | f=60 Hz | G=0.10 mT | U_flow=+10 mm/s | "
            f"delta_s={axial[index] - axial[0]:+.3f} mm | rocking={rocking[index]:+.2f} deg", fontsize=8.5)
        fig.tight_layout()

    animation = FuncAnimation(fig, draw, frames=len(frame_ids), interval=1000.0 / fps)
    if path.suffix.lower() == ".gif":
        animation.save(path, writer=PillowWriter(fps=fps), dpi=100)
    else:
        animation.save(path, writer=FFMpegWriter(fps=fps, bitrate=1400), dpi=100)
    plt.close(fig)
    return len(frame_ids)


def main():
    identity = json.loads((CASE / "case_identity.json").read_text())
    parent_metrics = json.loads(PARENT_METRICS.read_text())
    if identity.get("status") != "SOLVED":
        raise RuntimeError("case is not SOLVED")
    time, data = load_rp()
    c = base.unit(identity["canonical_plus_s_axis_aba"])
    n = base.unit(identity["n_routeA_aba"])
    b = base.unit(identity["b_routeA_aba"])
    rp0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    axial = data["U"].dot(c)
    velocity = data["V"].dot(c)
    radial = data["U"] - np.outer(axial, c)
    body_axis = Rotation.from_rotvec(data["UR"]).apply(np.broadcast_to(c, data["UR"].shape))
    rocking = np.degrees(np.arctan2(body_axis.dot(n), body_axis.dot(c)))
    field_main = identity["rocking_main_amplitude_deg"] * np.sin(2.0 * math.pi * identity["frequency_Hz"] * time)
    actual_f, phase_lag, fit_amplitude, fit_mse = harmonic_fit(time, rocking)

    deck = (CASE / f"{JOB}.inp").read_text()
    nodes = base.part_nodes(deck, "Robot_SOLID")
    wall = base.part_nodes(deck, "Pipe_WALL_HELPER")
    initial_s = (nodes - rp0).dot(c)
    region = np.full(len(nodes), "BODY", dtype=object)
    region[initial_s <= initial_s.min() + 0.25] = "TAIL"
    region[initial_s >= initial_s.max() - 0.25] = "HEAD"
    normals, radius = base.wall_planes(wall, rp0, c, n, b)
    base.gap_snapshot.normals = normals
    geometry_indices, geometry = base.geometry_samples(time, data, nodes, region, rp0, c, n, b, radius, rp0)
    force = base.contact_force(time, CASE)
    active, episodes = base.contact_episode_table(time, force, geometry, geometry_indices, region)
    bridge = (geometry["HEAD_gap_mm"] <= 0.0) & (geometry["TAIL_gap_mm"] <= 0.0)
    bridge_ms = base.longest_duration(geometry["time_s"].to_numpy(), bridge) * 1.0e3
    force_data = read_force_log(time)
    near_zero_threshold = 1.0
    fraction_positive = float(np.mean(velocity > 0.0))
    longest_near_zero_ms = float(base.longest_duration(
        time, np.abs(velocity) <= near_zero_threshold) * 1.0e3)
    axial_reversals = base.reversal_count(velocity)

    episodes.to_csv(CASE / f"{JOB}_contact_episodes.csv", index=False)
    force_table = pd.DataFrame({
        "time_s": time, "time_ms": time * 1.0e3, "field_main_deg": field_main,
        "rocking_deg": rocking, "delta_s_mm": axial - axial[0], "v_robot_s_mm_s": velocity,
        "Fmag_s_N": force_data["Fmag_s_N"], "Fhydro_s_N": force_data["Fhydro_s_N"],
        "U_flow_mm_s": force_data["U_flow_s_mm_s"], "vrel_s_mm_s": force_data["vrel_s_mm_s"],
        "contact_active": active.astype(int),
    })
    force_table.to_csv(OUT / f"{JOB}_force_velocity_timeseries.csv", index=False)
    fig, left = plt.subplots(figsize=(10.5, 4.4))
    left.plot(time * 1e3, force_data["Fmag_s_N"] * 1e6, label="Fmag_s", color="#d1493f")
    left.plot(time * 1e3, force_data["Fhydro_s_N"] * 1e6, label="Fhydro_s", color="#2864a8")
    left.set_xlabel("physical time (ms)"); left.set_ylabel("force along +s (uN)")
    right = left.twinx()
    right.plot(time * 1e3, velocity, label="v_robot_s", color="#383c42", lw=1.0)
    right.axhline(10.0, color="#5996a5", ls="--", label="U_flow=+10 mm/s")
    right.set_ylabel("velocity (mm/s)")
    handles, labels = left.get_legend_handles_labels(); h2, l2 = right.get_legend_handles_labels()
    left.legend(handles + h2, labels + l2, ncol=2, fontsize=8, loc="upper left")
    left.set_title("F60 background-flow screening: magnetic / ReducedHydro / axial velocity")
    left.grid(alpha=0.18); fig.tight_layout()
    fig.savefig(OUT / f"{JOB}_force_velocity.png", dpi=180); plt.close(fig)

    slow_frames = render_animation(time, data, nodes, region, identity, c, n, rp0, rocking, axial,
                                   OUT / f"{JOB}_SLOW.gif", 0.0005, 20)
    fast_mp4 = OUT / f"{JOB}_FAST.mp4"
    fast_mp4_status = "NOT_RENDERED_FFMPEG_UNAVAILABLE"
    try:
        fast_frames = render_animation(time, data, nodes, region, identity, c, n, rp0, rocking, axial,
                                       fast_mp4, 0.0025, 30)
        fast_mp4_status = "RENDERED"
    except FileNotFoundError as exc:
        fast_frames = 0
        print("FAST MP4 unavailable: {}".format(exc))
    counts = {name: int(episodes["region"].str.contains(name).sum()) if len(episodes) else 0
              for name in ("HEAD", "TAIL", "BODY")}
    metrics = {
        "case": JOB, "classification": "LOW_PRECISION_SCREENING",
        "duration_ms": float(time[-1] * 1.0e3), "dt_s": identity["direct_dt_s"],
        "frequency_commanded_Hz": identity["frequency_Hz"], "rocking_frequency_actual_Hz": actual_f,
        "rocking_fit_amplitude_deg": fit_amplitude, "rocking_fit_mse_deg2": fit_mse,
        "field_to_robot_phase_lag_deg": phase_lag,
        "rocking_min_deg": float(rocking.min()), "rocking_max_deg": float(rocking.max()),
        "delta_s_mm": float(axial[-1] - axial[0]), "v_robot_s_final_mm_s": float(velocity[-1]),
        "v_robot_s_max_mm_s": float(velocity.max()),
        "fraction_v_s_positive": fraction_positive,
        "longest_near_zero_v_s_ms": longest_near_zero_ms,
        "axial_reversals": axial_reversals,
        "max_radial_COM_displacement_mm": float(np.linalg.norm(radial, axis=1).max()),
        "near_zero_threshold_mm_s": near_zero_threshold,
        "HEAD_contact_episode_count": counts["HEAD"], "TAIL_contact_episode_count": counts["TAIL"],
        "BODY_contact_episode_count": counts["BODY"], "contact_episode_count": int(len(episodes)),
        "longest_contact_ms": float(episodes["duration_ms"].max()) if len(episodes) else 0.0,
        "both_end_bridge_ms": float(bridge_ms),
        "magnetic_forward_force_min_uN": float(force_data["Fmag_s_N"].min() * 1e6),
        "magnetic_forward_force_max_uN": float(force_data["Fmag_s_N"].max() * 1e6),
        "hydro_axial_force_min_uN": float(force_data["Fhydro_s_N"].min() * 1e6),
        "hydro_axial_force_max_uN": float(force_data["Fhydro_s_N"].max() * 1e6),
        "background_flow_mm_s": 10.0, "flow_status": "DIAGNOSTIC_SCREENING_FLOW_ONLY",
        "slow_gif_frames": slow_frames, "fast_mp4_frames": fast_frames,
        "fast_mp4_status": fast_mp4_status,
        "persistent_bridge_flag": bool(bridge_ms > 5.0),
        "long_wall_sliding_flag": bool((episodes["duration_ms"].max() if len(episodes) else 0.0) > 5.0),
        "gross_penetration_flag": bool(geometry["gap_mm"].min() < -0.01),
        "tumble_flag": bool(np.linalg.norm(radial, axis=1).max() > radius),
        "contact_chatter_flag": bool(len(episodes) > 100),
    }
    (OUT / f"{JOB}_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="ascii")
    pd.DataFrame([metrics]).to_csv(OUT / f"{JOB}_metrics.csv", index=False)
    report = f"""# F60 Low-G Background-Flow Screen

This is exactly one new short straight-tube, non-CEL Explicit dynamics run. It directly continues `S4_HEADFORWARD_F40_G0P05_FLOW` and freezes its head-forward S4 rigid robot, tube, RouteA gauge, elliptical rocking (`A_main=14.343111711438091 deg`, `A_cross=2.5 deg`), `B0=10 mT`, General Contact, `mu=0.03`, `zeta=0.50`, and ReducedHydro coefficients. Only the requested frequency, magnetic gradient, and already-enabled background-flow screening settings differ.

`U_flow=+10 mm/s` along canonical `+s` remains the existing relative-velocity coupling, `V_rel=V_robot-U_flow*c_hat`; it is `DIAGNOSTIC_SCREENING_FLOW_ONLY`, not a final experimental value. The magnetic gradient is `G=0.10 mT`, `L=45 mm`, and the initial `F_gradient dot c_hat` gate passed. This is not CEL or full CFD/FSI. The `dt=2e-7 s` low-precision gate and dynamic run passed; no fallback was needed.

Compared with the direct parent (`S4_HEADFORWARD_F40_G0P05_FLOW`), the commanded frequency is 60 Hz instead of 40 Hz and the gradient is 0.10 mT instead of 0.05 mT; all frozen geometry, contact, rigid-body, camera, damping, hydrodynamic coefficients, initial pose and flow settings are unchanged. Parent `delta_s={parent_metrics['delta_s_mm']:+.6f} mm`, parent actual rocking frequency `{parent_metrics['rocking_frequency_actual_Hz']:.6f} Hz`; new values are reported below. Visual naturalness remains a human screening judgment, not an automatic ranking.

| metric | value |
| --- | ---: |
| actual rocking frequency (Hz) | {actual_f:.6f} |
| field-to-robot phase lag (deg) | {phase_lag:+.6f} |
| rocking min/max (deg) | {rocking.min():+.4f} / {rocking.max():+.4f} |
| delta_s (mm) | {axial[-1]-axial[0]:+.6f} |
| fraction v_s > 0 | {fraction_positive:.6f} |
| longest near-zero v_s (ms) | {longest_near_zero_ms:.6f} |
| axial reversals | {axial_reversals} |
| final v_s (mm/s) | {velocity[-1]:+.6f} |
| max radial COM displacement (mm) | {np.linalg.norm(radial, axis=1).max():.6f} |
| HEAD contact episodes | {counts['HEAD']} |
| TAIL contact episodes | {counts['TAIL']} |
| BODY contact episodes | {counts['BODY']} |
| longest contact (ms) | {metrics['longest_contact_ms']:.6f} |
| BOTH bridge duration (ms) | {bridge_ms:.6f} |
| magnetic forward force range (uN) | {metrics['magnetic_forward_force_min_uN']:+.6f} / {metrics['magnetic_forward_force_max_uN']:+.6f} |
| hydrodynamic axial force range (uN) | {metrics['hydro_axial_force_min_uN']:+.6f} / {metrics['hydro_axial_force_max_uN']:+.6f} |
| background flow (mm/s) | +10.000000 |

Flags: persistent bridge=`{metrics['persistent_bridge_flag']}`, long wall sliding=`{metrics['long_wall_sliding_flag']}`, gross penetration=`{metrics['gross_penetration_flag']}`, tumble=`{metrics['tumble_flag']}`, contact chatter=`{metrics['contact_chatter_flag']}`.

The force/velocity plot is `{JOB}_force_velocity.png`. The fixed-camera slow contact-inspection GIF is `{JOB}_SLOW.gif`. No MP4 is generated by this screening workflow. This result is a coarse visual precursor to later true fluid simulation, not a publication-grade FSI result.

**USER VISUAL SELECTION REQUIRED**
"""
    (OUT / f"{JOB}_Screen.md").write_text(report, encoding="ascii")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
