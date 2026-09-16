"""Field-only regression for ROBOT_LOCAL_PRECESSING_ROCKING."""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
SERVER = REPO / "calibration_analysis" / "ProductionLocalFrameValidation" / "production" / "magpylib_socket_server.py"
CASE = REPO / "calibration_analysis" / "StraightPipeControl" / "case" / "PROD_LOCAL30_G0_STRAIGHT_CTRL"
DXF = CASE / "straight_control_centerline.dxf"
TRANSFORM = REPO.parent / "abaqus_magpylib_frame_transform.json"
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
MOMENT = 0.0010876227522174417
AMPLITUDE_DEG = 14.343111711438091
ROCK_HZ = 10.0
PREC_HZ = 2.5
ROUTE_A_DEG = -61.37284757596327
PSI0_DEG = -83.87284757596327
DURATION_S = 0.2


def load_server():
    spec = importlib.util.spec_from_file_location("precessing_rocking_server", str(SERVER))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def unit(vector):
    vector = np.asarray(vector, dtype=float)
    return vector / np.linalg.norm(vector)


def wrap_deg(angle):
    return (float(angle) + 180.0) % 360.0 - 180.0


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    module = load_server()
    model = module.MagneticCouplingModel(
        str(DXF), drive_type="analytic", robot_diameter_mm=.815, robot_height_mm=2.4,
        robot_br_t=1.46, robot_moment_Am2=MOMENT, robot_mass_mg=9.207793514166587,
        analytic_b_t=.010, analytic_follow_robot=True, analytic_gradient_b_t=0.0,
        spin_hz=ROCK_HZ, cone_half_angle_deg=30.0, cone_axis_bias_deg=0.0,
        field_frame_mode="ROBOT_LOCAL_PRECESSING_ROCKING",
        rocking_amplitude_deg=AMPLITUDE_DEG,
        precession_frequency_hz=PREC_HZ, precession_phase_deg=PSI0_DEG,
        ramp_time_s=0.0, endpoint_taper_mm=1.0)
    model.robot_polarity = -1.0
    identity = json.loads((CASE / "case_identity.json").read_text())
    target_c_aba = np.asarray(identity["straight_tangent_aba"], dtype=float)
    transform = module.RigidFrameTransform(str(TRANSFORM))
    model.robot_axis_global = transform.vector_to_mag(target_c_aba) / model.robot_polarity
    model.phase_offset_rad = 0.0
    wrapped = module.FrameMappedMagneticModel(model, transform)
    model.reset_robot_arc_continuity()
    first = wrapped.evaluate(0.0, RP0, np.zeros(3))
    robot_arc_mm = float(first["robot_arc_mm"])
    c_mag, e1_mag, e2_mag = model._robot_local_frame(robot_arc_mm)
    c, e1, e2 = [unit(transform.vector_to_aba(v)) for v in (c_mag, e1_mag, e2_mag)]

    time = np.linspace(0.0, DURATION_S, 2001)
    rows = []
    for t in time:
        result = wrapped.evaluate(t, RP0, np.zeros(3))
        field = np.asarray(result["B_aba_vec_T"], dtype=float)
        alpha_command = AMPLITUDE_DEG * math.sin(2.0 * math.pi * ROCK_HZ * t)
        psi_unwrapped = PSI0_DEG + 360.0 * PREC_HZ * t
        rows.append({
            "time_s": t, "alpha_command_deg": alpha_command,
            "psi_unwrapped_deg": psi_unwrapped, "psi_wrapped_deg": wrap_deg(psi_unwrapped),
            "B_T": np.linalg.norm(field), "Bc_T": np.dot(field, c),
            "qB1_T": np.dot(field, e1), "qB2_T": np.dot(field, e2),
            "Bx_aba_T": field[0], "By_aba_T": field[1], "Bz_aba_T": field[2],
        })
    data = pd.DataFrame(rows)
    csv_path = OUT / "PRECESSING_ROCKING_FieldPath.csv"
    data.to_csv(csv_path, index=False)

    transverse = data[["qB1_T", "qB2_T"]].to_numpy()
    centered = transverse - transverse.mean(axis=0)
    singular = np.linalg.svd(centered, compute_uv=False)
    radii = np.linalg.norm(transverse, axis=1)
    field_step = np.linalg.norm(np.diff(data[["Bx_aba_T", "By_aba_T", "Bz_aba_T"]], axis=0), axis=1)
    frame_dots = []
    # The straight validation run remains local to its initial COM arc.  Gate
    # the reachable +/-1 mm neighborhood rather than unrelated remote DXF
    # pieces that are never sampled by this no-gradient 200 ms case.
    s_min = max(0.0, robot_arc_mm - 1.0)
    s_max = min(float(model.arc_mm[-1]), robot_arc_mm + 1.0)
    for s in np.linspace(s_min, s_max, 1001):
        _, f1, f2 = model._robot_local_frame(s)
        if frame_dots:
            frame_dots[-1].extend([float(np.dot(previous1, f1)), float(np.dot(previous2, f2))])
        frame_dots.append([])
        previous1, previous2 = f1, f2
    adjacent_dots = [v for pair in frame_dots for v in pair]

    peak_times = np.array([.025, .075, .125, .175])
    peak_rows = []
    for t in peak_times:
        alpha = AMPLITUDE_DEG * math.sin(2.0 * math.pi * ROCK_HZ * t)
        psi = PSI0_DEG + 360.0 * PREC_HZ * t
        effective = psi if alpha > 0.0 else psi + 180.0
        peak_rows.append({
            "time_ms": 1000.0 * t, "alpha_deg": alpha,
            "plane_azimuth_from_e1_deg": wrap_deg(psi),
            "effective_signed_transverse_azimuth_from_e1_deg": wrap_deg(effective),
            "effective_signed_transverse_azimuth_from_RouteA_deg": wrap_deg(effective - ROUTE_A_DEG),
        })

    metrics = {
        "mode": "ROBOT_LOCAL_PRECESSING_ROCKING",
        "duration_s": DURATION_S, "B_target_T": .010,
        "reused_LIGHTCONTACT_amplitude_deg": AMPLITUDE_DEG,
        "rocking_frequency_Hz": ROCK_HZ, "precession_frequency_Hz": PREC_HZ,
        "RouteA_gauge_deg": ROUTE_A_DEG, "psi0_deg": PSI0_DEG,
        "first_positive_peak_alignment_error_deg": wrap_deg(PSI0_DEG + 360.0 * PREC_HZ * .025 - ROUTE_A_DEG),
        "B_magnitude_max_abs_error_T": float(np.max(np.abs(data.B_T - .010))),
        "alpha_min_deg": float(data.alpha_command_deg.min()),
        "alpha_max_deg": float(data.alpha_command_deg.max()),
        "precession_advance_deg": float(360.0 * PREC_HZ * DURATION_S),
        "field_max_adjacent_step_T": float(field_step.max()),
        "transport_frame_min_adjacent_dot": float(min(adjacent_dots)),
        "transverse_path_second_to_first_singular_ratio": float(singular[1] / singular[0]),
        "transverse_radius_min_T": float(radii.min()),
        "transverse_radius_max_T": float(radii.max()),
        "peak_effective_directions": peak_rows,
    }
    metrics["path_not_straight"] = metrics["transverse_path_second_to_first_singular_ratio"] > 0.1
    metrics["path_not_constant_radius"] = metrics["transverse_radius_min_T"] < 1e-12 and metrics["transverse_radius_max_T"] > 0.002
    metrics["pass"] = bool(
        metrics["B_magnitude_max_abs_error_T"] < 1e-14 and
        abs(metrics["alpha_min_deg"] + AMPLITUDE_DEG) < 1e-10 and
        abs(metrics["alpha_max_deg"] - AMPLITUDE_DEG) < 1e-10 and
        abs(metrics["first_positive_peak_alignment_error_deg"]) < 1e-12 and
        abs(metrics["precession_advance_deg"] - 180.0) < 1e-12 and
        metrics["transport_frame_min_adjacent_dot"] > 0.999 and
        metrics["field_max_adjacent_step_T"] < 2e-4 and
        metrics["path_not_straight"] and metrics["path_not_constant_radius"])
    (OUT / "PRECESSING_ROCKING_FieldPath_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="ascii")
    if not metrics["pass"]:
        raise SystemExit("PRECESSING_ROCKING_FIELD_REGRESSION_FAILED: " + json.dumps(metrics))

    colors = plt.cm.viridis(data.time_s / DURATION_S)
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4))
    axes[0].scatter(data.qB1_T * 1e3, data.qB2_T * 1e3, c=colors, s=5)
    axes[0].axhline(0, color="#777777", lw=.6); axes[0].axvline(0, color="#777777", lw=.6)
    axes[0].set_aspect("equal"); axes[0].set_xlabel("qB1 (mT)"); axes[0].set_ylabel("qB2 (mT)")
    axes[0].set_title("Rotating rocking field path")
    axes[1].plot(time * 1e3, data.alpha_command_deg, color="#c43c35", label="alpha")
    axes[1].plot(time * 1e3, data.psi_unwrapped_deg, color="#315a8a", label="psi")
    axes[1].set_xlabel("Time (ms)"); axes[1].set_ylabel("Angle (deg)"); axes[1].grid(alpha=.2)
    axes[1].legend(frameon=False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"PRECESSING_ROCKING_FieldPath.{suffix}", dpi=220)
    plt.close(fig)

    frame_indices = np.linspace(0, len(data) - 1, 101).astype(int)
    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    limit = 1.08 * AMPLITUDE_DEG * math.pi / 180.0 * 10.0
    def draw(frame):
        ax.cla()
        idx = frame_indices[frame]
        ax.plot(data.qB1_T[:idx + 1] * 1e3, data.qB2_T[:idx + 1] * 1e3, color="#237a57", lw=1.8)
        ax.scatter([data.qB1_T.iloc[idx] * 1e3], [data.qB2_T.iloc[idx] * 1e3], color="#c43c35", s=45)
        ax.axhline(0, color="#777777", lw=.6); ax.axvline(0, color="#777777", lw=.6)
        ax.set_aspect("equal"); ax.set(xlim=(-limit, limit), ylim=(-limit, limit), xlabel="qB1 (mT)", ylabel="qB2 (mT)")
        ax.set_title(f"t={data.time_s.iloc[idx]*1000:5.1f} ms  alpha={data.alpha_command_deg.iloc[idx]:+6.2f} deg\npsi={data.psi_wrapped_deg.iloc[idx]:+7.2f} deg")
    animation = FuncAnimation(fig, draw, frames=len(frame_indices), interval=50)
    animation.save(OUT / "PRECESSING_ROCKING_FieldPath.gif", writer=PillowWriter(fps=20), dpi=100)
    plt.close(fig)
    print(json.dumps(metrics, indent=2))
    print("PASS PRECESSING_ROCKING_FIELD_REGRESSION")


if __name__ == "__main__":
    main()
