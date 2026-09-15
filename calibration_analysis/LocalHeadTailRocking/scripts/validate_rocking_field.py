"""One-cycle offline regression and visualization for ROBOT_LOCAL_ROCKING."""
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


def load_server():
    spec = importlib.util.spec_from_file_location("rocking_server", str(SERVER))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def unit(v):
    return np.asarray(v, dtype=float) / np.linalg.norm(v)


def main():
    module = load_server()
    model = module.MagneticCouplingModel(
        str(DXF), drive_type="analytic", robot_diameter_mm=.815, robot_height_mm=2.4,
        robot_br_t=1.46, robot_moment_Am2=MOMENT, robot_mass_mg=9.207793514166587,
        analytic_b_t=.010, analytic_follow_robot=True, analytic_gradient_b_t=0.0,
        spin_hz=5.0, cone_half_angle_deg=30.0, cone_axis_bias_deg=0.0,
        field_frame_mode="ROBOT_LOCAL_ROCKING", rocking_amplitude_deg=10.0,
        rocking_frame_azimuth_deg=-61.37284757596327,
        ramp_time_s=.001, endpoint_taper_mm=1.0)
    model.robot_polarity = -1.0
    target_c_aba = np.asarray(json.loads((CASE / "case_identity.json").read_text())["straight_tangent_aba"])
    frame_transform = module.RigidFrameTransform(str(TRANSFORM))
    model.robot_axis_global = frame_transform.vector_to_mag(target_c_aba) / model.robot_polarity
    model.phase_offset_rad = 0.0
    wrapped = module.FrameMappedMagneticModel(model, frame_transform)
    model.reset_robot_arc_continuity()
    first = wrapped.evaluate(0.0, RP0, np.zeros(3))
    s = first["robot_arc_mm"]
    c_mag, n_mag, b_mag = model._robot_rocking_frame(s)
    transform = wrapped.frame_transform
    c, n, b = map(unit, (transform.vector_to_aba(c_mag), transform.vector_to_aba(n_mag), transform.vector_to_aba(b_mag)))

    time = np.linspace(0.0, .2, 1001)
    rows = []
    for t in time:
        result = wrapped.evaluate(t, RP0, np.zeros(3))
        field = np.asarray(result["B_aba_vec_T"])
        torque = np.asarray(result["torque_Nmm"])
        alpha = math.degrees(math.atan2(np.dot(field, n), np.dot(field, c)))
        rows.append({"time_s": t, "normalized_phase": t / .2, "alpha_B_deg": alpha,
                     "B_T": np.linalg.norm(field), "B_c_T": np.dot(field, c),
                     "B_n_rock_T": np.dot(field, n), "B_b_rock_T": np.dot(field, b),
                     "T_rock_Nmm": np.dot(torque, b), "Bx_T": field[0], "By_T": field[1], "Bz_T": field[2]})
    data = pd.DataFrame(rows)
    data.to_csv(OUT / "ROBOT_LOCAL_ROCKING_Field_OneCycle.csv", index=False)
    alpha = data.alpha_B_deg.to_numpy()
    field_angle_unwrapped = np.unwrap(np.arctan2(data.B_n_rock_T, data.B_c_T))
    active = np.abs(alpha) > .05
    sign_agreement = float(np.mean(np.sign(alpha[active]) == np.sign(data.T_rock_Nmm.to_numpy()[active])))
    metrics = {
        "B_target_T": .010,
        "B_magnitude_max_abs_error_T": float(np.max(np.abs(data.B_T - .010))),
        "alpha_B_min_deg": float(alpha.min()), "alpha_B_max_deg": float(alpha.max()),
        "B_cross_plane_max_abs_T": float(np.max(np.abs(data.B_b_rock_T))),
        "start_end_B_max_abs_error_T": float(np.max(np.abs(data.iloc[0][["Bx_T","By_T","Bz_T"]].to_numpy(float) - data.iloc[-1][["Bx_T","By_T","Bz_T"]].to_numpy(float)))),
        "azimuthal_winding_turn": float((field_angle_unwrapped[-1] - field_angle_unwrapped[0]) / (2 * np.pi)),
        "Trock_alpha_sign_agreement": sign_agreement,
        "c_hat_aba": c.tolist(), "n_rock_aba": n.tolist(), "b_rock_aba": b.tolist(),
        "right_handed_error": float(np.linalg.norm(np.cross(c, n) - b)),
    }
    model.analytic_rotation_sense = -1.0
    model.cone_half_angle_deg = 71.0
    invariant = np.asarray(wrapped.evaluate(.037, RP0, np.zeros(3))["B_aba_vec_T"])
    model.analytic_rotation_sense = 1.0
    model.cone_half_angle_deg = 30.0
    baseline = np.asarray(wrapped.evaluate(.037, RP0, np.zeros(3))["B_aba_vec_T"])
    metrics["sense_and_cone_parameter_invariance_T"] = float(np.max(np.abs(invariant-baseline)))
    metrics["pass"] = bool(metrics["B_magnitude_max_abs_error_T"] < 1e-14 and
                           abs(metrics["alpha_B_min_deg"] + 10) < 1e-9 and
                           abs(metrics["alpha_B_max_deg"] - 10) < 1e-9 and
                           metrics["B_cross_plane_max_abs_T"] < 1e-14 and
                           metrics["start_end_B_max_abs_error_T"] < 1e-14 and
                           abs(metrics["azimuthal_winding_turn"]) < 1e-12 and
                           metrics["Trock_alpha_sign_agreement"] == 1.0 and
                           metrics["sense_and_cone_parameter_invariance_T"] < 1e-14 and
                           metrics["right_handed_error"] < 1e-12)
    (OUT / "rocking_field_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    if not metrics["pass"]:
        raise SystemExit("ROBOT_LOCAL_ROCKING_FIELD_REGRESSION_FAILED: " + json.dumps(metrics))

    fig, axes = plt.subplots(3, 1, figsize=(8.0, 7.2), sharex=True)
    axes[0].plot(time * 1e3, alpha, color="#c43c35", label=r"$\alpha_B$")
    axes[0].axhline(10, color="#777", lw=.7); axes[0].axhline(-10, color="#777", lw=.7)
    axes[0].set_ylabel("Field angle (deg)"); axes[0].legend(frameon=False)
    axes[1].plot(time * 1e3, data.B_T * 1e3, color="#237a57")
    axes[1].set_ylabel("|B| (mT)")
    axes[2].plot(time * 1e3, data.T_rock_Nmm * 1e3, color="#315a8a")
    axes[2].axhline(0, color="#777", lw=.7); axes[2].set_ylabel(r"$T\cdot b_{rock}$ ($\mu$N mm)"); axes[2].set_xlabel("Time (ms)")
    for ax in axes: ax.grid(alpha=.2)
    fig.suptitle("ROBOT_LOCAL_ROCKING: one physical 5 Hz field cycle")
    fig.tight_layout()
    fig.savefig(OUT / "ROBOT_LOCAL_ROCKING_Field_OneCycle.png", dpi=220)
    plt.close(fig)

    sample = np.linspace(0, len(data) - 1, 81).astype(int)
    fig = plt.figure(figsize=(9.6, 4.8))
    ax3 = fig.add_subplot(121, projection="3d")
    ax2 = fig.add_subplot(122)
    colors = {"c": "#222222", "n": "#c43c35", "b": "#315a8a", "B": "#237a57"}
    def draw(k):
        ax3.cla(); ax2.cla(); row = data.iloc[sample[k]]
        B_local = np.array([row.B_c_T, row.B_n_rock_T, row.B_b_rock_T]) / .010
        for vec, color, label in [(np.eye(3)[0], colors["c"], "c_hat"), (np.eye(3)[1], colors["n"], "n_rock"), (np.eye(3)[2], colors["b"], "b_rock")]:
            ax3.quiver(0,0,0,*vec,color=color,length=.8,arrow_length_ratio=.12); ax3.text(*(vec*.88),label,color=color)
        ax3.quiver(0,0,0,*B_local,color=colors["B"],linewidth=3,arrow_length_ratio=.12)
        ax3.set(xlim=(-.2,1.05),ylim=(-.35,.35),zlim=(-.35,.35),xlabel="c",ylabel="n",zlabel="b",title="Local transported frame")
        arc = np.radians(np.linspace(-10,10,101)); ax2.plot(np.cos(arc),np.sin(arc),color="#bbbbbb",lw=2)
        ax2.arrow(0,0,B_local[0],B_local[1],color=colors["B"],width=.008,length_includes_head=True)
        ax2.axhline(0,color="#777",lw=.7); ax2.set_aspect("equal"); ax2.set(xlim=(.94,1.02),ylim=(-.2,.2),xlabel="c_hat",ylabel="n_rock",title="View along b_rock")
        fig.suptitle(f"t={row.time_s*1000:6.1f} ms   alpha_B={row.alpha_B_deg:+6.2f} deg   |B|={row.B_T*1000:.3f} mT")
    animation = FuncAnimation(fig, draw, frames=len(sample), interval=60)
    animation.save(OUT / "ROBOT_LOCAL_ROCKING_Field_OneCycle.gif", writer=PillowWriter(fps=16.667), dpi=100)
    plt.close(fig)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
