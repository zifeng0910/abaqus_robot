"""Analyze and render the single strict-CEL F60/G0.15 screening run."""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image
from scipy.optimize import minimize_scalar
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
JOB = "S4_HEADFORWARD_F60_G0P15_TRUECEL"
CASE = OUT / "case" / JOB
BASE_SCRIPTS = REPO / "calibration_analysis" / "S4HeadForwardLowGScreen" / "scripts"
PARENT_METRICS = REPO / "calibration_analysis" / "F60LowGBackgroundFlowScreen" / "S4_HEADFORWARD_F60_G0P10_FLOW_metrics.json"
sys.path.insert(0, str(BASE_SCRIPTS))
import analyze_and_render_lowg as base


BLUE = "#2864a8"
RED = "#d1493f"
NEAR_ZERO_MM_S = 1.0


def load_rp():
    archive = np.load(CASE / "private" / "rp_history_private.npz")
    time = base.get_array(archive, "V1")[:, 0]
    data = {prefix: np.column_stack([
        np.interp(time, base.get_array(archive, prefix + str(axis))[:, 0],
                  base.get_array(archive, prefix + str(axis))[:, 1])
        for axis in (1, 2, 3)]) for prefix in ("U", "UR", "V")}
    archive.close()
    return time, data


def harmonic_fit(time, value, lower=40.0, upper=80.0):
    def solve(freq):
        phase = 2.0 * math.pi * freq * time
        design = np.column_stack((np.ones_like(time), np.sin(phase), np.cos(phase)))
        coef, _, _, _ = np.linalg.lstsq(design, value, rcond=None)
        residual = design.dot(coef) - value
        return float(np.mean(residual * residual)), coef
    result = minimize_scalar(lambda f: solve(float(f))[0], bounds=(lower, upper), method="bounded")
    mse, coef = solve(float(result.x))
    return float(result.x), math.degrees(math.atan2(float(coef[2]), float(coef[1]))), float(math.hypot(coef[1], coef[2])), mse


def parse_elements(deck, part_name):
    part = re.search(rf"(?ms)^\*Part, name={re.escape(part_name)}\s*$.*?^\*End Part\s*$", deck).group(0)
    block = re.search(r"(?ms)^\*Element, type=EC3D8R\s*$\n(.*?)(?=^\*)", part).group(1)
    rows = [[int(value.strip()) for value in line.split(",") if value.strip()]
            for line in block.splitlines() if line.strip()]
    return {row[0]: row[1:] for row in rows}


def telemetry_force(time, c):
    table = pd.read_csv(CASE / f"{JOB}_telemetry.csv")
    table = table.drop_duplicates("t_s").sort_values("t_s")
    force = table[["fx_aba_N", "fy_aba_N", "fz_aba_N"]].to_numpy(float)
    fmag_s = force.dot(c)
    return np.interp(time, table.t_s.to_numpy(float), fmag_s), table


def contact_components(time):
    archive = np.load(CASE / "private" / "contact_history_private.npz")
    result = {}
    for surface_name in ("ROBOT", "PIPE_WALL_HELPER"):
        vectors = []
        for prefix in ("CFN", "CFS"):
            axes = []
            for axis in (1, 2, 3):
                keys = [key for key in archive.files if "|" + prefix + str(axis) in key and surface_name in key]
                if keys:
                    array = max((archive[key] for key in keys), key=len)
                    axes.append(np.interp(time, array[:, 0], array[:, 1]))
                else:
                    axes.append(np.zeros_like(time))
            vectors.append(np.column_stack(axes))
        result[surface_name] = vectors[0] + vectors[1]
    archive.close()
    return result


def solver_stats():
    text = (CASE / f"{JOB}.sta").read_text(encoding="latin1")
    rows = re.findall(r"(?m)^\s*(\d+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)\s+\S+\s+([0-9.E+-]+)\s+\d+", text)
    increments = np.asarray([int(row[0]) for row in rows], dtype=int)
    dt = np.asarray([float(row[3]) for row in rows], dtype=float)
    return {"total_increments": int(increments.max()), "minimum_dt_s": float(dt.min()),
            "typical_dt_s": float(np.median(dt[dt > 0.0]))}


def fixed_regions(nodes, rp0, c):
    s = (nodes - rp0).dot(c)
    region = np.full(len(nodes), "BODY", dtype=object)
    region[s <= s.min() + 0.25] = "TAIL"
    region[s >= s.max() - 0.25] = "HEAD"
    return region


def render_robot(time, data, nodes, region, identity, c, n, rp0, rocking, axial, output, fps):
    render_duration = min(float(identity["duration_s"]), float(time[-1]))
    frame_times = np.arange(0.0, render_duration + 0.00025, 0.0005)
    frame_ids = np.clip(np.searchsorted(time, frame_times), 0, len(time) - 1)
    colors = np.where(region == "HEAD", RED, np.where(region == "TAIL", BLUE, "#383c42"))
    tail_ref = nodes[region == "TAIL"].mean(axis=0) - rp0
    head_ref = nodes[region == "HEAD"].mean(axis=0) - rp0
    fig, ax = plt.subplots(figsize=(10.0, 2.8))

    def draw(frame):
        ax.clear()
        index = int(frame_ids[frame])
        rotation = Rotation.from_rotvec(data["UR"][index])
        points = rp0 + data["U"][index] + rotation.apply(nodes - rp0)
        tail = rp0 + data["U"][index] + rotation.apply(tail_ref)
        head = rp0 + data["U"][index] + rotation.apply(head_ref)
        local = np.column_stack(((points - rp0).dot(c), (points - rp0).dot(n)))
        tail_xy = ((tail - rp0).dot(c), (tail - rp0).dot(n))
        head_xy = ((head - rp0).dot(c), (head - rp0).dot(n))
        ax.axhspan(0.4075, 2.5, color="#d8edf2", alpha=0.55)
        ax.axhspan(-2.5, -0.4075, color="#d8edf2", alpha=0.55)
        ax.axhline(0.4075, color="#5996a5", lw=1.1)
        ax.axhline(-0.4075, color="#5996a5", lw=1.1)
        ax.scatter(local[::2, 0], local[::2, 1], c=colors[::2], s=8, linewidths=0)
        ax.text(*tail_xy, "TAIL", color=BLUE, weight="bold", fontsize=8, ha="right", va="bottom",
                bbox={"facecolor": "white", "alpha": 0.72, "pad": 1.2})
        ax.text(*head_xy, "HEAD", color=RED, weight="bold", fontsize=8, ha="left", va="bottom",
                bbox={"facecolor": "white", "alpha": 0.72, "pad": 1.2})
        ax.set_xlim(-3.0, 12.0); ax.set_ylim(-2.5, 2.5); ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("canonical +s (mm), fixed orthographic camera"); ax.set_ylabel("n_routeA (mm)")
        ax.grid(axis="x", alpha=0.16)
        ax.text(0.025, 0.95, "TAIL", color=BLUE, weight="bold", transform=ax.transAxes, va="top")
        ax.text(0.145, 0.95, "HEAD", color=RED, weight="bold", transform=ax.transAxes, va="top")
        ax.set_title(f"{JOB} | t={time[index]*1e3:5.1f} ms | f=60 Hz | G=0.15 mT | TRUE CEL | "
                     f"U0=+10 mm/s | delta_s={axial[index]-axial[0]:+.3f} mm | rocking={rocking[index]:+.2f} deg",
                     fontsize=8.1)
        fig.tight_layout()

    animation = FuncAnimation(fig, draw, frames=len(frame_ids), interval=1000.0 / fps)
    animation.save(output, writer=PillowWriter(fps=fps), dpi=100)
    plt.close(fig)
    return len(frame_ids)


def fluid_snapshots(field, deck, c, n, rp0, time, data, nodes, region):
    element_map = parse_elements(deck, "FLUID_EULERIAN")
    node_labels = field["fluid_node_labels"].astype(int)
    label_to_index = {label: index for index, label in enumerate(node_labels)}
    element_labels = field["fluid_element_labels"].astype(int)
    connectivity = np.asarray([[label_to_index[label] for label in element_map[int(eid)]] for eid in element_labels], dtype=int)
    coords = field["fluid_node_coordinates_mm"]
    centroids = coords[connectivity].mean(axis=1)
    local_centroids = np.column_stack(((centroids - rp0).dot(c), (centroids - rp0).dot(n)))
    indices = np.linspace(0, len(field["time"]) - 1, 4).round().astype(int)
    speeds = np.linalg.norm(field["fluid_velocity_mm_s"], axis=2)
    vmax = max(10.0, float(np.nanpercentile(speeds, 99.0)))
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 5.5), sharex=True, sharey=True)
    scatter = None
    for ax, frame in zip(axes.flat, indices):
        evf = field["fluid_evf"][frame]
        node_speed = speeds[frame]
        element_speed = np.nanmean(node_speed[connectivity], axis=1)
        active = np.isfinite(evf) & (evf > 0.02)
        scatter = ax.scatter(local_centroids[active, 0], local_centroids[active, 1],
                             c=element_speed[active], s=10, marker="s", cmap="viridis", vmin=0.0, vmax=vmax)
        rp_index = int(np.clip(np.searchsorted(time, field["time"][frame]), 0, len(time) - 1))
        rotation = Rotation.from_rotvec(data["UR"][rp_index])
        points = rp0 + data["U"][rp_index] + rotation.apply(nodes - rp0)
        local_robot = np.column_stack(((points - rp0).dot(c), (points - rp0).dot(n)))
        colors = np.where(region == "HEAD", RED, np.where(region == "TAIL", BLUE, "#383c42"))
        ax.scatter(local_robot[::3, 0], local_robot[::3, 1], c=colors[::3], s=5, linewidths=0)
        ax.axhline(0.4075, color="#5996a5", lw=1.0); ax.axhline(-0.4075, color="#5996a5", lw=1.0)
        ax.set_title(f"t={field['time'][frame]*1e3:.2f} ms | active EVF cells={int(active.sum())}")
        ax.set_xlim(-6.25, 6.25); ax.set_ylim(-0.62, 0.62); ax.grid(alpha=0.12)
    for ax in axes[1]: ax.set_xlabel("canonical +s (mm)")
    for ax in axes[:, 0]: ax.set_ylabel("n_routeA (mm)")
    fig.colorbar(scatter, ax=axes.ravel().tolist(), label="actual CEL nodal velocity magnitude (mm/s)", shrink=0.85)
    fig.suptitle("TRUE CEL fluid velocity and material volume fraction around the robot")
    fig.savefig(OUT / "TRUECEL_velocity_snapshots.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    return connectivity, centroids


def main():
    identity = json.loads((CASE / "case_identity.json").read_text())
    if identity.get("status") not in ("SOLVED", "ABORTED_PARTIAL") or identity.get("classification") != "STRICT_CEL_FLUID_PRESENT":
        raise RuntimeError("strict CEL case is not available for analysis")
    parent = json.loads(PARENT_METRICS.read_text())
    time, data = load_rp()
    c = base.unit(identity["canonical_plus_s_axis_aba"]); n = base.unit(identity["n_routeA_aba"]); b = base.unit(identity["b_routeA_aba"])
    rp0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    axial = data["U"].dot(c); velocity = data["V"].dot(c)
    acceleration = np.gradient(data["V"], time, axis=0)
    radial = data["U"] - np.outer(axial, c)
    body_axis = Rotation.from_rotvec(data["UR"]).apply(np.broadcast_to(c, data["UR"].shape))
    rocking = np.degrees(np.arctan2(body_axis.dot(n), body_axis.dot(c)))
    actual_f, phase_lag, fit_amplitude, fit_mse = harmonic_fit(time, rocking)
    deck = (CASE / f"{JOB}.inp").read_text(encoding="latin1")
    nodes = base.part_nodes(deck, "Robot_SOLID"); wall = base.part_nodes(deck, "Pipe_WALL_HELPER")
    region = fixed_regions(nodes, rp0, c)
    normals, radius = base.wall_planes(wall, rp0, c, n, b); base.gap_snapshot.normals = normals
    geom_indices, geometry = base.geometry_samples(time, data, nodes, region, rp0, c, n, b, radius, rp0)
    contact = contact_components(time); total_contact = contact["ROBOT"]
    active, episodes = base.contact_episode_table(time, total_contact, geometry, geom_indices, region)
    bridge = (geometry["HEAD_gap_mm"].to_numpy() <= 0.0) & (geometry["TAIL_gap_mm"].to_numpy() <= 0.0)
    bridge_ms = base.longest_duration(geometry["time_s"].to_numpy(), bridge) * 1e3
    episodes.to_csv(CASE / f"{JOB}_contact_episodes.csv", index=False)
    counts = {name: int(episodes["region"].str.contains(name).sum()) if len(episodes) else 0 for name in ("HEAD", "TAIL", "BODY")}

    fmag_s, telemetry = telemetry_force(time, c)
    inertial_s = identity["robot_mass_mg"] * 1.0e-9 * acceleration.dot(c)
    nonmagnetic_force_s = inertial_s - fmag_s
    field = np.load(CASE / "private" / "truecel_field_private.npz")
    connectivity, centroids = fluid_snapshots(field, deck, c, n, rp0, time, data, nodes, region)
    fluid_velocity = field["fluid_velocity_mm_s"]
    initial_flow = np.asarray(identity["flow_vector_aba_mm_s"])
    node_local_s = (field["fluid_node_coordinates_mm"] - rp0).dot(c)
    node_radial = field["fluid_node_coordinates_mm"] - rp0 - np.outer(node_local_s, c)
    near = (np.abs(node_local_s) <= 2.0) & (np.linalg.norm(node_radial, axis=1) <= 0.5)
    near_velocity = fluid_velocity[:, near, :]
    near_speed = np.linalg.norm(near_velocity, axis=2)
    disturbance = np.linalg.norm(near_velocity - initial_flow, axis=2)
    early = field["time"] <= 0.002
    early_disturbance = float(np.nanmax(disturbance[early]))
    initial_filled = int(np.sum(field["fluid_evf"][0] > 0.5))
    coupling_active = bool(initial_filled > 0 and early_disturbance > 1.0e-3 and np.nanmax(np.linalg.norm(total_contact, axis=1)) > 0.0)
    coupling_gate = {
        "classification": "TRUE_CEL_COUPLING_ACTIVE" if coupling_active else "CEL_COUPLING_NOT_ACTIVE",
        "passed": coupling_active, "initial_EVF_cells_gt_0p5": initial_filled,
        "early_near_robot_velocity_disturbance_max_mm_s": early_disturbance,
        "max_robot_contact_resultant_N": float(np.nanmax(np.linalg.norm(total_contact, axis=1))),
        "evidence": "actual ODB EVF and V fields plus nonzero general-contact resultant on robot",
    }
    (OUT / "truecel_coupling_gate.json").write_text(json.dumps(coupling_gate, indent=2) + "\n", encoding="ascii")
    if not coupling_active:
        raise RuntimeError("CEL_COUPLING_NOT_ACTIVE")

    primary = OUT / f"{JOB}_2D.gif"; slow = OUT / f"{JOB}_SLOW.gif"
    frame_count = render_robot(time, data, nodes, region, identity, c, n, rp0, rocking, axial, primary, 20)
    render_robot(time, data, nodes, region, identity, c, n, rp0, rocking, axial, slow, 10)
    with Image.open(primary) as gif:
        gif.seek(0); gif.convert("RGB").save(OUT / f"{JOB}_FRAME0.png")

    solver = solver_stats()
    mean_v = float((axial[-1] - axial[0]) / (time[-1] - time[0]))
    metrics = {
        "case": JOB, "run_status": identity.get("status"), "classification": "STRICT_CEL_FLUID_PRESENT", "flow_classification": "INITIALIZED_UNIFORM_FLOW_TRUE_CEL",
        "simulated_time_s": float(time[-1]), "target_duration_s": float(identity["duration_s"]),
        "frequency_commanded_Hz": 60.0, "rocking_frequency_actual_Hz": actual_f,
        "rocking_fit_amplitude_deg": fit_amplitude, "rocking_fit_mse_deg2": fit_mse,
        "field_to_robot_phase_lag_deg": phase_lag, "rocking_min_deg": float(rocking.min()), "rocking_max_deg": float(rocking.max()),
        "delta_s_mm": float(axial[-1] - axial[0]), "mean_axial_velocity_mm_s": mean_v,
        "final_axial_velocity_mm_s": float(velocity[-1]), "maximum_axial_velocity_mm_s": float(velocity.max()),
        "fraction_v_s_positive": float(np.mean(velocity > 0.0)), "axial_reversals": base.reversal_count(velocity),
        "max_radial_COM_displacement_mm": float(np.linalg.norm(radial, axis=1).max()),
        "HEAD_contact_episode_count": counts["HEAD"], "TAIL_contact_episode_count": counts["TAIL"],
        "BODY_contact_episode_count": counts["BODY"], "contact_episode_count": int(len(episodes)),
        "longest_contact_ms": float(episodes["duration_ms"].max()) if len(episodes) else 0.0,
        "both_end_bridge_ms": float(bridge_ms), "fluid_density_tonne_mm3": identity["fluid_density_tonne_mm3"],
        "fluid_EOS": {"type": "USUP", "c0_mm_s": identity["fluid_EOS_c0_mm_s"], "s": 0.0, "Gamma0": 0.0},
        "fluid_viscosity_N_s_mm2": identity["fluid_viscosity_N_s_mm2"],
        "fluid_filled_volume_mm3": identity["fluid_filled_volume_mm3"], "total_initial_fluid_mass_mg": identity["fluid_mass_mg"],
        "initial_fluid_velocity_mm_s": identity["flow_vector_aba_mm_s"],
        "near_robot_fluid_velocity_min_mm_s": float(np.nanmin(near_speed)), "near_robot_fluid_velocity_max_mm_s": float(np.nanmax(near_speed)),
        "CEL_axial_force_on_robot_estimate_min_uN": float(np.nanmin(nonmagnetic_force_s) * 1e6),
        "CEL_axial_force_on_robot_estimate_max_uN": float(np.nanmax(nonmagnetic_force_s) * 1e6),
        "CEL_transverse_force_on_robot_estimate_max_uN": float(np.nanmax(np.linalg.norm(acceleration - np.outer(acceleration.dot(c), c), axis=1) * identity["robot_mass_mg"] * 1e-9) * 1e6),
        "CEL_force_estimate_method": "robot momentum balance after subtracting Magpylib; includes any simultaneous robot-wall contact",
        "magnetic_axial_force_min_uN": float(np.min(fmag_s) * 1e6), "magnetic_axial_force_max_uN": float(np.max(fmag_s) * 1e6),
        "coupling_gate": coupling_gate, "tumble_flag": bool(np.any(body_axis.dot(c) < 0.0)),
        "gross_penetration_flag": bool(geometry["gap_mm"].min() < -0.01), "gif_frame_count": frame_count,
        "parent_case": parent["case"], "parent_delta_s_mm": parent["delta_s_mm"],
        "wallclock_s": identity.get("wallclock_s"), **solver,
    }
    (OUT / f"{JOB}_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="ascii")
    pd.DataFrame([{key: value for key, value in metrics.items() if not isinstance(value, (dict, list))}]).to_csv(OUT / f"{JOB}_metrics.csv", index=False)
    pd.DataFrame({"time_s": time, "Fmag_s_N": fmag_s, "Ffluid_plus_wall_est_s_N": nonmagnetic_force_s,
                  "v_robot_s_mm_s": velocity}).to_csv(OUT / f"{JOB}_force_velocity_timeseries.csv", index=False)
    fig, ax = plt.subplots(figsize=(10.5, 4.4)); ax.plot(time*1e3, fmag_s*1e6, label="Fmag_s", color=RED)
    ax.plot(time*1e3, nonmagnetic_force_s*1e6, label="Ffluid+wall estimate_s", color=BLUE)
    ax.set_xlabel("time (ms)"); ax.set_ylabel("force along +s (uN)"); right = ax.twinx()
    right.plot(time*1e3, velocity, color="#383c42", label="v_robot_s"); right.set_ylabel("robot axial velocity (mm/s)")
    h1,l1=ax.get_legend_handles_labels(); h2,l2=right.get_legend_handles_labels(); ax.legend(h1+h2,l1+l2,fontsize=8)
    ax.grid(alpha=0.18); fig.tight_layout(); fig.savefig(OUT / f"{JOB}_forces.png", dpi=180); plt.close(fig)

    comparison = {
        "parent": {key: parent.get(key) for key in ("rocking_frequency_actual_Hz", "rocking_fit_amplitude_deg", "field_to_robot_phase_lag_deg", "delta_s_mm", "v_robot_s_final_mm_s", "max_radial_COM_displacement_mm", "HEAD_contact_episode_count", "TAIL_contact_episode_count", "BODY_contact_episode_count", "both_end_bridge_ms")},
        "truecel": {key: metrics.get(key) for key in ("rocking_frequency_actual_Hz", "rocking_fit_amplitude_deg", "field_to_robot_phase_lag_deg", "delta_s_mm", "final_axial_velocity_mm_s", "max_radial_COM_displacement_mm", "HEAD_contact_episode_count", "TAIL_contact_episode_count", "BODY_contact_episode_count", "both_end_bridge_ms")},
    }
    (OUT / f"{JOB}_comparison.json").write_text(json.dumps(comparison, indent=2) + "\n", encoding="ascii")
    report = f"""# F60/G0.15 Strict CEL Screen

`{JOB}` is the only new case and is classified `STRICT_CEL_FLUID_PRESENT` / `INITIALIZED_UNIFORM_FLOW_TRUE_CEL`. It uses the authoritative F60 parent magnetic path unchanged except `G=0.15 mT`; ReducedHydro is absent. Actual ODB `EVF` and `V` fields show {initial_filled} initially filled cells and an early near-robot velocity disturbance of {early_disturbance:.6f} mm/s, so the CEL coupling gate passed.

The fixed view is canonical `+s` left-to-right, TAIL left/blue and HEAD right/red. The fluid snapshot is derived from actual ODB volume fraction and velocity fields, not a schematic.

| metric | value |
| --- | ---: |
| actual rocking frequency (Hz) | {actual_f:.6f} |
| fitted rocking amplitude (deg) | {fit_amplitude:.6f} |
| rocking min/max (deg) | {rocking.min():+.6f} / {rocking.max():+.6f} |
| phase lag (deg) | {phase_lag:+.6f} |
| delta_s (mm) | {metrics['delta_s_mm']:+.6f} |
| mean/final/max axial velocity (mm/s) | {mean_v:+.6f} / {velocity[-1]:+.6f} / {velocity.max():+.6f} |
| fraction v_s > 0 | {metrics['fraction_v_s_positive']:.6f} |
| axial reversals | {metrics['axial_reversals']} |
| max radial COM displacement (mm) | {metrics['max_radial_COM_displacement_mm']:.6f} |
| HEAD/TAIL/BODY contact episodes | {counts['HEAD']} / {counts['TAIL']} / {counts['BODY']} |
| longest contact / BOTH bridge (ms) | {metrics['longest_contact_ms']:.6f} / {bridge_ms:.6f} |
| fluid-filled volume / mass | {identity['fluid_filled_volume_mm3']:.6f} mm3 / {identity['fluid_mass_mg']:.6f} mg |
| CEL dt min/typical (s) | {solver['minimum_dt_s']:.6e} / {solver['typical_dt_s']:.6e} |

The CEL force series is a robot momentum-balance estimate after subtracting Magpylib and therefore includes simultaneous robot-wall contact; it is not falsely presented as a pair-isolated Abaqus contact output. Comparison is only against `S4_HEADFORWARD_F60_G0P10_FLOW`. This run was stopped at {time[-1]*1e3:.3f} ms before the requested 33.333 ms endpoint, so all metrics and GIF frames are partial-run screening evidence.
"""
    (OUT / f"{JOB}_Screen.md").write_text(report, encoding="ascii")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
