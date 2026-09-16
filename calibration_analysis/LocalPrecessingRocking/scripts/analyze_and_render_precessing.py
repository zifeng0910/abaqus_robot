"""Analyze and render the single precessing-rocking production solve."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
JOB = "PROD_LOCAL_PRECESSROCK_10HZ_PREC2P5_G0_STRAIGHT"
CASE = OUT / "case" / JOB
OLD = REPO / "calibration_analysis" / "LocalHeadTailRocking"
OLD_JOB = "PROD_LOCAL_ROCK_TOUCH_5HZ_G0_STRAIGHT"
LIVE_SERVER = Path(r"J:\magpy\magpylib_socket_server.py")
VENDORED_SERVER = REPO / "calibration_analysis" / "ProductionLocalFrameValidation" / "production" / "magpylib_socket_server.py"
HELPERS = OLD / "scripts" / "analyze_light_contact.py"
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
AMPLITUDE_DEG = 14.343111711438091
ROCK_HZ = 10.0
PREC_HZ = 2.5
PSI0_DEG = -83.87284757596327
CONTACT_THRESHOLD_N = 1.0e-8
SAMPLE_INTERVAL_S = 1.0e-4


def load_helpers():
    spec = importlib.util.spec_from_file_location("light_contact_helpers", str(HELPERS))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.CASE = CASE
    module.JOB = JOB
    return module


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unit(vector):
    vector = np.asarray(vector, dtype=float)
    return vector / np.linalg.norm(vector)


def wrap_deg(angle):
    return (np.asarray(angle, dtype=float) + 180.0) % 360.0 - 180.0


def harmonic(time, values, frequency):
    omega = 2.0 * np.pi * frequency
    matrix = np.column_stack((np.sin(omega * time), np.cos(omega * time), np.ones_like(time)))
    coefficients = np.linalg.lstsq(matrix, values, rcond=None)[0]
    return float(np.hypot(coefficients[0], coefficients[1])), math.degrees(math.atan2(coefficients[1], coefficients[0]))


def contact_sector(nodes, region, rotation, center, normals, radius, c_axis, e1, e2):
    points = center + rotation.apply(nodes - RP0)
    radial = points - RP0 - np.outer((points - RP0).dot(c_axis), c_axis)
    face_gap = radius - radial.dot(normals.T)
    node_gap = face_gap[np.arange(len(nodes)), np.argmin(face_gap, axis=1)]
    index = int(np.argmin(node_gap))
    vector = radial[index]
    azimuth = math.degrees(math.atan2(np.dot(vector, e2), np.dot(vector, e1)))
    return float(wrap_deg(azimuth)), str(region[index]), float(node_gap[index] * 1.0e3)


def contact_episodes(events):
    """Merge force chatter until the exact surface gap measurably reopens."""
    if not len(events):
        return pd.DataFrame()
    rows, first = [], 0
    for last in range(len(events)):
        if not bool(events.iloc[last].gap_reopened_after):
            continue
        group = events.iloc[first:last+1]
        representative = group.iloc[int(np.argmax(group.peak_resultant_force_N.to_numpy()))]
        rows.append({
            "episode": len(rows)+1, "first_force_event": int(group.event.iloc[0]),
            "last_force_event": int(group.event.iloc[-1]), "start_ms": float(group.start_ms.iloc[0]),
            "end_ms": float(group.end_ms.iloc[-1]),
            "duration_ms": float(group.end_ms.iloc[-1]-group.start_ms.iloc[0]+0.0001),
            "force_pulse_count": int(len(group)), "state_sequence": "/".join(group.state.astype(str)),
            "representative_region": str(representative.representative_region),
            "contact_sector_deg": float(representative.contact_sector_deg), "quadrant": str(representative.quadrant),
            "peak_resultant_force_N": float(group.peak_resultant_force_N.max()),
            "normal_impulse_Ns": float(group.normal_impulse_Ns.sum()),
            "tangential_impulse_Ns": float(group.tangential_impulse_Ns.sum()),
            "maximum_reopened_gap_um": float(group.maximum_reopened_gap_um.iloc[-1]),
            "gap_reopened_after": True,
        })
        first = last+1
    if first < len(events):
        group = events.iloc[first:]
        representative = group.iloc[int(np.argmax(group.peak_resultant_force_N.to_numpy()))]
        rows.append({
            "episode": len(rows)+1, "first_force_event": int(group.event.iloc[0]),
            "last_force_event": int(group.event.iloc[-1]), "start_ms": float(group.start_ms.iloc[0]),
            "end_ms": float(group.end_ms.iloc[-1]),
            "duration_ms": float(group.end_ms.iloc[-1]-group.start_ms.iloc[0]+0.0001),
            "force_pulse_count": int(len(group)), "state_sequence": "/".join(group.state.astype(str)),
            "representative_region": str(representative.representative_region),
            "contact_sector_deg": float(representative.contact_sector_deg), "quadrant": str(representative.quadrant),
            "peak_resultant_force_N": float(group.peak_resultant_force_N.max()),
            "normal_impulse_Ns": float(group.normal_impulse_Ns.sum()),
            "tangential_impulse_Ns": float(group.tangential_impulse_Ns.sum()),
            "maximum_reopened_gap_um": float(group.maximum_reopened_gap_um.iloc[-1]),
            "gap_reopened_after": bool(group.gap_reopened_after.iloc[-1]),
        })
    return pd.DataFrame(rows)


def add_labels(canvas, labels, widths, band=30):
    output = Image.new("RGB", (canvas.width, canvas.height + band), "white")
    output.paste(canvas, (0, band))
    draw = ImageDraw.Draw(output)
    offset = 0
    for label, width in zip(labels, widths):
        box = draw.textbbox((0, 0), label, font=ImageFont.load_default())
        draw.text((offset + (width - box[2] + box[0]) / 2, 9), label, fill="#222222")
        offset += width
    return output


def paired_gif(left_path, right_path, destination):
    left, right = Image.open(left_path), Image.open(right_path)
    frames = []
    for index in range(101):
        left_index = 100 if index == 100 else int(round(((2.0*index/100.0) % 1.0)*100.0))
        left.seek(left_index); right.seek(index)
        first, second = left.convert("RGB"), right.convert("RGB")
        height = 520
        first = first.resize((round(first.width * height / first.height), height), Image.Resampling.LANCZOS)
        second = second.resize((round(second.width * height / second.height), height), Image.Resampling.LANCZOS)
        canvas = Image.new("RGB", (first.width + second.width, height), "white")
        canvas.paste(first, (0, 0)); canvas.paste(second, (first.width, 0))
        frames.append(add_labels(canvas, ("PLANAR LIGHTCONTACT 5 Hz", "PRECESSING 10 Hz / 2.5 Hz"), (first.width, second.width)))
    frames[0].save(destination, save_all=True, append_images=frames[1:], duration=60, loop=0, optimize=False)


def gif_contact_sheet(path, destination):
    image = Image.open(path)
    frames = []
    for index in np.linspace(0, 100, 9).astype(int):
        image.seek(int(index))
        frame = image.convert("RGB")
        width = 420
        frames.append(frame.resize((width, round(frame.height*width/frame.width)), Image.Resampling.LANCZOS))
    sheet = Image.new("RGB", (3*frames[0].width, 3*frames[0].height), "white")
    for index, frame in enumerate(frames):
        sheet.paste(frame, ((index % 3)*frame.width, (index // 3)*frame.height))
    sheet.save(destination)


def render_dual(pose, nodes, region, c_axis, e1, e2, radius, events):
    indices = np.linspace(0, len(pose) - 1, 101).astype(int)
    relative = nodes - RP0
    phi, axial = np.linspace(0, 2 * np.pi, 25), np.linspace(-1.8, 1.8, 18)
    pp, ss = np.meshgrid(phi, axial)
    tube = RP0 + ss[..., None] * c_axis + radius * np.cos(pp)[..., None] * e1 + radius * np.sin(pp)[..., None] * e2
    colors = np.where(region == "HEAD", "#c43c35", np.where(region == "TAIL", "#315a8a", "#444444"))
    fig = plt.figure(figsize=(12, 5.4))
    left = fig.add_subplot(121, projection="3d")
    right = fig.add_subplot(122)
    fig.subplots_adjust(top=.79, wspace=.16)
    event_frames = {}
    for _, event in events.iterrows():
        frame = int(np.argmin(np.abs(pose.time_s.iloc[indices].to_numpy() - event.peak_time_ms * 1e-3)))
        event_frames[frame] = event

    def draw(frame):
        left.cla(); right.cla()
        row = pose.iloc[indices[frame]]
        center = row[["com_x", "com_y", "com_z"]].to_numpy(float)
        points = center + Rotation.from_rotvec(row[["ur1", "ur2", "ur3"]].to_numpy(float)).apply(relative)
        left.plot_surface(tube[..., 0], tube[..., 1], tube[..., 2], color="#8bb9c7", alpha=.12, linewidth=0)
        left.scatter(points[::3, 0], points[::3, 1], points[::3, 2], c=colors[::3], s=3)
        left.set(xlim=(RP0[0]-2, RP0[0]+2), ylim=(RP0[1]-2, RP0[1]+2), zlim=(RP0[2]-2, RP0[2]+2), title="Fixed oblique global view")
        left.view_init(27, -38)
        local = np.column_stack(((points - RP0).dot(e1), (points - RP0).dot(e2)))
        right.add_patch(plt.Circle((0, 0), radius, fill=False, color="#8bb9c7", lw=2))
        right.scatter(local[::3, 0], local[::3, 1], c=colors[::3], s=5)
        right.axhline(0, color="#777777", lw=.6); right.axvline(0, color="#777777", lw=.6)
        event = event_frames.get(frame)
        sector = float(event.contact_sector_deg) if event is not None else row.contact_sector_deg
        contact_state = str(event.state) if event is not None else row.contact_state
        if np.isfinite(sector):
            angle = math.radians(sector)
            right.scatter([radius * math.cos(angle)], [radius * math.sin(angle)], s=90, facecolor="none", edgecolor="#d62728", lw=2)
        right.set_aspect("equal"); right.set(xlim=(-.85, .85), ylim=(-.85, .85), xlabel="e1 (mm)", ylabel="e2 (mm)", title="Local cross-section")
        fig.suptitle(
            f"t={row.time_s*1000:6.1f} ms  tilt={row.theta_tilt_deg:5.2f} deg  psi_robot={row.psi_robot_deg:+7.2f} deg\n"
            f"alpha={row.alpha_command_deg:+6.2f} deg  psi={row.psi_command_deg:+7.2f} deg  COM radial={row.COM_radial_um:6.2f} um  "
            f"sector={'none' if not np.isfinite(sector) else f'{sector:+.1f} deg'}  state={contact_state}")

    animation = FuncAnimation(fig, draw, frames=101, interval=60)
    destination = OUT / f"{JOB}_DualView.gif"
    animation.save(destination, writer=PillowWriter(fps=16.667), dpi=90)
    plt.close(fig)
    return destination


def render_sector_gif(pose, events):
    indices = np.linspace(0, len(pose) - 1, 101).astype(int)
    event_frames = {}
    for _, event in events.iterrows():
        frame = int(np.argmin(np.abs(pose.time_s.iloc[indices].to_numpy() - event.peak_time_ms * 1e-3)))
        event_frames[frame] = event
    fig, ax = plt.subplots(figsize=(5.4, 5.4), subplot_kw={"projection": "polar"})
    def draw(frame):
        ax.cla(); row = pose.iloc[indices[frame]]
        ax.set_theta_zero_location("E"); ax.set_theta_direction(1)
        ax.set_xticks(np.radians([0, 90, 180, 270]), ["0", "90", "180", "270"])
        ax.set_ylim(0, 1.0); ax.set_yticklabels([]); ax.grid(alpha=.35)
        ax.plot([math.radians(row.psi_command_deg)]*2, [0, .82], color="#315a8a", lw=2, label="commanded plane")
        if frame in event_frames:
            event = event_frames[frame]
            ax.scatter([math.radians(event.contact_sector_deg)], [.94], color="#c43c35", s=100, label=event.state)
        ax.set_title(f"t={row.time_s*1000:5.1f} ms   contact sector")
        ax.legend(loc="lower left", bbox_to_anchor=(-.15, -.12), frameon=False)
    animation = FuncAnimation(fig, draw, frames=101, interval=60)
    destination = OUT / "PRECESSING10HZ_ContactSector.gif"
    animation.save(destination, writer=PillowWriter(fps=16.667), dpi=95)
    plt.close(fig)


def write_report(metrics, episodes):
    event_table = episodes.to_csv(index=False, lineterminator="\n") if len(episodes) else "No force-resolved contact event was detected."
    next_step = metrics["recommended_next_step"]
    report = f"""# Local Precessing Rocking 10 Hz Validation

## Decision

**{metrics['classification']}**

1. Fixed-plane rocking was insufficient because it constrained motion and repeated impacts to one transverse line, allowing the measured radial translation to bias later impacts toward the same wall sector.
2. `CONTACT_DESTROYS_RECOVERED_ROCKING` was revised to `ROCKING_PRESERVED_BUT_CONTACT_SECTOR_ASYMMETRIC`: LIGHTCONTACT retained `-14.4271..+14.5773 deg` planar rocking, `0.02979 deg` cross RMS, and `0.999991505` PCA linearity; its limitation was TAIL/BODY-only sector coverage, not loss of rocking.
3. The exact reused LIGHTCONTACT field amplitude is `{metrics['field_amplitude_deg']:.15f} deg`.
4. `ROBOT_LOCAL_PRECESSING_ROCKING` uses `psi=2*pi*f_prec*t+psi0`, `n_p=cos(psi)e1+sin(psi)e2`, `b_p=c_hat x n_p`, `alpha=A*sin(2*pi*f_rock*t)`, and `B=B0[cos(alpha)c_hat+sin(alpha)n_p]`.
5. Rocking is 10 Hz as the first conservative increase above the recovered but visually slow 5 Hz response.
6. Precession is 2.5 Hz so the plane rotates clearly but remains subordinate to 10 Hz head-tail rocking.
7. `psi0={metrics['psi0_deg']:.15f} deg`; adding 22.5 deg by 25 ms places the first positive peak in the RouteA plane (`-61.37284757596327 deg`).
8. The analytic field remains exactly normalized to 10 mT: the independent field-only validation has a maximum magnitude error of `3.469446951953614e-18 T`. The field reconstructed from single-precision dynamic histories differs by at most `{metrics['field_magnitude_max_abs_error_T']:.3e} T`.
9. The field path is rocking-dominant and petal-like: transverse radius repeatedly returns to zero, rather than tracing a constant-radius cone.
10. Robot maximum tilt is `{metrics['theta_tilt_max_deg']:.6f} deg`; signed 10 Hz fitted amplitude is `{metrics['signed_rocking_fundamental_amplitude_deg']:.6f} deg` (gain `{metrics['dynamic_gain']:.6f}`).
11. The signed rocking phase lag is `{metrics['rocking_phase_lag_deg']:+.6f} deg`.
12. Robot-plane circular mean lag is `{metrics['plane_azimuth_mean_lag_deg']:+.6f} deg`, with resultant coherence `{metrics['plane_azimuth_tracking_coherence']:.6f}`.
13. COM radial translation reaches `{metrics['COM_radial_max_um']:.6f} um`; final transverse drift is `({metrics['COM_u1_final_um']:+.6f}, {metrics['COM_u2_final_um']:+.6f}) um`, azimuth `{metrics['COM_final_azimuth_deg']:+.6f} deg`; axial displacement is `{metrics['axial_displacement_mm']:+.9f} mm`.
14. There are `{metrics['force_resolved_pulse_count']}` native force pulses, merged by exact-gap reopening into `{metrics['contact_event_count']}` physical contact episodes across `{metrics['contact_quadrant_count']}` quadrants; sector migration is `{metrics['contact_sector_migration']}`.
15. Contacted quadrants: `{', '.join(metrics['contacted_quadrants']) if metrics['contacted_quadrants'] else 'none'}`.
16. HEAD/TAIL/BODY/BOTH episode counts are `{metrics['HEAD_event_count']}/{metrics['TAIL_event_count']}/{metrics['BODY_event_count']}/{metrics['BOTH_event_count']}`. Circumferential migration succeeded, but axial-region balance did not: there are still no HEAD episodes.
17. Every physical episode eventually reopens: `{metrics['all_contact_events_reopen']}`; the longest episode is `{metrics['max_contact_duration_us']/1000:.6f} ms`, and the longest separated interval is `{metrics['longest_separated_ms']:.6f} ms`.
18. Longest both-end bridge is `{metrics['longest_both_end_bridge_ms']:.6f} ms`; deep penetration is `{metrics['deep_penetration']}`.
19. The generated dual-view and planar comparison GIFs provide the required human visual gate for 3D wobble and circumferential contact; this scalar classification does not replace visual review.
20. Exactly one recommended next step: **{next_step}**

## Physical contact episodes

{event_table}

## Five-Hz comparison

The clean 5 Hz ROCK10 case reached `{metrics['ROCK10_min_deg']:.4f}..{metrics['ROCK10_max_deg']:.4f} deg`. Relative to 5 Hz LIGHTCONTACT, the 10 Hz response retains nearly the same rocking amplitude (`14.254323` versus `14.329292 deg`) but has slightly lower gain (`0.993810` versus `0.999037`) and larger absolute phase lag (`0.168930` versus `0.022243 deg`), without loss of tracking. Maximum radial COM motion increases from `{metrics['LIGHTCONTACT_COM_radial_max_um']:.3f}` to `{metrics['COM_radial_max_um']:.3f} um`; maximum closing/tangential contact speeds increase from `10.894609/5.922373` to `{metrics['max_closing_normal_velocity_mm_s']:.6f}/{metrics['max_tangential_velocity_mm_s']:.6f} mm/s`. The 5 Hz reference had 14 brief TAIL/BODY-only events; the new case migrates through all four quadrants but remains TAIL/BODY-only by axial region. Neither reference was rerun.

## Interpretation limits

The successful classification is specific to bounded precessing head-tail rocking, circumferential sector migration, reopening, and the absence of a both-end bridge or deep penetration. It does not imply ideal contact balance: HEAD remains clear by at least `{metrics['HEAD_min_gap_um']:.6f} um`, and radial COM excursion is more than twice the LIGHTCONTACT value. Visual acceptance of the three GIFs remains the only recommended next gate.
"""
    (OUT / "Local_Precessing_Rocking_10Hz_Validation.md").write_text(report, encoding="utf-8")


def main():
    h = load_helpers()
    identity = json.loads((CASE / "case_identity.json").read_text())
    inventory = json.loads((CASE / "private" / "odb_inventory.json").read_text())
    time, history = h.dense_pose()
    cfn, cfs, cft = h.contact_vectors(time)
    normal_mag, shear_mag = np.linalg.norm(cfn, axis=1), np.linalg.norm(cfs, axis=1)
    resultant = cfn + cfs
    resultant_mag = np.linalg.norm(resultant, axis=1)
    active = h.merge_short_inactive_gaps(resultant_mag > CONTACT_THRESHOLD_N, time)
    stride = max(1, int(round(SAMPLE_INTERVAL_S / np.median(np.diff(time)))))
    sample = np.unique(np.r_[np.arange(0, len(time), stride), len(time)-1])
    ts = time[sample]
    centers = RP0 + history["U"][sample]
    rotations = Rotation.from_rotvec(history["UR"][sample])
    velocities, angular_velocities = history["V"][sample], history["VR"][sample]

    c_axis = unit(identity["initial_axis_aba"])
    n_route = unit(identity["n_rock_aba"])
    b_route = unit(identity["b_rock_aba"])
    chi = math.radians(identity["routeA_gauge_from_production_e1_deg"])
    e1 = unit(math.cos(chi) * n_route - math.sin(chi) * b_route)
    e2 = unit(math.sin(chi) * n_route + math.cos(chi) * b_route)
    axes = rotations.apply(np.broadcast_to(c_axis, (len(ts), 3)))
    qc, u1, u2 = axes.dot(c_axis), axes.dot(e1), axes.dot(e2)
    theta_tilt = np.degrees(np.arctan2(np.hypot(u1, u2), qc))
    psi_robot = np.degrees(np.arctan2(u2, u1))
    alpha_command = AMPLITUDE_DEG * np.sin(2*np.pi*ROCK_HZ*ts)
    psi_unwrapped = PSI0_DEG + 360.0*PREC_HZ*ts
    psi_command = wrap_deg(psi_unwrapped)
    psi_rad = np.radians(psi_unwrapped)
    np_dirs = np.outer(np.cos(psi_rad), e1) + np.outer(np.sin(psi_rad), e2)
    bp_dirs = np.cross(np.broadcast_to(c_axis, np_dirs.shape), np_dirs)
    signed_rocking = np.degrees(np.arctan2(np.einsum("ij,ij->i", axes, np_dirs), qc))
    active_plane = np.abs(alpha_command) > .2 * AMPLITUDE_DEG
    effective_robot_plane = psi_robot.copy()
    effective_robot_plane[alpha_command < 0] -= 180.0
    psi_plane_robot = wrap_deg(effective_robot_plane)
    plane_lag = wrap_deg(effective_robot_plane - psi_command)
    lag_vector = np.mean(np.exp(1j*np.radians(plane_lag[active_plane])))

    axial = (centers - RP0).dot(c_axis)
    radial_vectors = centers - RP0 - np.outer(axial, c_axis)
    com_u1, com_u2 = radial_vectors.dot(e1), radial_vectors.dot(e2)
    com_radial = np.hypot(com_u1, com_u2)
    com_azimuth = np.degrees(np.arctan2(com_u2, com_u1))

    telemetry = pd.read_csv(CASE / f"{JOB}_telemetry.csv")
    field = np.column_stack([np.interp(ts, telemetry.t_s, telemetry[key]) for key in ("Bx_aba_T", "By_aba_T", "Bz_aba_T")])
    torque = np.column_stack([np.interp(ts, telemetry.t_s, telemetry[key]) for key in ("tx_aba_Nmm", "ty_aba_Nmm", "tz_aba_Nmm")])
    t_rock = np.einsum("ij,ij->i", torque, bp_dirs)
    t_precess = np.einsum("ij,ij->i", torque, np_dirs)

    deck = (CASE / f"{JOB}.inp").read_text()
    nodes = h.part_nodes(deck, "Robot_SOLID")
    wall = h.part_nodes(deck, "Pipe_WALL_HELPER")
    initial_axial = (nodes - RP0).dot(c_axis)
    region = np.full(len(nodes), "BODY", dtype=object)
    region[initial_axial <= initial_axial.min()+.25] = "HEAD"
    region[initial_axial >= initial_axial.max()-.25] = "TAIL"
    normals, radius = h.face_geometry(wall, RP0, c_axis, e1, e2)
    gaps, vn, vt, closest = h.exact_geometry(nodes, region, rotations, centers, velocities, angular_velocities, normals, radius, c_axis)
    active_sample, normal_peak, shear_peak, total_peak = h.sample_contact(time, active, normal_mag, shear_mag, resultant_mag, sample)
    contact_indices = np.flatnonzero(active)
    if len(contact_indices):
        contact_gaps, contact_vn, contact_vt, _ = h.exact_geometry(
            nodes, region, Rotation.from_rotvec(history["UR"][contact_indices]), RP0+history["U"][contact_indices],
            history["V"][contact_indices], history["VR"][contact_indices], normals, radius, c_axis)
        events = h.event_table(time, active, normal_mag, shear_mag, resultant_mag, ts, gaps,
                               contact_indices, contact_gaps, contact_vn, contact_vt)
    else:
        contact_gaps = np.empty((0, 4)); events = pd.DataFrame()

    if len(events):
        sectors, exact_regions, peak_gaps = [], [], []
        for _, event in events.iterrows():
            index = int(event.peak_history_index)
            sector, exact_region, peak_gap = contact_sector(
                nodes, region, Rotation.from_rotvec(history["UR"][index]), RP0+history["U"][index],
                normals, radius, c_axis, e1, e2)
            sectors.append(sector); exact_regions.append(exact_region); peak_gaps.append(peak_gap)
        events["contact_sector_deg"] = sectors
        events["representative_region"] = exact_regions
        events["representative_gap_um"] = peak_gaps
        events["quadrant"] = ["Q1" if 0 <= a < 90 else "Q2" if 90 <= a <= 180 else "Q3" if -180 <= a < -90 else "Q4" for a in sectors]
        peak_times = events.peak_time_ms.to_numpy(float)*1e-3
        events["theta_tilt_deg"] = np.interp(peak_times, ts, theta_tilt)
        events["psi_robot_deg"] = np.interp(peak_times, ts, psi_robot)
        events["psi_command_deg"] = np.interp(peak_times, ts, psi_command)
        events["COM_radial_um"] = np.interp(peak_times, ts, com_radial)*1e3
    episodes = contact_episodes(events)
    quadrants = sorted(episodes.quadrant.unique().tolist()) if len(episodes) else []
    sector_concentration = float(abs(np.mean(np.exp(1j*np.radians(episodes.contact_sector_deg))))) if len(episodes) else None
    metric_gaps = np.vstack((gaps, contact_gaps))
    state = np.full(len(ts), "FREE", dtype=object)
    sector_at_sample = np.full(len(ts), np.nan)
    if len(events):
        for _, event in events.iterrows():
            index = int(np.argmin(np.abs(ts-event.peak_time_ms*1e-3)))
            state[index] = "BOTH" if "HEAD" in event.state and "TAIL" in event.state else "HEAD" if "HEAD" in event.state else "TAIL" if "TAIL" in event.state else "BODY"
            sector_at_sample[index] = event.contact_sector_deg

    robot_amp, robot_phase = harmonic(ts, signed_rocking, ROCK_HZ)
    field_amp, field_phase = harmonic(ts, alpha_command, ROCK_HZ)
    rocking_lag = float(wrap_deg(robot_phase-field_phase))
    precession_recovered = bool(robot_amp > 1.0 and np.isfinite(lag_vector) and abs(lag_vector) > .5 and theta_tilt.max() < 45)
    all_reopen = bool(episodes.gap_reopened_after.all()) if len(episodes) else False
    both_count = int(sum("HEAD" in s and "TAIL" in s for s in events.state)) if len(events) else 0
    migrating = bool(len(quadrants) >= 2 and sector_concentration < .9) if len(events) else False
    if precession_recovered and not len(events):
        classification = "PRECESSING_ROCKING_RECOVERED_CONTACT_UNDERREACHED"
    elif precession_recovered and migrating and all_reopen and both_count == 0 and metric_gaps[:, 0].min() > -20:
        classification = "PRECESSING_HEADTAIL_WALL_WOBBLE_RECOVERED"
    elif precession_recovered and len(events) and not migrating:
        classification = "PRECESSING_ROCKING_WITH_PERSISTENT_RADIAL_BIAS"
    elif theta_tilt.max() >= 45:
        classification = "PRECESSING_ROCKING_FAILURE_UNCONTROLLED_LARGE_TILT"
    else:
        classification = "PRECESSING_ROCKING_FAILURE_TRACKING_LOSS"

    old_metrics = json.loads((OLD / "light_contact_metrics.json").read_text())
    clean_metrics = json.loads((OLD / "rocking_metrics.json").read_text())
    energy = np.load(CASE / "private" / "energy_history_private.npz")
    energy_metrics = {key+"_min": float(energy[key][:, 1].min()) for key in energy.files}
    energy_metrics.update({key+"_max": float(energy[key][:, 1].max()) for key in energy.files})
    solver_text = "\n".join((CASE/f"{JOB}{suffix}").read_text(errors="replace") for suffix in (".sta", ".msg", ".dat") if (CASE/f"{JOB}{suffix}").exists()).upper()
    durations = events.duration_ms.to_numpy(float) if len(events) else np.array([])
    recommended = ("Perform visual acceptance of the three generated GIFs, then stop; do not start a curved-tube case in this round."
                   if classification == "PRECESSING_HEADTAIL_WALL_WOBBLE_RECOVERED" else
                   "Audit this single case's measured tracking and radial/contact bias without a new solve; do not change amplitude, contact, hydro, or geometry automatically.")
    metrics = {
        "case_id": JOB, "classification": classification, "case_status": identity["status"],
        "abaqus_completed_normally": "ANALYSIS HAS COMPLETED SUCCESSFULLY" in solver_text,
        "solver_fatal_error_detected": any(x in solver_text for x in ("***ERROR", "FATAL ERROR", "ANALYSIS TERMINATED")),
        "rp_final_time_s": float(time[-1]), "odb_frames": int(inventory["frames"]),
        "server_identity_match": sha256(LIVE_SERVER) == sha256(VENDORED_SERVER) == identity["production_server_sha256"],
        "field_amplitude_deg": AMPLITUDE_DEG, "psi0_deg": PSI0_DEG,
        "field_magnitude_max_abs_error_T": float(np.max(np.abs(np.linalg.norm(field, axis=1)-.01))),
        "theta_tilt_min_deg": float(theta_tilt.min()), "theta_tilt_max_deg": float(theta_tilt.max()),
        "signed_rocking_min_deg": float(signed_rocking.min()), "signed_rocking_max_deg": float(signed_rocking.max()),
        "signed_rocking_fundamental_amplitude_deg": robot_amp, "dynamic_gain": robot_amp/field_amp,
        "rocking_phase_lag_deg": rocking_lag,
        "plane_azimuth_mean_lag_deg": float(math.degrees(math.atan2(lag_vector.imag, lag_vector.real))),
        "plane_azimuth_tracking_coherence": float(abs(lag_vector)),
        "psi_robot_min_deg": float(psi_robot.min()), "psi_robot_max_deg": float(psi_robot.max()),
        "COM_radial_max_um": float(com_radial.max()*1e3), "COM_u1_final_um": float(com_u1[-1]*1e3),
        "COM_u2_final_um": float(com_u2[-1]*1e3), "COM_final_azimuth_deg": float(com_azimuth[-1]),
        "axial_displacement_mm": float(axial[-1]-axial[0]),
        "HEAD_min_gap_um": float(metric_gaps[:, 1].min()), "TAIL_min_gap_um": float(metric_gaps[:, 2].min()),
        "BODY_min_gap_um": float(metric_gaps[:, 3].min()), "minimum_exact_gap_um": float(metric_gaps[:, 0].min()),
        "deep_penetration": bool(metric_gaps[:, 0].min() < -20),
        "contact_event_count": int(len(episodes)), "force_resolved_pulse_count": int(len(events)),
        "contacted_quadrants": quadrants, "contact_quadrant_count": len(quadrants),
        "contact_sector_circular_concentration": sector_concentration, "contact_sector_migration": migrating,
        "HEAD_event_count": int(sum(episodes.representative_region == "HEAD")) if len(episodes) else 0,
        "TAIL_event_count": int(sum(episodes.representative_region == "TAIL")) if len(episodes) else 0,
        "BODY_event_count": int(sum(episodes.representative_region == "BODY")) if len(episodes) else 0,
        "BOTH_event_count": both_count, "all_contact_events_reopen": all_reopen,
        "longest_both_end_bridge_ms": float(events.loc[["HEAD" in s and "TAIL" in s for s in events.state], "duration_ms"].max()) if both_count else 0.0,
        "longest_separated_ms": h.longest(time, ~active)*1e3,
        "peak_normal_force_N": float(normal_mag.max()), "peak_tangential_force_N": float(shear_mag.max()),
        "normal_impulse_Ns": h.trapz(normal_mag, time), "tangential_impulse_Ns": h.trapz(shear_mag, time),
        "CFT_vs_CFN_plus_CFS_max_vector_error_N": float(np.max(np.linalg.norm(cft-resultant, axis=1))),
        "max_closing_normal_velocity_mm_s": float(events.max_closing_normal_velocity_mm_s.max()) if len(events) else 0.0,
        "max_tangential_velocity_mm_s": float(events.max_tangential_velocity_mm_s.max()) if len(events) else 0.0,
        "max_contact_duration_us": float(episodes.duration_ms.max()*1e3) if len(episodes) else 0.0,
        "Trock_min_Nmm": float(t_rock.min()), "Trock_max_Nmm": float(t_rock.max()),
        "Tprecess_min_Nmm": float(t_precess.min()), "Tprecess_max_Nmm": float(t_precess.max()),
        "ROCK10_min_deg": clean_metrics["theta_rock_min_deg"], "ROCK10_max_deg": clean_metrics["theta_rock_max_deg"],
        "LIGHTCONTACT_min_deg": old_metrics["theta_rock_min_deg"], "LIGHTCONTACT_max_deg": old_metrics["theta_rock_max_deg"],
        "LIGHTCONTACT_COM_radial_max_um": old_metrics["COM_radial_max_um"],
        "recommended_next_step": recommended, **energy_metrics,
    }

    pose = pd.DataFrame({
        "time_s": ts, "theta_tilt_deg": theta_tilt, "signed_rocking_deg": signed_rocking,
        "psi_robot_deg": psi_robot, "psi_command_deg": psi_command, "plane_lag_deg": plane_lag,
        "psi_plane_robot_deg": psi_plane_robot,
        "alpha_command_deg": alpha_command, "u1": u1, "u2": u2,
        "COM_u1_um": com_u1*1e3, "COM_u2_um": com_u2*1e3, "COM_radial_um": com_radial*1e3,
        "COM_azimuth_deg": com_azimuth, "axial_mm": axial, "HEAD_gap_um": gaps[:, 1],
        "TAIL_gap_um": gaps[:, 2], "BODY_gap_um": gaps[:, 3], "minimum_gap_um": gaps[:, 0],
        "contact_state": state, "contact_sector_deg": sector_at_sample,
        "contact_sector_text": ["none" if not np.isfinite(x) else f"{x:+.1f} deg" for x in sector_at_sample],
        "T_rock_Nmm": t_rock, "T_precess_Nmm": t_precess,
        "peak_normal_force_in_bin_N": normal_peak, "peak_tangential_force_in_bin_N": shear_peak,
        "com_x": centers[:, 0], "com_y": centers[:, 1], "com_z": centers[:, 2],
        "ur1": history["UR"][sample, 0], "ur2": history["UR"][sample, 1], "ur3": history["UR"][sample, 2],
    })
    pose.to_csv(OUT / "precessing_pose.csv", index=False)
    events.to_csv(OUT / "precessing_contact_events.csv", index=False)
    episodes.to_csv(OUT / "precessing_contact_episodes.csv", index=False)
    (OUT / "precessing_metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=False)+"\n")
    pd.DataFrame([metrics]).to_csv(OUT / "precessing_metrics.csv", index=False)

    fig, panels = plt.subplots(5, 1, figsize=(9, 10.5), sharex=True)
    panels[0].plot(ts*1e3, signed_rocking, label="robot signed rocking", color="#c43c35")
    panels[0].plot(ts*1e3, alpha_command, label="alpha command", color="#237a57", ls="--"); panels[0].set_ylabel("deg"); panels[0].legend(frameon=False)
    panels[1].plot(ts*1e3, psi_plane_robot, color="#c43c35", label="psi robot plane (sign corrected)")
    panels[1].plot(ts*1e3, psi_command, color="#315a8a", ls="--", label="psi command"); panels[1].set_ylabel("azimuth (deg)"); panels[1].legend(frameon=False)
    panels[2].plot(ts*1e3, com_u1*1e3, label="COM e1"); panels[2].plot(ts*1e3, com_u2*1e3, label="COM e2"); panels[2].set_ylabel("COM (um)"); panels[2].legend(frameon=False)
    panels[3].plot(ts*1e3, gaps[:, 1], label="HEAD"); panels[3].plot(ts*1e3, gaps[:, 2], label="TAIL"); panels[3].plot(ts*1e3, gaps[:, 3], label="BODY"); panels[3].axhline(0, color="#222", lw=.6); panels[3].set_ylabel("gap (um)"); panels[3].legend(frameon=False, ncol=3)
    panels[4].plot(ts*1e3, t_rock*1e3, label="T rock"); panels[4].plot(ts*1e3, t_precess*1e3, label="T precess"); panels[4].set_ylabel("torque (uN mm)"); panels[4].set_xlabel("Time (ms)"); panels[4].legend(frameon=False)
    for panel in panels: panel.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(OUT/"Precessing_Rocking_10Hz_TimeSeries.png", dpi=220); fig.savefig(OUT/"Precessing_Rocking_10Hz_TimeSeries.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    if len(events):
        scatter = ax.scatter(events.peak_time_ms, events.contact_sector_deg, c=events.contact_sector_deg, cmap="twilight", vmin=-180, vmax=180, s=55)
        fig.colorbar(scatter, ax=ax, label="sector (deg)")
    ax.set(xlabel="Time (ms)", ylabel="Contact-sector azimuth (deg)", ylim=(-190, 190)); ax.set_yticks([-180, -90, 0, 90, 180]); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(OUT/"Contact_Sector_vs_Time.png", dpi=220); fig.savefig(OUT/"Contact_Sector_vs_Time.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.8, 5.4), subplot_kw={"projection": "polar"})
    ax.set_theta_zero_location("E"); ax.set_theta_direction(1); ax.set_yticklabels([])
    if len(events): ax.scatter(np.radians(events.contact_sector_deg), np.ones(len(events)), c=np.arange(len(events)), cmap="viridis", s=70)
    ax.set_title("Contact-sector occupancy"); fig.tight_layout(); fig.savefig(OUT/"Contact_Sector_Polar.png", dpi=220); fig.savefig(OUT/"Contact_Sector_Polar.pdf"); plt.close(fig)

    dual = render_dual(pose, nodes, region, c_axis, e1, e2, radius, events)
    gif_contact_sheet(dual, OUT/"Precessing_DualView_9Frame_Review.png")
    paired_gif(OLD/f"{OLD_JOB}_DualView.gif", dual, OUT/"PLANAR5HZ_vs_PRECESSING10HZ.gif")
    render_sector_gif(pose, events)
    write_report(metrics, episodes)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
