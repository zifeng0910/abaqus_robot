"""Analyze and render the single amplitude-derived light-contact solve."""
from __future__ import annotations

import hashlib
import json
import math
import re
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
JOB = "PROD_LOCAL_ROCK_TOUCH_5HZ_G0_STRAIGHT"
CASE = OUT / "case" / JOB
ROCK10_JOB = "PROD_LOCAL_ROCK10_5HZ_G0_STRAIGHT"
ROCK10_CASE = OUT / "case" / ROCK10_JOB
LEGACY = REPO / "calibration_analysis" / "Legacy_Wobble_GIFs" / "Job_RouteA_CEL_SOLID_headtail_rock_probe.gif"
LIVE_SERVER = Path(r"J:\magpy\magpylib_socket_server.py")
VENDORED_SERVER = REPO / "calibration_analysis" / "ProductionLocalFrameValidation" / "production" / "magpylib_socket_server.py"
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
CONTACT_THRESHOLD_N = 1.0e-8
SAMPLE_INTERVAL_S = 1.0e-4


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def get_array(archive, name):
    keys = [key for key in archive.files if key == name or key.startswith(name + " (Repeated:")]
    if not keys:
        raise KeyError(name)
    return max((archive[key] for key in keys), key=len)


def dense_pose():
    archive = np.load(CASE / "private" / "rp_history_private.npz")
    time = get_array(archive, "V1")[:, 0]
    data = {}
    for prefix in ("U", "UR", "V", "VR"):
        columns = []
        for axis in (1, 2, 3):
            array = get_array(archive, prefix + str(axis))
            columns.append(np.interp(time, array[:, 0], array[:, 1]))
        data[prefix] = np.column_stack(columns)
    return time, data


def contact_vectors(time):
    archive = np.load(CASE / "private" / "contact_history_private.npz")

    def vector(prefix):
        columns = []
        for axis in (1, 2, 3):
            keys = [
                key for key in archive.files
                if "|" + prefix + str(axis) + " on surface " in key and "ASSEMBLY_ROBOT" in key
            ]
            array = max((archive[key] for key in keys), key=len)
            columns.append(np.interp(time, array[:, 0], array[:, 1]))
        return np.column_stack(columns)

    return vector("CFN"), vector("CFS"), vector("CFT")


def part_nodes(deck, name):
    part = re.search(rf"(?ms)^\*Part, name={re.escape(name)}\s*$.*?^\*End Part\s*$", deck).group(0)
    block = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", part).group(1)
    return np.asarray([[float(x) for x in row.split(",")[1:4]] for row in block.splitlines() if row.strip()])


def unit(vector):
    vector = np.asarray(vector, dtype=float)
    return vector / np.linalg.norm(vector)


def face_geometry(wall, center, c_axis, n_axis, b_axis):
    axial = (wall - center).dot(c_axis)
    ring = wall[np.abs(axial - axial.min()) < 1.0e-7]
    points = np.column_stack(((ring - center).dot(n_axis), (ring - center).dot(b_axis)))
    points = points[np.argsort(np.arctan2(points[:, 1], points[:, 0]))]
    normals, radii = [], []
    for first, second in zip(points, np.roll(points, -1, axis=0)):
        edge = second - first
        normal2 = unit([edge[1], -edge[0]])
        if np.dot(normal2, first + second) < 0:
            normal2 *= -1
        normals.append(normal2[0] * n_axis + normal2[1] * b_axis)
        radii.append(np.dot(first, normal2))
    return np.asarray(normals), float(np.mean(radii))


def harmonic(time, values, frequency=5.0):
    omega = 2.0 * np.pi * frequency
    matrix = np.column_stack((np.sin(omega * time), np.cos(omega * time), np.ones_like(time)))
    coef = np.linalg.lstsq(matrix, values, rcond=None)[0]
    amplitude = float(np.hypot(coef[0], coef[1]))
    phase = math.degrees(math.atan2(coef[1], coef[0]))
    linear_drift = float(np.polyfit(time, values - matrix.dot(coef), 1)[0])
    return amplitude, phase, float(coef[2]), linear_drift


def intervals(flag):
    edges = np.diff(np.r_[False, flag, False].astype(int))
    return list(zip(np.where(edges == 1)[0], np.where(edges == -1)[0] - 1))


def duration(time, start, end):
    dt = float(time[1] - time[0])
    return float(time[end] - time[start] + dt)


def longest(time, flag):
    return max((duration(time, first, last) for first, last in intervals(flag)), default=0.0)


def trapz(values, time):
    return float(np.trapezoid(values, time)) if hasattr(np, "trapezoid") else float(np.trapz(values, time))


def merge_short_inactive_gaps(active, time, maximum_gap_s=1.0e-5):
    merged = active.copy()
    for first, last in intervals(~active):
        if first == 0 or last == len(active) - 1:
            continue
        if duration(time, first, last) <= maximum_gap_s:
            merged[first:last + 1] = True
    return merged


def exact_geometry(nodes, region, rotations, centers, velocities, angular_velocities, normals, radius, c_axis):
    gaps = np.empty((len(centers), 4))
    normal_velocity = np.empty((len(centers), 3))
    tangential_velocity = np.empty((len(centers), 3))
    closest_labels = np.empty(len(centers), dtype=object)
    region_names = ("HEAD", "TAIL", "BODY")
    relative = nodes - RP0
    for row, (center, rotation, velocity, omega) in enumerate(zip(centers, rotations, velocities, angular_velocities)):
        rotated = rotation.apply(relative)
        points = center + rotated
        radial = points - RP0 - np.outer((points - RP0).dot(c_axis), c_axis)
        face_gap = radius - radial.dot(normals.T)
        face_index = np.argmin(face_gap, axis=1)
        node_gap = face_gap[np.arange(len(nodes)), face_index]
        gaps[row, 0] = node_gap.min()
        for column, name in enumerate(region_names, start=1):
            indices = np.where(region == name)[0]
            local = indices[np.argmin(node_gap[indices])]
            gaps[row, column] = node_gap[local]
            wall_normal = normals[face_index[local]]
            point_velocity = velocity + np.cross(omega, rotated[local])
            outward_speed = float(np.dot(point_velocity, wall_normal))
            normal_velocity[row, column - 1] = outward_speed
            tangential_velocity[row, column - 1] = np.linalg.norm(point_velocity - outward_speed * wall_normal)
        closest_labels[row] = region_names[int(np.argmin(gaps[row, 1:]))]
    return gaps * 1.0e3, normal_velocity, tangential_velocity, closest_labels


def sample_contact(time, active, normal_mag, shear_mag, total_mag, sample):
    edges = np.r_[0, ((sample[:-1] + sample[1:]) // 2), len(time) - 1]
    active_bin = np.zeros(len(sample), dtype=bool)
    normal_peak = np.zeros(len(sample))
    shear_peak = np.zeros(len(sample))
    total_peak = np.zeros(len(sample))
    for index in range(len(sample)):
        first, last = int(edges[index]), int(edges[index + 1]) + 1
        active_bin[index] = active[first:last].any()
        normal_peak[index] = normal_mag[first:last].max(initial=0.0)
        shear_peak[index] = shear_mag[first:last].max(initial=0.0)
        total_peak[index] = total_mag[first:last].max(initial=0.0)
    return active_bin, normal_peak, shear_peak, total_peak


def support_states(active, gaps):
    head = active & (gaps[:, 1] <= 20.0)
    tail = active & (gaps[:, 2] <= 20.0)
    body = active & ~(head | tail)
    both = head & tail
    head_only = head & ~tail
    tail_only = tail & ~head
    free = ~active
    state = np.where(both, "BOTH", np.where(head_only, "HEAD", np.where(tail_only, "TAIL", np.where(body, "BODY", "FREE"))))
    return state, head_only, tail_only, both, body, free


def event_table(time, active, normal_mag, shear_mag, total_mag, ts, sampled_gaps,
                contact_indices, contact_gaps, contact_vn, contact_vt):
    rows = []
    for number, (first, last) in enumerate(intervals(active), start=1):
        geometry_rows = np.where((contact_indices >= first) & (contact_indices <= last))[0]
        event_gaps = contact_gaps[geometry_rows]
        minimum_by_region = event_gaps[:, 1:].min(axis=0)
        present = minimum_by_region <= 1.0
        if present[0] and present[1]:
            label = "BOTH" + ("+BODY" if present[2] else "")
        elif present[0]:
            label = "HEAD" + ("+BODY" if present[2] else "")
        elif present[1]:
            label = "TAIL" + ("+BODY" if present[2] else "")
        else:
            label = "BODY"
        closest_column = np.argmin(event_gaps[:, 1:], axis=1)
        event_vn = contact_vn[geometry_rows, closest_column]
        event_vt = contact_vt[geometry_rows, closest_column]
        segment = slice(first, last + 1)
        time_segment = time[segment]
        peak_index = first + int(np.argmax(total_mag[segment]))
        rows.append({
            "event": number,
            "start_ms": time[first] * 1.0e3,
            "end_ms": time[last] * 1.0e3,
            "peak_time_ms": time[peak_index] * 1.0e3,
            "peak_history_index": peak_index,
            "duration_ms": duration(time, first, last) * 1.0e3,
            "state": label,
            "peak_normal_force_N": float(normal_mag[segment].max()),
            "peak_tangential_force_N": float(shear_mag[segment].max()),
            "peak_resultant_force_N": float(total_mag[segment].max()),
            "normal_impulse_Ns": trapz(normal_mag[segment], time_segment),
            "tangential_impulse_Ns": trapz(shear_mag[segment], time_segment),
            "HEAD_min_gap_um": float(minimum_by_region[0]),
            "TAIL_min_gap_um": float(minimum_by_region[1]),
            "BODY_min_gap_um": float(minimum_by_region[2]),
            "max_closing_normal_velocity_mm_s": float(event_vn.max()),
            "max_tangential_velocity_mm_s": float(event_vt.max()),
            "separation_before_ms": (time[first] - time[rows[-1]["_last_index"]]) * 1.0e3 if rows else time[first] * 1.0e3,
            "_first_index": first,
            "_last_index": last,
        })
    dt = float(time[1] - time[0])
    for index, row in enumerate(rows):
        next_first = rows[index + 1]["_first_index"] if index + 1 < len(rows) else len(time) - 1
        reopen_start = int(np.searchsorted(ts, time[row["_last_index"]], side="left"))
        reopen_end = int(np.searchsorted(ts, time[next_first], side="right"))
        reopen_start = min(reopen_start, len(ts) - 1)
        reopen_end = max(reopen_start + 1, min(reopen_end, len(ts)))
        reopened_gap = float(sampled_gaps[reopen_start:reopen_end, 0].max())
        row["separation_after_ms"] = max(0.0, time[next_first] - time[row["_last_index"]] - dt) * 1.0e3
        row["maximum_reopened_gap_um"] = reopened_gap
        row["gap_reopened_after"] = reopened_gap > 0.05
        row.pop("_first_index")
        row.pop("_last_index")
    return pd.DataFrame(rows)


def classify(theta, theta_cross, linearity, gaps, both_fraction, active, events):
    clean = (
        theta.min() < -8.0 and theta.max() > 8.0
        and np.max(np.abs(theta_cross)) < 1.0
        and linearity > 0.99
        and gaps[:, 0].min() > -20.0
        and both_fraction < 0.25
    )
    event_states = events["state"].tolist() if len(events) else []
    alternating = any("HEAD" in event_states[first] and "TAIL" in event_states[second]
                      or "TAIL" in event_states[first] and "HEAD" in event_states[second]
                      for first in range(len(event_states)) for second in range(first + 1, len(event_states)))
    has_separation = len(intervals(active)) >= 2
    if clean and alternating and has_separation:
        return "HEADTAIL_WALL_ROCKING_RECOVERED", bool(clean), bool(alternating), bool(has_separation)
    if clean and not active.any():
        return "ROCKING_PRESERVED_BUT_LIGHT_CONTACT_THRESHOLD_NOT_REACHED", bool(clean), bool(alternating), bool(has_separation)
    return "CONTACT_DESTROYS_RECOVERED_ROCKING", bool(clean), bool(alternating), bool(has_separation)


def add_labels(canvas, labels, widths, band=32):
    output = Image.new("RGB", (canvas.width, canvas.height + band), "white")
    output.paste(canvas, (0, band))
    draw = ImageDraw.Draw(output)
    offset = 0
    for label, width in zip(labels, widths):
        box = draw.textbbox((0, 0), label, font=ImageFont.load_default())
        draw.text((offset + (width - (box[2] - box[0])) / 2, 10), label, fill="#222222")
        offset += width
    return output


def paired_gif(left_path, right_path, destination, left_label, right_label, left_indices=None):
    left = Image.open(left_path)
    right = Image.open(right_path)
    frames = []
    for frame in range(101):
        left.seek(left_indices[frame] if left_indices is not None else frame)
        right.seek(frame)
        first, second = left.convert("RGB"), right.convert("RGB")
        height = 540
        first = first.resize((round(first.width * height / first.height), height), Image.Resampling.LANCZOS)
        second = second.resize((round(second.width * height / second.height), height), Image.Resampling.LANCZOS)
        canvas = Image.new("RGB", (first.width + second.width, height), "white")
        canvas.paste(first, (0, 0))
        canvas.paste(second, (first.width, 0))
        frames.append(add_labels(canvas, (left_label, right_label), (first.width, second.width)))
    frames[0].save(destination, save_all=True, append_images=frames[1:], duration=60, loop=0, optimize=False)


def dataframe_markdown(frame):
    columns = list(frame.columns)
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for _, row in frame.iterrows():
        values = []
        for column in columns:
            value = row[column]
            values.append(f"{value:.6g}" if isinstance(value, (float, np.floating)) else str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def render_dual_view(pose, nodes, region, c_axis, n_axis, b_axis, radius):
    ids = np.linspace(0, len(pose) - 1, 101).astype(int)
    relative = nodes - RP0
    phi, axial = np.linspace(0, 2 * np.pi, 25), np.linspace(-1.8, 1.8, 18)
    pp, ss = np.meshgrid(phi, axial)
    tube = RP0 + ss[..., None] * c_axis + radius * np.cos(pp)[..., None] * n_axis + radius * np.sin(pp)[..., None] * b_axis
    colors = np.where(region == "HEAD", "#c43c35", np.where(region == "TAIL", "#315a8a", "#444444"))
    fig = plt.figure(figsize=(12, 5.4))
    left = fig.add_subplot(121, projection="3d")
    right = fig.add_subplot(122)
    fig.subplots_adjust(top=0.78, wspace=0.18)

    def draw(frame):
        left.cla()
        right.cla()
        row = pose.iloc[ids[frame]]
        rotation = Rotation.from_rotvec(row[["ur1", "ur2", "ur3"]].to_numpy(float))
        center = row[["com_x", "com_y", "com_z"]].to_numpy(float)
        points = center + rotation.apply(relative)
        left.plot_surface(tube[..., 0], tube[..., 1], tube[..., 2], color="#8bb9c7", alpha=0.12, linewidth=0)
        left.scatter(points[::3, 0], points[::3, 1], points[::3, 2], c=colors[::3], s=3)
        left.plot([RP0[0] - 1.8 * c_axis[0], RP0[0] + 1.8 * c_axis[0]],
                  [RP0[1] - 1.8 * c_axis[1], RP0[1] + 1.8 * c_axis[1]],
                  [RP0[2] - 1.8 * c_axis[2], RP0[2] + 1.8 * c_axis[2]], color="#237a57", lw=1)
        left.set(xlim=(RP0[0] - 2, RP0[0] + 2), ylim=(RP0[1] - 2, RP0[1] + 2),
                 zlim=(RP0[2] - 2, RP0[2] + 2), title="Fixed global camera")
        left.view_init(24, -58)
        local = np.column_stack(((points - RP0).dot(c_axis), (points - RP0).dot(n_axis)))
        right.scatter(local[::3, 0], local[::3, 1], c=colors[::3], s=5)
        for name, color in (("HEAD", "#c43c35"), ("TAIL", "#315a8a")):
            center2 = local[region == name].mean(axis=0)
            right.text(center2[0], center2[1], name, color=color, weight="bold", fontsize=8)
        right.axhline(radius, color="#8bb9c7")
        right.axhline(-radius, color="#8bb9c7")
        right.axhline(0, color="#237a57", lw=0.8)
        right.set_aspect("equal")
        right.set(xlim=(-1.8, 1.8), ylim=(-0.8, 0.8), xlabel="c_hat (mm)", ylabel="n_rock (mm)",
                  title="Rocking plane; view along b_rock")
        fig.suptitle(
            f"t={row.time_s*1000:6.1f} ms  alpha_B={row.alpha_B_deg:+6.2f} deg  "
            f"theta_rock={row.theta_rock_deg:+6.2f} deg  theta_cross={row.theta_cross_deg:+5.2f} deg\n"
            f"HEAD gap={row.HEAD_gap_um:+7.2f} um  TAIL gap={row.TAIL_gap_um:+7.2f} um  "
            f"state={row.contact_state}  T_rock={row.T_rock_Nmm*1e3:+7.3f} uN mm"
        )

    animation = FuncAnimation(fig, draw, frames=len(ids), interval=60)
    destination = OUT / f"{JOB}_DualView.gif"
    animation.save(destination, writer=PillowWriter(fps=16.667), dpi=90)
    plt.close(fig)
    return destination


def write_report(metrics, events):
    selection = json.loads((OUT / "Selected_LightContact_Amplitude.json").read_text())
    event_text = dataframe_markdown(events) if len(events) else "No force-resolved contact event was detected."
    if metrics["classification"] == "HEADTAIL_WALL_ROCKING_RECOVERED":
        recommendation = "Transfer this identical rocking plane, amplitude, 5 Hz/10 mT actuation, contact, and Reduced-Hydro model unchanged into the curved tube for one turn-section validation."
    elif metrics["classification"] == "ROCKING_PRESERVED_BUT_LIGHT_CONTACT_THRESHOLD_NOT_REACHED":
        recommendation = "Before any later solve, reconcile the static touch prediction with the measured dynamic COM shift and derive one corrected amplitude; do not change contact or geometry."
    else:
        recommendation = "Before any later solve, revise the zero-solve contact target model to include the measured dynamic COM shift and the observed BODY/TAIL event sequence; do not change contact, hydro, or geometry."
    taxonomy_note = ""
    if metrics["classification"] == "CONTACT_DESTROYS_RECOVERED_ROCKING" and metrics["line_like_rocking_preserved"]:
        taxonomy_note = (
            "The required three-label taxonomy has no label for the observed fourth outcome: line-like rocking is preserved, "
            "but contact occurs at BODY/TAIL rather than alternating HEAD/TAIL. The failure label therefore denotes destruction "
            "of the requested HEAD-TAIL contact topology, not destruction of planar rocking."
        )
    report = f"""# Local Head-Tail Light-Contact Validation

## Decision

**{metrics['classification']}**

{taxonomy_note}

1. The completed ROCK10 reference genuinely recovered RouteA-style planar rocking: robot range `{metrics['ROCK10_min_deg']:.6f}..{metrics['ROCK10_max_deg']:.6f} deg`, fitted amplitude `{metrics['ROCK10_amplitude_deg']:.6f} deg`, and PCA line-likeness `{metrics['ROCK10_pca_linearity']:.12f}`.
2. The light-contact robot reached `{metrics['positive_peak_deg']:+.6f} deg` and `{metrics['negative_peak_deg']:+.6f} deg`; these are robot angles, not field commands.
3. Light-contact `theta_cross` RMS/max-absolute is `{metrics['theta_cross_rms_deg']:.6g}/{metrics['theta_cross_max_abs_deg']:.6g} deg`; PCA line-likeness is `{metrics['pca_rocking_linearity']:.12f}`.
4. Exact positive first touch is `{selection['theta_touch_pos_deg']:+.9f} deg`.
5. Exact negative first touch is `{selection['theta_touch_neg_deg']:+.9f} deg`.
6. Static exact geometry says `{selection['first_touch_pos_region']}` controls positive first touch and `{selection['first_touch_neg_region']}` controls negative first touch. This asymmetry was retained rather than relabelled to fit the desired motion.
7. The ROCK10 no-drift dynamic gain used for design was `{selection['G_theta']:.12f}`, with phase lag `{selection['phase_lag_robot_vs_field_deg']:.9f} deg`.
8. The target was the more conservative of a 10 um free-pose overtravel and a 0.5 deg cap on each side, then the minimum symmetric field amplitude satisfying both biased, gain-corrected peak inequalities.
9. The selected field amplitude was `{selection['selected_A_B_required_deg']:.12f} deg`; no other dynamic amplitude was run.
10. Force-resolved HEAD-only/TAIL-only fractions are `{metrics['HEAD_only_contact_fraction']:.6f}/{metrics['TAIL_only_contact_fraction']:.6f}`. Alternating HEAD/TAIL events: `{metrics['alternating_head_tail_events']}`.
11. There are `{metrics['contact_event_count']}` force-resolved events and `{metrics['separation_interval_count']}` separate intervals; longest separated interval is `{metrics['longest_separated_ms']:.6f} ms`. Every event measurably reopens by more than 0.05 um before the next event/end: `{metrics['all_contact_events_reopen']}`.
12. Both-end fraction is `{metrics['both_end_contact_fraction']:.6f}` and longest bridge is `{metrics['longest_both_end_bridge_ms']:.6f} ms`.
13. Line-like rocking preserved: `{metrics['line_like_rocking_preserved']}`. Minimum exact gap is `{metrics['minimum_exact_gap_um']:.6f} um`; deep penetration is `{metrics['deep_penetration']}`.
14. Final classification: **{metrics['classification']}**.
15. Exactly one next recommendation: **{recommendation}**

## Contact Event Audit

`CFN` and `CFS` are the Abaqus whole-robot surface history components. Event duration, impulses, and HEAD/TAIL/BODY attribution use native `1e-7 s` contact increments; attribution recomputes the exact faceted-wall gap for every active increment. The full-cycle pose/GIF grid is sampled at 0.1 ms. Relative normal and tangential velocities are rigid-surface kinematic estimates at each region's closest solver node because Abaqus did not write a direct relative-contact-velocity history.

{event_text}

The first force-resolved event occurs at `{metrics['first_contact_peak_time_ms']:.6f} ms`: `theta_rock={metrics['first_contact_theta_rock_deg']:+.6f} deg`, `alpha_B={metrics['first_contact_alpha_B_deg']:+.6f} deg`, `T_rock={metrics['first_contact_T_rock_Nmm']:.6g} N mm`, normal/tangential impulse `{metrics['first_contact_normal_impulse_Ns']:.6g}/{metrics['first_contact_tangential_impulse_Ns']:.6g} N s`, closest-point closing/tangential speed `{metrics['first_contact_closing_velocity_mm_s']:.6g}/{metrics['first_contact_tangential_velocity_mm_s']:.6g} mm/s`, and COM radial shift `{metrics['first_contact_COM_radial_um']:.6f} um`. Magnetic torque pushes toward larger absolute rocking angle: `{metrics['first_contact_magnetic_torque_pushes_deeper']}`. The event is `{metrics['first_contact_state']}`, then opens to `{metrics['first_contact_reopened_gap_um']:.6f} um` before the next event.

## Solver And Identity Gate

- Case status: `{metrics['case_status']}`; Abaqus normal completion: `{metrics['abaqus_completed_normally']}`; fatal-error scan: `{metrics['solver_fatal_error_detected']}`.
- Final RP/contact/energy history times: `{metrics['rp_final_time_s']:.9f}`, `{metrics['contact_final_time_s']:.9f}`, `{metrics['energy_final_time_s']:.9f} s`.
- Explicit increment: `{metrics['explicit_dt_s']:.9g} s`; ODB frames: `{metrics['odb_frames']}`; socket stderr bytes: `{metrics['socket_stderr_bytes']}`.
- Live/vendored magnetic-server hashes match the frozen identity: `{metrics['server_identity_match']}`.
- Field magnitude maximum absolute error from 10 mT: `{metrics['field_magnitude_max_abs_error_T']:.6g} T`.
- Final robot rocking angle at 200 ms: `{metrics['final_theta_rock_deg']:+.6f} deg`; zero crossings: `{metrics['theta_rock_zero_crossings']}`.
- COM radial excursion: `{metrics['COM_radial_max_um']:.6f} um`; net axial displacement: `{metrics['axial_displacement_mm']:.9g} mm`.
- Normal/tangential impulses: `{metrics['normal_impulse_Ns']:.6g}/{metrics['tangential_impulse_Ns']:.6g} N s`; peak normal/tangential/resultant force: `{metrics['peak_normal_force_N']:.6g}/{metrics['peak_tangential_force_N']:.6g}/{metrics['peak_resultant_force_N']:.6g} N`.
- Energy ranges (N mm): ALLKE `{metrics['ALLKE_min']:.6g}..{metrics['ALLKE_max']:.6g}`, ALLIE `{metrics['ALLIE_min']:.6g}..{metrics['ALLIE_max']:.6g}`, ETOTAL `{metrics['ETOTAL_min']:.6g}..{metrics['ETOTAL_max']:.6g}`.

The complete fixed-camera GIFs, not scalar metrics alone, are the final topology gate. The legacy RouteA animation is used only for normalized rocking-phase comparison.
"""
    (OUT / "Local_HeadTail_LightContact_Validation.md").write_text(report, encoding="utf-8")


def main():
    identity = json.loads((CASE / "case_identity.json").read_text())
    inventory = json.loads((CASE / "private" / "odb_inventory.json").read_text())
    time, history = dense_pose()
    cfn, cfs, cft = contact_vectors(time)
    normal_mag, shear_mag = np.linalg.norm(cfn, axis=1), np.linalg.norm(cfs, axis=1)
    resultant = cfn + cfs
    resultant_mag = np.linalg.norm(resultant, axis=1)
    active_raw = resultant_mag > CONTACT_THRESHOLD_N
    active = merge_short_inactive_gaps(active_raw, time)

    stride = max(1, int(round(SAMPLE_INTERVAL_S / np.median(np.diff(time)))))
    sample = np.unique(np.r_[np.arange(0, len(time), stride), len(time) - 1])
    ts = time[sample]
    centers = RP0 + history["U"][sample]
    rotations = Rotation.from_rotvec(history["UR"][sample])
    velocities = history["V"][sample]
    angular_velocities = history["VR"][sample]
    global c_axis
    c_axis, n_axis, b_axis = map(unit, (identity["initial_axis_aba"], identity["n_rock_aba"], identity["b_rock_aba"]))
    axes = rotations.apply(np.broadcast_to(c_axis, (len(sample), 3)))
    qc, qrock, qcross = axes.dot(c_axis), axes.dot(n_axis), axes.dot(b_axis)
    theta = np.degrees(np.arctan2(qrock, qc))
    theta_cross = np.degrees(np.arctan2(qcross, qc))
    axial = (centers - RP0).dot(c_axis)
    radial_vectors = centers - RP0 - np.outer(axial, c_axis)
    radial = np.linalg.norm(radial_vectors, axis=1)

    telemetry = pd.read_csv(CASE / f"{JOB}_telemetry.csv")
    field = np.column_stack([np.interp(ts, telemetry.t_s, telemetry[key]) for key in ("Bx_aba_T", "By_aba_T", "Bz_aba_T")])
    torque = np.column_stack([np.interp(ts, telemetry.t_s, telemetry[key]) for key in ("tx_aba_Nmm", "ty_aba_Nmm", "tz_aba_Nmm")])
    alpha = np.degrees(np.arctan2(field.dot(n_axis), field.dot(c_axis)))
    rocking_torque = torque.dot(b_axis)

    deck = (CASE / f"{JOB}.inp").read_text()
    nodes = part_nodes(deck, "Robot_SOLID")
    wall = part_nodes(deck, "Pipe_WALL_HELPER")
    initial_axial = (nodes - RP0).dot(c_axis)
    region = np.full(len(nodes), "BODY", dtype=object)
    region[initial_axial <= initial_axial.min() + 0.25] = "HEAD"
    region[initial_axial >= initial_axial.max() - 0.25] = "TAIL"
    normals, radius = face_geometry(wall, RP0, c_axis, n_axis, b_axis)
    gaps, normal_velocity, tangential_velocity, closest = exact_geometry(
        nodes, region, rotations, centers, velocities, angular_velocities, normals, radius, c_axis
    )
    active_sample, normal_peak, shear_peak, total_peak = sample_contact(
        time, active, normal_mag, shear_mag, resultant_mag, sample
    )
    region_column = np.choose(np.argmin(gaps[:, 1:], axis=1), [0, 1, 2])
    vn_closest = normal_velocity[np.arange(len(ts)), region_column]
    vt_closest = tangential_velocity[np.arange(len(ts)), region_column]
    contact_indices = np.flatnonzero(active)
    contact_gaps, contact_vn, contact_vt, _ = exact_geometry(
        nodes, region,
        Rotation.from_rotvec(history["UR"][contact_indices]),
        RP0 + history["U"][contact_indices],
        history["V"][contact_indices], history["VR"][contact_indices],
        normals, radius, c_axis,
    )
    events = event_table(
        time, active, normal_mag, shear_mag, resultant_mag, ts, gaps,
        contact_indices, contact_gaps, contact_vn, contact_vt,
    )
    peak_indices = events["peak_history_index"].to_numpy(dtype=int)
    peak_times = time[peak_indices]
    peak_rotations = Rotation.from_rotvec(history["UR"][peak_indices])
    peak_axes = peak_rotations.apply(np.broadcast_to(c_axis, (len(events), 3)))
    events["theta_rock_deg"] = np.degrees(np.arctan2(peak_axes.dot(n_axis), peak_axes.dot(c_axis)))
    events["theta_cross_deg"] = np.degrees(np.arctan2(peak_axes.dot(b_axis), peak_axes.dot(c_axis)))
    events["alpha_B_deg"] = np.interp(peak_times, ts, alpha)
    events["T_rock_Nmm"] = np.interp(peak_times, ts, rocking_torque)
    peak_centers = RP0 + history["U"][peak_indices]
    peak_axial = (peak_centers - RP0).dot(c_axis)
    peak_radial = peak_centers - RP0 - np.outer(peak_axial, c_axis)
    events["COM_radial_um"] = np.linalg.norm(peak_radial, axis=1) * 1.0e3
    events["magnetic_torque_pushes_deeper"] = events["theta_rock_deg"] * events["T_rock_Nmm"] > 0.0
    state = np.full(len(ts), "FREE", dtype=object)
    for _, event in events.iterrows():
        event_state = "BOTH" if "HEAD" in event.state and "TAIL" in event.state else (
            "HEAD" if "HEAD" in event.state else "TAIL" if "TAIL" in event.state else "BODY"
        )
        state[int(np.argmin(np.abs(ts - event.peak_time_ms * 1.0e-3)))] = event_state
    total_duration_s = float(time[-1] - time[0] + (time[1] - time[0]))
    event_duration_s = events["duration_ms"].to_numpy(float) * 1.0e-3
    event_states = events["state"].to_numpy(dtype=str)
    head_only_event = np.char.find(event_states, "HEAD") >= 0
    tail_only_event = np.char.find(event_states, "TAIL") >= 0
    both_event = head_only_event & tail_only_event
    head_only_event &= ~both_event
    tail_only_event &= ~both_event
    body_only_event = ~(head_only_event | tail_only_event | both_event)
    head_fraction = float(event_duration_s[head_only_event].sum() / total_duration_s)
    tail_fraction = float(event_duration_s[tail_only_event].sum() / total_duration_s)
    both_fraction = float(event_duration_s[both_event].sum() / total_duration_s)
    body_fraction = float(event_duration_s[body_only_event].sum() / total_duration_s)
    separated_fraction = 1.0 - float(event_duration_s.sum() / total_duration_s)
    metric_gaps = np.vstack((gaps, contact_gaps))

    covariance = np.cov(np.column_stack((qrock, qcross)), rowvar=False)
    eigenvalues = np.linalg.eigvalsh(covariance)
    linearity = float(eigenvalues[-1] / max(eigenvalues.sum(), 1.0e-30))
    valid_sign = np.abs(theta) > 0.1
    signs = np.sign(theta[valid_sign])
    zero_crossings = int(np.sum(signs[1:] != signs[:-1]))
    robot_amplitude, robot_phase, mean_offset, drift = harmonic(ts, theta)
    field_amplitude, field_phase, _, _ = harmonic(ts, alpha)
    phase_lag = ((robot_phase - field_phase + 180.0) % 360.0) - 180.0
    classification, clean, alternating, separated_events = classify(
        theta, theta_cross, linearity, metric_gaps, both_fraction, active, events
    )

    energy = np.load(CASE / "private" / "energy_history_private.npz")
    energy_metrics = {key + "_min": float(energy[key][:, 1].min()) for key in energy.files}
    energy_metrics.update({key + "_max": float(energy[key][:, 1].max()) for key in energy.files})
    rock10 = json.loads((OUT / "rocking_metrics.json").read_text())
    selection = json.loads((OUT / "Selected_LightContact_Amplitude.json").read_text())
    # The interpolated vectors no longer carry time, so use a native contact array here.
    contact_archive = np.load(CASE / "private" / "contact_history_private.npz")
    contact_final = max(float(contact_archive[key][-1, 0]) for key in contact_archive.files)
    energy_final = max(float(energy[key][-1, 0]) for key in energy.files)
    separation_intervals = intervals(~active)
    solver_text = "\n".join(
        (CASE / f"{JOB}{suffix}").read_text(errors="replace")
        for suffix in (".sta", ".msg", ".dat")
        if (CASE / f"{JOB}{suffix}").exists()
    )
    solver_text_upper = solver_text.upper()

    metrics = {
        "case_id": JOB,
        "classification": classification,
        "case_status": identity["status"],
        "abaqus_completed_normally": "ANALYSIS HAS COMPLETED SUCCESSFULLY" in solver_text_upper,
        "solver_fatal_error_detected": any(token in solver_text_upper for token in ("***ERROR", "FATAL ERROR", "ANALYSIS TERMINATED")),
        "rp_final_time_s": float(time[-1]),
        "contact_final_time_s": contact_final,
        "energy_final_time_s": energy_final,
        "explicit_dt_s": float(identity["direct_dt_s"]),
        "history_time_spacing_median_s": float(np.median(np.diff(time))),
        "odb_frames": int(inventory["frames"]),
        "socket_stderr_bytes": (CASE / f"{JOB}_socket_stderr.log").stat().st_size,
        "live_server_sha256": sha256(LIVE_SERVER),
        "vendored_server_sha256": sha256(VENDORED_SERVER),
        "server_identity_match": sha256(LIVE_SERVER) == sha256(VENDORED_SERVER) == identity["production_server_sha256"],
        "theta_rock_min_deg": float(theta.min()),
        "theta_rock_max_deg": float(theta.max()),
        "positive_peak_deg": float(theta.max()),
        "negative_peak_deg": float(theta.min()),
        "theta_rock_mean_offset_deg": mean_offset,
        "theta_rock_sample_mean_deg": float(theta.mean()),
        "theta_rock_linear_drift_deg_per_s": drift,
        "theta_rock_fundamental_amplitude_deg": robot_amplitude,
        "field_fundamental_amplitude_deg": field_amplitude,
        "field_to_robot_amplitude_ratio": robot_amplitude / field_amplitude,
        "phase_lag_robot_vs_field_deg": phase_lag,
        "final_theta_rock_deg": float(theta[-1]),
        "theta_rock_zero_crossings": zero_crossings,
        "theta_cross_rms_deg": float(np.sqrt(np.mean(theta_cross ** 2))),
        "theta_cross_max_abs_deg": float(np.max(np.abs(theta_cross))),
        "pca_rocking_linearity": linearity,
        "line_like_rocking_preserved": clean,
        "COM_radial_max_um": float(radial.max() * 1.0e3),
        "axial_displacement_mm": float(axial[-1] - axial[0]),
        "HEAD_min_gap_um": float(metric_gaps[:, 1].min()),
        "TAIL_min_gap_um": float(metric_gaps[:, 2].min()),
        "BODY_min_gap_um": float(metric_gaps[:, 3].min()),
        "minimum_exact_gap_um": float(metric_gaps[:, 0].min()),
        "deep_penetration": bool(metric_gaps[:, 0].min() < -20.0),
        "HEAD_only_contact_fraction": head_fraction,
        "TAIL_only_contact_fraction": tail_fraction,
        "BODY_only_contact_fraction": body_fraction,
        "both_end_contact_fraction": both_fraction,
        "separated_fraction": separated_fraction,
        "longest_HEAD_only_ms": float(events.loc[head_only_event, "duration_ms"].max()) if head_only_event.any() else 0.0,
        "longest_TAIL_only_ms": float(events.loc[tail_only_event, "duration_ms"].max()) if tail_only_event.any() else 0.0,
        "longest_both_end_bridge_ms": float(events.loc[both_event, "duration_ms"].max()) if both_event.any() else 0.0,
        "longest_separated_ms": longest(time, ~active) * 1.0e3,
        "contact_event_count": int(len(events)),
        "separation_interval_count": int(len(separation_intervals)),
        "all_contact_events_reopen": bool(events["gap_reopened_after"].all()) if len(events) else False,
        "alternating_head_tail_events": alternating,
        "separated_contact_events": separated_events,
        "first_contact_peak_time_ms": float(events.iloc[0].peak_time_ms),
        "first_contact_state": str(events.iloc[0].state),
        "first_contact_theta_rock_deg": float(events.iloc[0].theta_rock_deg),
        "first_contact_alpha_B_deg": float(events.iloc[0].alpha_B_deg),
        "first_contact_T_rock_Nmm": float(events.iloc[0].T_rock_Nmm),
        "first_contact_normal_impulse_Ns": float(events.iloc[0].normal_impulse_Ns),
        "first_contact_tangential_impulse_Ns": float(events.iloc[0].tangential_impulse_Ns),
        "first_contact_closing_velocity_mm_s": float(events.iloc[0].max_closing_normal_velocity_mm_s),
        "first_contact_tangential_velocity_mm_s": float(events.iloc[0].max_tangential_velocity_mm_s),
        "first_contact_COM_radial_um": float(events.iloc[0].COM_radial_um),
        "first_contact_magnetic_torque_pushes_deeper": bool(events.iloc[0].magnetic_torque_pushes_deeper),
        "first_contact_reopened_gap_um": float(events.iloc[0].maximum_reopened_gap_um),
        "peak_normal_force_N": float(normal_mag.max()),
        "peak_tangential_force_N": float(shear_mag.max()),
        "peak_resultant_force_N": float(resultant_mag.max()),
        "normal_impulse_Ns": trapz(normal_mag, time),
        "tangential_impulse_Ns": trapz(shear_mag, time),
        "resultant_impulse_Ns": trapz(resultant_mag, time),
        "max_closing_normal_velocity_mm_s": float(events["max_closing_normal_velocity_mm_s"].max()) if len(events) else 0.0,
        "max_tangential_velocity_mm_s": float(events["max_tangential_velocity_mm_s"].max()) if len(events) else 0.0,
        "Trock_min_Nmm": float(rocking_torque.min()),
        "Trock_max_Nmm": float(rocking_torque.max()),
        "field_magnitude_max_abs_error_T": float(np.max(np.abs(np.linalg.norm(field, axis=1) - 0.01))),
        "CFT_vs_CFN_plus_CFS_max_vector_error_N": float(np.max(np.linalg.norm(cft - resultant, axis=1))),
        "ROCK10_min_deg": rock10["theta_rock_min_deg"],
        "ROCK10_max_deg": rock10["theta_rock_max_deg"],
        "ROCK10_amplitude_deg": selection["A_robot_current_deg"],
        "ROCK10_pca_linearity": rock10["pca_rocking_linearity"],
        **energy_metrics,
    }

    pose = pd.DataFrame({
        "time_s": ts,
        "normalized_phase": ts / 0.2,
        "theta_rock_deg": theta,
        "theta_cross_deg": theta_cross,
        "alpha_B_deg": alpha,
        "q_rock": qrock,
        "q_cross": qcross,
        "HEAD_gap_um": gaps[:, 1],
        "TAIL_gap_um": gaps[:, 2],
        "BODY_gap_um": gaps[:, 3],
        "minimum_gap_um": gaps[:, 0],
        "closest_region": closest,
        "T_rock_Nmm": rocking_torque,
        "contact_state": state,
        "peak_normal_force_in_bin_N": normal_peak,
        "peak_tangential_force_in_bin_N": shear_peak,
        "peak_resultant_force_in_bin_N": total_peak,
        "CFN1_N": cfn[sample, 0], "CFN2_N": cfn[sample, 1], "CFN3_N": cfn[sample, 2],
        "CFS1_N": cfs[sample, 0], "CFS2_N": cfs[sample, 1], "CFS3_N": cfs[sample, 2],
        "CFT1_N": cft[sample, 0], "CFT2_N": cft[sample, 1], "CFT3_N": cft[sample, 2],
        "closest_normal_closing_velocity_mm_s": vn_closest,
        "closest_tangential_velocity_mm_s": vt_closest,
        "com_x": centers[:, 0], "com_y": centers[:, 1], "com_z": centers[:, 2],
        "axis_x": axes[:, 0], "axis_y": axes[:, 1], "axis_z": axes[:, 2],
        "ur1": history["UR"][sample, 0], "ur2": history["UR"][sample, 1], "ur3": history["UR"][sample, 2],
        "s_mm": axial, "COM_radial_mm": radial,
    })
    pose.to_csv(OUT / "light_contact_pose.csv", index=False)
    events.to_csv(OUT / "light_contact_events.csv", index=False)
    (OUT / "light_contact_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    pd.DataFrame([metrics]).to_csv(OUT / "light_contact_metrics.csv", index=False)

    fig, panels = plt.subplots(6, 1, figsize=(9, 12), sharex=True)
    panels[0].plot(ts * 1.0e3, theta, label="robot theta_rock", color="#c43c35")
    panels[0].plot(ts * 1.0e3, alpha, label="field alpha_B", color="#237a57", ls="--")
    panels[0].set_ylabel("deg"); panels[0].legend(ncol=2, frameon=False)
    panels[1].plot(ts * 1.0e3, theta_cross, color="#315a8a"); panels[1].set_ylabel("theta_cross\n(deg)")
    panels[2].plot(ts * 1.0e3, gaps[:, 1], label="HEAD", color="#c43c35")
    panels[2].plot(ts * 1.0e3, gaps[:, 2], label="TAIL", color="#315a8a")
    panels[2].plot(ts * 1.0e3, gaps[:, 3], label="BODY", color="#777777", alpha=0.7)
    panels[2].axhline(0, color="#222222", lw=0.7); panels[2].set_ylabel("exact gap\n(um)"); panels[2].legend(ncol=3, frameon=False)
    codes = pd.Categorical(state, categories=["FREE", "BODY", "HEAD", "TAIL", "BOTH"]).codes
    panels[3].step(ts * 1.0e3, codes, where="post", color="#555555"); panels[3].set_ylabel("contact")
    panels[3].set_yticks(range(5), ["FREE", "BODY", "HEAD", "TAIL", "BOTH"])
    panels[4].plot(ts * 1.0e3, normal_peak, color="#c43c35", label="normal")
    panels[4].plot(ts * 1.0e3, shear_peak, color="#315a8a", label="tangential")
    panels[4].set_ylabel("bin peak force\n(N)"); panels[4].legend(frameon=False)
    panels[5].plot(ts * 1.0e3, rocking_torque * 1.0e3, color="#237a57")
    panels[5].set_ylabel("T_rock\n(uN mm)"); panels[5].set_xlabel("Time (ms)")
    for marker in (0, 50, 100, 150, 200):
        for panel in panels:
            panel.axvline(marker, color="#777777", lw=0.55, alpha=0.5)
    for panel in panels:
        panel.grid(alpha=0.2)
    fig.tight_layout(); fig.savefig(OUT / "LightContact_5Hz_TimeSeries.png", dpi=220); fig.savefig(OUT / "LightContact_5Hz_TimeSeries.pdf"); plt.close(fig)

    fig, panel = plt.subplots(figsize=(6.4, 5))
    panel.plot(qrock, qcross, color="#c43c35", lw=1.6, label="Light-contact response")
    clean_pose = pd.read_csv(OUT / "rocking_pose.csv")
    panel.plot(clean_pose.q_rock, clean_pose.q_cross, color="#777777", lw=2.4, alpha=0.75, label="ROCK10 response")
    panel.set(xlabel="q_rock", ylabel="q_cross", title=f"Rocking orbit; light-contact PCA linearity={linearity:.6f}")
    panel.grid(alpha=0.2); panel.legend(frameon=False); fig.tight_layout()
    fig.savefig(OUT / "LightContact_Rocking_Orbit.png", dpi=220); fig.savefig(OUT / "LightContact_Rocking_Orbit.pdf"); plt.close(fig)

    dual = render_dual_view(pose, nodes, region, c_axis, n_axis, b_axis, radius)
    paired_gif(OUT / f"{ROCK10_JOB}_DualView.gif", dual, OUT / "ROCK10_vs_LIGHTCONTACT.gif", "ROCK10 clean rocking", "LIGHTCONTACT")
    legacy_indices = [24 + int(round(6 * frame / 100)) for frame in range(101)]
    paired_gif(LEGACY, dual, OUT / "RouteA_vs_LightContactMagneticRocking.gif", "Legacy RouteA normalized phase", "Magnetic LIGHTCONTACT", legacy_indices)

    review_ids = [int(np.argmin(np.abs(ts - target))) for target in (0.0, 0.05, 0.1, 0.15, 0.2)]
    dual_image = Image.open(dual)
    review_frames = []
    for review_id in review_ids:
        gif_id = int(round(review_id * 100 / (len(pose) - 1)))
        dual_image.seek(gif_id)
        review_frames.append(dual_image.convert("RGB"))
    review_width = 420
    resized = [frame.resize((review_width, round(frame.height * review_width / frame.width)), Image.Resampling.LANCZOS) for frame in review_frames]
    sheet = Image.new("RGB", (review_width, sum(frame.height for frame in resized)), "white")
    offset = 0
    for frame in resized:
        sheet.paste(frame, (0, offset)); offset += frame.height
    sheet.save(OUT / "LightContact_VisualGate_0_50_100_150_200ms.png")

    full_cycle = []
    for gif_id in np.linspace(0, 100, 21).astype(int):
        dual_image.seek(int(gif_id))
        frame = dual_image.convert("RGB")
        width = 300
        full_cycle.append(frame.resize((width, round(frame.height * width / frame.width)), Image.Resampling.LANCZOS))
    columns = 3
    rows = math.ceil(len(full_cycle) / columns)
    review = Image.new("RGB", (columns * full_cycle[0].width, rows * full_cycle[0].height), "white")
    for index, frame in enumerate(full_cycle):
        review.paste(frame, ((index % columns) * frame.width, (index // columns) * frame.height))
    review.save(OUT / "LightContact_FullCycle_21Frame_Review.png")

    write_report(metrics, events)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
