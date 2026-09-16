"""Extract screening metrics and render the four fixed-camera low-G GIFs."""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
CASES = (
    ("G0", "S4_HEADFORWARD_G0", "S4_HEADFORWARD_G0_2D.gif", 0.0),
    ("G0P25", "S4_HEADFORWARD_G0P25", "S4_HEADFORWARD_G0P25_2D.gif", 0.25),
    ("G0P5", "S4_HEADFORWARD_G0P5", "S4_HEADFORWARD_G0P5_2D.gif", 0.50),
    ("G1P0", "S4_HEADFORWARD_G1P0", "S4_HEADFORWARD_G1P0_2D.gif", 1.00),
)
OLD_REF_JOB = "FAST_ELLIPTIC_20HZ_A14P343_X2P5_FORWARD"
NEAR_ZERO_MM_S = 1.0
CONTACT_THRESHOLD_N = 1.0e-8
CONTACT_MERGE_US = 20.0
GEOMETRY_STEP_S = 20.0e-6
GIF_STEP_S = 1.0e-3
GIF_FPS = 20.0
XLIM = (-3.0, 17.0)
YLIM = (-2.5, 2.5)


def unit(value):
    value = np.asarray(value, dtype=float)
    return value / max(np.linalg.norm(value), 1.0e-30)


def get_array(archive, name):
    keys = [key for key in archive.files if key == name or key.startswith(name + " (Repeated:")]
    if not keys:
        raise KeyError(name)
    return max((archive[key] for key in keys), key=len)


def intervals(flag):
    edge = np.diff(np.r_[False, flag, False].astype(int))
    return list(zip(np.where(edge == 1)[0], np.where(edge == -1)[0] - 1))


def longest_duration(time, flag):
    spans = intervals(flag)
    if not spans:
        return 0.0
    dt = float(np.median(np.diff(time)))
    return max(float(time[j] - time[i] + dt) for i, j in spans)


def part_nodes(deck, name):
    part = re.search(rf"(?ms)^\*Part, name={re.escape(name)}\s*$.*?^\*End Part\s*$", deck).group(0)
    block = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", part).group(1)
    return np.asarray([[float(value) for value in row.split(",")[1:4]]
                       for row in block.splitlines() if row.strip()])


def wall_planes(wall, center, c, n, b):
    axial = (wall - center).dot(c)
    ring = wall[np.abs(axial - axial.min()) < 1.0e-7]
    xy = np.column_stack(((ring - center).dot(n), (ring - center).dot(b)))
    xy = xy[np.argsort(np.arctan2(xy[:, 1], xy[:, 0]))]
    normals, radii = [], []
    for p, q in zip(xy, np.roll(xy, -1, axis=0)):
        edge = q - p
        normal2 = unit([edge[1], -edge[0]])
        if np.dot(normal2, p + q) < 0.0:
            normal2 *= -1.0
        normals.append(normal2[0] * n + normal2[1] * b)
        radii.append(np.dot(p, normal2))
    return np.asarray(normals), float(np.mean(radii))


def contact_force(time, case):
    archive = np.load(case / "private" / "contact_history_private.npz")
    columns = []
    for prefix in ("CFN", "CFS"):
        vector = []
        for axis in (1, 2, 3):
            keys = [key for key in archive.files
                    if "|" + prefix + str(axis) + " on surface " in key and "ASSEMBLY_ROBOT" in key]
            if not keys:
                vector.append(np.zeros_like(time))
            else:
                array = max((archive[key] for key in keys), key=len)
                vector.append(np.interp(time, array[:, 0], array[:, 1]))
        columns.append(np.column_stack(vector))
    archive.close()
    return columns[0] + columns[1]


def load_rp(case):
    archive = np.load(case / "private" / "rp_history_private.npz")
    time = get_array(archive, "V1")[:, 0]
    data = {prefix: np.column_stack([
        np.interp(time, get_array(archive, prefix + str(axis))[:, 0],
                  get_array(archive, prefix + str(axis))[:, 1])
        for axis in (1, 2, 3)]) for prefix in ("U", "UR", "V")}
    archive.close()
    return time, data


def gap_snapshot(index, time, data, nodes, region, center, c, n, b, radius, rp0):
    points = rp0 + data["U"][index] + Rotation.from_rotvec(data["UR"][index]).apply(nodes - rp0)
    relative = points - center
    transverse = relative - np.outer(relative.dot(c), c)
    node_gap = radius - np.max(transverse.dot(gap_snapshot.normals.T), axis=1)
    result = {"all": float(node_gap.min())}
    for name in ("HEAD", "TAIL", "BODY"):
        result[name] = float(node_gap[region == name].min())
    result["points"] = points
    return result


def reversal_count(velocity):
    signs = np.sign(velocity[np.abs(velocity) > NEAR_ZERO_MM_S])
    return int(np.sum(signs[1:] != signs[:-1])) if len(signs) > 1 else 0


def fit_effective_acceleration(time, s):
    design = np.column_stack((np.ones_like(time), time, 0.5 * time * time))
    coefficients, _, _, _ = np.linalg.lstsq(design, s, rcond=None)
    return float(coefficients[2])


def geometry_samples(time, data, nodes, region, center, c, n, b, radius, rp0):
    step = max(1, int(round(GEOMETRY_STEP_S / np.median(np.diff(time)))))
    indices = np.unique(np.r_[np.arange(0, len(time), step), len(time) - 1])
    rows = []
    for index in indices:
        snap = gap_snapshot(int(index), time, data, nodes, region, center, c, n, b, radius, rp0)
        rows.append([time[index], index, snap["all"], snap["HEAD"], snap["TAIL"], snap["BODY"]])
    return indices, pd.DataFrame(rows, columns=["time_s", "index", "gap_mm", "HEAD_gap_mm", "TAIL_gap_mm", "BODY_gap_mm"])


def contact_episode_table(time, force, geom, geom_indices, region):
    active = np.linalg.norm(force, axis=1) > CONTACT_THRESHOLD_N
    raw_events = intervals(active)
    episodes = []
    merge_steps = max(1, int(round(CONTACT_MERGE_US * 1.0e-6 / np.median(np.diff(time)))))
    for number, (start, end) in enumerate(raw_events, 1):
        selected = geom[(geom.index >= 0) & (geom["index"] >= start) & (geom["index"] <= end)]
        if selected.empty:
            selected = geom.iloc[[int(np.argmin(np.abs(geom["index"].to_numpy() - (start + end) / 2.0)))] ]
        mins = {name: float(selected[f"{name}_gap_mm"].min()) for name in ("HEAD", "TAIL", "BODY")}
        close = [name for name, value in mins.items() if value <= min(mins.values()) + 0.02]
        region_label = "+".join(close)
        reopen_end = min(len(active), end + merge_steps + 1)
        reopened = not active[end + 1:reopen_end].any() if end + 1 < len(active) else True
        episodes.append({
            "episode": number,
            "start_s": float(time[start]), "end_s": float(time[end]),
            "duration_ms": float((time[end] - time[start] + np.median(np.diff(time))) * 1.0e3),
            "region": region_label,
            "minimum_gap_mm": float(selected["gap_mm"].min()),
            "HEAD_min_gap_mm": mins["HEAD"], "TAIL_min_gap_mm": mins["TAIL"], "BODY_min_gap_mm": mins["BODY"],
            "peak_contact_force_N": float(np.linalg.norm(force[start:end + 1], axis=1).max()),
            "reopen_after_episode": bool(reopened),
        })
    return active, pd.DataFrame(episodes)


def render_case(case, identity, time, data, nodes, region, c, n, b, rp0, radius, rocking, axial, label, gif_name):
    frame_times = np.arange(0.0, identity["duration_s"] + 0.5 * GIF_STEP_S, GIF_STEP_S)
    frame_ids = np.clip(np.searchsorted(time, frame_times), 0, len(time) - 1)
    colors = np.where(region == "HEAD", "#d1493f", np.where(region == "TAIL", "#2864a8", "#383c42"))
    fig, ax = plt.subplots(figsize=(10.0, 2.5))

    def draw(frame):
        ax.clear()
        index = int(frame_ids[frame])
        points = rp0 + data["U"][index] + Rotation.from_rotvec(data["UR"][index]).apply(nodes - rp0)
        local = np.column_stack(((points - rp0).dot(c), (points - rp0).dot(n)))
        ax.axhspan(radius, YLIM[1], color="#d8edf2", alpha=0.55)
        ax.axhspan(YLIM[0], -radius, color="#d8edf2", alpha=0.55)
        ax.axhline(radius, color="#5996a5", lw=1.1)
        ax.axhline(-radius, color="#5996a5", lw=1.1)
        ax.scatter(local[::2, 0], local[::2, 1], c=colors[::2], s=8, linewidths=0)
        ax.set_xlim(XLIM); ax.set_ylim(YLIM); ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("canonical +s (mm), fixed orthographic camera")
        ax.set_ylabel("n_routeA (mm)")
        ax.grid(axis="x", alpha=0.16)
        ax.text(0.025, 0.95, "HEAD", color="#d1493f", weight="bold", transform=ax.transAxes, va="top")
        ax.text(0.095, 0.95, "TAIL", color="#2864a8", weight="bold", transform=ax.transAxes, va="top")
        ax.set_title(
            f"{label}  |  t={time[index]*1e3:5.1f} ms  |  f=20 Hz  |  "
            f"delta_s={axial[index]-axial[0]:+.3f} mm  |  rocking={rocking[index]:+.2f} deg", fontsize=9)
        fig.tight_layout()

    animation = FuncAnimation(fig, draw, frames=len(frame_ids), interval=1000.0 / GIF_FPS)
    animation.save(OUT / gif_name, writer=PillowWriter(fps=GIF_FPS), dpi=100)
    plt.close(fig)


def render_contact_sheet(records):
    fig, axes = plt.subplots(2, 2, figsize=(14, 6.5), sharex=True, sharey=True)
    for ax, record in zip(axes.flat, records):
        index = int(record["contact_frame_index"])
        points = record["snapshot"]["points"]
        local = np.column_stack(((points - record["rp0"]).dot(record["c"]),
                                 (points - record["rp0"]).dot(record["n"])))
        colors = np.where(record["region"] == "HEAD", "#d1493f",
                          np.where(record["region"] == "TAIL", "#2864a8", "#383c42"))
        radius = record["radius"]
        ax.axhspan(radius, YLIM[1], color="#d8edf2", alpha=0.55)
        ax.axhspan(YLIM[0], -radius, color="#d8edf2", alpha=0.55)
        ax.axhline(radius, color="#5996a5", lw=1.1); ax.axhline(-radius, color="#5996a5", lw=1.1)
        ax.scatter(local[::2, 0], local[::2, 1], c=colors[::2], s=8, linewidths=0)
        ax.set_xlim(XLIM); ax.set_ylim(YLIM); ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"{record['label']}  t={record['time_s']*1e3:.1f} ms  {record['gradient_mT']:g} mT\n"
                     f"contact={record['contact_now']}  nearest={record['contact_region']}", fontsize=9)
        ax.grid(axis="x", alpha=0.16)
    for ax in axes[1, :]: ax.set_xlabel("canonical +s (mm), fixed camera")
    for ax in axes[:, 0]: ax.set_ylabel("n_routeA (mm)")
    fig.tight_layout()
    fig.savefig(OUT / "S4_HEADFORWARD_LowG_ContactSheet.png", dpi=180)
    plt.close(fig)


def analyze_case(label, job, gif_name, gradient_mT):
    case = OUT / "cases" / job
    identity = json.loads((case / "case_identity.json").read_text())
    if identity["status"] != "SOLVED":
        raise RuntimeError(f"{job} is not solved")
    time, data = load_rp(case)
    c = unit(identity["canonical_plus_s_axis_aba"]); n = unit(identity["n_routeA_aba"]); b = unit(identity["b_routeA_aba"])
    rp0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    axial = data["U"].dot(c)
    velocity = data["V"].dot(c)
    acceleration = np.gradient(velocity, time)
    radial = data["U"] - np.outer(axial, c)
    body_axis = Rotation.from_rotvec(data["UR"]).apply(np.broadcast_to(c, data["UR"].shape))
    rocking = np.degrees(np.arctan2(body_axis.dot(n), body_axis.dot(c)))
    deck = (case / f"{job}.inp").read_text()
    nodes = part_nodes(deck, "Robot_SOLID"); wall = part_nodes(deck, "Pipe_WALL_HELPER")
    initial_s = (nodes - rp0).dot(c)
    region = np.full(len(nodes), "BODY", dtype=object)
    region[initial_s <= initial_s.min() + 0.25] = "TAIL"
    region[initial_s >= initial_s.max() - 0.25] = "HEAD"
    normals, radius = wall_planes(wall, rp0, c, n, b)
    gap_snapshot.normals = normals
    geom_indices, geom = geometry_samples(time, data, nodes, region, rp0, c, n, b, radius, rp0)
    force = contact_force(time, case)
    active, episodes = contact_episode_table(time, force, geom, geom_indices, region)
    bridge = (geom["HEAD_gap_mm"].to_numpy() <= 0.0) & (geom["TAIL_gap_mm"].to_numpy() <= 0.0)
    bridge_ms = longest_duration(geom["time_s"].to_numpy(), bridge) * 1.0e3
    s = rp0.dot(c) + axial
    ts = pd.DataFrame({"time_s": time, "s_mm": s, "delta_s_mm": axial - axial[0],
                       "v_s_mm_s": velocity, "a_s_mm_s2": acceleration,
                       "rocking_deg": rocking, "contact_active": active.astype(int)})
    ts.to_csv(case / f"{job}_v_s_a_s_timeseries.csv", index=False)
    episodes.to_csv(case / f"{job}_contact_episodes.csv", index=False)
    if episodes.empty:
        region_counts = {f"{name}_contact_episodes": 0 for name in ("HEAD", "TAIL", "BODY")}
    else:
        region_counts = {f"{name}_contact_episodes": int(episodes["region"].str.contains(name).sum()) for name in ("HEAD", "TAIL", "BODY")}
    metrics = {
        "case": label, "job": job, "gradient_mT": gradient_mT,
        "delta_s_mm": float(axial[-1] - axial[0]),
        "fraction_v_s_positive": float(np.mean(velocity > 0.0)),
        "longest_near_zero_v_s_ms": float(longest_duration(time, np.abs(velocity) <= NEAR_ZERO_MM_S) * 1.0e3),
        "axial_reversals": reversal_count(velocity),
        "max_forward_velocity_mm_s": float(velocity.max()),
        "final_forward_velocity_mm_s": float(velocity[-1]),
        "fitted_a_eff_mm_s2": fit_effective_acceleration(time, s),
        "rocking_min_deg": float(rocking.min()), "rocking_max_deg": float(rocking.max()),
        "max_radial_COM_drift_mm": float(np.linalg.norm(radial, axis=1).max()),
        "contact_count": int(len(episodes)),
        "longest_contact_ms": float(episodes["duration_ms"].max()) if len(episodes) else 0.0,
        "both_end_bridge_ms": float(bridge_ms),
        "contact_force_threshold_N": CONTACT_THRESHOLD_N,
        "near_zero_threshold_mm_s": NEAR_ZERO_MM_S,
        "contact_episode_csv": f"cases/{job}/{job}_contact_episodes.csv",
        "v_s_a_s_timeseries_csv": f"cases/{job}/{job}_v_s_a_s_timeseries.csv",
        **region_counts,
        "contact_episodes": episodes.to_dict(orient="records"),
    }
    render_case(case, identity, time, data, nodes, region, c, n, b, rp0, radius, rocking, axial, label, gif_name)
    contact_idx = int(np.argmin(np.abs(geom["gap_mm"].to_numpy())))
    contact_time = float(geom.iloc[contact_idx]["time_s"])
    raw_contact = bool(np.any(active[max(0, int(geom.iloc[contact_idx]["index"]) - 1):
                                  min(len(active), int(geom.iloc[contact_idx]["index"]) + 2)]))
    contact_region = "none"
    if len(episodes):
        contact_region = str(episodes.iloc[int(np.argmin(np.abs(episodes["start_s"].to_numpy() - contact_time)))]["region"])
    record = {"label": label, "gradient_mT": gradient_mT, "time_s": contact_time,
              "contact_now": raw_contact, "contact_region": contact_region,
              "contact_frame_index": int(geom.iloc[contact_idx]["index"]),
              "snapshot": gap_snapshot(int(geom.iloc[contact_idx]["index"]), time, data, nodes, region, rp0, c, n, b, radius, rp0),
              "region": region, "rp0": rp0, "c": c, "n": n, "radius": radius}
    return metrics, record


def old_reference():
    path = REPO / "calibration_analysis" / "FastStraightDynamicScreen" / "fast_straight_screen_metrics.json"
    items = json.loads(path.read_text())
    return next(item for item in items if item["job"] == OLD_REF_JOB)


def write_report(metrics):
    screen_keys = (
        "case", "job", "gradient_mT", "delta_s_mm", "fraction_v_s_positive",
        "longest_near_zero_v_s_ms", "axial_reversals", "rocking_min_deg",
        "rocking_max_deg", "max_radial_COM_drift_mm", "contact_count",
        "longest_contact_ms", "both_end_bridge_ms",
    )
    screened = [{key: item[key] for key in screen_keys} for item in metrics]
    (OUT / "lowg_screen_metrics.json").write_text(json.dumps(screened, indent=2) + "\n", encoding="ascii")
    pd.DataFrame(screened).to_csv(OUT / "lowg_screen_metrics.csv", index=False)
    rows = []
    for item in metrics:
        rows.append(f"| {item['case']} | {item['gradient_mT']:.2f} | {item['delta_s_mm']:+.6f} | {item['fraction_v_s_positive']:.4f} | {item['longest_near_zero_v_s_ms']:.4f} | {item['axial_reversals']} | {item['rocking_min_deg']:+.3f}/{item['rocking_max_deg']:+.3f} | {item['max_radial_COM_drift_mm']:.6f} | {item['contact_count']} | {item['longest_contact_ms']:.4f} | {item['both_end_bridge_ms']:.4f} |")
    ref = old_reference()
    report = """# S4 HEAD-FORWARD Low-G Dynamic Screen

These are four new 75 ms Explicit dynamics runs. The only candidate variable is the analytic Gaussian axial gradient: `G=0`, `0.25`, `0.50`, and `1.00 mT`, all with `L=45 mm`. The robot mesh is a true rigid 180 degree flip about fixed `n_routeA`; COM, straight tube, rigid mesh topology, one-wall General Contact, `mu=0.03`, `zeta=0.50`, ReducedHydro coefficients, `B0=10 mT`, initial pose, RouteA gauge, and `dt=1e-7 s` are frozen. These are non-CEL cases.

The fixed orthographic camera is identical in every GIF: horizontal canonical `+s` left-to-right, vertical `n_routeA`, view direction `b_routeA`, and a fixed 20 mm by 5 mm display window. HEAD is the forward/right endpoint and TAIL is the rear/left endpoint. No camera following, auto-fit, camera reversal, or label exchange is used.

| case | G (mT) | delta_s (mm) | fraction v_s > 0 | longest near-zero (ms) | axial reversals | rocking min/max (deg) | max radial drift (mm) | contact count | longest contact (ms) | both-end bridge (ms) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
""" + "\n".join(rows) + f"\n| OLD_S4_REFERENCE_G6 | 6.00 | {ref['delta_s_mm']:+.6f} | {ref['fraction_v_s_positive']:.4f} | {ref['longest_near_zero_v_s_ms']:.4f} | {ref['axial_reversals']} | {ref['rocking_min_deg']:+.3f}/{ref['rocking_max_deg']:+.3f} | {ref['max_radial_COM_drift_mm']:.6f} | {ref['contact_count']} | {ref['longest_contact_ms']:.4f} | {ref['both_end_bridge_ms']:.4f} |\n\n"""
    report += "Per-case `v_s(t)` and derived `a_s(t)` are in each case's `*_v_s_a_s_timeseries.csv`. Contact episodes, episode region, minimum gap, peak force, and reopen flag are in each case's `*_contact_episodes.csv`. Both-end bridge is a 20 us screening-resolution exact-gap check; no publication-grade contact analysis or automatic ranking is assigned.\n\n**USER VISUAL SELECTION REQUIRED**\n"
    (OUT / "S4_HEADFORWARD_LowG_Screen.md").write_text(report, encoding="ascii")


def main():
    metrics, records = zip(*(analyze_case(*item) for item in CASES))
    images = [Image.open(OUT / item[2]) for item in CASES]
    frames = []
    for index in range(max(image.n_frames for image in images)):
        tiles = []
        for image in images:
            image.seek(min(index, image.n_frames - 1)); tiles.append(image.convert("RGB"))
        width = max(tile.width for tile in tiles); height = max(tile.height for tile in tiles)
        canvas = Image.new("RGB", (2 * width, 2 * height), "white")
        for position, tile in enumerate(tiles): canvas.paste(tile, ((position % 2) * width, (position // 2) * height))
        frames.append(canvas)
    frames[0].save(OUT / "S4_HEADFORWARD_LowG_4Way.gif", save_all=True,
                   append_images=frames[1:], duration=round(1000.0 / GIF_FPS), loop=0, optimize=False)
    render_contact_sheet(list(records))
    write_report(list(metrics))
    print(json.dumps(list(metrics), indent=2))
    print("USER VISUAL SELECTION REQUIRED")


if __name__ == "__main__":
    main()
