"""Compute screening-only metrics and render fixed 2D fast-straight GIFs."""
from __future__ import annotations

import json
import math
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
CASES = (
    ("S1", "FAST_PLANAR_15HZ_A12_FORWARD", "FAST_S1_2D.gif"),
    ("S2", "FAST_PLANAR_20HZ_A12_FORWARD", "FAST_S2_2D.gif"),
    ("S3", "FAST_PLANAR_20HZ_A14P343_FORWARD", "FAST_S3_2D.gif"),
    ("S4", "FAST_ELLIPTIC_20HZ_A14P343_X2P5_FORWARD", "FAST_S4_2D.gif"),
)
NEAR_ZERO_MM_S = 1.0
CONTACT_THRESHOLD_N = 1.0e-8
GIF_PHYSICAL_STEP_S = 0.001
GIF_FPS = 20.0
XLIM = (-2.0, 48.0)
YLIM = (-2.5, 2.5)


def unit(vector):
    vector = np.asarray(vector, dtype=float)
    return vector / np.linalg.norm(vector)


def get_array(archive, name):
    keys = [key for key in archive.files if key == name or key.startswith(name + " (Repeated:")]
    return max((archive[key] for key in keys), key=len)


def intervals(flag):
    edge = np.diff(np.r_[False, flag, False].astype(int))
    return list(zip(np.where(edge == 1)[0], np.where(edge == -1)[0] - 1))


def longest(time, flag):
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
                continue
            array = max((archive[key] for key in keys), key=len)
            vector.append(np.interp(time, array[:, 0], array[:, 1]))
        columns.append(np.column_stack(vector))
    archive.close()
    return columns[0] + columns[1]


def reversal_count(velocity):
    signs = np.sign(velocity[np.abs(velocity) > NEAR_ZERO_MM_S])
    return int(np.sum(signs[1:] != signs[:-1])) if len(signs) > 1 else 0


def analyze_case(label, job, gif_name):
    case = OUT / "cases" / job
    identity = json.loads((case / "case_identity.json").read_text())
    if identity["status"] != "SOLVED":
        raise RuntimeError(f"{job} is not solved")
    rp = np.load(case / "private" / "rp_history_private.npz")
    time = get_array(rp, "V1")[:, 0]
    data = {prefix: np.column_stack([
        np.interp(time, get_array(rp, prefix + str(axis))[:, 0], get_array(rp, prefix + str(axis))[:, 1])
        for axis in (1, 2, 3)]) for prefix in ("U", "UR", "V")}
    rp.close()

    c = unit(identity["initial_axis_aba"])
    n = unit(identity["n_routeA_aba"])
    b = unit(identity["b_routeA_aba"])
    rp0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    axial = data["U"].dot(c)
    velocity = data["V"].dot(c)
    radial = data["U"] - np.outer(axial, c)
    axis = Rotation.from_rotvec(data["UR"]).apply(np.broadcast_to(c, data["UR"].shape))
    rocking = np.degrees(np.arctan2(axis.dot(n), axis.dot(c)))

    force = contact_force(time, case)
    active = np.linalg.norm(force, axis=1) > CONTACT_THRESHOLD_N
    events = intervals(active)
    contact_longest = longest(time, active)

    deck = (case / f"{job}.inp").read_text()
    nodes = part_nodes(deck, "Robot_SOLID")
    wall = part_nodes(deck, "Pipe_WALL_HELPER")
    initial_node_s = (nodes - rp0).dot(c)
    region = np.full(len(nodes), "BODY", dtype=object)
    region[initial_node_s <= initial_node_s.min() + 0.25] = "HEAD"
    region[initial_node_s >= initial_node_s.max() - 0.25] = "TAIL"
    normals, radius = wall_planes(wall, rp0, c, n, b)

    gap_times = np.arange(0.0, identity["duration_s"] + 5.0e-5, 1.0e-4)
    gap_ids = np.clip(np.searchsorted(time, gap_times), 0, len(time) - 1)
    both = []
    rel0 = nodes - rp0
    for index in gap_ids:
        points = rp0 + data["U"][index] + Rotation.from_rotvec(data["UR"][index]).apply(rel0)
        relative = points - rp0
        transverse = relative - np.outer(relative.dot(c), c)
        node_gap = radius - np.max(transverse.dot(normals.T), axis=1)
        both.append(node_gap[region == "HEAD"].min() <= 0.0 and
                    node_gap[region == "TAIL"].min() <= 0.0)
    both_bridge = longest(gap_times, np.asarray(both, dtype=bool))

    metrics = {
        "case": label,
        "job": job,
        "delta_s_mm": float(axial[-1] - axial[0]),
        "fraction_v_s_positive": float(np.mean(velocity > 0.0)),
        "longest_near_zero_v_s_ms": 1.0e3 * longest(time, np.abs(velocity) <= NEAR_ZERO_MM_S),
        "axial_reversals": reversal_count(velocity),
        "rocking_min_deg": float(rocking.min()),
        "rocking_max_deg": float(rocking.max()),
        "max_radial_COM_drift_mm": float(np.linalg.norm(radial, axis=1).max()),
        "contact_count": len(events),
        "longest_contact_ms": 1.0e3 * contact_longest,
        "both_end_bridge_ms": 1.0e3 * both_bridge,
        "near_zero_threshold_mm_s": NEAR_ZERO_MM_S,
        "contact_force_threshold_N": CONTACT_THRESHOLD_N,
    }
    render_case_gif(case, identity, time, data, nodes, region, c, n, rp0, radius,
                    rocking, axial, label, gif_name)
    return metrics


def render_case_gif(case, identity, time, data, nodes, region, c, n, rp0, radius,
                    rocking, axial, label, gif_name):
    frame_times = np.arange(0.0, identity["duration_s"] + 0.5 * GIF_PHYSICAL_STEP_S,
                            GIF_PHYSICAL_STEP_S)
    frame_ids = np.clip(np.searchsorted(time, frame_times), 0, len(time) - 1)
    rel0 = nodes - rp0
    colors = np.where(region == "HEAD", "#d1493f", np.where(region == "TAIL", "#2864a8", "#383c42"))
    fig, ax = plt.subplots(figsize=(12.0, 2.2))

    def draw(frame):
        ax.clear()
        index = frame_ids[frame]
        points = rp0 + data["U"][index] + Rotation.from_rotvec(data["UR"][index]).apply(rel0)
        local = np.column_stack(((points - rp0).dot(c), (points - rp0).dot(n)))
        ax.axhspan(radius, YLIM[1], color="#d8edf2", alpha=0.55)
        ax.axhspan(YLIM[0], -radius, color="#d8edf2", alpha=0.55)
        ax.axhline(radius, color="#5996a5", lw=1.2)
        ax.axhline(-radius, color="#5996a5", lw=1.2)
        ax.scatter(local[::2, 0], local[::2, 1], c=colors[::2], s=8, linewidths=0)
        ax.set_xlim(XLIM)
        ax.set_ylim(YLIM)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("canonical +s (mm), fixed camera")
        ax.set_ylabel("n_routeA (mm)")
        ax.grid(axis="x", alpha=0.16)
        ax.text(0.012, 0.96, "HEAD", color="#d1493f", weight="bold", transform=ax.transAxes, va="top")
        ax.text(0.085, 0.96, "TAIL", color="#2864a8", weight="bold", transform=ax.transAxes, va="top")
        ax.set_title(
            f"{label}  |  t={time[index]*1e3:5.1f} ms  |  f={identity['frequency_Hz']:.0f} Hz  |  "
            f"delta_s={axial[index]-axial[0]:+.3f} mm  |  rocking={rocking[index]:+.2f} deg",
            fontsize=10)
        fig.tight_layout()

    animation = FuncAnimation(fig, draw, frames=len(frame_ids), interval=1000.0 / GIF_FPS)
    animation.save(OUT / gif_name, writer=PillowWriter(fps=GIF_FPS), dpi=100)
    plt.close(fig)


def four_way():
    images = [Image.open(OUT / gif_name) for _, _, gif_name in CASES]
    counts = [image.n_frames for image in images]
    frames = []
    for frame_index in range(max(counts)):
        tiles = []
        for image, count in zip(images, counts):
            image.seek(min(frame_index, count - 1))
            tiles.append(image.convert("RGB"))
        width = max(tile.width for tile in tiles)
        height = max(tile.height for tile in tiles)
        canvas = Image.new("RGB", (2 * width, 2 * height), "white")
        for index, tile in enumerate(tiles):
            canvas.paste(tile, ((index % 2) * width, (index // 2) * height))
        frames.append(canvas)
    frames[0].save(OUT / "FastStraightScreen_4Way.gif", save_all=True,
                   append_images=frames[1:], duration=round(1000.0 / GIF_FPS),
                   loop=0, optimize=False)


def write_report(metrics):
    frame = pd.DataFrame(metrics)
    frame.to_csv(OUT / "fast_straight_screen_metrics.csv", index=False)
    (OUT / "fast_straight_screen_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="ascii")
    rows = []
    for item in metrics:
        rows.append(
            f"| {item['case']} | {item['delta_s_mm']:+.6f} | {item['fraction_v_s_positive']:.4f} | "
            f"{item['longest_near_zero_v_s_ms']:.4f} | {item['axial_reversals']} | "
            f"{item['rocking_min_deg']:+.3f} / {item['rocking_max_deg']:+.3f} | "
            f"{item['max_radial_COM_drift_mm']:.6f} | {item['contact_count']} | "
            f"{item['longest_contact_ms']:.4f} | {item['both_end_bridge_ms']:.4f} |")
    report = """# Fast Straight Dynamic Screen

All four cases are new short dynamics runs. They use the same non-CEL rigid model, one-wall General Contact, ReducedHydro, `G=6 mT`, `L=45 mm`, `B0=10 mT`, `dt=1e-7 s`, initial pose, and RouteA gauge.

The GIF camera is fixed and identical: horizontal canonical `+s`, vertical `n_routeA`, view direction `b_routeA`, 50 mm displayed tube length, 1 ms physical frame spacing, and 20 fps playback. Shorter 75 ms cases hold their last frame while S1 completes its final 25 ms in the four-way comparison.

| case | delta_s (mm) | fraction v_s > 0 | longest near-zero v_s (ms) | axial reversals | rocking min/max (deg) | max radial COM drift (mm) | contact count | longest contact (ms) | both-end bridge (ms) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
""" + "\n".join(rows) + """

Near-zero means `|v_s| <= 1 mm/s`; contact uses whole-robot resultant force above `1e-8 N`. Both-end bridge is a screening-resolution exact-gap check at 0.1 ms, not a publication-grade contact audit.

No automatic ranking is assigned.

**USER VISUAL SELECTION REQUIRED**
"""
    (OUT / "Fast_Straight_Dynamic_Screen.md").write_text(report, encoding="ascii")


def main():
    metrics = [analyze_case(*case) for case in CASES]
    four_way()
    write_report(metrics)
    print(json.dumps(metrics, indent=2))
    print("USER VISUAL SELECTION REQUIRED")


if __name__ == "__main__":
    main()
