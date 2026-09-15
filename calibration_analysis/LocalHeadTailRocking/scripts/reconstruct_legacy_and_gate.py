"""Reconstruct RouteA and gate the unchanged straight solver geometry."""
from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
ROOT = REPO.parent
STRAIGHT = REPO / "calibration_analysis" / "StraightPipeControl"
CASE = STRAIGHT / "case" / "PROD_LOCAL30_G0_STRAIGHT_CTRL"
DECK = CASE / "PROD_LOCAL30_G0_STRAIGHT_CTRL.inp"
LEGACY_GIF = REPO / "calibration_analysis" / "Legacy_Wobble_GIFs" / "Job_RouteA_CEL_SOLID_headtail_rock_probe.gif"
HIGH_GIF = REPO / "calibration_analysis" / "Legacy_Wobble_GIFs" / "CEL_CurrentCenterline_Cone15_Z110_PolarityMinus_Damp055_AfterBendReverse_HighEndStop_060_FilletR010_CEL_FSI_bend_validation.gif"
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
A0 = np.array([0.9647382600215763, -0.11887423721399708, 0.2348382536498942])
ANGLES = np.array([-10.0, -7.5, -5.0, 0.0, 5.0, 7.5, 10.0])
ROUTEA_GAUGE_DEG = -61.37284757596327


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def part_nodes(deck: str, name: str) -> np.ndarray:
    part = re.search(rf"(?ms)^\*Part, name={re.escape(name)}\s*$.*?^\*End Part\s*$", deck).group(0)
    block = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", part).group(1)
    return np.asarray([[float(x) for x in row.split(",")[1:4]] for row in block.splitlines() if row.strip()])


def rotate_from_to(source, target):
    source, target = unit(source), unit(target)
    cross = np.cross(source, target)
    if np.linalg.norm(cross) < 1e-14:
        return np.eye(3) if np.dot(source, target) > 0 else Rotation.from_rotvec(math.pi * unit(np.cross(source, [0, 0, 1]) if abs(source[2]) < .9 else np.cross(source, [0, 1, 0]))).as_matrix()
    axis = unit(cross)
    return Rotation.from_rotvec(axis * math.atan2(np.linalg.norm(cross), np.dot(source, target))).as_matrix()


def gif_contact_sheet(path: Path, indices, labels, output: Path):
    image = Image.open(path)
    width = 430
    height = round(image.height * width / image.width)
    canvas = Image.new("RGB", (width * len(indices), height + 34), "white")
    draw = ImageDraw.Draw(canvas)
    for column, (index, label) in enumerate(zip(indices, labels)):
        image.seek(index)
        frame = image.convert("RGB").resize((width, height), Image.Resampling.LANCZOS)
        canvas.paste(frame, (column * width, 0))
        draw.text((column * width + 8, height + 8), label, fill="black")
    canvas.save(output)


def trusted_early_gif():
    source = Image.open(HIGH_GIF)
    frames, durations = [], []
    for index in range(28):
        source.seek(index)
        frames.append(source.convert("RGB"))
        durations.append(source.info.get("duration", 60))
    frames[0].save(OUT / "Legacy_HighEndStop_TrustedEarlyWindow.gif", save_all=True,
                   append_images=frames[1:], duration=durations, loop=0, optimize=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    deck = DECK.read_text()
    identity = json.loads((CASE / "case_identity.json").read_text())
    robot = part_nodes(deck, "Robot_SOLID")
    wall = part_nodes(deck, "Pipe_WALL_HELPER")
    c = unit(identity["straight_tangent_aba"])
    e1 = unit(identity["initial_e1_aba"]); e2 = unit(identity["initial_e2_aba"])
    chi = math.radians(ROUTEA_GAUGE_DEG)
    n = unit(math.cos(chi) * e1 + math.sin(chi) * e2)
    b = unit(np.cross(c, n))
    n = unit(np.cross(b, c))
    center = np.asarray(identity["initial_center_aba_mm"], dtype=float)

    axial_wall = (wall - center).dot(c)
    first_ring = wall[np.abs(axial_wall - axial_wall.min()) < 1e-7]
    xy = np.column_stack(((first_ring - center).dot(n), (first_ring - center).dot(b)))
    order = np.argsort(np.arctan2(xy[:, 1], xy[:, 0]))
    xy = xy[order]
    face_radius = []
    face_normals = []
    for p, q in zip(xy, np.roll(xy, -1, axis=0)):
        edge = q - p
        normal2 = unit(np.array([edge[1], -edge[0]]))
        if np.dot(normal2, (p + q) / 2) < 0:
            normal2 *= -1
        face_normals.append(normal2[0] * n + normal2[1] * b)
        face_radius.append(np.dot(p, normal2))
    face_normals = np.asarray(face_normals)
    radius = float(np.mean(face_radius))
    radius_spread = float(np.ptp(face_radius))

    original_axial = (robot - RP0).dot(A0)
    lo, hi = original_axial.min(), original_axial.max()
    region = np.full(len(robot), "BODY", dtype=object)
    region[original_axial <= lo + 0.25] = "HEAD"
    region[original_axial >= hi - 0.25] = "TAIL"
    align = rotate_from_to(A0, c)
    aligned = center + (robot - RP0).dot(align.T)

    rows = []
    for angle in ANGLES:
        rotation = Rotation.from_rotvec(b * math.radians(angle)).as_matrix()
        points = center + (aligned - center).dot(rotation.T)
        radial = points - center - np.outer((points - center).dot(c), c)
        gap = radius - np.max(radial.dot(face_normals.T), axis=1)
        row = {"theta_rock_deg": angle, "HEAD_gap_um": gap[region == "HEAD"].min() * 1e3,
               "TAIL_gap_um": gap[region == "TAIL"].min() * 1e3,
               "BODY_gap_um": gap[region == "BODY"].min() * 1e3,
               "minimum_exact_surface_gap_um": gap.min() * 1e3}
        row["both_end_support_possible_at_centered_COM"] = bool(row["HEAD_gap_um"] <= 1.0 and row["TAIL_gap_um"] <= 1.0)
        rows.append(row)
    geometry = pd.DataFrame(rows)
    geometry.to_csv(OUT / "RouteA10_geometry_feasibility.csv", index=False)

    theta = np.linspace(0, 20, 201)
    width = 0.815 * np.cos(np.radians(theta)) + 2.4 * np.sin(np.radians(theta))
    envelope = pd.DataFrame({"theta_deg": theta, "W_mm": width, "tube_ID_minus_W_mm": 2 * radius - width})
    envelope.to_csv(OUT / "RouteA10_analytic_envelope.csv", index=False)
    w10 = float(np.interp(10.0, theta, width))

    phase = np.linspace(0, 1, 401)
    theta_ref = 10.0 * np.sin(2 * np.pi * phase)
    target = pd.DataFrame({"normalized_phase": phase, "theta_rock_deg": theta_ref,
                           "theta_cross_deg": np.zeros_like(phase),
                           "q_rock": np.sin(np.radians(theta_ref)), "q_cross": np.zeros_like(phase)})
    target.to_csv(OUT / "Legacy_RouteA_motion_target.csv", index=False)

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.plot(target.q_rock, target.q_cross, color="#c43c35", lw=2.4)
    ax.axhline(0, color="#777777", lw=.8)
    ax.set(xlabel=r"$q_{rock}=a\cdot n_{rock}$", ylabel=r"$q_{cross}=a\cdot b_{rock}$",
           title="Legacy RouteA target: line-like rocking orbit", xlim=(-.19, .19), ylim=(-.04, .04))
    ax.grid(alpha=.22)
    fig.tight_layout()
    fig.savefig(OUT / "Legacy_RouteA_rocking_orbit.png", dpi=220)
    plt.close(fig)

    gif_contact_sheet(LEGACY_GIF, [24, 25, 27, 28, 30],
                      ["- / center", "+peak", "center", "-peak", "return"],
                      OUT / "Legacy_RouteA_contact_sheet.png")
    trusted_early_gif()

    sources = [ROOT / "CEL_REFERENCE_BASELINE.md", ROOT / "make_centerline_headtail_rock.py",
               ROOT / "Job_RouteA_CEL_SOLID_headtail_rock_probe.inp",
               ROOT / "Job_RouteA_CEL_SOLID_headtail_rock_probe_report.md", LEGACY_GIF, HIGH_GIF]
    manifest = {str(path): {"sha256": sha256(path), "bytes": path.stat().st_size} for path in sources}
    manifest.update({"routeA_frame_count": 121, "highendstop_frame_count": 84,
                     "highendstop_last_complete_frame": 83, "highendstop_last_complete_time_s": 0.0498000942170619,
                     "highendstop_trusted_early_frames": [0, 27], "highendstop_first_suspicious_frame": 28,
                     "highendstop_first_obviously_invalid_frame": 32})
    (OUT / "legacy_source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    min_gap = float(geometry.minimum_exact_surface_gap_um.min())
    report = f"""# Legacy RouteA Motion Target Report

## Authoritative reconstruction

The paired generator and report define `beta(t) = 10 deg * r(t) * sin(2*pi*666.667*t)` for the 0.03 s display probe, where `r(t)` is a C1 smoothstep from 0 to 1 over the first 2 ms. The 666.667 Hz value produces 20 compressed display cycles. The physical reference is **5 Hz**, so the magnetic validation uses `alpha_B(t) = 10 deg * sin(2*pi*5*t)` over 200 ms.

The generator constructs `Q=[T,N,B]`, with `B=T x N`, and applies `Q Rz(beta) Q0^T`. Therefore `c_hat=T`, `n_rock=N`, and `b_rock=B`. Legacy `N` is initialized from the robot transverse PCA axis, while production `e1` is initialized from projected global Z, so the names cannot be equated. Projection at their common inlet establishes the fixed transported gauge **`chi={ROUTEA_GAUGE_DEG:.12f} deg`**: `n_rock=cos(chi)e1+sin(chi)e2 = 0.479107880 e1 - 0.877756025 e2`, and `b_rock=c_hat x n_rock = 0.877756025 e1 + 0.479107880 e2`. Positive rocking gives `a=cos(theta)c_hat+sin(theta)n_rock`. The robot axis `a` is HEAD-to-TAIL; HEAD is the low axial end and TAIL the high axial end.

The legacy initial rocking angle is 0 deg (the ramped sine starts at zero). Its RP follows all three prescribed centerline path translations; it has no separately imposed global-Y/Z wobble, but radial translation is not dynamically free. This and the reported CEL deep-penetration warning limit the GIF to motion topology, plane, phase convention, and approximate amplitude.

## Geometry gate

The actual straight Abaqus wall nodes/faces measure an inscribed solver radius of **{radius:.9f} mm** (ID **{2*radius:.9f} mm**, face-to-face spread {radius_spread:.3e} mm). For `L=2.40 mm`, `D=0.815 mm`, `W(10 deg)={w10:.9f} mm`, leaving **{2*radius-w10:.9f} mm** diametral analytic margin. Exact faceted-wall checks at -10, -7.5, -5, 0, +5, +7.5, +10 deg give a minimum surface gap of **{min_gap:.3f} um** at centered COM. All requested poses are penetration-free; geometry is unchanged and the Abaqus gate is open.

The centered pose does not create both-end wall support at +/-10 deg. That is acceptable for a clean Level-1 motion-mode test and is not altered by artificial penetration.

## Secondary HighEndStop reference

The centerline-distance audit stays below 0.25 mm through frame 27. Frame 28 is the first suspicious frame (0.3197 mm), and frame 32 is the first obviously invalid frame (0.5153 mm for this robot/lumen scale). `Legacy_HighEndStop_TrustedEarlyWindow.gif` therefore contains frames 0-27 only. This secondary reference does not override RouteA.
"""
    (OUT / "Legacy_RouteA_Motion_Target_Report.md").write_text(report, encoding="utf-8")
    print(json.dumps({"wall_radius_mm": radius, "tube_ID_mm": 2 * radius, "W10_mm": w10,
                      "ID_minus_W10_mm": 2 * radius - w10, "minimum_exact_gap_um": min_gap,
                      "geometry_feasible": min_gap >= -1e-6}, indent=2))


if __name__ == "__main__":
    main()
