"""Bounded exact-mesh refinement around GEO-C; no dynamics are run here."""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
SOURCE = REPO / "calibration_analysis" / "F60G015TrueCEL" / "case" / "S4_HEADFORWARD_F60_G0P15_TRUECEL"
BASE_SCRIPTS = REPO / "calibration_analysis" / "S4HeadForwardLowGScreen" / "scripts"
sys.path.insert(0, str(BASE_SCRIPTS))
import analyze_and_render_lowg as base


BASE_D = 0.815
COMMAND_ANGLE = 14.343111711438091
LENGTHS = (2.60, 2.65, 2.70)
DIAMETERS = (0.815,)
OFFSETS_N = (0.0, 0.005, 0.010, 0.015)
# The exact mesh has a smaller HEAD radial envelope than TAIL.  Apply only a
# smooth, local radial profile correction over the final 25% of axial length.
HEAD_RADIAL_CORRECTIONS = tuple(np.linspace(0.035, 0.045, 9))
HEAD_CORRECTION_DIRECTIONS = (-1.0, 1.0)


def parse_part_nodes(deck, part_name):
    part = re.search(rf"(?ms)^\*Part, name={re.escape(part_name)}\s*$.*?^\*End Part\s*$", deck).group(0)
    block = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", part).group(1)
    labels, coords = [], []
    for line in block.splitlines():
        if not line.strip():
            continue
        values = [value.strip() for value in line.split(",")]
        labels.append(int(values[0])); coords.append([float(value) for value in values[1:4]])
    return np.asarray(labels, dtype=int), np.asarray(coords, dtype=float)


def transform_nodes(nodes, rp0, c, n, length, diameter, offset_n, head_radial_correction,
                    head_correction_direction):
    relative = nodes - rp0
    axial = relative.dot(c)
    radial = relative - np.outer(axial, c)
    base_length = float(axial.max() - axial.min())
    axial_normalized = (axial - axial.min()) / base_length
    blend = np.clip((axial_normalized - 0.75) / 0.25, 0.0, 1.0)
    blend = blend * blend * (3.0 - 2.0 * blend)
    scaled_radial = radial * (diameter / BASE_D)
    directed_coordinate = head_correction_direction * scaled_radial.dot(n)
    side_weight = np.clip(directed_coordinate / (0.5 * diameter), 0.0, 1.0)
    local_profile = np.outer(blend * side_weight * head_radial_correction * head_correction_direction, n)
    return rp0 + offset_n * n + np.outer(axial * length / base_length, c) + scaled_radial + local_profile


def gap_state(nodes, regions, robot_center, pipe_center, c, n, b, normals, radius, angle_deg):
    rotated = robot_center + Rotation.from_rotvec(b * math.radians(angle_deg)).apply(nodes - robot_center)
    relative = rotated - pipe_center
    transverse = relative - np.outer(relative.dot(c), c)
    support = transverse.dot(normals.T)
    gap = radius - support.max(axis=1)
    return {name: float(gap[regions == name].min()) for name in ("HEAD", "TAIL", "BODY")}


def touch_angle(nodes, regions, robot_center, pipe_center, c, n, b, normals, radius, region, sign):
    low, high = 0.0, 20.0
    if gap_state(nodes, regions, robot_center, pipe_center, c, n, b, normals, radius, sign * high)[region] > 0.0:
        return None
    if gap_state(nodes, regions, robot_center, pipe_center, c, n, b, normals, radius, 0.0)[region] <= 0.0:
        return 0.0
    for _ in range(32):
        middle = 0.5 * (low + high)
        if gap_state(nodes, regions, robot_center, pipe_center, c, n, b, normals, radius,
                     sign * middle)[region] <= 0.0:
            high = middle
        else:
            low = middle
    return 0.5 * (low + high)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    identity = json.loads((SOURCE / "case_identity.json").read_text())
    deck = (SOURCE / "S4_HEADFORWARD_F60_G0P15_TRUECEL.inp").read_text(encoding="latin1")
    c = base.unit(identity["canonical_plus_s_axis_aba"])
    n = base.unit(identity["n_routeA_aba"])
    b = base.unit(identity["b_routeA_aba"])
    rp0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    _, robot = parse_part_nodes(deck, "Robot_SOLID")
    _, wall = parse_part_nodes(deck, "Pipe_WALL_HELPER")
    normals, wall_radius = base.wall_planes(wall, rp0, c, n, b)
    axial = (robot - rp0).dot(c)
    regions = np.full(len(robot), "BODY", dtype=object)
    regions[axial <= axial.min() + 0.25] = "TAIL"
    regions[axial >= axial.max() - 0.25] = "HEAD"

    rows = []
    for length in LENGTHS:
        for diameter in DIAMETERS:
            for offset_n in OFFSETS_N:
              for head_correction in HEAD_RADIAL_CORRECTIONS:
               for correction_direction in HEAD_CORRECTION_DIRECTIONS:
                center = rp0 + offset_n * n
                candidate = transform_nodes(robot, rp0, c, n, length, diameter, offset_n,
                                            head_correction, correction_direction)
                initial = gap_state(candidate, regions, center, rp0, c, n, b, normals, wall_radius, 0.0)
                plus = gap_state(candidate, regions, center, rp0, c, n, b, normals, wall_radius, COMMAND_ANGLE)
                minus = gap_state(candidate, regions, center, rp0, c, n, b, normals, wall_radius, -COMMAND_ANGLE)
                phase_options = []
                for head_sign, tail_sign in ((1.0, -1.0), (-1.0, 1.0)):
                    head_touch = touch_angle(candidate, regions, center, rp0, c, n, b, normals, wall_radius, "HEAD", head_sign)
                    tail_touch = touch_angle(candidate, regions, center, rp0, c, n, b, normals, wall_radius, "TAIL", tail_sign)
                    if head_touch is None or tail_touch is None:
                        continue
                    phase_options.append((abs(head_touch - tail_touch), head_sign, tail_sign, head_touch, tail_touch))
                if not phase_options:
                    continue
                _, head_sign, tail_sign, head_touch, tail_touch = min(phase_options)
                head_phase = plus if head_sign > 0 else minus
                tail_phase = plus if tail_sign > 0 else minus
                opposite_head_gap = head_phase["TAIL"]
                opposite_tail_gap = tail_phase["HEAD"]
                initial_gap = min(initial.values())
                same_phase_bridge = ((plus["HEAD"] <= 0 and plus["TAIL"] <= 0) or
                                     (minus["HEAD"] <= 0 and minus["TAIL"] <= 0))
                window_penalty = max(13.5 - head_touch, 0.0) + max(head_touch - 14.3, 0.0)
                window_penalty += max(13.0 - tail_touch, 0.0) + max(tail_touch - 14.3, 0.0)
                score = (20.0 * window_penalty + 10.0 * max(abs(head_touch - tail_touch) - 1.0, 0.0) +
                         100.0 * max(-initial_gap, 0.0) + 20.0 * float(same_phase_bridge) +
                         0.15 * abs(length - 2.6) / 0.05 + 0.1 * abs(diameter - BASE_D) / 0.015 +
                         0.05 * abs(offset_n) / 0.015)
                rows.append({
                    "length_mm": length, "diameter_mm": diameter, "L_over_D": length / diameter,
                    "offset_n_mm": offset_n, "head_radial_correction_mm": head_correction,
                    "head_correction_direction_n": int(correction_direction),
                    "head_contact_phase_sign": int(head_sign),
                    "tail_contact_phase_sign": int(tail_sign), "HEAD_touch_deg": head_touch,
                    "TAIL_touch_deg": tail_touch, "touch_difference_deg": abs(head_touch - tail_touch),
                    "initial_min_gap_mm": initial_gap, "HEAD_phase_opposite_end_gap_mm": opposite_head_gap,
                    "TAIL_phase_opposite_end_gap_mm": opposite_tail_gap,
                    "same_phase_STATIC_BOTH_bridge": same_phase_bridge, "score": score,
                })
    table = pd.DataFrame(rows)
    table["passes_geometry_gate"] = (
        table.HEAD_touch_deg.between(13.5, 14.3) &
        table.TAIL_touch_deg.between(13.0, 14.3) &
        (table.touch_difference_deg <= 1.0) &
        (table.initial_min_gap_mm > 0.0) &
        ~table.same_phase_STATIC_BOTH_bridge
    )
    table = table.sort_values(["passes_geometry_gate", "score"], ascending=[False, True]).reset_index(drop=True)
    shortlist = table.head(60).copy()
    shortlist.to_csv(OUT / "micro_geometry_shortlist.csv", index=False)
    passed = table[table.passes_geometry_gate]
    if passed.empty:
        best = table.iloc[0].to_dict()
        failure = {
            "selection_status": "NO_FEASIBLE_GEOMETRY_WITHIN_BOUNDED_REFINEMENT",
            "reason": ("No candidate simultaneously puts HEAD and TAIL touch thresholds in the requested "
                       "windows with <=1 deg difference and avoids a same-phase static BOTH bridge."),
            "best_near_candidate_not_selected": best,
            "asymmetry_source": "shape-driven: centered pipe/COM; HEAD max radius 0.3891 mm vs TAIL 0.4075 mm",
            "candidate_count": int(len(table)),
            "abaqus_run_authorized": False,
        }
        (OUT / "selected_refined_geometry.json").write_text(
            json.dumps(failure, indent=2) + "\n", encoding="ascii")
        print(json.dumps(failure, indent=2))
        return
    selected = passed.iloc[0].to_dict()
    selected["selection_status"] = "SELECTED_FOR_SINGLE_10MS_TRUECEL_GATE"
    selected["asymmetry_source"] = "shape-driven: centered pipe/COM; HEAD max radius 0.3891 mm vs TAIL 0.4075 mm"
    selected["search_bounds"] = {"L_mm": [2.55, 2.75], "D_mm": [0.815, 0.86],
                                 "offset_n_mm_tested": list(OFFSETS_N),
                                 "head_radial_correction_mm": [0.035, 0.045]}
    selected["candidate_count"] = int(len(table))
    (OUT / "selected_refined_geometry.json").write_text(json.dumps(selected, indent=2) + "\n", encoding="ascii")

    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    shown = shortlist.head(8)
    labels = [f"L={row.length_mm:.2f}, D={row.diameter_mm:.3f}, dn={row.offset_n_mm:+.3f}, dRh={row.head_radial_correction_mm:.4f} ({row.head_correction_direction_n:+.0f}n)" for _, row in shown.iterrows()]
    y = np.arange(len(shown))
    ax.scatter(shown.HEAD_touch_deg, y, label="HEAD touch", marker="o", color="#d1493f")
    ax.scatter(shown.TAIL_touch_deg, y, label="TAIL touch", marker="s", color="#2864a8")
    ax.axvspan(13.5, 14.3, color="#d1493f", alpha=0.08)
    ax.axvspan(13.0, 14.3, color="#2864a8", alpha=0.08)
    ax.set_yticks(y, labels); ax.invert_yaxis(); ax.set_xlabel("exact touch angle (deg)")
    ax.set_title("Bounded GEO-C micro-refinement shortlist"); ax.grid(axis="x", alpha=0.2); ax.legend()
    fig.tight_layout(); fig.savefig(OUT / "micro_geometry_shortlist.png", dpi=180); plt.close(fig)
    print(json.dumps(selected, indent=2))


if __name__ == "__main__":
    main()
