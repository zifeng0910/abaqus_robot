"""Native-history geometric HEAD/TAIL wall exposure; no contact-force inference."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from analyze_f100_g2p20_final import kinematics


ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "case" / "TRUECEL_B0P11_G2P20_A14P5_F100_FAST"


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def part_nodes(deck, name):
    block = re.search(r"(?ms)^\*Part, name=" + re.escape(name) + r"\s*$.*?^\*End Part\s*$", deck).group()
    body = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", block).group(1)
    return np.array([[float(x) for x in line.split(",")[1:4]] for line in body.splitlines() if line.strip()])


def geometry(case, ident):
    deck = (case / f"{case.name}.inp").read_text(encoding="latin1")
    robot, wall = part_nodes(deck, "Robot_SOLID"), part_nodes(deck, "Pipe_WALL_HELPER")
    rp0 = np.asarray(ident["initial_center_aba_mm"], float)
    c, n, rock = [unit(ident[k]) for k in ("canonical_plus_s_axis_aba", "n_routeA_aba", "b_routeA_aba")]
    pipe = rp0 - float(ident["initial_axial_shift_mm"]) * c - float(ident["radial_offset_n_mm"]) * n
    axial = (wall - pipe) @ c
    ring = wall[np.abs(axial - axial.min()) < 1e-7]
    xy = np.column_stack(((ring - pipe) @ n, (ring - pipe) @ rock))
    xy = xy[np.argsort(np.arctan2(xy[:, 1], xy[:, 0]))]
    normals, radii = [], []
    for p, q in zip(xy, np.roll(xy, -1, axis=0)):
        edge = q - p
        nv = unit([edge[1], -edge[0]])
        if nv @ (p + q) < 0:
            nv = -nv
        normals.append(nv[0] * n + nv[1] * rock)
        radii.append(p @ nv)
    material = (robot - rp0) @ c
    points = {"TAIL": robot[material <= material.min() + .25] - rp0,
              "HEAD": robot[material >= material.max() - .25] - rp0}
    return rp0, pipe, c, n, np.asarray(normals), float(np.mean(radii)), points


def exposure(case, continuation):
    ident, t, u, ur, _, _ = kinematics(case, continuation)
    rp0, pipe, c, n, normals, radius, points = geometry(case, ident)
    gaps = {side: np.empty(len(t)) for side in points}
    angle = np.empty(len(t))
    for begin in range(0, len(t), 512):
        end = min(begin + 512, len(t))
        rot = Rotation.from_rotvec(ur[begin:end]).as_matrix()
        axis = np.einsum("tij,j->ti", rot, c)
        angle[begin:end] = np.rad2deg(np.arctan2(axis @ n, axis @ c))
        center = rp0 + u[begin:end] - pipe
        center_projection = center @ normals.T
        for side, relative in points.items():
            rotated = np.einsum("tij,pj->tpi", rot, relative)
            projection = np.einsum("tpi,ji->tpj", rotated, normals)
            gaps[side][begin:end] = radius - np.max(projection + center_projection[:, None, :], axis=(1, 2))
    return t, angle, gaps


def dwell(t, gap, threshold):
    # Linear crossing interpolation on every native-history interval.
    a, b = gap[:-1] - threshold, gap[1:] - threshold
    dt = np.diff(t)
    both = (a < 0) & (b < 0)
    enter = (a >= 0) & (b < 0)
    leave = (a < 0) & (b >= 0)
    fraction = np.divide(np.abs(a), np.abs(a) + np.abs(b), out=np.zeros_like(a), where=(np.abs(a) + np.abs(b)) > 0)
    return float(np.sum(dt[both]) + np.sum(dt[enter] * (1 - fraction[enter])) + np.sum(dt[leave] * fraction[leave]))


def summarize(t, angle, gaps):
    rows = []
    for cycle in range(1, min(3, int(np.ceil((t[-1] - 1e-10) / .01))) + 1):
        lo, hi = (cycle - 1) * .01, min(cycle * .01, t[-1])
        mask = (t >= lo - 1e-12) & (t <= hi + 1e-12)
        tc, ac = t[mask], angle[mask]
        if len(tc) < 2:
            continue
        for side, full in gaps.items():
            g = full[mask]
            idx = int(np.argmin(g))
            rows.append({"cycle": cycle, "complete_cycle": bool(hi >= cycle * .01 - 1e-8),
                         "side": side, "min_signed_gap_mm": float(g[idx]),
                         "closest_time_ms": float(tc[idx] * 1000), "rocking_at_closest_deg": float(ac[idx]),
                         "dwell_gap_negative_ms": 1000 * dwell(tc, g, 0),
                         "dwell_within_5um_ms": 1000 * dwell(tc, g, .005),
                         "dwell_within_10um_ms": 1000 * dwell(tc, g, .010)})
    return rows


def main():
    t, angle, gaps = exposure(CASE, True)
    rows = summarize(t, angle, gaps)
    output = ROOT / "F100_G2P20_WALL_EXPOSURE_NATIVE.csv"
    with output.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"samples": len(t), "duration_ms": float(t[-1] * 1000),
                      "angle_min_deg": float(angle.min()), "angle_max_deg": float(angle.max()),
                      "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
