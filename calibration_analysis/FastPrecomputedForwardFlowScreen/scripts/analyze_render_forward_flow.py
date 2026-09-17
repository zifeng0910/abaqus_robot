"""Analyze and render the one-cycle G0.30 FAST_SURROGATE forward screen."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
JOB = sys.argv[1] if len(sys.argv) > 1 else "FAST_PRECOMP_F120_G0P30_FLOW_1CYCLE"
CASE = OUT / "case" / JOB
PRIVATE = CASE / "private"
PARENT = REPO / "calibration_analysis" / "FastPrecomputedNoReboundScreen" / "case" / "FAST_PRECOMP_F120_DUALEND_NOREBOUND"
BLUE, RED, BODY = "#2864a8", "#d1493f", "#383c42"


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def vector(archive, prefix, time):
    return np.column_stack([np.interp(time, archive[prefix + str(i)][:, 0], archive[prefix + str(i)][:, 1]) for i in (1, 2, 3)])


def part_nodes(deck, name):
    part = re.search(r"(?ms)^\*Part, name=" + re.escape(name) + r"\s*$.*?^\*End Part\s*$", deck).group(0)
    block = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", part).group(1)
    rows = [[float(x) for x in line.split(",")[:4]] for line in block.splitlines() if line.strip()]
    return np.asarray([row[0] for row in rows], int), np.asarray([row[1:] for row in rows], float)


def wall_planes(wall, center, c, n, b):
    s = (wall - center).dot(c)
    ring = wall[np.abs(s - s.min()) < 1e-7]
    xy = np.column_stack(((ring - center).dot(n), (ring - center).dot(b)))
    xy = xy[np.argsort(np.arctan2(xy[:, 1], xy[:, 0]))]
    normals, radii = [], []
    for p, q in zip(xy, np.roll(xy, -1, axis=0)):
        edge = q - p
        normal2 = unit([edge[1], -edge[0]])
        if np.dot(normal2, p + q) < 0:
            normal2 *= -1
        normals.append(normal2[0] * n + normal2[1] * b)
        radii.append(np.dot(p, normal2))
    return np.asarray(normals), float(np.mean(radii))


def runs(active):
    edge = np.diff(np.r_[False, active, False].astype(int))
    return list(zip(np.where(edge == 1)[0], np.where(edge == -1)[0] - 1))


def parent_first_cycle(c):
    archive = np.load(PARENT / "private" / "rp_history_private.npz")
    time = archive["V1"][:, 0]
    U, V = vector(archive, "U", time), vector(archive, "V", time)
    archive.close()
    keep = time <= 1.0 / 120.0 + 1e-9
    return {"delta_s_mm": float(U[keep][-1].dot(c) - U[keep][0].dot(c)),
            "fraction_v_s_positive": float(np.mean(V[keep].dot(c) > 0))}


def main():
    identity = json.loads((CASE / "case_identity.json").read_text())
    c, n, b = (unit(identity[key]) for key in ("canonical_plus_s_axis_aba", "n_routeA_aba", "b_routeA_aba"))
    rp0 = np.asarray(identity["initial_center_aba_mm"], float)
    pipe_center = rp0 - identity["radial_offset_n_mm"] * n
    deck = (CASE / (JOB + ".inp")).read_text(encoding="latin1")
    labels, nodes = part_nodes(deck, "Robot_SOLID")
    _, wall = part_nodes(deck, "Pipe_WALL_HELPER")
    normals, radius = wall_planes(wall, pipe_center, c, n, b)
    s0 = (nodes - rp0).dot(c)
    region = np.full(len(nodes), "BODY", dtype=object)
    region[s0 <= s0.min() + .25] = "TAIL"
    region[s0 >= s0.max() - .25] = "HEAD"

    rp = np.load(PRIVATE / "rp_history_private.npz")
    time = rp["V1"][:, 0]
    U, UR, V = (vector(rp, prefix, time) for prefix in ("U", "UR", "V"))
    rp.close()
    axial = U.dot(c); vs = V.dot(c)
    body_axis = Rotation.from_rotvec(UR).apply(np.broadcast_to(c, UR.shape))
    rocking = np.degrees(np.arctan2(body_axis.dot(n), body_axis.dot(c)))

    field = np.load(PRIVATE / "robot_contact_fields_private.npz")
    ctime = field["time"].astype(float)
    if not np.array_equal(field["node_labels"], labels):
        raise RuntimeError("node ordering mismatch")
    normal = field["CNORMF"].astype(float)
    field.close()
    Ui = np.column_stack([np.interp(ctime, time, U[:, i]) for i in range(3)])
    URi = np.column_stack([np.interp(ctime, time, UR[:, i]) for i in range(3)])
    positions = np.asarray([rp0 + Ui[k] + Rotation.from_rotvec(URi[k]).apply(nodes - rp0) for k in range(len(ctime))])
    active, gaps, peak_force, episode_count = {}, {}, {}, {}
    for name in ("HEAD", "TAIL"):
        mask = region == name
        values = []
        for pts in positions:
            rel = pts[mask] - pipe_center
            transverse = rel - np.outer(rel.dot(c), c)
            values.append(radius - transverse.dot(normals.T).max(axis=1).max())
        gaps[name] = np.asarray(values)
        force = np.linalg.norm(normal[:, mask, :].sum(axis=1), axis=1)
        active[name] = (force > 1e-10) & (gaps[name] <= .015)
        peak_force[name] = float(force.max())
        episode_count[name] = len(runs(active[name]))

    parent = parent_first_cycle(c)
    metrics = {
        "case": JOB, "duration_ms": float(time[-1] * 1e3), "frequency_Hz": 120.0,
        "gradient_mT": identity["gradient_mT"], "U_flow_mm_s": 10.0, "fluid_mode": "FAST_SURROGATE_FLUID",
        "delta_s_mm": float(axial[-1] - axial[0]), "fraction_v_s_positive": float(np.mean(vs > 0)),
        "mean_v_s_mm_s": float(np.mean(vs)), "final_v_s_mm_s": float(vs[-1]),
        "rocking_min_deg": float(rocking.min()), "rocking_max_deg": float(rocking.max()),
        "tumble": bool(np.any(body_axis.dot(c) < 0)),
        "HEAD_contact_episodes": episode_count["HEAD"], "TAIL_contact_episodes": episode_count["TAIL"],
        "HEAD_peak_normal_force_N": peak_force["HEAD"], "TAIL_peak_normal_force_N": peak_force["TAIL"],
        "parent_G0p15_first_cycle": parent,
        "delta_s_improvement_mm": float((axial[-1] - axial[0]) - parent["delta_s_mm"]),
        "forward_gate_passed": bool(axial[-1] - axial[0] > 0 and np.mean(vs > 0) > .5),
        "socket_calls": 0, "wallclock_s": identity["wallclock_s"],
    }
    (OUT / (JOB + "_metrics.json")).write_text(json.dumps(metrics, indent=2) + "\n", encoding="ascii")
    pd.DataFrame({"time_s": time, "delta_s_mm": axial - axial[0], "v_s_mm_s": vs,
                  "rocking_deg": rocking}).to_csv(OUT / (JOB + "_timeseries.csv"), index=False)

    frame_ids = np.unique(np.linspace(0, len(ctime) - 1, 60).astype(int))
    colors = np.where(region == "HEAD", RED, np.where(region == "TAIL", BLUE, BODY))
    rock_i = np.interp(ctime, time, rocking); axial_i = np.interp(ctime, time, axial)
    fig, ax = plt.subplots(figsize=(10, 3))
    def draw(frame):
        ax.clear(); k = int(frame_ids[frame]); pts = positions[k]
        x, y = (pts - pipe_center).dot(c), (pts - pipe_center).dot(n)
        ax.axhspan(radius, .9, color="#d8edf2", alpha=.55); ax.axhspan(-.9, -radius, color="#d8edf2", alpha=.55)
        ax.axhline(radius, color="#5996a5"); ax.axhline(-radius, color="#5996a5")
        ax.scatter(x[::2], y[::2], c=colors[::2], s=7, linewidths=0)
        tail = [x[region == "TAIL"].mean(), y[region == "TAIL"].mean()]
        head = [x[region == "HEAD"].mean(), y[region == "HEAD"].mean()]
        ax.text(tail[0]-.06, tail[1], "TAIL", color=BLUE, weight="bold", ha="right", va="center")
        ax.text(head[0]+.06, head[1], "HEAD", color=RED, weight="bold", ha="left", va="center")
        ax.text(.02, .93, "TAIL:{}  HEAD:{}".format("CONTACT" if active["TAIL"][k] else "free", "CONTACT" if active["HEAD"][k] else "free"),
                transform=ax.transAxes, va="top", weight="bold", fontsize=9)
        ax.set(xlim=(-3, 3), ylim=(-.9, .9), aspect="equal", xlabel="canonical +s (mm), left to right", ylabel="n_routeA (mm)")
        ax.grid(axis="x", alpha=.15)
        ax.set_title("{} | t={:.3f} ms | G={:.2f} mT | flow=+10 mm/s | rock={:+.2f} deg | delta_s={:+.4f} mm".format(
            JOB, ctime[k]*1e3, identity["gradient_mT"], rock_i[k], axial_i[k]-axial_i[0]), fontsize=8)
        fig.tight_layout()
    FuncAnimation(fig, draw, frames=len(frame_ids), interval=1000/30).save(OUT / (JOB + ".gif"), writer=PillowWriter(fps=30), dpi=100)
    plt.close(fig)

    conclusion = "PASSED_FORWARD_SHORT_SCREEN" if metrics["forward_gate_passed"] else "FAILED_FORWARD_SHORT_SCREEN"
    report = """# FAST_SURROGATE forward short screen

Classification: `{classification}`

- delta_s: {delta:+.8f} mm
- fraction v_s > 0: {fraction:.6f}
- rocking: {rmin:+.4f} to {rmax:+.4f} deg
- HEAD/TAIL contact episodes: {head}/{tail}
- HEAD/TAIL peak normal force: {hf:.8g}/{tf:.8g} N
- parent G0.15 first-cycle delta_s: {parent:+.8f} mm
- improvement: {improvement:+.8f} mm
- fluid: FAST_SURROGATE relative-velocity model, U_flow=+10 mm/s, previous coefficients unchanged
- socket calls: 0
""".format(classification=conclusion, delta=metrics["delta_s_mm"], fraction=metrics["fraction_v_s_positive"],
           rmin=metrics["rocking_min_deg"], rmax=metrics["rocking_max_deg"], head=episode_count["HEAD"], tail=episode_count["TAIL"],
           hf=peak_force["HEAD"], tf=peak_force["TAIL"], parent=parent["delta_s_mm"], improvement=metrics["delta_s_improvement_mm"])
    (OUT / (JOB + "_Report.md")).write_text(report, encoding="ascii")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
