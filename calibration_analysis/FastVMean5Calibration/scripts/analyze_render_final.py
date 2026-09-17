"""Write calibration results and render the selected five-cycle verification."""
from __future__ import annotations

import json
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
sys.path.insert(0, str(HERE))
import analyze_initial_candidates as calibration

JOB = "FAST_PRECOMP_F120_VMEAN5_SELECTED"
CASE = OUT / "case" / JOB
PRIVATE = CASE / "private"
GIF = OUT / "FAST_PRECOMP_F120_VMEAN5_5CYCLES.gif"
BLUE, RED, BODY = "#2864a8", "#d1493f", "#383c42"
PERIOD = 1.0 / 120.0


def main():
    initial = json.loads((OUT / "initial_candidate_metrics.json").read_text())
    refinement = calibration.candidate("FAST_VMEAN5_G1P184_UFLOW0P5_2CYCLES")
    final = calibration.candidate(JOB)
    (OUT / "refinement_metrics.json").write_text(json.dumps(refinement, indent=2) + "\n", encoding="ascii")

    identity = json.loads((CASE / "case_identity.json").read_text())
    c, n, b = (calibration.base.unit(identity[key]) for key in
               ("canonical_plus_s_axis_aba", "n_routeA_aba", "b_routeA_aba"))
    rp0 = np.asarray(identity["initial_center_aba_mm"], float)
    pipe_center = rp0 - identity["radial_offset_n_mm"] * n
    deck = (CASE / (JOB + ".inp")).read_text(encoding="latin1")
    labels, nodes = calibration.base.part_nodes(deck, "Robot_SOLID")
    _, wall = calibration.base.part_nodes(deck, "Pipe_WALL_HELPER")
    normals, radius = calibration.base.wall_planes(wall, pipe_center, c, n, b)
    s0 = (nodes - rp0).dot(c)
    region = np.full(len(nodes), "BODY", dtype=object)
    region[s0 <= s0.min() + .25] = "TAIL"
    region[s0 >= s0.max() - .25] = "HEAD"

    rp = np.load(PRIVATE / "rp_history_private.npz")
    time = rp["V1"][:, 0]
    U, UR, V = (calibration.base.vector(rp, prefix, time) for prefix in ("U", "UR", "V"))
    rp.close()
    axial, vs = U.dot(c), V.dot(c)
    body_axis = Rotation.from_rotvec(UR).apply(np.broadcast_to(c, UR.shape))
    rocking = np.degrees(np.arctan2(body_axis.dot(n), body_axis.dot(c)))
    cycle_edges = np.arange(6, dtype=float) * PERIOD
    edge_s = np.interp(cycle_edges, time, axial)
    cycle_means = np.diff(edge_s) / PERIOD
    cycles_2_5_mean = float((edge_s[5] - edge_s[1]) / (4 * PERIOD))
    final.update({
        "overall_5cycle_mean_mm_s": float((edge_s[5] - edge_s[0]) / (5 * PERIOD)),
        "cycles_2_to_5_mean_mm_s": cycles_2_5_mean,
        "per_cycle_mean_mm_s": [float(value) for value in cycle_means],
        "five_cycle_target_passed": bool(4.5 <= cycles_2_5_mean <= 5.5),
        "two_cycle_calibration_target_passed": bool(4.5 <= refinement["v_mean_mm_s"] <= 5.5),
        "interpretation": "two-cycle calibration matched, but five-cycle verification accelerated and did not reach a steady 5 mm/s mean",
    })
    (OUT / "five_cycle_metrics.json").write_text(json.dumps(final, indent=2) + "\n", encoding="ascii")

    comparison = initial + [refinement]
    pd.DataFrame(comparison).to_csv(OUT / "calibration_candidate_metrics.csv", index=False)
    pd.DataFrame({"cycle": np.arange(1, 6), "mean_v_s_mm_s": cycle_means}).to_csv(
        OUT / "five_cycle_per_cycle_speed.csv", index=False)

    field = np.load(PRIVATE / "robot_contact_fields_private.npz")
    ctime = field["time"].astype(float)
    if not np.array_equal(field["node_labels"], labels):
        raise RuntimeError("node ordering mismatch")
    normal = field["CNORMF"].astype(float)
    shear = field["CSHEARF"].astype(float)
    field.close()
    Ui = np.column_stack([np.interp(ctime, time, U[:, i]) for i in range(3)])
    URi = np.column_stack([np.interp(ctime, time, UR[:, i]) for i in range(3)])
    positions = np.asarray([rp0 + Ui[k] + Rotation.from_rotvec(URi[k]).apply(nodes - rp0)
                            for k in range(len(ctime))])
    active = {}
    for name in ("HEAD", "TAIL"):
        mask = region == name
        gap = []
        for pts in positions:
            rel = pts[mask] - pipe_center
            transverse = rel - np.outer(rel.dot(c), c)
            gap.append(radius - transverse.dot(normals.T).max(axis=1).max())
        force = np.linalg.norm((normal[:, mask, :] + shear[:, mask, :]).sum(axis=1), axis=1)
        active[name] = (force > calibration.FORCE_THRESHOLD) & (np.asarray(gap) <= .015)

    frame_ids = np.unique(np.linspace(0, len(ctime) - 1, 150).astype(int))
    colors = np.where(region == "HEAD", RED, np.where(region == "TAIL", BLUE, BODY))
    rock_i = np.interp(ctime, time, rocking)
    axial_i = np.interp(ctime, time, axial)
    vs_i = np.interp(ctime, time, vs)
    fig, ax = plt.subplots(figsize=(11, 3.4))

    def draw(frame):
        ax.clear()
        k = int(frame_ids[frame])
        pts = positions[k]
        x, y = (pts - pipe_center).dot(c), (pts - pipe_center).dot(n)
        ax.axhspan(radius, .9, color="#d8edf2", alpha=.55)
        ax.axhspan(-.9, -radius, color="#d8edf2", alpha=.55)
        ax.axhline(radius, color="#5996a5")
        ax.axhline(-radius, color="#5996a5")
        ax.scatter(x[::2], y[::2], c=colors[::2], s=7, linewidths=0)
        tail = [x[region == "TAIL"].mean(), y[region == "TAIL"].mean()]
        head = [x[region == "HEAD"].mean(), y[region == "HEAD"].mean()]
        ax.text(tail[0] - .06, tail[1], "TAIL", color=BLUE, weight="bold", ha="right", va="center")
        ax.text(head[0] + .06, head[1], "HEAD", color=RED, weight="bold", ha="left", va="center")
        contact = "BOTH" if active["TAIL"][k] and active["HEAD"][k] else (
            "TAIL" if active["TAIL"][k] else ("HEAD" if active["HEAD"][k] else "free"))
        elapsed = ctime[k]
        running_mean = (axial_i[k] - axial_i[0]) / elapsed if elapsed > 0 else 0.0
        cycle = min(5.0, elapsed / PERIOD)
        ax.text(.02, .93, "cycle={:.2f}/5 | contact={}".format(cycle, contact),
                transform=ax.transAxes, va="top", weight="bold", fontsize=9)
        ax.set(xlim=(-3, 3), ylim=(-.9, .9), aspect="equal",
               xlabel="canonical +s (mm), left to right", ylabel="n_routeA (mm)")
        ax.grid(axis="x", alpha=.15)
        ax.set_title(
            "{} | t={:.2f} ms | G=1.184 mT | U_flow=+0.5 mm/s\n"
            "delta_s={:+.4f} mm | mean={:+.2f} mm/s | v_s={:+.2f} mm/s | rock={:+.2f} deg".format(
                JOB, elapsed * 1e3, axial_i[k] - axial_i[0], running_mean, vs_i[k], rock_i[k]), fontsize=8)
        fig.tight_layout()

    FuncAnimation(fig, draw, frames=len(frame_ids), interval=1000 / 30).save(
        GIF, writer=PillowWriter(fps=30), dpi=100)
    draw(0)
    fig.savefig(OUT / "FAST_PRECOMP_F120_VMEAN5_5CYCLES_FRAME0.png", dpi=100)
    plt.close(fig)

    rows = []
    for item in comparison:
        rows.append("| {G:.3f} | {d:+.6f} | {v:+.4f} | {vf:+.3f} | {r0:+.2f}..{r1:+.2f} | {tc} | {hc} | {pf:.4f} | {reb} | {rt:.1f} |".format(
            G=item["G_mT"], d=item["delta_s_mm"], v=item["v_mean_mm_s"], vf=item["v_final_mm_s"],
            r0=item["rocking_min_deg"], r1=item["rocking_max_deg"], tc=item["TAIL_contact_episodes"],
            hc=item["HEAD_contact_episodes"], pf=item["peak_contact_resultant_N"], reb=item["obvious_rebound"], rt=item["runtime_s"]))
    report = """# FAST precomputed magnetic v_mean calibration

`SELECTED_G = 1.184 mT`

| G (mT) | delta_s (mm) | v_mean (mm/s) | v_final (mm/s) | rocking (deg) | TAIL contacts | HEAD contacts | peak contact (N) | rebound | runtime (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|:---:|---:|
{rows}

## Five-cycle verification

- Overall five-cycle mean: `{overall:+.6f} mm/s`
- Cycles 2-5 mean: `{steady:+.6f} mm/s`
- Per-cycle means: `{per_cycle}` mm/s
- Five-cycle target gate (4.5-5.5 mm/s): **FAILED**
- Rocking: `{rmin:+.4f} to {rmax:+.4f} deg`, dominant frequency `{freq:.4f} Hz`
- Contact safety hard gate: `{safety}`; no tumble, no gross penetration, no BOTH-wall bridge
- Rebound diagnostic: `{rebound}`; TAIL same-end recontact remains present
- Backend: `PRECOMPUTED_TABLE`; fluid: `FAST_SURROGATE_FLUID`; socket calls: `0`; CEL: `false`

The fitted two-cycle case matched the requested mean (`{short:+.6f} mm/s`, `{ratio:.3f}x` the background flow), but the five-cycle run continued accelerating. It therefore is the bounded calibration selection, not a validated steady 5 mm/s propulsion case. No additional refinement was run because the calibration limit allowed only one fitted follow-up.
""".format(rows="\n".join(rows), overall=final["overall_5cycle_mean_mm_s"], steady=cycles_2_5_mean,
           per_cycle=", ".join("{:+.3f}".format(value) for value in cycle_means),
           rmin=final["rocking_min_deg"], rmax=final["rocking_max_deg"],
           freq=final["dominant_rocking_frequency_Hz"], safety=final["contact_safety_passed"],
           rebound=final["obvious_rebound"], short=refinement["v_mean_mm_s"],
           ratio=refinement["v_mean_over_U_flow"])
    (OUT / "FastVMean5Calibration_Report.md").write_text(report, encoding="ascii")
    print(json.dumps({"selected_G_mT": 1.184, "five_cycle": final, "gif": str(GIF)}, indent=2))


if __name__ == "__main__":
    main()
