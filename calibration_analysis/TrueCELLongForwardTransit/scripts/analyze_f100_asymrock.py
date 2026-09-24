"""Compare the sole asymmetric F100 candidate with the symmetric parent."""
from __future__ import annotations

import csv
import json
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
from scipy.spatial.transform import Rotation

from analyze_f100_g2p20_final import geometry, kinematics, pose
from audit_f100_wall_exposure import exposure, summarize, ROOT, CASE, unit

NAME = "TRUECEL_B0P11_G2P20_F100_ASYMROCK_FAST"
NEW = ROOT / "case" / NAME
PERIOD = .01


def longest_negative(t, v):
    starts = np.flatnonzero((v < 0) & np.r_[True, v[:-1] >= 0])
    ends = np.flatnonzero((v < 0) & np.r_[v[1:] >= 0, True])
    return float(max((t[b] - t[a] for a, b in zip(starts, ends)), default=0))


def cycles(t, s, v):
    rows = []
    count = int(np.floor((t[-1] + 1e-8) / PERIOD))
    for k in range(1, count + 1):
        lo, hi = (k - 1) * PERIOD, k * PERIOD
        tt = np.r_[lo, t[(t > lo) & (t < hi)], hi]
        ss, vv = np.interp(tt, t, s), np.interp(tt, t, v)
        rows.append({"cycle": k, "delta_s_mm": float(ss[-1] - ss[0]),
                     "mean_v_s_mm_s": float((ss[-1] - ss[0]) / PERIOD),
                     "end_v_s_mm_s": float(vv[-1]), "min_v_s_mm_s": float(vv.min()),
                     "max_backtrack_mm": float(np.max(np.maximum.accumulate(ss) - ss)),
                     "longest_negative_v_ms": 1000 * longest_negative(tt, vv)})
    return rows


def actual_angle(ident, ur):
    c, n = unit(ident["canonical_plus_s_axis_aba"]), unit(ident["n_routeA_aba"])
    axis = Rotation.from_rotvec(ur).apply(np.broadcast_to(c, ur.shape))
    return np.rad2deg(np.arctan2(axis @ n, axis @ c))


def render_single(data, path):
    ident, t, u, ur, s, v = data
    rp, c, n, rel, colors, pipe = geometry(NEW, ident)
    angle = actual_angle(ident, ur)
    final = t[-1]
    frames = np.linspace(0, final, min(161, max(33, round(final / .03 * 97))))
    fig, ax = plt.subplots(figsize=(10, 4.4), constrained_layout=True)
    r = float(ident["lumen_radius_mm"])
    ax.axhline(r, color="#466d8a"); ax.axhline(-r, color="#466d8a")
    ax.set(xlim=tuple(ident["axial_CEL_extent_s_mm"]), ylim=(-.9, .9), aspect="equal",
           xlabel="canonical +s (mm), TAIL left -> HEAD right", ylabel="n (mm)", title=NAME)
    p = pose(0, t, u, ur, rp, rel)
    sc = ax.scatter(float(ident["s_start_mm"]) + (p - rp) @ c, (p - pipe) @ n, s=5, c=colors, linewidths=0)
    tx = ax.text(.99, .98, "", transform=ax.transAxes, ha="right", va="top", family="monospace",
                 bbox={"boxstyle": "round", "facecolor": "white", "alpha": .9})
    def update(i):
        ti = frames[i]
        p = pose(ti, t, u, ur, rp, rel)
        sc.set_offsets(np.c_[float(ident["s_start_mm"]) + (p - rp) @ c, (p - pipe) @ n])
        mask = t <= ti
        back = float(np.max(np.maximum.accumulate(s[mask]) - s[mask])) if mask.any() else 0
        tx.set_text(f"t={ti*1000:6.2f} ms  C{min(5, int(ti/PERIOD)+1)}  angle={np.interp(ti,t,angle):+.2f} deg\n"
                    f"delta_s={np.interp(ti,t,s)-s[0]:+.5f} mm  v_s={np.interp(ti,t,v):+.2f} mm/s\n"
                    f"MAX_BACKTRACK={back:.5f} mm")
        return sc, tx
    FuncAnimation(fig, update, frames=len(frames), interval=62).save(path, PillowWriter(fps=16), dpi=96)
    plt.close(fig)


def render_sync(new, old, path):
    fig, axes = plt.subplots(2, 1, figsize=(10, 7.3), constrained_layout=True)
    cases = ((CASE, old, "SYMMETRIC F100 / G2.20"), (NEW, new, "ASYMMETRIC F100 / G2.20"))
    artists = []
    final = min(new[1][-1], old[1][-1], .03)
    frames = np.linspace(0, final, 97)
    for ax, (case, data, label) in zip(axes, cases):
        ident, t, u, ur, s, v = data
        rp, c, n, rel, colors, pipe = geometry(case, ident)
        r = float(ident["lumen_radius_mm"])
        ax.axhline(r, color="#466d8a"); ax.axhline(-r, color="#466d8a")
        ax.set(xlim=tuple(ident["axial_CEL_extent_s_mm"]), ylim=(-.9, .9), aspect="equal",
               xlabel="canonical +s (mm), TAIL left -> HEAD right", ylabel="n (mm)", title=label)
        p = pose(0, t, u, ur, rp, rel)
        sc = ax.scatter(float(ident["s_start_mm"]) + (p - rp) @ c, (p - pipe) @ n, s=5, c=colors, linewidths=0)
        tx = ax.text(.99, .98, "", transform=ax.transAxes, ha="right", va="top", family="monospace",
                     bbox={"boxstyle": "round", "facecolor": "white", "alpha": .9})
        artists.append((sc, tx, ident, t, u, ur, s, v, actual_angle(ident, ur), rp, c, n, rel, pipe))
    def update(i):
        ti = frames[i]
        out = []
        for sc, tx, ident, t, u, ur, s, v, angle, rp, c, n, rel, pipe in artists:
            p = pose(ti, t, u, ur, rp, rel)
            sc.set_offsets(np.c_[float(ident["s_start_mm"]) + (p - rp) @ c, (p - pipe) @ n])
            mask = t <= ti
            back = float(np.max(np.maximum.accumulate(s[mask]) - s[mask]))
            tx.set_text(f"t={ti*1000:5.2f} ms  C{min(3,int(ti/PERIOD)+1)}  angle={np.interp(ti,t,angle):+.2f} deg\n"
                        f"delta_s={np.interp(ti,t,s)-s[0]:+.5f} mm  v_s={np.interp(ti,t,v):+.2f} mm/s\n"
                        f"MAX_BACKTRACK={back:.5f} mm")
            out.extend((sc, tx))
        return out
    FuncAnimation(fig, update, frames=len(frames), interval=62).save(path, PillowWriter(fps=16), dpi=96)
    plt.close(fig)


def main():
    new, old = kinematics(NEW, False), kinematics(CASE, True)
    ident, t, u, ur, s, v = new
    rows = cycles(t, s, v)
    with (ROOT / f"{NAME}_CYCLE_SUMMARY.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    gt, ga, gg = exposure(NEW, False)
    gap_rows = summarize(gt, ga, gg)
    with (ROOT / f"{NAME}_WALL_EXPOSURE_NATIVE.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(gap_rows[0]))
        writer.writeheader(); writer.writerows(gap_rows)
    back = np.maximum.accumulate(s) - s
    partial = None
    start = len(rows) * PERIOD
    if t[-1] > start + 1e-8:
        tt = np.r_[start, t[(t > start) & (t < t[-1])], t[-1]]
        ss, vv = np.interp(tt, t, s), np.interp(tt, t, v)
        partial = {"cycle": len(rows)+1, "complete_cycle": False, "start_ms": start*1000,
                   "end_ms": float(t[-1]*1000), "duration_ms": float((t[-1]-start)*1000),
                   "delta_s_mm": float(ss[-1]-ss[0]),
                   "mean_v_s_mm_s": float((ss[-1]-ss[0])/(t[-1]-start)),
                   "end_v_s_mm_s": float(vv[-1]), "min_v_s_mm_s": float(vv.min()),
                   "max_backtrack_mm": float(np.max(np.maximum.accumulate(ss)-ss)),
                   "longest_negative_v_ms": 1000*longest_negative(tt,vv)}
    bad = [row["cycle"] for row in rows if row["delta_s_mm"] <= 0]
    result = {"case": NAME, "duration_ms": float(t[-1] * 1000), "cycle_metrics": rows,
              "partial_cycle_metrics": partial, "classification": "ASYMROCK_RECOIL_FAIL",
              "follow_on_assessment": "SCALAR_AND_SIMPLE_WAVEFORM_TUNING_INSUFFICIENT",
              "max_backtrack_mm": float(back.max()), "mean_speed_mm_s": float((s[-1]-s[0])/t[-1]),
              "first_net_negative_cycle": bad[0] if bad else None,
              "first_negative_v_cycle": int(np.floor(t[np.flatnonzero(v<0)[0]]/PERIOD))+1 if np.any(v<0) else None,
              "wall_exposure": gap_rows,
              "three_cycle_preferred_gate": bool(len(rows)>=3 and all(row["delta_s_mm"]>0 for row in rows[:3]) and back.max()<.01)}
    (ROOT / f"{NAME}_METRICS.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    render_single(new, ROOT / f"{NAME}.gif")
    render_sync(new, old, ROOT / "F100_SYMMETRIC_VS_ASYMROCK_PHASE_SYNC.gif")
    with zipfile.ZipFile(ROOT / f"{NAME}_GIFS.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for filename in (f"{NAME}.gif", "F100_SYMMETRIC_VS_ASYMROCK_PHASE_SYNC.gif",
                         f"{NAME}_CYCLE_SUMMARY.csv", f"{NAME}_WALL_EXPOSURE_NATIVE.csv", f"{NAME}_METRICS.json"):
            archive.write(ROOT / filename, filename)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
