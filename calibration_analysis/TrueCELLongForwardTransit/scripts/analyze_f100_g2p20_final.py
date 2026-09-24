"""Final three-cycle analysis and synchronized GIFs for the sole G=2.20 run."""
from __future__ import annotations

import csv
import json
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.spatial.transform import Rotation

from analyze_f100_candidate import read_robot_nodes, unit, longest_negative


OUT = Path(__file__).resolve().parents[1]
SEL_JOB = "TRUECEL_B0P11_G2P20_A14P5_F100_FAST"
BASE_JOB = "TRUECEL_B0P11_A14P5_F100_FAST"
SEL = OUT / "case" / SEL_JOB
BASE = OUT / "case" / BASE_JOB
PERIOD = 0.01


def load_history(case: Path, continuation: bool):
    z = np.load(case / "private" / "rp_history_private.npz")
    keys = ("U1", "U2", "U3", "UR1", "UR2", "UR3", "V1", "V2", "V3")
    data = {k: z[k].astype(float) for k in keys}
    z.close()
    if continuation:
        z = np.load(case / "F100_CYCLE3_RESTART_private" / "rp_history_private.npz")
        offset = data["U1"][-1, 0]
        for k in keys:
            a = z[k].astype(float).copy()
            a[:, 0] += offset
            data[k] = np.vstack((data[k], a[a[:, 0] > offset + 1e-12]))
        z.close()
    return data


def kinematics(case: Path, continuation: bool):
    ident = json.loads((case / "case_identity.json").read_text(encoding="utf-8"))
    h = load_history(case, continuation)
    t = h["U1"][:, 0]
    u = np.column_stack([h[f"U{i}"][:, 1] for i in (1, 2, 3)])
    ur = np.column_stack([h[f"UR{i}"][:, 1] for i in (1, 2, 3)])
    v = np.column_stack([h[f"V{i}"][:, 1] for i in (1, 2, 3)])
    c = unit(ident["canonical_plus_s_axis_aba"])
    s = u @ c
    vs = v @ c
    return ident, t, u, ur, s, vs


def cycle_metrics(t, s, vs):
    rows = []
    for k in range(1, 4):
        t0, t1 = (k - 1) * PERIOD, k * PERIOD
        tc = np.r_[t0, t[(t > t0) & (t < t1)], t1]
        sc, vc = np.interp(tc, t, s), np.interp(tc, t, vs)
        local_back = np.maximum.accumulate(sc) - sc
        rows.append({
            "cycle": k, "t_start_s": t0, "t_end_s": t1,
            "delta_s_mm": float(sc[-1] - sc[0]),
            "mean_v_s_mm_s": float((sc[-1] - sc[0]) / PERIOD),
            "end_v_s_mm_s": float(vc[-1]),
            "minimum_v_s_mm_s": float(vc.min()),
            "maximum_backward_excursion_mm": float(local_back.max()),
            "negative_velocity_duration_s": float(longest_negative(tc, vc)),
        })
    return rows


def geometry(case, ident):
    rp0 = np.asarray(ident["initial_center_aba_mm"], float)
    c = unit(ident["canonical_plus_s_axis_aba"])
    n = unit(ident["n_routeA_aba"])
    nodes = read_robot_nodes((case / f"{case.name}.inp").read_text(encoding="latin1"))
    rel = nodes - rp0
    axial = rel @ c
    colors = np.full(len(nodes), "#3e434b", dtype=object)
    colors[axial <= axial.min() + 0.25] = "#2864a8"
    colors[axial >= axial.max() - 0.25] = "#d1493f"
    pipe_center = rp0 - float(ident["s_start_mm"]) * c - float(ident["radial_offset_n_mm"]) * n
    return rp0, c, n, rel, colors, pipe_center


def pose(ti, t, u, ur, rp0, rel):
    ui = np.array([np.interp(ti, t, u[:, j]) for j in range(3)])
    uri = np.array([np.interp(ti, t, ur[:, j]) for j in range(3)])
    return rp0 + ui + Rotation.from_rotvec(uri).apply(rel)


def setup_axis(ax, ident, title):
    r = float(ident["lumen_radius_mm"])
    ax.axhline(r, color="#5d8794", lw=2); ax.axhline(-r, color="#5d8794", lw=2)
    ax.axvline(float(ident["s_start_mm"]), color="#777", ls="--", lw=1)
    ax.set(xlim=tuple(ident["axial_CEL_extent_s_mm"]), ylim=(-0.9, 0.9), aspect="equal",
           xlabel="canonical +s (mm), LEFT -> RIGHT", ylabel="n_routeA (mm)", title=title)


def render_selected(ident, t, u, ur, s, vs, path):
    rp0, c, n, rel, colors, pipe_center = geometry(SEL, ident)
    frames = np.linspace(0, 0.03, 97)
    fig, ax = plt.subplots(figsize=(11, 4.3), constrained_layout=True)
    setup_axis(ax, ident, "TRUE-CEL F100: selected G=2.20 mT")
    p0 = pose(0, t, u, ur, rp0, rel)
    scat = ax.scatter(float(ident["s_start_mm"]) + (p0-rp0)@c, (p0-pipe_center)@n,
                      s=5, c=colors, linewidths=0)
    txt = ax.text(.985, .985, "", transform=ax.transAxes, ha="right", va="top", family="monospace",
                  bbox={"boxstyle":"round", "facecolor":"white", "alpha":.9})
    def update(i):
        ti = frames[i]; p = pose(ti, t, u, ur, rp0, rel)
        scat.set_offsets(np.c_[float(ident["s_start_mm"]) + (p-rp0)@c, (p-pipe_center)@n])
        mask = t <= ti
        b = float((np.maximum.accumulate(s[mask])-s[mask]).max()) if mask.any() else 0.0
        txt.set_text(f"t={ti*1e3:6.2f} ms  cycle={min(3,int(ti/PERIOD)+1)}\n"
                     f"delta_s={np.interp(ti,t,s):+.5f} mm  v_s={np.interp(ti,t,vs):+.3f} mm/s\n"
                     f"MAX_BACKTRACK={b:.5f} mm")
        return scat, txt
    FuncAnimation(fig, update, frames=len(frames), interval=62).save(path, PillowWriter(fps=16), dpi=96)
    plt.close(fig)


def render_comparison(sel_data, base_data, path):
    si, st, su, sur, ss, sv = sel_data
    bi, bt, bu, bur, bs, bv = base_data
    sg = geometry(SEL, si); bg = geometry(BASE, bi)
    frames = np.linspace(0, 0.02, 65)
    fig, axes = plt.subplots(2, 1, figsize=(11, 7.2), constrained_layout=True)
    setup_axis(axes[0], bi, "Baseline G=2.00 mT")
    setup_axis(axes[1], si, "Selected G=2.20 mT")
    artists = []
    for ax, ident, tt, uu, rr, spos, vel, geo in (
        (axes[0], bi, bt, bu, bur, bs, bv, bg), (axes[1], si, st, su, sur, ss, sv, sg)):
        rp0,c,n,rel,colors,pc = geo; p = pose(0,tt,uu,rr,rp0,rel)
        sc = ax.scatter(float(ident["s_start_mm"])+(p-rp0)@c,(p-pc)@n,s=5,c=colors,linewidths=0)
        tx = ax.text(.985,.985,"",transform=ax.transAxes,ha="right",va="top",family="monospace",
                     bbox={"boxstyle":"round","facecolor":"white","alpha":.9})
        artists.append((sc,tx,ident,tt,uu,rr,spos,vel,geo))
    def update(i):
        ti=frames[i]; out=[]
        for sc,tx,ident,tt,uu,rr,spos,vel,geo in artists:
            rp0,c,n,rel,colors,pc=geo; p=pose(ti,tt,uu,rr,rp0,rel)
            sc.set_offsets(np.c_[float(ident["s_start_mm"])+(p-rp0)@c,(p-pc)@n])
            mask=tt<=ti; back=float((np.maximum.accumulate(spos[mask])-spos[mask]).max())
            tx.set_text(f"t={ti*1e3:6.2f} ms  delta_s={np.interp(ti,tt,spos):+.5f} mm\n"
                        f"v_s={np.interp(ti,tt,vel):+.3f} mm/s  back={back:.5f} mm")
            out += [sc,tx]
        return out
    fig.suptitle("Synchronized F100 comparison — identical time/phase scale")
    FuncAnimation(fig,update,frames=len(frames),interval=62).save(path,PillowWriter(fps=16),dpi=96)
    plt.close(fig)


def main():
    selected = kinematics(SEL, True); baseline = kinematics(BASE, False)
    ident,t,u,ur,s,vs = selected
    rows = cycle_metrics(t,s,vs)
    back = np.maximum.accumulate(s)-s
    classification = "TWO_CYCLE_IMPROVED_BUT_RESIDUAL_RECOIL"
    metrics = {
        "candidate": SEL_JOB, "classification": classification, "frequency_Hz": 100,
        "B0_mT": 11.0, "G_mT": 2.2, "total_time_s": float(t[-1]),
        "cycle_metrics": rows, "total_delta_s_mm": float(s[-1]-s[0]),
        "mean_speed_mm_s": float((s[-1]-s[0])/(t[-1]-t[0])),
        "max_backtrack_mm": float(back.max()), "minimum_v_s_mm_s": float(vs.min()),
        "three_cycle_gate_pass": bool(all(r["delta_s_mm"]>0 for r in rows) and back.max()<0.01),
        "stage1_two_cycle_gate_pass": bool(rows[0]["delta_s_mm"]>0 and rows[1]["delta_s_mm"]>0),
        "selection_basis": "G=2.15 first open-loop <0.005 mm plus one allowed +0.05 robustness increment",
    }
    csv_path=OUT/f"{SEL_JOB}_CYCLE_SUMMARY.csv"
    with csv_path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    gif=OUT/f"{SEL_JOB}.gif"; comp=OUT/"F100_G2P0_VS_G2P20.gif"
    render_selected(ident,t,u,ur,s,vs,gif); render_comparison(selected,baseline,comp)
    zip_path=OUT/f"{SEL_JOB}_GIFS.zip"
    with zipfile.ZipFile(zip_path,"w",zipfile.ZIP_DEFLATED) as z:
        z.write(gif,gif.name); z.write(comp,comp.name); z.write(csv_path,csv_path.name)
    metrics.update({"selected_gif":gif.name,"comparison_gif":comp.name,"gif_zip":zip_path.name})
    (OUT/f"{SEL_JOB}_METRICS.json").write_text(json.dumps(metrics,indent=2)+"\n",encoding="utf-8")
    report=(f"# F100 G=2.20 final\n\nClassification: **{classification}**\n\n"
            f"Cycles Δs: {rows[0]['delta_s_mm']:+.6f}, {rows[1]['delta_s_mm']:+.6f}, "
            f"{rows[2]['delta_s_mm']:+.6f} mm. Overall max backtrack: {back.max():.6f} mm. "
            f"The 20-ms gate passed, but Cycle 3 became negative and the recoil exceeded 0.010 mm.\n")
    (OUT/f"{SEL_JOB}_REPORT.md").write_text(report,encoding="utf-8")
    print(json.dumps(metrics,indent=2))


if __name__ == "__main__":
    main()
