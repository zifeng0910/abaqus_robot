"""Validate a clean 0-30 ms run before evaluating its third cycle."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.spatial.transform import Rotation

from analyze_f100_candidate import longest_negative, unit
from analyze_f100_g2p20_final import geometry, pose, setup_axis


ROOT = Path(__file__).resolve().parents[1]
OLD_NAME = "TRUECEL_B0P11_G2P20_A14P5_F100_FAST"
NEW_NAME = "TRUECEL_B0P11_G2P20_A14P5_F100_CLEAN30"
OLD = ROOT / "case" / OLD_NAME
NEW = ROOT / "case" / NEW_NAME
PERIOD = .01
KEYS = tuple(p + str(i) for p in ("U", "UR", "V", "VR") for i in (1, 2, 3))


def load_rp(folder: Path, continuation: bool = False):
    with np.load(folder / "private" / "rp_history_private.npz") as z:
        h = {k: z[k].astype(float) for k in KEYS}
    if continuation:
        with np.load(folder / "F100_CYCLE3_RESTART_private" / "rp_history_private.npz") as z:
            offset = h["U1"][-1, 0]
            for k in KEYS:
                a = z[k].astype(float)
                a[:, 0] += offset
                h[k] = np.vstack((h[k], a[a[:, 0] > offset + 1e-12]))
    return h


def motion(h, ident):
    t = h["U1"][:, 0]
    u = np.column_stack([np.interp(t, h[f"U{i}"][:, 0], h[f"U{i}"][:, 1]) for i in (1, 2, 3)])
    ur = np.column_stack([np.interp(t, h[f"UR{i}"][:, 0], h[f"UR{i}"][:, 1]) for i in (1, 2, 3)])
    v = np.column_stack([np.interp(t, h[f"V{i}"][:, 0], h[f"V{i}"][:, 1]) for i in (1, 2, 3)])
    vr = np.column_stack([np.interp(t, h[f"VR{i}"][:, 0], h[f"VR{i}"][:, 1]) for i in (1, 2, 3)])
    c = unit(ident["canonical_plus_s_axis_aba"])
    n = unit(ident["n_routeA_aba"])
    rock = unit(ident["b_routeA_aba"])
    axis0 = unit(ident["head_tail_axis_aba"])
    body = Rotation.from_rotvec(ur).apply(np.broadcast_to(axis0, ur.shape))
    alpha = np.rad2deg(np.arctan2(body @ n, body @ c))
    return {"t": t, "u": u, "ur": ur, "s": u @ c, "v_s": v @ c,
            "rocking_angle": alpha, "omega_rock": vr @ rock}


def sample(m, key, t):
    return np.interp(t, m["t"], m[key])


def stats(a):
    d = np.abs(a)
    return {"rms": float(np.sqrt(np.mean(a*a))),
            "p99": float(np.percentile(d, 99)), "max": float(d.max())}


def finish_identity(ident, classification):
    ident.update({"status":"SOLVED", "dynamics_run_count":1,
                  "classification":classification, "cpus":1,
                  "precision":"BOTH"})
    (NEW / "case_identity.json").write_text(
        json.dumps(ident,indent=2)+"\n",encoding="ascii")


def compare_history(name, file_name, t):
    with np.load(OLD / "private" / file_name) as a, np.load(NEW / "private" / file_name) as b:
        keys = sorted(set(a.files) & set(b.files))
        if name == "whole_general_contact":
            keys = [k for k in keys if "ROBOT_SOLID" in k and
                    any(f"CFN{i} " in k for i in (1, 2, 3)) or
                    ("ROBOT_SOLID" in k and "CFNM " in k)]
        else:
            keys = [k for k in keys if k in ("ALLKE", "ALLPW", "ETOTAL")]
        return {k: stats(np.interp(t, b[k][:, 0], b[k][:, 1]) -
                         np.interp(t, a[k][:, 0], a[k][:, 1])) for k in keys}


def compare_magnetic(t, c):
    a = np.loadtxt(OLD / "magnetic_increment_g2p20_f100.csv", delimiter=",", skiprows=1)
    b = np.loadtxt(NEW / "magnetic_increment_g2p20_f100.csv", delimiter=",", skiprows=1)
    out = {}
    for col, key in enumerate(("Fmag1", "Fmag2", "Fmag3", "Mmag1", "Mmag2", "Mmag3"), 3):
        out[key] = stats(np.interp(t, b[:, 0], b[:, col]) -
                         np.interp(t, a[:, 0], a[:, col]))
    fa = np.column_stack([np.interp(t, a[:, 0], a[:, i]) for i in (3, 4, 5)])
    fb = np.column_stack([np.interp(t, b[:, 0], b[:, i]) for i in (3, 4, 5)])
    out["Fmag_s"] = stats((fb-fa) @ c)
    return out


def cycle_metrics(m):
    rows = []
    t, s, vs = m["t"], m["s"], m["v_s"]
    for cycle in (1, 2, 3):
        t0, t1 = (cycle-1)*PERIOD, cycle*PERIOD
        tc = np.r_[t0, t[(t > t0) & (t < t1)], t1]
        sc, vc = np.interp(tc, t, s), np.interp(tc, t, vs)
        back = np.maximum.accumulate(sc)-sc
        rows.append({"cycle": cycle, "t_start_s": t0, "t_end_s": t1,
                     "delta_s_mm": float(sc[-1]-sc[0]),
                     "mean_v_s_mm_s": float((sc[-1]-sc[0])/PERIOD),
                     "end_v_s_mm_s": float(vc[-1]),
                     "minimum_v_s_mm_s": float(vc.min()),
                     "maximum_backtrack_mm": float(back.max()),
                     "longest_negative_v_duration_ms": float(longest_negative(tc, vc)*1000)})
    return rows


def draw_motion(m, old, output):
    fig, axes = plt.subplots(4, 1, figsize=(10, 9), sharex=True, constrained_layout=True)
    for ax, key, ylabel in zip(axes,
                               ("s", "v_s", "rocking_angle", "omega_rock"),
                               ("s (mm)", "v_s (mm/s)", "rocking angle (deg)", "omega (rad/s)")):
        ax.plot(m["t"]*1000, m[key], color="#136f63", lw=1.1, label="clean 0-30 ms")
        mask = old["t"] >= .02
        ax.plot(old["t"][mask]*1000, old[key][mask], color="#b24c3b", lw=.9,
                alpha=.75, label="old restart C3")
        for x in (10, 20):
            ax.axvline(x, color="#999", lw=.7, ls="--")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=.15)
    axes[0].legend(loc="upper left")
    axes[-1].set_xlabel("absolute time (ms)")
    fig.savefig(output, dpi=160)
    plt.close(fig)


def render_full(m, ident, output, validated=True):
    rp0, c, n, rel, colors, pipe = geometry(NEW, ident)
    fig, ax = plt.subplots(figsize=(11, 4.3), constrained_layout=True)
    title = "F100 / G2.20 uninterrupted 0-30 ms"
    if not validated:
        title += " | 0-20 ms NOT REPRODUCIBLE"
    setup_axis(ax, ident, title)
    p0 = pose(0, m["t"], m["u"], m["ur"], rp0, rel)
    sc = ax.scatter(float(ident["s_start_mm"])+(p0-rp0)@c,
                    (p0-pipe)@n, s=5, c=colors, linewidths=0)
    tx = ax.text(.99, .99, "", transform=ax.transAxes, ha="right", va="top",
                 family="monospace", bbox={"facecolor":"white", "alpha":.92})
    frames = np.linspace(0, .03, 121)
    def update(i):
        ti = frames[i]
        p = pose(ti, m["t"], m["u"], m["ur"], rp0, rel)
        sc.set_offsets(np.c_[float(ident["s_start_mm"])+(p-rp0)@c, (p-pipe)@n])
        j = min(3, int(ti/PERIOD)+1)
        mask = m["t"] <= ti
        sv = m["s"][mask]
        back = float((np.maximum.accumulate(sv)-sv).max()) if len(sv) else 0.0
        ds = float(np.interp(ti,m["t"],m["s"])-np.interp((j-1)*PERIOD,m["t"],m["s"]))
        tx.set_text(f"t={ti*1000:5.2f} ms   cycle={j}\n"
                    f"delta_s={ds:+.5f} mm   v_s={sample(m,'v_s',ti):+.3f} mm/s\n"
                    f"MAX_BACKTRACK={back:.5f} mm")
        return sc, tx
    FuncAnimation(fig, update, frames=len(frames), interval=62).save(
        output, PillowWriter(fps=16), dpi=96)
    plt.close(fig)


def render_c3(old, clean, ident, output, validated=True):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.4), constrained_layout=True)
    frames = np.linspace(0, PERIOD, 65)
    artists = []
    for ax, folder, m, title in ((axes[0], OLD, old, "Old restart-derived C3"),
                                  (axes[1], NEW, clean, "Uninterrupted C3")):
        setup_axis(ax, ident, title)
        rp0,c,n,rel,colors,pipe = geometry(folder,ident)
        p0 = pose(.02,m["t"],m["u"],m["ur"],rp0,rel)
        sc = ax.scatter(float(ident["s_start_mm"])+(p0-rp0)@c,
                        (p0-pipe)@n,s=4,c=colors,linewidths=0)
        tx = ax.text(.99,.99,"",transform=ax.transAxes,ha="right",va="top",
                     family="monospace",bbox={"facecolor":"white","alpha":.92})
        artists.append((sc,tx,m,rp0,c,n,rel,pipe))
    if not validated:
        fig.suptitle("0-20 ms NOT REPRODUCIBLE | visual comparison only")
    def update(i):
        ti = .02+frames[i]
        out = []
        for sc,tx,m,rp0,c,n,rel,pipe in artists:
            p = pose(ti,m["t"],m["u"],m["ur"],rp0,rel)
            sc.set_offsets(np.c_[float(ident["s_start_mm"])+(p-rp0)@c,(p-pipe)@n])
            s0 = sample(m,"s",.02)
            ds = sample(m,"s",ti)-s0
            mask = (m["t"] >= .02) & (m["t"] <= ti)
            ss = np.r_[s0,m["s"][mask],sample(m,"s",ti)]
            back = float((np.maximum.accumulate(ss)-ss).max())
            tx.set_text(f"phase={frames[i]*1000:5.2f} ms\n"
                        f"delta_s={ds:+.5f} mm   v_s={sample(m,'v_s',ti):+.3f} mm/s\n"
                        f"MAX_BACKTRACK={back:.5f} mm")
            out += [sc,tx]
        return out
    FuncAnimation(fig,update,frames=len(frames),interval=62).save(
        output,PillowWriter(fps=16),dpi=96)
    plt.close(fig)


def main():
    ident = json.loads((NEW / "case_identity.json").read_text(encoding="utf-8"))
    old_ident = json.loads((OLD / "case_identity.json").read_text(encoding="utf-8"))
    if ident["magnetic_table_sha256"] != old_ident["magnetic_table_sha256"]:
        raise RuntimeError("Magnetic table identity differs")
    new = motion(load_rp(NEW), ident)
    old_clean = motion(load_rp(OLD), old_ident)
    if new["t"][0] > 1e-10 or new["t"][-1] < .03-1e-8:
        raise RuntimeError("CLEAN30 does not span 0-30 ms")
    t = old_clean["t"][old_clean["t"] <= .02]
    check = {key: stats(sample(new,key,t)-sample(old_clean,key,t))
             for key in ("s","v_s","rocking_angle","omega_rock")}
    first_floors = {"s":1e-6,"v_s":.01,"rocking_angle":.001,"omega_rock":.05}
    first = {}
    for key,floor in first_floors.items():
        delta = sample(new,key,t)-sample(old_clean,key,t)
        ix = np.flatnonzero(np.abs(delta)>floor)
        first[key] = ({"time_ms":float(t[ix[0]]*1000),"difference":float(delta[ix[0]]),
                       "threshold":floor} if len(ix) else None)
    check["magnetic"] = compare_magnetic(t, unit(ident["canonical_plus_s_axis_aba"]))
    check["whole_general_contact"] = compare_history(
        "whole_general_contact", "contact_history_private.npz", t)
    check["energy"] = compare_history("energy", "energy_history_private.npz", t)
    # Clean-start histories should clear the stricter trajectory gates used in
    # the restart audit. A load/contact mismatch is reported even if motion passes.
    gates = {"s_max_mm":check["s"]["max"] < .001,
             "v_s_p99_mm_s":check["v_s"]["p99"] < .1,
             "rocking_angle_max_deg":check["rocking_angle"]["max"] < .02,
             "omega_rock_p99_rad_s":check["omega_rock"]["p99"] < .5,
             "magnetic_force_max_N":max(check["magnetic"][k]["max"] for k in ("Fmag1","Fmag2","Fmag3")) < 1e-8,
             "magnetic_torque_max_Nmm":max(check["magnetic"][k]["max"] for k in ("Mmag1","Mmag2","Mmag3")) < 1e-5,
             "contact_p99_N":max(x["p99"] for x in check["whole_general_contact"].values()) < 1e-5}
    report = {"job":NEW_NAME,"baseline":OLD_NAME,
              "old_c3_evidence_status":"RESTART_CONTAMINATED_C3_RESULT",
              "clean20_comparison":check,"clean20_gates":gates,
              "clean20_first_divergence":first,
              "restart_write_times_ms":{"original_clean20":[5,10,15,20],
                                         "clean30":[7.5,15,22.5,30]},
              "clean20_reproducible":all(gates.values())}
    if not report["clean20_reproducible"]:
        report["classification"] = "CLEAN_RERUN_NOT_REPRODUCIBLE"
        old = motion(load_rp(OLD,True),old_ident)
        render_full(new,ident,ROOT/(NEW_NAME+".gif"),validated=False)
        render_c3(old,new,ident,ROOT/"F100_G2P20_RESTART_C3_VS_CLEAN_C3.gif",
                  validated=False)
        (ROOT / "F100_G2P20_CLEAN30_METRICS.json").write_text(json.dumps(report,indent=2)+"\n")
        lines=["# F100/G2.20 CLEAN30 restart sensitivity", "",
               "**CLEAN_RERUN_NOT_REPRODUCIBLE**", "",
               "CLEAN30 used one uninterrupted 0-30 ms Explicit step, B0=11 mT, "
               "G=2.20 mT, f=100 Hz, A_main=14.5 deg, A_cross=2.5 deg, "
               "c0=100000 mm/s, a 44x20x20 CEL mesh, scale factor 0.4, "
               "no mass scaling, one CPU and double precision. The original "
               "field interval of 8.333333 ms and per-increment RP/contact "
               "histories were retained.", "",
               "The 10-ms restart audit demonstrated non-reproducible "
               "contact/CEL continuation dynamics. The old C3 "
               "delta_s=-0.023527 mm came from a 20-30 ms restart "
               "and is retained as `RESTART_CONTAMINATED_C3_RESULT`. The new "
               "single-step 0-30 ms "
               "run completed, but its 0-20 ms trajectory failed the clean "
               "baseline reproduction gate. Its 20-30 ms history is recorded "
               "but C3 metrics and physical interpretation are withheld.", "",
               "The original 20 ms deck wrote restart states at 5/10/15/20 ms. "
               "Keeping `number interval=4` in the 30 ms deck moved those "
               "writes to 7.5/15/22.5/30 ms. This output-timing difference "
               "coincides with the first observed online departure near 5 ms; "
               "causality is not established by timing alone.", "",
               "| Clean 0-20 ms quantity (s: mm; v_s: mm/s; angle: deg; omega: rad/s) | RMS | p99 absolute | Max absolute | First threshold crossing (ms) |",
               "|---|---:|---:|---:|---:|"]
        for k in ("s","v_s","rocking_angle","omega_rock"):
            x=check[k];f=first[k]
            lines.append(f"| {k} | {x['rms']:.6g} | {x['p99']:.6g} | {x['max']:.6g} | "
                         f"{f['time_ms']:.6f} |" if f else
                         f"| {k} | {x['rms']:.6g} | {x['p99']:.6g} | {x['max']:.6g} | none |")
        lines.extend(["", "| Magnetic / whole General Contact quantity (force: N; torque: N mm) | RMS | p99 absolute | Max absolute |",
                      "|---|---:|---:|---:|"])
        for group in ("magnetic","whole_general_contact"):
            for k,x in check[group].items():
                label=k.split("|")[-1] if group=="whole_general_contact" else k
                lines.append(f"| {label} | {x['rms']:.6g} | {x['p99']:.6g} | {x['max']:.6g} |")
        lines.extend(["", "Magnetic loads depend on the evolving pose and position; "
                      "their later differences do not independently establish an "
                      "excitation mismatch. The General Contact values cover the robot surface's "
                      "whole contact domain, not an isolated wall-pair force. "
                      "The metrics JSON records each gate. The requested GIFs "
                      "show the recorded motion with a non-authoritative label; "
                      "no clean C3 cycle metrics or physical interpretation "
                      "were produced after this failed gate."])
        (ROOT / "F100_G2P20_CLEAN30_RESTART_SENSITIVITY_REPORT.md").write_text(
            "\n".join(lines)+"\n",encoding="utf-8")
        finish_identity(ident,report["classification"])
        print(json.dumps(report,indent=2))
        return
    old = motion(load_rp(OLD,True),old_ident)
    rows = cycle_metrics(new)
    old_rows = cycle_metrics(old)
    phase_t = np.linspace(.02, .03, 1001)
    c3_differences = {key: stats(sample(new,key,phase_t)-sample(old,key,phase_t))
                      for key in ("s","v_s","rocking_angle","omega_rock")}
    c3 = rows[2]
    if c3["delta_s_mm"] > 0:
        classification = ("CLEAN30_C3_FORWARD_SUSTAINED" if c3["maximum_backtrack_mm"] < .01
                          else "CLEAN30_C3_FORWARD_WITH_MINOR_RECOIL")
    else:
        classification = "CLEAN30_C3_RECOIL_CONFIRMED"
    report.update({"classification":classification,"clean_cycles":rows,
                   "old_restart_cycles":old_rows,
                   "cycle3_clean_minus_restart":c3_differences,
                   "old_c3_failure_physically_credible":classification=="CLEAN30_C3_RECOIL_CONFIRMED"})
    with (ROOT / "F100_G2P20_CLEAN30_CYCLE_METRICS.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    draw_motion(new,old,ROOT/"F100_G2P20_CLEAN30_C3_MOTION.png")
    render_full(new,ident,ROOT/(NEW_NAME+".gif"))
    render_c3(old,new,ident,ROOT/"F100_G2P20_RESTART_C3_VS_CLEAN_C3.gif")
    (ROOT / "F100_G2P20_CLEAN30_METRICS.json").write_text(json.dumps(report,indent=2)+"\n")
    lines=["# F100/G2.20 CLEAN30 restart sensitivity", "",
           f"**{classification}**", "",
           "CLEAN30 used one uninterrupted 0-30 ms Explicit step, B0=11 mT, "
           "G=2.20 mT, f=100 Hz, A_main=14.5 deg, A_cross=2.5 deg, "
           "c0=100000 mm/s, a 44x20x20 CEL mesh, scale factor 0.4, "
           "no mass scaling, one CPU and double precision. The original "
           "field interval of 8.333333 ms and per-increment RP/contact "
           "histories were retained.", "",
           "The 10-ms restart audit showed non-reproducible contact/CEL dynamics. "
           "The old Cycle-3 result was produced by a 20-30 ms restart and is "
           "labelled `RESTART_CONTAMINATED_C3_RESULT`; it is retained for comparison, "
           "not used as authoritative physical evidence.", "",
           "## Clean 0-20 ms validation", "",
           "The single-step CLEAN30 solve reproduces the original clean 0-20 ms "
           "baseline under all listed gates. Differences are CLEAN30 minus baseline.", "",
           "| Quantity | RMS | p99 absolute | Max absolute |",
           "|---|---:|---:|---:|"]
    for k in ("s","v_s","rocking_angle","omega_rock"):
        x=check[k];lines.append(f"| {k} | {x['rms']:.6g} | {x['p99']:.6g} | {x['max']:.6g} |")
    lines.extend(["", "| Load / contact quantity | RMS | p99 absolute | Max absolute |",
                  "|---|---:|---:|---:|"])
    for group in ("magnetic","whole_general_contact"):
        for k,x in check[group].items():
            label=k.split("|")[-1] if group=="whole_general_contact" else k
            lines.append(f"| {label} | {x['rms']:.6g} | {x['p99']:.6g} | {x['max']:.6g} |")
    lines.extend(["", "Contact entries are whole General Contact resultants on the "
                  "robot surface, not isolated robot-wall pair forces.", "",
                  "## Cycle metrics", "",
                  "| Source | Cycle | delta_s (mm) | mean v_s (mm/s) | end v_s (mm/s) | "
                  "min v_s (mm/s) | max backtrack (mm) | longest negative v (ms) |",
                  "|---|---:|---:|---:|---:|---:|---:|---:|"])
    for source,r in [("clean",x) for x in rows]+[("old restart",old_rows[2])]:
        lines.append(f"| {source} | {r['cycle']} | {r['delta_s_mm']:+.6f} | "
                     f"{r['mean_v_s_mm_s']:+.4f} | {r['end_v_s_mm_s']:+.4f} | "
                     f"{r['minimum_v_s_mm_s']:+.4f} | {r['maximum_backtrack_mm']:.6f} | "
                     f"{r['longest_negative_v_duration_ms']:.4f} |")
    lines.extend(["", "## C3 history difference", "",
                  "New uninterrupted C3 minus old restart C3 on 1,001 equal-phase "
                  "timestamps from 20 to 30 ms.", "",
                  "| Quantity | RMS | p99 absolute | Max absolute |",
                  "|---|---:|---:|---:|"])
    for k,x in c3_differences.items():
        lines.append(f"| {k} | {x['rms']:.6g} | {x['p99']:.6g} | {x['max']:.6g} |")
    if classification=="CLEAN30_C3_RECOIL_CONFIRMED":
        lines.extend(["", "The uninterrupted C3 also develops net recoil, so the "
                      "late-cycle failure remains physically credible without restart. "
                      "A dense-output clean diagnostic would require separate authorization."])
    else:
        lines.extend(["", "The uninterrupted C3 does not reproduce the old net "
                      "negative displacement. The old F100/G2.20 Cycle-3 failure "
                      "claim is withdrawn and labelled `OLD_C3_FAILURE_WAS_RESTART_SENSITIVE`."])
    (ROOT / "F100_G2P20_CLEAN30_RESTART_SENSITIVITY_REPORT.md").write_text(
        "\n".join(lines)+"\n",encoding="utf-8")
    finish_identity(ident,classification)
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
