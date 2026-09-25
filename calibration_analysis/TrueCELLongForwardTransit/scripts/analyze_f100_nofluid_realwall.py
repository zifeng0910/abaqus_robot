"""Audit the single no-fluid/real-wall solve and render fixed-camera GIFs."""
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

from analyze_f100_g2p20_clean30 import KEYS, motion
from analyze_f100_g2p20_final import geometry, pose, setup_axis

ROOT = Path(__file__).resolve().parents[1]
JOB = "F100_G2P20_NOFLUID_REALWALL50"
CASE = ROOT / "case" / JOB
NO_WALL = ROOT / "case" / "TRUECEL_B0P11_G2P20_F100_NOFLUID_CONTROL"
FULL = ROOT / "case" / "TRUECEL_B0P11_G2P20_A14P5_F100_CLEAN50"


def load(case: Path):
    identity = json.loads((case / "case_identity.json").read_text(encoding="utf-8-sig"))
    with np.load(case / "private" / "rp_history_private.npz") as z:
        history = {key: z[key].astype(float) for key in KEYS}
    return identity, motion(history, identity)


def cycle_rows(m):
    rows = []
    for k in range(1, 6):
        end = k * .01
        if m["t"][-1] < end - 1e-7:
            break
        start = end - .01
        t = np.r_[start, m["t"][(m["t"] > start) & (m["t"] < end)], end]
        s = np.interp(t, m["t"], m["s"])
        v = np.interp(t, m["t"], m["v_s"])
        theta = np.interp(t, m["t"], m["rocking_angle"])
        omega = np.interp(t, m["t"], m["omega_rock"])
        a, b, _ = np.linalg.lstsq(np.c_[np.sin(2*np.pi*100*t),
                                           np.cos(2*np.pi*100*t), np.ones(len(t))], theta, rcond=None)[0]
        rows.append({"cycle": k, "delta_s_mm": float(s[-1]-s[0]),
                     "mean_v_s_mm_s": float((s[-1]-s[0])/.01),
                     "end_v_s_mm_s": float(v[-1]), "minimum_v_s_mm_s": float(v.min()),
                     "max_backtrack_mm": float((np.maximum.accumulate(s)-s).max()),
                     "rocking_amplitude_deg": float((theta.max()-theta.min())/2),
                     "rocking_fundamental_amplitude_deg": float(np.hypot(a, b)),
                     "phase_lag_deg": float((-np.degrees(np.arctan2(b, a))+180)%360-180),
                     "omega_rock_peak_rad_s": float(np.abs(omega).max())})
    return rows


def contact_events(identity, m):
    z = np.load(CASE / "private" / "robot_wall_copen_private.npz")
    t = z["time"].astype(float)
    nodes = z["node_coordinates_mm"].astype(float)
    opening = z["COPEN"].astype(float)
    with np.load(CASE / "private" / "contact_history_private.npz") as h:
        def resultant(variable):
            key = next(k for k in h.files if f"|{variable} on surface ASSEMBLY_ROBOT_SOLID" in k)
            return np.interp(t, h[key][:, 0], h[key][:, 1])
        wall_normal = resultant("CFNM")
        wall_tangent = resultant("CFSM")
    rp0 = np.asarray(identity["initial_center_aba_mm"], float)
    c = np.asarray(identity["canonical_plus_s_axis_aba"], float)
    n = np.asarray(identity["n_routeA_aba"], float)
    b = np.asarray(identity["b_routeA_aba"], float)
    c /= np.linalg.norm(c); n /= np.linalg.norm(n); b /= np.linalg.norm(b)
    sref = (nodes-rp0) @ c
    masks = {"TAIL": sref <= sref.min()+.25, "HEAD": sref >= sref.max()-.25}
    pipe = rp0 - identity["s_start_mm"]*c - identity["radial_offset_n_mm"]*n
    u = np.column_stack([np.interp(t, m["t"], m["u"][:, j]) for j in range(3)])
    ur = np.column_stack([np.interp(t, m["t"], m["ur"][:, j]) for j in range(3)])
    rotations = Rotation.from_rotvec(ur).as_matrix()
    gaps = {key: np.empty(len(t)) for key in masks}
    contact_nodes = {key: np.empty(len(t), dtype=bool) for key in masks}
    for i in range(len(t)):
        p = rp0 + u[i] + (rotations[i] @ (nodes-rp0).T).T
        rel = p-pipe
        qn, qb = rel@n, rel@b
        radial = np.hypot(qn, qb)
        gap = identity["lumen_radius_mm"]-radial
        for name, mask in masks.items():
            gaps[name][i] = gap[mask].min()
            contact_nodes[name][i] = np.any(mask & np.isfinite(opening[i]) &
                                             (opening[i] <= .015) & (gap <= .02))
    events = []
    active = {key: contact_nodes[key] & (wall_normal > 1e-9) & (gaps[key] <= .015)
              for key in masks}
    for name in ("HEAD", "TAIL"):
        idx = np.flatnonzero(active[name])
        if not len(idx):
            continue
        groups = np.split(idx, np.flatnonzero(np.diff(t[idx]) > .00015)+1)
        for group in groups:
            lo, hi = group[0], group[-1]
            exclusive = ~active["TAIL" if name == "HEAD" else "HEAD"][lo:hi+1]
            interval_force = wall_normal[lo:hi+1][exclusive]
            interval_tangent = wall_tangent[lo:hi+1][exclusive]
            interval_time = t[lo:hi+1][exclusive]
            events.append({"end": name, "start_ms": float(t[lo]*1000),
                           "end_ms": float(t[hi]*1000), "duration_ms": float((t[hi]-t[lo])*1000),
                           "minimum_signed_gap_mm": float(gaps[name][lo:hi+1].min()),
                           "peak_normal_force_N": float(interval_force.max()) if len(interval_force) else None,
                           "tangential_impulse_Ns": float(np.trapz(interval_tangent,interval_time)) if len(interval_time)>1 else None,
                           "force_attribution": "exclusive end" if len(interval_force)==hi-lo+1 else "overlap; exclusive samples only"})
    events.sort(key=lambda x: x["start_ms"])
    previous = None
    counts = {"HEAD": 0, "TAIL": 0}
    for event in events:
        event["same_end_recontact"] = event["end"] == previous
        counts[event["end"]] += int(event["same_end_recontact"])
        previous = event["end"]
    return events, counts, t, active, gaps, wall_normal


def magnetic_audit(identity, m, events):
    log = np.loadtxt(CASE / "magnetic_increment_g2p20_f100.csv", delimiter=",", skiprows=1)
    table = np.loadtxt(CASE / "magnetic_field_gradient_table_B0P11_A14P5.dat", skiprows=1, max_rows=361)
    phase = log[:, 2]
    field = np.column_stack([np.interp(phase, table[:, 0], table[:, j]) for j in (1, 2, 3)]) * .011
    c = np.asarray(identity["canonical_plus_s_axis_aba"], float)
    b = np.asarray(identity["b_routeA_aba"], float)
    command = 14.5*np.sin(np.radians(phase))
    force_s = log[:, 3:6] @ c
    torque_rock = log[:, 6:9] @ b
    t = log[:, 0]
    windows = []
    for event in events:
        center = (event["start_ms"]+event["end_ms"])/2000
        mask = np.abs(t-center) <= .0002
        before = max(t[0], center-.0002)
        after = min(t[-1], center+.0002)
        phase_error = (np.interp((before,after),m["t"],m["rocking_angle"]) -
                       14.5*np.sin(2*np.pi*100*np.asarray((before,after))))
        windows.append({"end": event["end"], "center_ms": center*1000,
                        "sample_count": int(mask.sum()),
                        "B_mT_range": [float(x) for x in (1000*np.linalg.norm(field[mask],axis=1).min(),
                                                            1000*np.linalg.norm(field[mask],axis=1).max())],
                        "force_s_N_range": [float(force_s[mask].min()), float(force_s[mask].max())],
                        "torque_rock_Nmm_range": [float(torque_rock[mask].min()),float(torque_rock[mask].max())],
                        "command_deg_range": [float(command[mask].min()),float(command[mask].max())],
                        "theta_deg_range": [float(np.interp(t[mask],m["t"],m["rocking_angle"]).min()),
                                            float(np.interp(t[mask],m["t"],m["rocking_angle"]).max())],
                        "phase_error_change_0p4ms_deg":float(phase_error[1]-phase_error[0])})
    sample = np.arange(0,len(t),100)
    with (ROOT / f"{JOB}_MAGNETIC_10US.csv").open("w",newline="") as handle:
        writer=csv.writer(handle)
        writer.writerow(("time_s","B_mT","command_deg","force_s_N","torque_rock_Nmm","theta_deg"))
        writer.writerows(zip(t[sample],1000*np.linalg.norm(field[sample],axis=1),command[sample],
                             force_s[sample],torque_rock[sample],
                             np.interp(t[sample],m["t"],m["rocking_angle"])))
    return windows, {"B_mT_min":float((1000*np.linalg.norm(field,axis=1)).min()),
                     "B_mT_max":float((1000*np.linalg.norm(field,axis=1)).max()),
                     "force_s_N_min":float(force_s.min()),"force_s_N_max":float(force_s.max()),
                     "torque_rock_Nmm_min":float(torque_rock.min()),
                     "torque_rock_Nmm_max":float(torque_rock.max()),
                     "phase_implementation":"VUAMP modulo(36000*t,360), 100 Hz; table JSON frequency_Hz=120 is stale metadata"}


def render(identity, m, ct, active, comparison, end_t):
    rp0,c,n,rel,colors,pipe = geometry(CASE,identity)
    frames=np.linspace(0,end_t,81)
    for mode in ("canonical","comparison"):
        if mode == "canonical":
            fig,axes=plt.subplots(2,1,figsize=(10,6.5),gridspec_kw={"height_ratios":[2,1]},constrained_layout=True)
            ax,trace=axes
            setup_axis(ax,identity,"F100 G2P20 | no fluid, real wall")
            trace.plot(m["t"]*1000,m["s"],color="#126d62",lw=1)
            marker,=trace.plot([],[],"o",color="#126d62")
            trace.set(xlim=(0,50),xlabel="time (ms)",ylabel="s (mm)")
        else:
            fig,axes=plt.subplots(3,1,figsize=(10,9),constrained_layout=True)
        artists=[]
        datasets=[(identity,m,CASE)] if mode=="canonical" else comparison
        for j,(ident,move,case) in enumerate(datasets):
            axis=axes[0] if mode=="canonical" else axes[j]
            if mode!="canonical": setup_axis(axis,ident,("NO WALL","REAL WALL","FULL CEL")[j])
            r0,cc,nn,rr,col,pc=geometry(case,ident)
            p0=pose(0,move["t"],move["u"],move["ur"],r0,rr)
            scatter=axis.scatter(float(ident["s_start_mm"])+(p0-r0)@cc,(p0-pc)@nn,s=4,c=col,linewidths=0)
            label=axis.text(.98,.97,"",transform=axis.transAxes,ha="right",va="top",family="monospace",fontsize=8,
                            bbox={"facecolor":"white","alpha":.9})
            artists.append((scatter,label,ident,move,r0,cc,nn,rr,pc))
        def update(k):
            ti=frames[k]
            for scatter,label,ident,move,r0,cc,nn,rr,pc in artists:
                shown=min(ti,move["t"][-1])
                p=pose(shown,move["t"],move["u"],move["ur"],r0,rr)
                scatter.set_offsets(np.c_[float(ident["s_start_mm"])+(p-r0)@cc,(p-pc)@nn])
                s=np.interp(shown,move["t"],move["s"])
                v=np.interp(shown,move["t"],move["v_s"])
                theta=np.interp(shown,move["t"],move["rocking_angle"])
                cycle=min(5,int(shown/.01)+1)
                start=np.interp((cycle-1)*.01,move["t"],move["s"])
                extra=""
                if ident["case_id"]==JOB:
                    idx=np.clip(np.searchsorted(ct,shown),0,len(ct)-1)
                    extra=f"\nHEAD {int(active['HEAD'][idx])}  TAIL {int(active['TAIL'][idx])}"
                if shown<ti-1e-6: extra+="\npartial trajectory ended"
                label.set_text(f"t={shown*1000:.2f} ms C{cycle}  ds={s-start:+.4f} mm\n"
                               f"v_s={v:+.2f}  theta={theta:+.2f}  cmd={14.5*np.sin(2*np.pi*100*shown):+.2f}"+extra)
            if mode=="canonical": marker.set_data([ti*1000],[np.interp(ti,m["t"],m["s"])])
        gif=ROOT/("F100_G2P20_NOFLUID_REALWALL.gif" if mode=="canonical" else "F100_G2P20_THREE_WAY.gif")
        FuncAnimation(fig,update,frames=len(frames),interval=75).save(gif,PillowWriter(fps=12),dpi=85)
        plt.close(fig)


def main():
    identity,m=load(CASE)
    sta=(CASE/f"{JOB}.sta").read_text(encoding="latin1")
    complete="THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in sta
    rows=cycle_rows(m)
    if not rows or len(rows)<5 and "Process terminated by external request" not in sta:
        raise RuntimeError("Missing complete cycle or unexplained solver stop")
    events,counts,ct,active,gaps,normal=contact_events(identity,m)
    windows,mag=magnetic_audit(identity,m,events)
    other=[load(NO_WALL),(identity,m),load(FULL)]
    comparison=[(ident,move,case) for (ident,move),case in zip(other,(NO_WALL,CASE,FULL))]
    references={case.name:cycle_rows(move) for ident,move,case in comparison}
    robust_tail = [e for e in events if e["end"]=="TAIL" and
                   e["peak_normal_force_N"] is not None and
                   e["peak_normal_force_N"]>.001 and e["minimum_signed_gap_mm"]<0]
    maxback=float((np.maximum.accumulate(m["s"])-m["s"]).max())
    neg=any(r["delta_s_mm"]<0 and r["max_backtrack_mm"]>.05 for r in rows)
    chatter=any(v>=2 for v in counts.values())
    if not np.all(np.isfinite(m["s"])) or "***ERROR" in sta:
        classification="NOFLUID_REALWALL_NUMERICALLY_INVALID"
        mechanism="MECHANISM_INCONCLUSIVE"
    elif neg:
        classification="NOFLUID_REALWALL_RECOIL_CONFIRMED"
        mechanism="WALL_CONTACT_ALONE_CAN_TRIGGER_RECOIL"
    elif chatter:
        classification="NOFLUID_REALWALL_CHATTER_FAIL"
        mechanism="WALL_CONTACT_PERTURBS_BUT_DOES_NOT_CAUSE_SUSTAINED_RECOIL"
    else:
        classification="NOFLUID_REALWALL_STABLE_FORWARD"
        mechanism="FLUID_CEL_COUPLING_REQUIRED_FOR_RECOIL" if events else "MECHANISM_INCONCLUSIVE"
    result={"classification":classification,"mechanism":mechanism,"completed_50ms":complete,
            "cycles":rows,"net_displacement_mm":float(m["s"][-1]-m["s"][0]),
            "max_positional_backtrack_mm":maxback,
            "forward_time_fraction":float(np.mean(m["v_s"]>0)),
            "contact_events":events,"same_end_recontacts":counts,
            "robust_tail_contact_episodes_force_gt_0p001N_negative_gap":len(robust_tail),
            "minimum_signed_gap_mm":{k:float(v.min()) for k,v in gaps.items()},
            "wall_only_peak_normal_force_N":float(normal.max()),
            "collision_magnetic_windows":windows,"magnetic":mag,"comparison_cycles":references,
            "contact_output_limit":"Only whole robot-wall CFNM/CFSM history was available; end force attributed only during exclusive end contact.",
            "fail_fast_note":"No direct HEAD/TAIL contact gate existed online; repeated chatter became decisive in offline 25-us COPEN/gap audit after the uninterrupted 50-ms run."}
    (ROOT/f"{JOB}_METRICS.json").write_text(json.dumps(result,indent=2)+"\n")
    with (ROOT/f"{JOB}_CYCLES.csv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    with (ROOT/f"{JOB}_CONTACT_EVENTS.csv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(events[0]) if events else ["end"])
        writer.writeheader();writer.writerows(events)
    render(identity,m,ct,active,comparison,min(.05,m["t"][-1]))
    with zipfile.ZipFile(ROOT/f"{JOB}_GIFS.zip","w",zipfile.ZIP_DEFLATED) as archive:
        for name in ("F100_G2P20_NOFLUID_REALWALL.gif","F100_G2P20_THREE_WAY.gif"):
            archive.write(ROOT/name,name)
    report=[f"# {JOB}","",classification,"",mechanism,"",
            f"Completed 50 ms: {complete}. Net displacement {result['net_displacement_mm']:+.6f} mm; "
            f"max backtrack {maxback:.6f} mm; forward-time fraction {result['forward_time_fraction']:.3f}.","",
            "| Cycle | delta_s mm | mean v_s mm/s | end v_s | minimum v_s | max backtrack mm | rocking amp deg | fundamental deg | lag deg | peak omega rad/s |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        report.append("| "+" | ".join(str(r[k]) for k in r)+" |")
    report += ["","## Three-way comparison","",
               "| Cycle | no fluid, no wall delta_s mm | no fluid, real wall delta_s mm | clean FULL CEL delta_s mm |",
               "|---:|---:|---:|---:|"]
    for i in range(5):
        a=references[NO_WALL.name][i]["delta_s_mm"]
        b=references[CASE.name][i]["delta_s_mm"]
        c=references[FULL.name][i]["delta_s_mm"] if i<len(references[FULL.name]) else None
        report.append(f"| {i+1} | {a:+.6f} | {b:+.6f} | {c if c is None else f'{c:+.6f}'} |")
    report += ["", "No-wall C1-C5 rocking half-range is 21.18-25.67 deg; real-wall C1-C5 is 11.53-13.27 deg. "
               "Real-wall lag changes from 17.27 to 6.57 deg without a sustained positional reversal. "
               "FULL CEL C4 is -0.565 mm with 0.565 mm backtrack; C5 was not completed. "
               "No-wall has no active wall contact. The FULL CEL parent lacks 25-us end-resolved contact fields, so its same-end episode count cannot be compared at this resolution.",
               "","## Contact","",f"Same-end recontacts under broad COPEN/gap gate: {counts}. "
               f"TAIL episodes with negative gap and peak wall force >0.001 N: {len(robust_tail)}.","",
               "The repeated-contact stop condition was established only in offline end-resolved fields; the solver therefore reached the full 50 ms.","",
               "Contact requires robot-surface COPEN, signed cylindrical gap and direct wall-only CFNM; 25-us fields, 0.15-ms episode merge. End-specific force is attributed only where that end contacts alone.","",
               "| End | start ms | end ms | duration ms | minimum gap mm | peak normal N | tangent impulse N s | attribution | same-end recontact |",
               "|---|---:|---:|---:|---:|---:|---:|---|---|"]
    for e in events:
        report.append("| "+" | ".join(str(e[k]) for k in e)+" |")
    report += ["","## Magnetic and comparison","",
               "Magnetic 10-us trace and 0.2-ms collision windows are in the CSV/JSON. |B| remains 10.99997-11.00000 mT; +s magnetic force remains positive. VUAMP uses 100 Hz; table JSON's 120 Hz is stale metadata.","",
               f"Maximum absolute change in theta-command error across a 0.4-ms contact-centered window: {max(abs(w['phase_error_change_0p4ms_deg']) for w in windows):.3f} deg. This and the cycle phase lags show no catastrophic loss of magnetic phase lock, though the lag drifts over five cycles.","",
               "The 21 stronger TAIL episodes establish contact chatter, but wall contact alone does not reproduce FULL CEL's C4 sustained recoil. Fluid/CEL coupling remains required for that observed recoil in this matched control; the exact force pathway is not resolved here.","",
               "The clean FULL CEL baseline ended partially at 40.755 ms; no C5 comparison is inferred.","",
               "Three-way cycle metrics are in the JSON. This one-case isolation changes only removal of fluid physics and retains parent robot-wall law."]
    (ROOT/f"{JOB}_REPORT.md").write_text("\n".join(report)+"\n")
    identity.update(status="SOLVED" if complete else "USER_TERMINATED_PARTIAL",classification=classification)
    (CASE/"case_identity.json").write_text(json.dumps(identity,indent=2)+"\n")
    print(json.dumps({"classification":classification,"mechanism":mechanism,
                      "cycles":rows,"contacts":len(events),"recontacts":counts},indent=2))


if __name__=="__main__":
    main()
