"""Analyze forces only after the fixed V/VR CEL replay passes its kinematic gate."""
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
SOURCE = ROOT / "case" / "F100_G2P20_NOFLUID_REALWALL50"
JOB = "F100_G2P20_REALWALL_FIXEDVVR_FULLCEL40"
CASE = ROOT / "case" / JOB


def history_vector(z, prefix, t):
    return np.column_stack([np.interp(t, z[f"{prefix}{i}"][:, 0], z[f"{prefix}{i}"][:, 1])
                            for i in (1, 2, 3)])


def contact_vector(z, kind, pair, t):
    values = []
    for axis in (1, 2, 3):
        matches = [key for key in z.files if f"|{kind}{axis} on surface ASSEMBLY_ROBOT_SOLID-1_ROBOT_SOLID_SURF" in key
                   and ("/ASSEMBLY_PIPE_WALL_HELPER-1_PIPE_WALL_HELPER_SURF" in key) == pair]
        if len(matches) != 1:
            raise RuntimeError(f"Expected one {kind}{axis} pair={pair}, found {matches}")
        data = z[matches[0]]
        values.append(np.interp(t, data[:, 0], data[:, 1]))
    return np.column_stack(values)


def integrate(t, y):
    return np.r_[0., np.cumsum((y[1:]+y[:-1])*.5*np.diff(t))]


def window_impulse(t, y, lo, hi):
    mask = (t > lo) & (t < hi)
    x = np.r_[lo, t[mask], hi]
    v = np.r_[np.interp(lo, t, y), y[mask], np.interp(hi, t, y)]
    return float(np.trapz(v, x))


def events(identity, source_motion, replay_motion, t, source_normal, replay_normal,
           source_wall_s, replay_wall_s, source_tangent, replay_tangent):
    with np.load(SOURCE / "private" / "robot_wall_copen_private.npz") as z:
        frame_t = z["time"].astype(float)
        nodes = z["node_coordinates_mm"].astype(float)
        source_open = z["COPEN"].astype(float)
    frame_t = frame_t[frame_t <= .04 + 1e-8]
    source_open = source_open[:len(frame_t)]
    rp0 = np.asarray(identity["initial_center_aba_mm"], float)
    c = np.asarray(identity["canonical_plus_s_axis_aba"], float)
    n = np.asarray(identity["n_routeA_aba"], float)
    b = np.asarray(identity["b_routeA_aba"], float)
    axial = (nodes-rp0)@c
    masks = {"HEAD":axial>=axial.max()-.25, "TAIL":axial<=axial.min()+.25}
    pipe = rp0 - identity["s_start_mm"]*c - identity["radial_offset_n_mm"]*n
    def gap_series(m):
        u = np.column_stack([np.interp(frame_t,m["t"],m["u"][:,j]) for j in range(3)])
        ur = np.column_stack([np.interp(frame_t,m["t"],m["ur"][:,j]) for j in range(3)])
        rot = Rotation.from_rotvec(ur).as_matrix()
        gaps = {name:np.empty(len(frame_t)) for name in masks}
        for i in range(len(frame_t)):
            xyz = rp0+u[i]+(rot[i]@(nodes-rp0).T).T
            rel = xyz-pipe
            rad = np.hypot(rel@n,rel@b)
            gap = identity["lumen_radius_mm"]-rad
            for name, mask in masks.items():
                gaps[name][i] = gap[mask].min()
        return gaps

    source_gaps = gap_series(source_motion)
    replay_gaps = gap_series(replay_motion)
    rows = []
    active_by_case = {}
    for label,gaps,normal,wall_s,tangent in (
            ("NO_FLUID_REAL_WALL",source_gaps,source_normal,source_wall_s,source_tangent),
            ("FULL_CEL_FIXED_REPLAY",replay_gaps,replay_normal,replay_wall_s,replay_tangent)):
        sampled = np.interp(frame_t,t,normal)
        active = {}
        for name,mask in masks.items():
            if label=="NO_FLUID_REAL_WALL":
                near = np.any(mask[None,:] & np.isfinite(source_open) &
                              (source_open<=.015),axis=1)
            else:
                # CEL COPEN can include fluid contact, so wall force plus signed gap are decisive.
                near = np.ones(len(frame_t),dtype=bool)
            active[name] = (gaps[name]<=.015)&near&(sampled>1e-9)
        active_by_case[label] = active
        previous = None
        for name in ("HEAD","TAIL"):
            idx = np.flatnonzero(active[name])
            if not len(idx):
                continue
            groups = np.split(idx,np.flatnonzero(np.diff(frame_t[idx])>.00015)+1)
            for group in groups:
                a,bidx = group[0],group[-1]
                other = "TAIL" if name=="HEAD" else "HEAD"
                exclusive = ~active[other][a:bidx+1]
                lo = max(0.,frame_t[a]-.0000125)
                hi = min(.04,frame_t[bidx]+.0000125)
                overlaps = not bool(np.all(exclusive))
                event = {"case":label,"end":name,"start_ms":float(frame_t[a]*1000),
                         "end_ms":float(frame_t[bidx]*1000),
                         "duration_ms_resolution_limited":float((hi-lo)*1000),
                         "minimum_signed_gap_mm":float(gaps[name][a:bidx+1].min()),
                         "normal_resultant_peak_N":float(sampled[a:bidx+1][exclusive].max()) if np.any(exclusive) else None,
                         "tangential_resultant_peak_N":float(np.interp(frame_t[a:bidx+1][exclusive],t,tangent).max()) if np.any(exclusive) else None,
                         "tangential_impulse_Ns":window_impulse(t,tangent,lo,hi) if not overlaps else None,
                         "wall_axial_impulse_Ns":window_impulse(t,wall_s,lo,hi) if not overlaps else None,
                         "force_attribution":"exclusive end" if not overlaps else "overlapping ends; force partition unresolved"}
                rows.append(event)
        subset = sorted((r for r in rows if r["case"]==label),key=lambda r:r["start_ms"])
        for event in subset:
            event["same_end_recontact"] = event["end"]==previous
            previous=event["end"]
    return rows, source_gaps, replay_gaps, active_by_case, frame_t


def plot_forces(t, series, impulse):
    labels = {"magnetic":"Magnetic applied", "wall":"Robot-wall direct",
              "fluid":"INFERRED ROBOT-CEL", "reaction":"Prescribed-motion reaction"}
    colors = {"magnetic":"#257c67","wall":"#bd563f","fluid":"#426ea8","reaction":"#37383d"}
    sample = np.unique(np.r_[np.arange(0,len(t),max(1,len(t)//4000)),len(t)-1])
    for name,data,ylabel in (("AXIAL_FORCE",series,"axial force (N)"),
                              ("CUMULATIVE_IMPULSE",impulse,"cumulative axial impulse (N s)")):
        fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
        for key in labels:
            ax.plot(t[sample]*1000,data[key][sample],label=labels[key],color=colors[key],lw=1)
        for boundary in (10,20,30,40): ax.axvline(boundary,color="#999",ls="--",lw=.6)
        ax.axhline(0,color="#666",lw=.7)
        ax.set(xlim=(0,40),xlabel="time (ms)",ylabel=ylabel,title=f"F100 fixed-trajectory FULL CEL | {name.replace('_',' ').lower()}")
        ax.grid(alpha=.2);ax.legend(loc="best")
        fig.savefig(ROOT/f"{JOB}_{name}.png",dpi=180)
        plt.close(fig)


def render_gif(identity, source, replay, time_end):
    fig,axes=plt.subplots(2,1,figsize=(10,6.5),constrained_layout=True)
    datasets=((SOURCE,source,"No fluid + real wall"),(CASE,replay,"FULL CEL fixed V/VR"))
    artists=[]
    for ax,(case,m,title) in zip(axes,datasets):
        setup_axis(ax,identity,title)
        rp0,c,n,rel,colors,pipe=geometry(case,identity)
        p=pose(0,m["t"],m["u"],m["ur"],rp0,rel)
        dots=ax.scatter(identity["s_start_mm"]+(p-rp0)@c,(p-pipe)@n,s=4,c=colors,linewidths=0)
        label=ax.text(.98,.97,"",transform=ax.transAxes,ha="right",va="top",family="monospace",fontsize=8,
                      bbox={"facecolor":"white","alpha":.9})
        artists.append((dots,label,m,rp0,c,n,rel,pipe))
    frames=np.linspace(0,time_end,81)
    def update(i):
        ti=frames[i]
        for dots,label,m,rp0,c,n,rel,pipe in artists:
            p=pose(ti,m["t"],m["u"],m["ur"],rp0,rel)
            dots.set_offsets(np.c_[identity["s_start_mm"]+(p-rp0)@c,(p-pipe)@n])
            label.set_text(f"t={ti*1000:.2f} ms C{min(4,int(ti/.01)+1)}\n"
                           f"s={np.interp(ti,m['t'],m['s']):+.4f} mm  "
                           f"theta={np.interp(ti,m['t'],m['rocking_angle']):+.2f} deg")
    path=ROOT/f"{JOB}_SYNC.gif"
    FuncAnimation(fig,update,frames=len(frames),interval=75).save(path,PillowWriter(fps=12),dpi=85)
    plt.close(fig)
    with zipfile.ZipFile(ROOT/f"{JOB}_GIF.zip","w",zipfile.ZIP_DEFLATED) as archive:
        archive.write(path,path.name)


def main():
    gate=json.loads((ROOT/f"{JOB}_KINEMATIC_GATE.json").read_text())
    if not gate["passed"]:
        raise RuntimeError("KINEMATIC_REPLAY_INVALID: force interpretation forbidden")
    identity=json.loads((CASE/"case_identity.json").read_text(encoding="utf-8-sig"))
    sta=(CASE/f"{JOB}.sta").read_text(encoding="latin1")
    if "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" not in sta:
        raise RuntimeError("Replay solver did not complete 40 ms")
    with np.load(SOURCE/"private"/"rp_history_private.npz") as source_rp, \
         np.load(CASE/"private"/"rp_replay_private.npz") as replay_rp:
        source_m=motion({k:source_rp[k].astype(float) for k in KEYS},identity)
        replay_m=motion({k:replay_rp[k].astype(float) for k in KEYS},identity)
        t=replay_rp["U1"][:,0].astype(float)
        reaction=history_vector(replay_rp,"RF",t)
        acceleration=history_vector(replay_rp,"A",t)
        angular_acceleration=history_vector(replay_rp,"AR",t)
        omega=history_vector(replay_rp,"VR",t)
    with np.load(SOURCE/"private"/"contact_history_private.npz") as src_contact, \
         np.load(CASE/"private"/"contact_replay_private.npz") as cel_contact:
        source_wall=contact_vector(src_contact,"CFT",False,t)
        wall=contact_vector(cel_contact,"CFT",True,t)
        whole=contact_vector(cel_contact,"CFT",False,t)
        source_normal=np.linalg.norm(contact_vector(src_contact,"CFN",False,t),axis=1)
        replay_normal=np.linalg.norm(contact_vector(cel_contact,"CFN",True,t),axis=1)
        source_tangent=np.linalg.norm(contact_vector(src_contact,"CFS",False,t),axis=1)
        replay_tangent=np.linalg.norm(contact_vector(cel_contact,"CFS",True,t),axis=1)
    fluid=whole-wall
    magnetic_log=np.loadtxt(CASE/"magnetic_increment_g2p20_f100.csv",delimiter=",",skiprows=1)
    magnetic=np.column_stack([np.interp(t,magnetic_log[:,0],magnetic_log[:,j],left=0)
                              for j in (3,4,5)])
    c=np.asarray(identity["canonical_plus_s_axis_aba"],float)
    series={"magnetic":magnetic@c,"wall":wall@c,"fluid":fluid@c,"reaction":reaction@c,
            "no_fluid_wall":source_wall@c}
    audit=json.loads((SOURCE/"setup_audit.json").read_text(encoding="utf-8-sig"))
    mass=audit["robot_mass_mg"]*1e-9
    offset=np.asarray(audit["robot_com_aba_mm"],float)-np.asarray(identity["initial_center_aba_mm"],float)
    rotated=Rotation.from_rotvec(replay_m["ur"]).apply(np.broadcast_to(offset,replay_m["ur"].shape))
    acom=acceleration+np.cross(angular_acceleration,rotated)+np.cross(omega,np.cross(omega,rotated))
    series["inertia"]=(mass*acom)@c
    series["closure_residual"]=series["inertia"]-(series["reaction"]+series["magnetic"]+
                                                  series["wall"]+series["fluid"])
    impulses={key:integrate(t,value) for key,value in series.items()}
    rows=[]
    for i in range(4):
        lo,hi=i*.01,(i+1)*.01
        rows.append({"cycle":i+1,
                     **{f"{key}_J_s_Ns":window_impulse(t,series[key],lo,hi)
                        for key in ("no_fluid_wall","wall","fluid","magnetic","reaction","inertia","closure_residual")},
                     "delta_wall_J_s_Ns":window_impulse(t,series["wall"]-series["no_fluid_wall"],lo,hi)})
    contact_rows,source_gaps,replay_gaps,active,frame_t=events(
        identity,source_m,replay_m,t,source_normal,replay_normal,
        series["no_fluid_wall"],series["wall"],source_tangent,replay_tangent)
    adverse_fluid=sum(min(0.,r["fluid_J_s_Ns"]) for r in rows)
    adverse_wall=sum(min(0.,r["delta_wall_J_s_Ns"]) for r in rows)
    unresolved=False
    peak=max(np.percentile(np.abs(series["inertia"]),99),1e-12)
    closure_ratio=float(np.percentile(np.abs(series["closure_residual"]),99)/peak)
    if not np.all(np.isfinite(np.column_stack(list(series.values())))) or closure_ratio>.1:
        unresolved=True
    if unresolved or adverse_fluid>=0 and adverse_wall>=0:
        classification="FORCE_CAUSAL_ROLE_UNRESOLVED"
    elif adverse_fluid<0 and adverse_wall<0 and min(abs(adverse_fluid),abs(adverse_wall)) >= .25*max(abs(adverse_fluid),abs(adverse_wall)):
        classification="COUPLED_FLUID_WALL_LOAD"
    elif abs(adverse_fluid)>3*abs(adverse_wall):
        classification="DIRECT_FLUID_AXIAL_LOAD_DOMINATES"
    elif abs(adverse_wall)>3*abs(adverse_fluid):
        classification="FLUID_MODIFIES_WALL_CONTACT_DOMINANTLY"
    else:
        classification="FORCE_CAUSAL_ROLE_UNRESOLVED"
    result={"classification":classification,"kinematic_gate":gate,
            "fluid_force_label":"INFERRED ROBOT-CEL FORCE",
            "fluid_force_method":"whole robot General Contact CFT minus direct robot-wall pair CFT",
            "force_sign":"positive along canonical +s", "cycles":rows,
            "contact_events":contact_rows,"closure_p99_to_inertia_p99":closure_ratio,
            "adverse_fluid_impulse_Ns":adverse_fluid,"adverse_delta_wall_impulse_Ns":adverse_wall,
            "force_partition_limitation":"Eulerian robot-fluid SECOND SURFACE output rejected in datacheck; no direct fluid traction",
            "end_force_limitation":"HEAD/TAIL force attribution only when geometric gap singles out one end; overlapping-end force is unresolved"}
    (ROOT/f"{JOB}_FORCE_METRICS.json").write_text(json.dumps(result,indent=2)+"\n")
    with (ROOT/f"{JOB}_CYCLE_IMPULSES.csv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    with (ROOT/f"{JOB}_CONTACT_EVENTS.csv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(contact_rows[0]) if contact_rows else ["case"])
        writer.writeheader();writer.writerows(contact_rows)
    sample=np.unique(np.r_[np.arange(0,len(t),max(1,len(t)//4000)),len(t)-1])
    with (ROOT/f"{JOB}_AXIAL_FORCE_10US.csv").open("w",newline="") as handle:
        writer=csv.writer(handle)
        keys=("magnetic","wall","fluid","reaction","no_fluid_wall","inertia","closure_residual")
        writer.writerow(("time_s",)+tuple(key+"_N" for key in keys)+tuple(key+"_J_Ns" for key in keys))
        writer.writerows(zip(t[sample],*(series[k][sample] for k in keys),*(impulses[k][sample] for k in keys)))
    plot_forces(t,series,impulses)
    render_gif(identity,source_m,replay_m,.04)
    lines=[f"# {JOB}","",classification,"", "Kinematic replay: PASS", "",
           f"Maximum RP position error {gate['results']['max_position_mm']:.6g} mm; "
           f"maximum orientation error {gate['results']['max_orientation_deg']:.6g} deg; "
           f"p99 |delta v_s| {gate['results']['p99_abs_delta_v_s_mm_s']:.6g} mm/s; "
           f"p99 |delta omega_rock| {gate['results']['p99_abs_delta_omega_rock_rad_s']:.6g} rad/s.","",
           "| Axial impulse (N s) | C1 | C2 | C3 | C4 |",
           "|---|---:|---:|---:|---:|"]
    for key,title in (("no_fluid_wall_J_s_Ns","wall J_s no fluid"),("wall_J_s_Ns","wall J_s with CEL"),
                      ("delta_wall_J_s_Ns","delta wall J_s"),("fluid_J_s_Ns","INFERRED ROBOT-CEL J_s"),
                      ("magnetic_J_s_Ns","magnetic J_s"),("reaction_J_s_Ns","reaction J_s")):
        lines.append("| "+title+" | "+" | ".join(f"{r[key]:+.4e}" for r in rows)+" |")
    lines += ["",f"Force closure p99 residual / inertia p99: {closure_ratio:.3g}.",
              "The robot-CEL term is INFERRED ROBOT-CEL FORCE, calculated as whole-robot General Contact minus the direct robot-wall pair. It is not direct fluid traction.",
              "HEAD/TAIL event forces are attributable only for exclusive-end contact; overlapping-end partitions remain unresolved.",
              "Contact events and 10-us force/impulse traces are in adjacent CSV files. No restart continuation or free dynamics candidate was run."]
    lines += ["", "## Contact events", ""]
    for label in ("NO_FLUID_REAL_WALL", "FULL_CEL_FIXED_REPLAY"):
        for end in ("HEAD", "TAIL"):
            subset = [event for event in contact_rows if event["case"] == label and event["end"] == end]
            exclusive = sum(event["force_attribution"] == "exclusive end" for event in subset)
            lines.append(f"- {label} {end}: {len(subset)} sampled events; {exclusive} exclusive-end events.")
    lines += ["", "The matched prescribed motion preserves the sampled contact-event times. The force differences describe the loads required along this trajectory; they do not establish a new free-dynamics trajectory or a change in collision timing."]
    (ROOT/f"{JOB}_FORCE_REPORT.md").write_text("\n".join(lines)+"\n")
    identity.update(status="SOLVED",classification=classification,kinematic_gate_passed=True)
    (CASE/"case_identity.json").write_text(json.dumps(identity,indent=2)+"\n")
    print(json.dumps({"classification":classification,"cycles":rows,
                      "closure_ratio":closure_ratio,"events":len(contact_rows)},indent=2))


if __name__=="__main__":
    main()
