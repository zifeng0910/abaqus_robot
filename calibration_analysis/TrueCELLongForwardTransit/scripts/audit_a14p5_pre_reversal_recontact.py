"""Read-only native-ledger and sampled-geometry audit of the solved A14P5 parent."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon
import numpy as np
from scipy.signal import savgol_filter
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation

import analyze_a14p5_tail_contact_refine as ledger
import analyze_a14p5_tailgap_refine as geometry

ROOT = Path(__file__).resolve().parent.parent
JOB = ledger.JOB
START, ZERO = 1 / 120, .010682020978413481
RELEASE, RECONTACT = .011550149880349636, .011675059795379639
COLORS = {"magnetic": "#b25335", "wall": "#26638d", "fluid_inferred": "#22846f",
          "momentum": "#242a31", "TAIL": "#2864a8", "HEAD": "#d1493f"}


def writetable(name, rows):
    with (ROOT / name).open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def spans(time, active, start=None, end=None):
    edges = np.diff(np.r_[False, active, False].astype(int))
    intervals = zip(np.where(edges == 1)[0], np.where(edges == -1)[0] - 1)
    return [(max(time[a], start or time[0]), min(time[b], end or time[-1]))
            for a, b in intervals if (start is None or time[b] >= start) and
            (end is None or time[a] <= end)]


def duration_start(t, flag, minimum_s):
    for lo, hi in spans(t, flag):
        if hi - lo >= minimum_s:
            return lo
    return None


def sampled_state(case, t):
    return {key: np.interp(t, case["t"], case[key]) for key in ("angle", "command", "omega", "vs")}


def pressure_history(case, field):
    times = field["time"].astype(float)
    centers = field["fluid_element_centroids_mm"].astype(float)
    s = (centers - case["identity"]["initial_center_aba_mm"]) @ case["s"]
    pressure = field["fluid_pressure"]
    evf = field["fluid_evf"]
    u_s = np.interp(times, case["t"], case["u"] @ case["s"])
    tail, head = [], []
    for k in range(len(times)):
        local_s = s - u_s[k]
        values = []
        for low, high in ((-1.8, -.6), (.6, 1.8)):
            mask = (local_s >= low) & (local_s < high) & (evf[k] > .01)
            w = evf[k, mask]
            values.append(float(np.average(pressure[k, mask], weights=w)) if w.size else np.nan)
        tail.append(values[0]); head.append(values[1])
    return times, np.array(tail), np.array(head)


def wall_geometry(case):
    inp = (case["case"] / f"{JOB}.inp").read_text(encoding="latin1")
    _, wall = geometry.base.part_nodes(inp, "Pipe_WALL_HELPER")
    center = np.array(case["identity"]["initial_center_aba_mm"])
    pipe = center - case["identity"]["initial_axial_shift_mm"] * case["s"] - case["identity"]["radial_offset_n_mm"] * case["n"]
    normals, radius = geometry.base.wall_planes(wall, pipe, case["s"], case["n"], case["b"])
    contact = np.load(case["case"] / "private" / "robot_contact_fields_private.npz")
    nodes = contact["node_coordinates_mm"].astype(float)
    material = (nodes - center) @ case["s"]
    tail = nodes[material <= material.min() + .25] - center
    return dict(center=center, pipe=pipe, normals=normals, radius=radius, tail=tail,
                nodes=nodes, contact=contact)


def gap_from_pose(geo, u, ur):
    points = geo["center"] + u + Rotation.from_rotvec(ur).apply(geo["tail"])
    support = (points - geo["pipe"]) @ geo["normals"].T
    return float(geo["radius"] - np.max(support))


def gap_series(case, geo, times):
    u = ledger.interp(case["t"], case["u"], times)
    ur = ledger.interp(case["t"], case["ur"], times)
    return np.array([gap_from_pose(geo, uk, rk) for uk, rk in zip(u, ur)]), u, ur


def decompose_gap(geo, u, ur):
    rot, trans = [], []
    for i in range(1, len(u)):
        g00 = gap_from_pose(geo, u[i-1], ur[i-1])
        g10 = gap_from_pose(geo, u[i], ur[i-1])
        g01 = gap_from_pose(geo, u[i-1], ur[i])
        g11 = gap_from_pose(geo, u[i], ur[i])
        trans.append(.5 * ((g10 - g00) + (g11 - g01)))
        rot.append(.5 * ((g01 - g00) + (g11 - g10)))
    return np.r_[0, np.cumsum(rot)], np.r_[0, np.cumsum(trans)]


def annotate(ax, times=(ZERO,)):
    for t, label in ((START, "Cycle 2"), (.00930870045, "wall onset"),
                     (ZERO, "v_s=0"), (RELEASE, "TAIL release"), (RECONTACT, "TAIL recontact")):
        if times == "all" or t in times:
            ax.axvline(t * 1000, color="#82888d", ls="--", lw=.8)


def main():
    case = ledger.load(JOB)
    rp_source = np.load(case["case"] / "private" / "rp_history_private.npz")
    ur_time, case["ur"] = ledger.vectors(rp_source, "UR")
    assert np.array_equal(ur_time, case["t"])
    assert abs(case["mass"] * 1e9 - 10.018809) < 1e-5
    ctime = case["contact_time"] if "contact_time" in case else case["geom"]["contact_time"]
    active = case["geom"]["active"]
    ftime = case["geom"]["contact_time"]
    wall = case["forces"][1]
    inferred = case["forces"][2]
    magnetic = case["forces"][0]
    t = case["t"]
    dt = float(np.median(np.diff(t)))
    vs = case["vs"]
    # Compare three finite, symmetric Savitzky-Golay differentiation windows;
    # do not equate an algorithmic departure to a physically unique trigger.
    widths = [max(7, int(w / dt) // 2 * 2 + 1) for w in (15e-6, 30e-6, 60e-6)]
    deriv = [savgol_filter(vs, width, 3, deriv=1, delta=dt, mode="interp") for width in widths]
    accel = deriv[1]
    first_negative = [duration_start(t, (a < 0) & (t >= START) & (t <= ZERO), 30e-6) for a in deriv]
    baseline = (t >= START) & (t <= .00850)
    slope, intercept = np.polyfit(t[baseline]-START, vs[baseline], 1)
    predicted = intercept + slope * (t-START)
    residual = vs-predicted
    sigma = 1.4826*np.median(abs(residual[baseline]-np.median(residual[baseline])))
    extra = [duration_start(t, (residual < -max(mult*sigma, floor)) & (t >= .0085) & (t <= ZERO),
                            hold) for mult, floor, hold in ((6,.3,15e-6),(3,.15,30e-6))]
    aidx = int(np.argmin(np.where((t >= START) & (t <= ZERO), accel, np.inf)))
    closure = ledger.impulse(case, False, START, ZERO)
    angular = ledger.impulse(case, True, RELEASE, RECONTACT)
    if closure["relative_residual"] >= .10:
        raise RuntimeError("Linear momentum closure failed; stop physical force attribution")

    boundaries = [START,.0088,.0093,.0098,.0103,ZERO]
    # Native wall onset and long TAIL support boundary are kept as distinct
    # diagnostic boundaries; they are not interchangeable time measurements.
    boundaries = sorted(set(boundaries + [.009250123053789139,.009308700449764729]))
    windows = []
    for lo, hi in zip(boundaries[:-1], boundaries[1:]):
        result = ledger.impulse(case, False, lo, hi)
        windows.append(dict(start_ms=lo*1000,end_ms=hi*1000,
                            delta_p_Ns=result["measured_momentum_change"],
                            Jmag_Ns=result["magnetic"],Jwall_Ns=result["wall"],
                            Jfluid_inferred_Ns=result["fluid_inferred"],
                            closure_residual_Ns=result["residual"]))
    writetable("A14P5_PRE_REVERSAL_NATIVE_IMPULSE_WINDOWS.csv", windows)

    field = np.load(case["case"] / "private" / "truecel_field_private.npz")
    pt, tail_p, head_p = pressure_history(case, field)
    asym = tail_p - head_p
    writetable("A14P5_PRE_REVERSAL_PRESSURE_ASYMMETRY.csv", [
        dict(time_ms=q*1000,pressure_tail_N_mm2=p1,pressure_head_N_mm2=p2,
             tail_minus_head_N_mm2=p1-p2,definition="EVF-weighted; material s [-1.8,-.6] vs [.6,1.8] mm")
        for q,p1,p2 in zip(pt,tail_p,head_p)])

    # Native force sources are never integrated on a resampled common grid.
    show = np.linspace(START,.0108,2200)
    ms = show*1000
    fig, axes = plt.subplots(4,1,figsize=(12,10),sharex=True,layout="constrained")
    axes[0].plot(ms,np.interp(show,t,vs),color=COLORS["momentum"],label="v_s (mm/s)")
    a2=axes[0].twinx(); a2.plot(ms,np.interp(show,t,accel),color=COLORS["magnetic"],alpha=.75,lw=.8)
    a2.set_ylabel("robust a_s (mm/s2)")
    axes[0].set_ylabel("v_s (mm/s)")
    for y,label,color,src in ((magnetic,"magnetic",COLORS["magnetic"],case["mt"]),
                              (wall,"wall",COLORS["wall"],case["ct"]),
                              (inferred,"inferred robot-fluid",COLORS["fluid_inferred"],case["ct"])):
        axes[1].plot(ms,np.interp(show,src,y),label=label,color=color,lw=.85)
    axes[1].set_ylabel("axial force (N)"); axes[1].legend(ncol=3,fontsize=8)
    for values,label,color,src in ((magnetic,"Jmag",COLORS["magnetic"],case["mt"]),
                                   (wall,"Jwall",COLORS["wall"],case["ct"]),
                                   (inferred,"Jfluid inferred",COLORS["fluid_inferred"],case["ct"])):
        axes[2].plot(ms,ledger.cumulative(src,values,START,show)*1e9,label=label,color=color)
    axes[2].plot(ms,case["mass"]*(np.interp(show,t,vs)-np.interp(START,t,vs))*1e9,
                 color=COLORS["momentum"],ls="--",label="measured delta p")
    axes[2].set_ylabel("cumulative impulse / delta p (nN s)");axes[2].legend(ncol=4,fontsize=8)
    for key in ("TAIL","HEAD"):
        for lo,hi in spans(ftime,active[key],START,.0108):
            axes[3].axvspan(lo*1000,hi*1000,alpha=.55,color=COLORS[key],label=key if not any(p.get_label()==key for p in axes[3].patches[:-1]) else None)
    axes[3].step(ftime*1000,active["TAIL"].astype(int)+active["HEAD"].astype(int),where="mid",color="#333",lw=.8)
    axes[3].set(xlim=(START*1000,10.8),ylim=(-.1,2.15),
                ylabel="whole-GC regional flags\nHEAD+TAIL",xlabel="time (ms)")
    for ax in axes:
        annotate(ax,(ZERO,));ax.grid(alpha=.17)
    fig.suptitle("A14P5 refined: source-native integrated axial ledger; whole-GC regional flags are not wall-pair contact")
    fig.savefig(ROOT/"A14P5_PRE_REVERSAL_FORCE_IMPULSE_BALANCE.png",dpi=160);plt.close(fig)

    fig,axes=plt.subplots(3,1,figsize=(11,8),sharex=True,layout="constrained")
    subset=(pt>=START)&(pt<=.0108)
    axes[0].plot(pt[subset]*1000,asym[subset],color="#b25335",marker=".",ms=2)
    axes[0].axhline(0,color="gray",lw=.6)
    axes[0].set_ylabel("TAIL minus HEAD pressure (N/mm2)")
    axes[1].plot(case["ct"]*1000,inferred,color=COLORS["fluid_inferred"],lw=.75)
    axes[1].set_ylabel("inferred robot-fluid F_s (N)")
    axes[2].plot(t*1000,vs,color=COLORS["momentum"],lw=.9)
    axes[2].set(xlim=(START*1000,10.8),ylim=(-2,17),ylabel="v_s (mm/s)",xlabel="time (ms)")
    for ax in axes:annotate(ax,(ZERO,));ax.grid(alpha=.17)
    fig.suptitle("EVF-weighted local CEL pressure asymmetry: positive TAIL-HEAD favors +s (proxy only)")
    fig.savefig(ROOT/"A14P5_PRE_REVERSAL_FLUID_ASYMMETRY.png",dpi=160);plt.close(fig)

    # Fixed-range maps: pressure, axial velocity and EVF at six requested field frames.
    fields=geometry.load_case(JOB)
    cent=field["fluid_element_centroids_mm"].astype(float)-fields["pipe"]
    fs,fn,fb=cent@case["s"],cent@case["n"],cent@case["b"]
    p=field["fluid_pressure"];evf=field["fluid_evf"]
    conn=field["fluid_connectivity_index"].astype(int)
    queries=np.array([8.5,9.,9.5,10.,10.4,ZERO*1000])
    idx=np.array([np.argmin(abs(pt-q/1000)) for q in queries])
    pressure_limit=max(float(np.percentile(abs(p[idx][evf[idx]>.01]),99)),1e-9)
    velocity=[]
    for k in idx:
        nodal=field["fluid_velocity_mm_s"][k].astype(float)
        velocity.append(np.nanmean(nodal[conn],axis=1)@case["s"])
    velocity=np.array(velocity)
    velocity_limit=max(float(np.percentile(abs(velocity[np.isfinite(velocity)]),99)),1)
    fig,axes=plt.subplots(3,6,figsize=(18,6.6),sharex=True,sharey=True)
    scatters=[]
    for j,k in enumerate(idx):
        xy=fields["center0"]+ledger.interp(fields["time"],fields["u"],pt[k])+Rotation.from_rotvec(
            ledger.interp(fields["time"],fields["ur"],pt[k])).apply(
            np.load(case["case"]/"private"/"robot_contact_fields_private.npz")["node_coordinates_mm"]-fields["center0"])
        projected=np.column_stack(((xy-fields["pipe"])@case["s"],(xy-fields["pipe"])@case["n"]))
        hull=ConvexHull(projected).vertices
        mask=(abs(fb)<.066)&(evf[k]>.01)&(fs>-5)&(fs<-1.1)
        for i,(values,vmin,vmax,cmap) in enumerate(((p[k],-pressure_limit,pressure_limit,"coolwarm"),
                                                     (velocity[j],-velocity_limit,velocity_limit,"coolwarm"),
                                                     (evf[k],0,1,"viridis"))):
            ax=axes[i,j]
            selected=(abs(fb)<.066)&(fs>-5)&(fs<-1.1) if i==2 else mask
            im=ax.scatter(fs[selected],fn[selected],c=values[selected],cmap=cmap,vmin=vmin,vmax=vmax,
                          s=10,marker="s",linewidths=0,rasterized=True)
            ax.add_patch(Polygon(projected[hull],facecolor="#e4e8e8",edgecolor="#242a31",lw=.8))
            ax.axhline(fields["identity"]["lumen_radius_mm"],color="#777",lw=.8)
            ax.axhline(-fields["identity"]["lumen_radius_mm"],color="#777",lw=.8)
            ax.set(xlim=(-5,-1.1),ylim=(-.75,.75),aspect="equal")
            if i==0:ax.set_title(f"{pt[k]*1000:.3f} ms")
            if j==0:ax.set_ylabel(("pressure N/mm2","fluid v_s mm/s","EVF (0-1)")[i]+"\nn (mm)")
            if i==2:ax.set_xlabel("+s (mm)")
            if j==5:scatters.append(im)
    fig.subplots_adjust(left=.065,right=.91,top=.88,bottom=.10,hspace=.21,wspace=.10)
    for i,im in enumerate(scatters):
        cax=fig.add_axes((.93,.72-i*.28,.012,.16))
        fig.colorbar(im,cax=cax)
    fig.suptitle("Fixed b=0 CEL side slice; each row has one range shared across all times",y=.97)
    fig.savefig(ROOT/"A14P5_PRE_REVERSAL_LOCAL_FLUID_MAPS.png",dpi=140);plt.close(fig)

    geo=wall_geometry(case)
    regional=(ftime>=.0113)&(ftime<=.0118)
    gt=ftime[regional]
    gaps,u,ur=gap_series(case,geo,gt)
    gap_rate=np.gradient(gaps,gt)
    # Include exact event boundaries by native-RP interpolation, independent of
    # the 25-us field output used for the contact ON/OFF flag.
    sample_t=np.sort(np.unique(np.r_[gt,RELEASE,RECONTACT]))
    sample_gap,su,sur=gap_series(case,geo,sample_t)
    rot_contrib,trans_contrib=decompose_gap(geo,su,sur)
    scope=(t>=.0113)&(t<=.0118)
    mag_rock=np.interp(t[scope],case["mt"],case["moments"][0])
    fig,axes=plt.subplots(5,1,figsize=(11,10),sharex=True,layout="constrained")
    axes[0].plot(t[scope]*1000,case["angle"][scope],label="actual",color="#242a31")
    axes[0].plot(t[scope]*1000,case["command"][scope],label="command",color="#b25335",ls="--")
    axes[0].set_ylabel("rocking angle (deg)");axes[0].legend()
    axes[1].plot(t[scope]*1000,case["omega"][scope]*180/np.pi,color="#242a31")
    axes[1].set_ylabel("omega_rock (deg/s)")
    axes[2].plot(t[scope]*1000,mag_rock,color=COLORS["magnetic"])
    axes[2].set_ylabel("magnetic torque (N mm)")
    axes[3].plot(gt*1000,gaps*1000,color=COLORS["TAIL"],marker=".")
    axes[3].set_ylabel("signed TAIL gap (um)")
    axes[4].plot(gt*1000,gap_rate,color=COLORS["TAIL"],label="dg/dt")
    axes[4].step(ftime*1000,active["TAIL"].astype(float)*max(10,np.nanmax(abs(gap_rate))*.7),
                 where="mid",color=COLORS["wall"],alpha=.6,label="TAIL whole-GC flag (scaled)")
    axes[4].set(xlim=(11.3,11.8),ylabel="gap rate (mm/s)",xlabel="time (ms)")
    axes[4].legend(fontsize=8)
    for ax in axes:annotate(ax,"all");ax.grid(alpha=.17)
    fig.suptitle("TAIL regional-flag toggles: signed gap to actual discrete wall planes; wall-pair force = 0")
    fig.savefig(ROOT/"A14P5_TAIL_RECONTACT_KINEMATICS.png",dpi=160);plt.close(fig)

    subset=(sample_t>=RELEASE)&(sample_t<=RECONTACT)
    g0=sample_gap[np.flatnonzero(subset)[0]]
    rs,ts=decompose_gap(geo,su[subset],sur[subset])
    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True,layout="constrained")
    axes[0].plot(sample_t[subset]*1000,(sample_gap[subset]-g0)*1000,color="#242a31",label="net gap change")
    axes[0].plot(sample_t[subset]*1000,rs*1000,color="#b25335",label="rotation contribution")
    axes[0].plot(sample_t[subset]*1000,ts*1000,color="#26638d",label="translation contribution")
    axes[0].set_ylabel("change since release (um)");axes[0].legend()
    axes[1].plot(sample_t[subset]*1000,sample_gap[subset]*1000,color=COLORS["TAIL"])
    axes[1].axhline(0,color="gray",lw=.7)
    axes[1].set(xlabel="time (ms)",ylabel="actual signed gap (um)")
    for ax in axes:
        annotate(ax,(RELEASE,RECONTACT))
        ax.set_xlim(RELEASE*1000,RECONTACT*1000)
        ax.grid(alpha=.17)
    fig.suptitle("Exact endpoint gap partition (two-order average); fixed straight wall curvature = 0")
    fig.savefig(ROOT/"A14P5_TAIL_RECONTACT_GAP_DECOMPOSITION.png",dpi=160);plt.close(fig)

    gif_times=np.linspace(.0113,.0118,101)
    fig,ax=plt.subplots(figsize=(10.7,4.5))
    nodes=geo["nodes"]
    material=(nodes-geo["center"])@case["s"]
    tail_nodes=material<=material.min()+.25
    def draw(frame):
        moment=gif_times[frame]
        ax.clear()
        pos=geo["center"]+ledger.interp(case["t"],case["u"],moment)+Rotation.from_rotvec(
            ledger.interp(case["t"],case["ur"],moment)).apply(nodes-geo["center"])
        projection=np.column_stack(((pos-geo["pipe"])@case["s"],(pos-geo["pipe"])@case["n"]))
        hull=ConvexHull(projection).vertices
        ax.add_patch(Polygon(projection[hull],facecolor="#e5e9e9",edgecolor="#242a31",lw=1.1))
        ax.scatter(projection[tail_nodes,0],projection[tail_nodes,1],color=COLORS["TAIL"],s=5,linewidths=0)
        ax.axhline(case["identity"]["lumen_radius_mm"],color="#71929b")
        ax.axhline(-case["identity"]["lumen_radius_mm"],color="#71929b")
        state=sampled_state(case,moment)
        gap=gap_from_pose(geo,ledger.interp(case["t"],case["u"],moment),
                              ledger.interp(case["t"],case["ur"],moment))
        contact=active["TAIL"][np.argmin(abs(ftime-moment))]
        ax.set(xlim=(-4.95,-2.8),ylim=(-.75,.75),aspect="equal",xlabel="canonical +s (mm)",ylabel="n (mm)")
        ax.set_title(f"TAIL regional-flag toggle | t={moment*1000:.3f} ms | "
                     f"whole-GC flag {'ON' if contact else 'OFF'} (25-us sampled); wall-pair F=0\n"
                     f"gap={gap*1000:+.2f} um | actual={state['angle']:+.2f} deg | "
                     f"cmd={state['command']:+.2f} deg | omega={state['omega']*180/np.pi:+.0f} deg/s | "
                     f"v_s={state['vs']:+.2f} mm/s",fontsize=10,loc="left")
        fig.tight_layout()
    FuncAnimation(fig,draw,frames=len(gif_times),interval=110).save(
        ROOT/"A14P5_TAIL_RELEASE_RECONTACT_ULTRASLOW.gif",writer=PillowWriter(fps=9),dpi=95)
    plt.close(fig)

    events=[]
    for label,q in (("TAIL release",RELEASE),("TAIL recontact",RECONTACT)):
        state=sampled_state(case,q)
        events.append(dict(event=label,time_ms=q*1000,actual_deg=state["angle"],command_deg=state["command"],
                           error_deg=state["angle"]-state["command"],omega_deg_s=state["omega"]*180/np.pi,
                           torque_Nmm=np.interp(q,case["mt"],case["moments"][0]),
                           gap_um=gap_from_pose(geo,ledger.interp(case["t"],case["u"],q),
                                                ledger.interp(case["t"],case["ur"],q))*1000))
    writetable("A14P5_TAIL_RECONTACT_NATIVE_EVENTS.csv",events)
    summary=dict(mass_mg=case["mass"]*1e9,native_RP_dt_ms=dt*1000,
                 native_field_dt_ms=float(np.median(np.diff(ftime)))*1000,
                 first_negative_acceleration_ms=[None if v is None else v*1000 for v in first_negative],
                 first_additional_deceleration_ms=[None if v is None else v*1000 for v in extra],
                 early_cycle2_baseline_acceleration_mm_s2=slope,
                 maximum_negative_acceleration_ms=t[aidx]*1000,
                 maximum_negative_acceleration_mm_s2=accel[aidx],
                 velocity_zero_ms=ZERO*1000,axial_closure=closure,
                 angular_release_to_recontact=angular,events=events,
                 gap_min_between_release_recontact_um=float(np.min(sample_gap[subset])*1000),
                 gap_max_between_release_recontact_um=float(np.max(sample_gap[subset])*1000),
                 gap_at_release_um=events[0]["gap_um"],gap_at_recontact_um=events[1]["gap_um"],
                 gap_rotation_change_um=float(rs[-1]*1000),gap_translation_change_um=float(ts[-1]*1000),
                 field_contact_state_at_release=bool(active["TAIL"][np.argmin(abs(ftime-RELEASE))]),
                 field_contact_state_at_recontact=bool(active["TAIL"][np.argmin(abs(ftime-RECONTACT))]),
                 max_direct_wall_pair_force_11p3_11p8_N=float(np.max(np.linalg.norm(
                     case["wall_force_vector"][(case["ct"]>=.0113)&(case["ct"]<=.0118)],axis=1))))
    (ROOT/"A14P5_PRE_REVERSAL_RECONTACT_AUDIT.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="ascii")
    print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__":
    main()
