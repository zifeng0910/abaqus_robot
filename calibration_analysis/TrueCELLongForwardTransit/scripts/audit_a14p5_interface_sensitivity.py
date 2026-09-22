"""Offline comparison of the two already-solved A14P5 robot interfaces.

No ODB access, Abaqus invocation, or model mutation. Force histories retain
their own native integration clocks; nodal fields are whole General Contact.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation

import analyze_a14p5_tail_contact_refine as ledger
from prepare_a14p5_tail_contact_refine import parse_mesh, exterior

ROOT = Path(__file__).resolve().parents[1]
START = 1 / 120
END = .010682020978413481
ENDS = (.009, .0095, .010, .0103, END)
NAMES = (ledger.BASE, ledger.JOB)
REGIONS = ("TAIL", "MID", "HEAD")
COLORS = ("#267b73", "#b54f39")


def region(x):
    return np.where(x < -.6, 0, np.where(x > .6, 2, 1))


def geometry(data):
    deck = (data["case"] / (data["job"] + ".inp")).read_text(encoding="latin1")
    _, _, _, nodes, tets, _ = parse_mesh(deck)
    faces = exterior(tets)
    tri = np.array([[nodes[i] for i in face[2]] for face in faces])
    area = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    center = (tri.mean(axis=1) - data["identity"]["initial_center_aba_mm"]) @ data["s"]
    bins = region(center)
    return np.bincount(bins, weights=area, minlength=3), np.bincount(bins, minlength=3), nodes, faces


def fields(data, geo):
    f = np.load(data["case"] / "private/robot_contact_fields_private.npz")
    t = f["time"].astype(float)
    x = (f["node_coordinates_mm"].astype(float) - data["identity"]["initial_center_aba_mm"]) @ data["s"]
    bins = region(x)
    force = f["CNORMF"].astype(float) + f["CSHEARF"].astype(float)
    projected = np.stack([force @ axis for axis in (data["s"], data["n"], data["b"])], axis=2)
    regional = np.stack([projected[:, bins == k].sum(axis=1) for k in range(3)], axis=1)
    # A nodal area is a geometric tributary, not a contact-patch area.
    nodal_area = {int(label): 0. for label in f["node_labels"]}
    for _, _, face in geo[3]:
        v = np.array([geo[2][i] for i in face])
        a = np.linalg.norm(np.cross(v[1]-v[0], v[2]-v[0])) / 6
        for label in face:
            nodal_area[label] = nodal_area.get(label, 0.) + a
    area_node = np.array([nodal_area[int(label)] for label in f["node_labels"]])
    active = np.linalg.norm(force,axis=2)>1e-8
    return dict(t=t, x=x, bins=bins, projected=projected, regional=regional,
                node_area=area_node, node_count=np.bincount(bins, minlength=3),
                active_count=np.stack([active[:,bins==k].sum(axis=1) for k in range(3)],axis=1),
                source=list(f["source_names"]))


def integral(t, y, lo, hi):
    interior = (t > lo) & (t < hi)
    tt = np.r_[lo, t[interior], hi]
    return float(np.trapz(np.interp(tt, t, y), tt))


def write_csv(name, rows):
    with (ROOT / name).open("w", newline="", encoding="ascii") as fp:
        writer = csv.DictWriter(fp, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def native_value(d, t):
    return [float(np.interp(t, clock, value)) for clock, value in
            ((d["t"], d["disp"]), (d["t"], d["vs"]),
             (d["t"], d["angle"]), (d["t"], d["omega"]),
             (d["mt"], d["forces"][0]), (d["ct"], d["forces"][1]),
             (d["ct"], d["forces"][2]))]


def onset(t, y, threshold, hold=30e-6):
    # Threshold is an operational sustained-load marker, not a unique physical onset.
    k = np.flatnonzero((t >= START) & (t <= END) & (y < threshold))
    for i in k:
        stop = np.searchsorted(t, t[i] + hold)
        if stop < len(t) and np.all(y[i:stop+1] < threshold):
            return float(t[i])
    return None


def rolling_force_onset(t, force, window=50e-6, threshold=-5e-5):
    cumulative=np.r_[0,np.cumsum((force[:-1]+force[1:])*np.diff(t)/2)]
    eligible=(t>=START+window)&(t<=END)
    average=(cumulative-np.interp(t-window,t,cumulative,left=0))/window
    ids=np.flatnonzero(eligible&(average<threshold))
    return None if not len(ids) else float(t[ids[0]])


def cel_metrics(d):
    z = np.load(d["case"] / "private/truecel_field_private.npz")
    t = z["time"].astype(float)
    center = z["fluid_element_centroids_mm"]
    s0 = (center - d["identity"]["initial_center_aba_mm"]) @ d["s"]
    pressure = z["fluid_pressure"]
    evf = z["fluid_evf"]
    values = np.full(len(t), np.nan)
    for i in np.flatnonzero((t >= START) & (t <= .0108)):
        local = s0 - np.interp(t[i], d["t"], d["disp"])
        mask = (np.abs(local) <= 1.8) & (evf[i] > .01)
        values[i] = np.average(np.abs(pressure[i, mask]), weights=evf[i, mask])
    return z, t, values


def main():
    coarse, fine = [ledger.load(n, i == 0) for i, n in enumerate(NAMES)]
    original = np.load(ROOT / "case/TRUECEL_A14P5_DIAG_2CYCLES/private/rp_history_private.npz")
    proxy = np.load(coarse["case"] / "private/rp_history_private.npz")
    proxy_error = {key: float(np.max(np.abs(original[key] - proxy[key]))) for key in
                   ("U1", "U2", "U3", "UR1", "UR2", "UR3", "V1", "V2", "V3", "VR1", "VR2", "VR3")}
    assert max(proxy_error.values()) == 0 and np.array_equal(coarse["t"], original["U1"][:, 0])
    data = (coarse, fine)
    geos = [geometry(d) for d in data]
    contact = [fields(d, g) for d, g in zip(data, geos)]
    cel = [cel_metrics(d) for d in data]
    rows = []
    for end in ENDS:
        pair = [ledger.impulse(d, False, START, end) for d in data]
        sensitivity = abs(pair[1]["fluid_inferred"] - pair[0]["fluid_inferred"]) / max(abs(pair[0]["fluid_inferred"]), 1e-15)
        for i, d in enumerate(data):
            q, g, f = pair[i], geos[i], contact[i]
            region_imp = [integral(f["t"], f["regional"][:, k, 0], START, end) for k in range(3)]
            rows.append(dict(case=d["job"], start_ms=START*1e3, end_ms=end*1e3,
                magnetic_phase_start_deg=360*120*START, magnetic_phase_end_deg=360*120*end,
                J_mag_s_Ns=q["magnetic"], J_wall_s_Ns=q["wall"],
                J_fluid_inferred_s_Ns=q["fluid_inferred"], Delta_p_s_Ns=q["measured_momentum_change"],
                closure_residual_Ns=q["residual"], closure_relative=q["relative_residual"],
                TAIL_whole_GC_nodal_impulse_Ns=region_imp[0],
                MID_whole_GC_nodal_impulse_Ns=region_imp[1],
                HEAD_whole_GC_nodal_impulse_Ns=region_imp[2],
                TAIL_area_mm2=g[0][0], TAIL_facet_count=g[1][0],
                TAIL_surface_node_count=f["node_count"][0],
                TAIL_peak_active_contact_node_count=int(f["active_count"][(f["t"]>=START)&(f["t"]<=end),0].max()),
                TAIL_whole_GC_nodal_impulse_per_area_Ns_mm2=region_imp[0]/g[0][0],
                v_s_end_mm_s=np.interp(end,d["t"],d["vs"]),
                delta_s_mm=np.interp(end,d["t"],d["disp"])-np.interp(START,d["t"],d["disp"]),
                S_J=sensitivity, local_force_semantics="whole_GC_nodal_not_pair_isolated"))
    write_csv("A14P5_ROBOT_CEL_INTERFACE_SENSITIVITY.csv", rows)

    phase_rows=[]
    for ms in (8.5, 9., 9.5, 10., 10.3, END*1e3):
        for d in data:
            phase_rows.append(dict(case=d["job"], time_ms=ms, magnetic_phase_deg=360*120*ms/1e3,
                **dict(zip(("s_mm","v_s_mm_s","rocking_deg","omega_rad_s","F_mag_N","F_wall_N","F_fluid_inferred_N"),native_value(d,ms/1e3)))))
    write_csv("A14P5_INTERFACE_COMMON_PHASE_STATES.csv", phase_rows)

    # Require closeness in all three independent coordinates, before wall-induced divergence.
    states=[]
    for tm in np.arange(.00835, .00930, 25e-6):
        cf=native_value(coarse,tm)
        candidates=fine["t"][(fine["t"]>=START)&(fine["t"]<=.00930)]
        cost=((np.interp(candidates,fine["t"],fine["angle"])-cf[2])/.15)**2 + \
             ((np.interp(candidates,fine["t"],fine["omega"])-cf[3])/5.)**2 + \
             ((candidates-tm)/25e-6)**2
        tf=candidates[np.argmin(cost)]
        ff=native_value(fine,tf)
        if abs(tf-tm)<=25e-6 and abs(ff[2]-cf[2])<=.15 and abs(ff[3]-cf[3])<=5.:
            states.append(dict(coarse_time_ms=tm*1e3, refined_time_ms=tf*1e3,
                phase_delta_deg=(tf-tm)*360*120, angle_delta_deg=ff[2]-cf[2],
                omega_delta_rad_s=ff[3]-cf[3], coarse_fluid_N=cf[6], refined_fluid_N=ff[6],
                fluid_force_delta_N=ff[6]-cf[6]))
    if states: write_csv("A14P5_INTERFACE_STATE_MATCHED.csv", states)

    show=np.linspace(START,.0108,1800)
    fig,ax=plt.subplots(4,1,figsize=(11,9),sharex=True,layout="constrained")
    for i,d in enumerate(data):
        for k,clock in enumerate((d["ct"],d["ct"],d["mt"])):
            source=(2,1,0)[k]
            ax[k].plot(show*1e3,np.interp(show,clock,d["forces"][source]),color=COLORS[i],lw=.85,label=("coarse","TAIL refined")[i])
        ax[3].plot(show*1e3,np.interp(show,d["t"],d["vs"]),color=COLORS[i],label=("coarse","TAIL refined")[i])
    for k,label in enumerate(("inferred robot-CEL F_s (N)","direct wall F_s (N)","magnetic F_s (N)","v_s (mm/s)")):
        ax[k].set_ylabel(label);ax[k].grid(alpha=.2);ax[k].legend(fontsize=8)
    ax[-1].set_xlabel("physical time (ms)")
    fig.suptitle("A14P5: native resultant histories on common physical time")
    fig.savefig(ROOT/"A14P5_COARSE_VS_REFINE_FLUID_FORCE_HISTORY.png",dpi=160);plt.close(fig)

    fig,ax=plt.subplots(4,1,figsize=(11,9),sharex=True,layout="constrained")
    for i,d in enumerate(data):
        for k,src in enumerate((2,1,0)):
            ax[k].plot(show*1e3,ledger.cumulative((d["ct"],d["ct"],d["mt"])[k],d["forces"][src],START,show)*1e9,
                       color=COLORS[i],label=("coarse","TAIL refined")[i])
        ax[3].plot(show*1e3,d["mass"]*(np.interp(show,d["t"],d["vs"])-np.interp(START,d["t"],d["vs"]))*1e9,
                   color=COLORS[i],label=("coarse","TAIL refined")[i])
    for k,label in enumerate(("inferred robot-CEL J_s", "direct wall J_s", "magnetic J_s", "measured delta p_s")):
        ax[k].set_ylabel(label+" (nN s)");ax[k].grid(alpha=.2);ax[k].legend(fontsize=8)
    ax[-1].set_xlabel("physical time (ms)")
    fig.suptitle("A14P5: separately native-integrated impulses, identical start 8.333333 ms")
    fig.savefig(ROOT/"A14P5_COARSE_VS_REFINE_FLUID_IMPULSE.png",dpi=160);plt.close(fig)

    # Equal body-coordinate bins and one symmetric color limit across cases.
    edges=np.linspace(-1.35,1.35,55)
    maps=[]
    for f in contact:
        image=np.zeros((len(f["t"]),len(edges)-1))
        b=np.digitize(f["x"],edges)-1
        for j in range(len(edges)-1):
            image[:,j]=f["projected"][:,b==j,0].sum(axis=1)/(edges[j+1]-edges[j])
        maps.append(image)
    limit=max(np.percentile(np.abs(m[(contact[i]["t"]>=START)&(contact[i]["t"]<=.0108)]),99.5) for i,m in enumerate(maps))
    fig,ax=plt.subplots(2,1,figsize=(11,6.5),sharex=True,layout="constrained")
    for i,(f,m) in enumerate(zip(contact,maps)):
        image=ax[i].imshow(m.T,origin="lower",aspect="auto",extent=(f["t"][0]*1e3,f["t"][-1]*1e3,edges[0],edges[-1]),
                           cmap="RdBu_r",vmin=-limit,vmax=limit,interpolation="nearest")
        ax[i].set(xlim=(START*1e3,10.8),ylabel=("coarse","TAIL refined")[i]+"\nbody s (mm)")
        for bnd in (-.6,.6):ax[i].axhline(bnd,color="black",ls="--",lw=.6)
    ax[-1].set_xlabel("physical time (ms)")
    fig.colorbar(image,ax=ax,label="whole-GC robot nodal axial force / body-s bin width (N/mm)")
    fig.suptitle("A14P5: whole-GC axial force localization; color limits shared")
    fig.savefig(ROOT/"A14P5_ROBOT_INTERFACE_AXIAL_FORCE_MAP.png",dpi=160);plt.close(fig)

    fig,ax=plt.subplots(3,1,figsize=(11,7.5),sharex=True,layout="constrained")
    for k,rname in enumerate(REGIONS):
        for i,f in enumerate(contact):
            ax[k].plot(f["t"]*1e3,f["regional"][:,k,0]/geos[i][0][k],color=COLORS[i],label=("coarse","TAIL refined")[i])
        ax[k].set_ylabel(rname+"\nwhole-GC F_s/A\n(N/mm2)");ax[k].set_xlim(START*1e3,10.8);ax[k].grid(alpha=.2);ax[k].legend(fontsize=8)
    ax[-1].set_xlabel("physical time (ms)")
    fig.suptitle("A14P5: whole-GC nodal force / physical exterior area (not fluid-only traction)")
    fig.savefig(ROOT/"A14P5_INTERFACE_FORCE_PER_AREA.png",dpi=160);plt.close(fig)

    # Identical world-coordinate b=0 section and shared pressure/velocity scales.
    times=(8.5,9.,9.5,10.,10.3,10.682)
    frames=[[int(np.argmin(abs(item[1]-t/1e3))) for t in times] for item in cel]
    entries=[]
    pipe=coarse["identity"]["initial_center_aba_mm"]-coarse["identity"]["initial_axial_shift_mm"]*coarse["s"]-coarse["identity"]["radial_offset_n_mm"]*coarse["n"]
    for i,(d,(z,tt,metric)) in enumerate(zip(data,cel)):
        cc=z["fluid_element_centroids_mm"]-pipe
        ss,nn,bb=cc@d["s"],cc@d["n"],cc@d["b"]
        mask=(abs(bb)<.066)&(ss>-5)&(ss<-1.1)
        conn=z["fluid_connectivity_index"]
        for j,k in enumerate(frames[i]):
            p=z["fluid_pressure"][k]
            evf=z["fluid_evf"][k]
            v=z["fluid_velocity_mm_s"][k].astype(float)[conn].mean(axis=1)@d["s"]
            entries.append((i,j,ss,nn,mask,p,v,evf,float(tt[k])))
    wet=np.concatenate([e[5][e[4]&(e[7]>.01)] for e in entries])
    vel=np.concatenate([e[6][e[4]&(e[7]>.01)] for e in entries])
    limits=(max(np.percentile(abs(wet),99),1e-8),max(np.percentile(abs(vel),99),1))
    fig,ax=plt.subplots(6,6,figsize=(17,12),sharex=True,sharey=True)
    outlines={}
    for i,d in enumerate(data):
        cf=np.load(d["case"]/"private/robot_contact_fields_private.npz")
        x0=cf["node_coordinates_mm"].astype(float)-d["identity"]["initial_center_aba_mm"]
        rp=np.load(d["case"]/"private/rp_history_private.npz")
        rt,ur=ledger.vectors(rp,"UR")
        for j,k in enumerate(frames[i]):
            t0=cel[i][1][k]
            x=d["identity"]["initial_center_aba_mm"]+ledger.interp(d["t"],d["u"],t0)+Rotation.from_rotvec(ledger.interp(rt,ur,t0)).apply(x0)
            x=x-pipe
            xy=np.column_stack((x@d["s"],x@d["n"]))
            outlines[i,j]=xy[ConvexHull(xy).vertices]
    for i,j,ss,nn,mask,p,v,evf,tt in entries:
        for k,(values,vmin,vmax,cmap) in enumerate(((p,-limits[0],limits[0],"RdBu_r"),(v,-limits[1],limits[1],"RdBu_r"),(evf,0,1,"viridis"))):
            row=3*i+k
            selected=mask if k==2 else mask&(evf>.01)
            im=ax[row,j].scatter(ss[selected],nn[selected],c=values[selected],vmin=vmin,vmax=vmax,cmap=cmap,
                                 marker="s",s=9,lw=0,rasterized=True)
            ax[row,j].add_patch(Polygon(outlines[i,j],facecolor="white",edgecolor="#333333",lw=.55,zorder=3))
            ax[row,j].set(xlim=(-5,-1.1),ylim=(-.75,.75))
            if j==0:ax[row,j].set_ylabel(("coarse","refined")[i]+" "+("p N/mm2","v_s mm/s","EVF")[k]+"\nn (mm)",fontsize=8)
            if row==0:ax[row,j].set_title(f"{tt*1e3:.3f} ms",fontsize=9)
            if row==5:ax[row,j].set_xlabel("+s (mm)")
        if j==5:
            for k in range(3):
                row=3*i+k
                # Row-specific legend, same bounds for both cases and all six times.
                fig.colorbar(ax[row,j].collections[0],ax=ax[row,:],shrink=.7,pad=.008)
    fig.suptitle("A14P5: fixed b=0 CEL coordinates and common row-pair scales; pressure is cell-centered water-phase proxy")
    fig.savefig(ROOT/"A14P5_COARSE_VS_REFINE_CEL_FIELDS.png",dpi=130);plt.close(fig)

    summary=[]
    for i,d in enumerate(data):
        f=contact[i];z,ft,pm=cel[i]; q=ledger.impulse(d,False,START,END)
        whole=np.interp(f["t"],d["ct"],d["forces"][1]+d["forces"][2])
        nodal=f["regional"][:,:,0].sum(axis=1)
        sl=(f["t"]>=START)&(f["t"]<=END)
        discrepancy=float(np.max(abs(whole[sl]-nodal[sl]))/max(np.max(abs(whole[sl])),1e-12))
        fluid=d["forces"][2]
        flu_mask=(d["ct"]>=START)&(d["ct"]<=END)
        field_mask=(ft>=START)&(ft<=END)
        peak_idx=np.flatnonzero(flu_mask)[np.argmin(fluid[flu_mask])]
        first=rolling_force_onset(d["ct"],fluid)
        # Output-cadence pressure metric departure is a descriptive threshold,
        # not a traction or causal onset. Reference is the first cycle-2 frame.
        first_field=np.flatnonzero(field_mask)[0]
        baseline_p=pm[first_field]
        p_depart=np.flatnonzero(field_mask & (abs(pm-baseline_p)>max(.5*abs(baseline_p),1e-4)))
        accel=np.gradient(d["vs"],d["t"])
        accel_first=onset(d["t"],accel,-1000.,30e-6)
        summary.append(dict(case=d["job"],impulse=q,areas_mm2=geos[i][0].tolist(),facets=geos[i][1].tolist(),
            node_count=f["node_count"].tolist(),whole_GC_nodal_vs_native_max_relative=discrepancy,
            trailing_50us_fluid_mean_below_minus_5e_minus_5N_ms=None if first is None else first*1e3,
            native_contact_sample_half_width_ms=float(np.median(np.diff(d["ct"]))*500),
            fluid_force_peak_negative_N=float(fluid[peak_idx]),fluid_force_peak_time_ms=float(d["ct"][peak_idx]*1e3),
            acceleration_below_minus_1000_mm_s2_onset_ms=None if accel_first is None else accel_first*1e3,
            pressure_metric_departure_50pct_or_1e_minus_4_ms=None if not len(p_depart) else float(ft[p_depart[0]]*1e3),
            field_sample_half_width_ms=float(np.median(np.diff(ft))*500),
            pressure_metric_at_9p5ms_N_mm2=float(np.interp(.0095,ft[np.isfinite(pm)],pm[np.isfinite(pm)])),
            velocity_at_end_mm_s=float(np.interp(END,d["t"],d["vs"])),
            motion_at_end_mm=float(np.interp(END,d["t"],d["disp"])-np.interp(START,d["t"],d["disp"]))))
    outcome={"proxy_RP_U_UR_V_VR_max_absolute_differences":proxy_error,"cases":summary,
        "state_match_count":len(states),"state_match_max_abs_fluid_force_delta_N":max((abs(r["fluid_force_delta_N"]) for r in states),default=None),
        "CEL_pressure_and_velocity_shared_color_limits":limits,
        "limitations":["Nodal force field is whole General Contact, not pair-isolated robot-CEL.",
            "Available CPRESS CSV contains only whole-GC frame maxima and counts; no location, pair identity or patch area.",
            "CEL pressure is cell-centered water-phase pressure; no validated interface mapping or surface shear traction.",
            "Local CV boundary viscous stress, moving solid cut-cell momentum transfer and face flux are absent; a pressure-only CV is not an independent fluid reaction."]}
    (ROOT/"A14P5_ROBOT_CEL_INTERFACE_SENSITIVITY.json").write_text(json.dumps(outcome,indent=2)+"\n",encoding="ascii")
    print(json.dumps(outcome,indent=2))


if __name__ == "__main__":
    main()
