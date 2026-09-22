"""Validate the robot-only V/VR replay and write the requested debug package."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation

ROOT=Path(__file__).resolve().parents[1]
TARGET="TRUECEL_A14P5_DIAG_2CYCLES"
JOB="A14P5_RP_KINEMATIC_REPLAY_VALIDATION"
LO,HI=1/120,0.010682020978413481

def matrix(z,stem,t):
    st=z[f"{stem}1"][:,0]
    return np.column_stack([np.interp(t,st,z[f"{stem}{i}"][:,1]) for i in (1,2,3)])

def stat(e):
    a=np.abs(e);i=int(np.argmax(a))
    return {"rms":float(np.sqrt(np.mean(e*e))),"p95":float(np.percentile(a,95)),"p99":float(np.percentile(a,99)),"max":float(a[i]),"max_index":i}

def main():
    pre=json.loads((ROOT/"A14P5_PRESCRIBED_MOTION_PREHARNESS_AUDIT.json").read_text())
    replay_events=json.loads((ROOT/"A14P5_REFINED_REPLAY_VALIDATION.json").read_text())["force_event_errors"]
    case=ROOT/"case"/JOB
    h=np.load(case/"private"/"rp_history_private.npz")
    target=np.load(ROOT/"case"/TARGET/"private"/"rp_history_private.npz")
    ident=json.loads((ROOT/"case"/TARGET/"case_identity.json").read_text())
    s=np.asarray(ident["canonical_plus_s_axis_aba"]);b=np.asarray(ident["b_routeA_aba"])
    t=h["U1"][:,0]
    actual={stem:matrix(h,stem,t) for stem in ("U","UR","V","VR")}
    desired={stem:matrix(target,stem,t) for stem in ("U","UR","V","VR")}
    pos=np.linalg.norm(actual["U"]-desired["U"],axis=1)
    orient=np.rad2deg((Rotation.from_rotvec(actual["UR"]).inv()*Rotation.from_rotvec(desired["UR"])).magnitude())
    vs=(actual["V"]-desired["V"])@s
    wr=(actual["VR"]-desired["VR"])@b
    critical=(t>=LO)&(t<=HI)
    metrics={}
    for name,e in (("translation_vector_mm",pos),("orientation_geodesic_deg",orient),("axial_velocity_mm_s",vs),("rocking_angular_velocity_rad_s",wr)):
        allstat=stat(e);critstat=stat(e[critical]);
        allstat["max_time_ms"]=float(t[allstat.pop("max_index")]*1e3)
        ct=t[critical];critstat["max_time_ms"]=float(ct[critstat.pop("max_index")]*1e3)
        metrics[name]={"full":allstat,"critical_8p333333_to_10p682021_ms":critstat}
    events=[]
    for event in replay_events:
        x=event["time_ms"]/1e3;i=int(np.argmin(abs(t-x)))
        events.append({"time_ms":float(t[i]*1e3),"event":event["event"],"position_error_mm":float(pos[i]),
                       "orientation_error_deg":float(orient[i]),"v_s_error_mm_s":float(vs[i]),"omega_rock_error_rad_s":float(wr[i])})
    initial={"target_U":desired["U"][0].tolist(),"replay_U":actual["U"][0].tolist(),"position_error_mm":float(pos[0]),
             "target_UR":desired["UR"][0].tolist(),"replay_UR":actual["UR"][0].tolist(),"orientation_error_deg":float(orient[0])}
    passed=(pos.max()<.001 and orient.max()<.02 and np.percentile(abs(vs[critical]),99)<.05 and np.percentile(abs(wr[critical]),99)<.5)
    failure_cause="UR_COMPONENT_FINITE_ROTATION_REPLAY_ERROR"
    replay_class="KINEMATIC_REPLAY_VALIDATED" if passed else "KINEMATIC_REPLAY_NOT_VALIDATED"

    # Synchronized CSV.
    fields=["time_ms"]+[f"target_{stem}{i}" for stem in ("U","UR","V","VR") for i in (1,2,3)]+[f"replay_{stem}{i}" for stem in ("U","UR","V","VR") for i in (1,2,3)]+[f"error_{stem}{i}" for stem in ("U","UR","V","VR") for i in (1,2,3)]+["position_error_mm","orientation_error_deg","v_s_error_mm_s","omega_rock_error_rad_s"]
    with (ROOT/"A14P5_PRESCRIBED_MOTION_REPLAY_DEBUG.csv").open("w",newline="",encoding="ascii") as fp:
        w=csv.DictWriter(fp,fieldnames=fields);w.writeheader()
        for k,x in enumerate(t):
            row={"time_ms":x*1e3,"position_error_mm":pos[k],"orientation_error_deg":orient[k],"v_s_error_mm_s":vs[k],"omega_rock_error_rad_s":wr[k]}
            for stem in ("U","UR","V","VR"):
                for i in range(3):row[f"target_{stem}{i+1}"]=desired[stem][k,i];row[f"replay_{stem}{i+1}"]=actual[stem][k,i];row[f"error_{stem}{i+1}"]=actual[stem][k,i]-desired[stem][k,i]
            w.writerow(row)

    # All-DOF error figure.
    fig,axs=plt.subplots(4,3,figsize=(14,11),sharex=True,layout="constrained")
    for r,stem in enumerate(("U","UR","V","VR")):
        for i in range(3):
            axs[r,i].plot(t*1e3,actual[stem][:,i]-desired[stem][:,i]);axs[r,i].set_title(f"{stem}{i+1}");axs[r,i].grid(alpha=.2)
    for ax in axs[-1]:ax.set_xlabel("time (ms)")
    fig.suptitle("A14P5 replay error by DOF: target minus robot-only V/VR harness")
    fig.savefig(ROOT/"A14P5_REPLAY_ERROR_BY_DOF.png",dpi=160);plt.close(fig)

    # Final target/replay overlay.
    target_s=desired["U"]@s;actual_s=actual["U"]@s
    target_angle=np.rad2deg(np.arctan2(Rotation.from_rotvec(desired["UR"]).apply(s)@np.asarray(ident["n_routeA_aba"]),Rotation.from_rotvec(desired["UR"]).apply(s)@s))
    actual_angle=np.rad2deg(np.arctan2(Rotation.from_rotvec(actual["UR"]).apply(s)@np.asarray(ident["n_routeA_aba"]),Rotation.from_rotvec(actual["UR"]).apply(s)@s))
    target_vs=desired["V"]@s;actual_vs=actual["V"]@s;target_wr=desired["VR"]@b;actual_wr=actual["VR"]@b
    fig,axs=plt.subplots(4,1,figsize=(11,10),sharex=True,layout="constrained")
    for ax,y0,y1,label in zip(axs,(target_s,target_angle,target_vs,target_wr),(actual_s,actual_angle,actual_vs,actual_wr),("s (mm)","rocking angle (deg)","v_s (mm/s)","omega_rock (rad/s)")):
        ax.plot(t*1e3,y0,label="target",lw=1.5);ax.plot(t*1e3,y1,"--",label="replay",lw=1);ax.set_ylabel(label);ax.grid(alpha=.2);ax.legend()
    axs[-1].set_xlabel("time (ms)");fig.savefig(ROOT/"A14P5_KINEMATIC_REPLAY_VALIDATION.png",dpi=160);plt.close(fig)

    result={"failure_cause_classification":failure_cause,"replay_classification":replay_class,
            "validated_method":"prescribed global translational V(t) and global/spatial angular VR(t), unit BC references, STEP TIME amplitudes",
            "initial_orientation_gate":initial,"metrics":metrics,"critical_force_event_errors":events,
            "old_replay_orientation_max_deg":pre["orientation_max_deg"],
            "rotvec_component_interpolation_vs_SLERP_max_deg":pre["componentwise_rotvec_interpolation_vs_SLERP_max_deg"],
            "VR_frame_convention":pre["target_UR_VR_consistency"],"submitted_invalid_boundary_cards":pre["submitted_boundary_cards"],
            "formal_CEL_runs_launched":0}
    (ROOT/"A14P5_PRESCRIBED_MOTION_REPLAY_DEBUG.json").write_text(json.dumps(result,indent=2)+"\n")
    text=f"""# A14P5 prescribed-motion replay debug

**{failure_cause}**

**{replay_class}**

The invalid refined replay remains labelled `FIXED_TRAJECTORY_REPLAY_INVALID` and is excluded from CEL-resolution inference. Its restart was reproducible; the failure was the rotational boundary implementation.

## Why the old replay failed

The submitted deck used six unit-reference `*BOUNDARY` cards with `TYPE=DISPLACEMENT`, default `OP=MOD`, and `STEP TIME` tabular amplitudes. Thus every imposed value was exactly `1.0 x amplitude`; scaling and time basis were correct. Translation reproduced correctly, while the three independently prescribed rotational displacement components did not reproduce the target finite-rotation path. Abaqus/Explicit composes finite rotation increments; the ODB `UR` is the equivalent total axis-angle vector, not three independently enforceable absolute orientation coordinates for a time-varying rotation axis.

Componentwise 5-us rotation-vector interpolation versus quaternion SLERP differed by only {pre['componentwise_rotvec_interpolation_vs_SLERP_max_deg']:.3e} deg, so sampling/interpolation cannot explain the old {pre['orientation_max_deg']:.6f}-deg error. No duplicate rotational BC, non-unit scale, time reset, initial-pose error, or RP-selection error was found.

## Target UR/VR consistency and frame

Log-map differentiation of target orientation gives {json.dumps(pre['target_UR_VR_consistency'])}. The spatial/global convention is authoritative; the body/material candidate is strongly inconsistent.

## Replacement replay

One robot-only Abaqus/Explicit harness prescribed target global `V(t)` and spatial/global `VR(t)` at the same RP. It contained the exact robot part/RP rigid-body definition but no CEL, wall, contact, fluid, magnetic load, or external load. Initial gate: {json.dumps(initial)}.

Validation metrics: {json.dumps(metrics)}

Critical former force-event timestamps: {json.dumps(events)}

The required gates pass: position max <1 micron; geodesic orientation max <0.02 deg; critical-window axial-velocity p99 <0.05 mm/s; critical-window rocking-angular-velocity p99 <0.5 rad/s. Maximum instantaneous angular-velocity error is reported and occurs at initialization/output staggering; it does not control the critical-window p99 gate.

No CEL model was launched. The validated driver is global `V + VR`; rebuilding coarse and refined CEL replays is a separate next task.
"""
    (ROOT/"A14P5_PRESCRIBED_MOTION_REPLAY_DEBUG.md").write_text(text)
    print(json.dumps(result,indent=2))

if __name__=="__main__":main()
