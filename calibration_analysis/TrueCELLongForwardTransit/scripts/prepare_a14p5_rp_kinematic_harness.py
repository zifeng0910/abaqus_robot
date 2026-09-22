"""Audit the invalid replay and prepare one robot-only V/VR kinematic harness."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation, Slerp

ROOT = Path(__file__).resolve().parents[1]
INVALID = "TRUECEL_A14P5_COARSETRAJ_CEL_HREFINE_REPLAY"
TARGET = "TRUECEL_A14P5_DIAG_2CYCLES"
JOB = "A14P5_RP_KINEMATIC_REPLAY_VALIDATION"
CASE = ROOT / "case" / JOB
STOP = 0.0108
CRIT_LO, CRIT_HI = 1/120, 0.010682020978413481


def vectors(z, stem):
    t = z[f"{stem}1"][:, 0]
    return t, np.column_stack([z[f"{stem}{i}"][:, 1] for i in (1,2,3)])


def interp(z, stem, t):
    st = z[f"{stem}1"][:, 0]
    return np.column_stack([np.interp(t, st, z[f"{stem}{i}"][:, 1]) for i in (1,2,3)])


def metric(t, e, meaningful):
    a = np.abs(e)
    above = np.flatnonzero(a > meaningful)
    imax = int(np.argmax(a))
    return {"rms": float(np.sqrt(np.mean(e*e))), "p95": float(np.percentile(a,95)),
            "p99": float(np.percentile(a,99)), "max": float(a[imax]),
            "max_time_ms": float(t[imax]*1e3),
            "first_meaningful_divergence_ms": None if not len(above) else float(t[above[0]]*1e3),
            "meaningful_threshold": meaningful}


def block(text, pattern):
    match = re.search(pattern, text)
    if not match: raise RuntimeError(f"Missing deck block: {pattern}")
    return match.group(0)


def main():
    if CASE.exists(): raise RuntimeError(f"Harness exists; refusing overwrite: {CASE}")
    invalid_case = ROOT/"case"/INVALID
    target_case = ROOT/"case"/TARGET
    invalid = np.load(invalid_case/"private"/"rp_history_authoritative_private.npz")
    target = np.load(target_case/"private"/"rp_history_private.npz")
    t = invalid["U1"][:,0]
    mask = t <= STOP
    t = t[mask]
    actual = {stem: interp(invalid,stem,t) for stem in ("U","UR","V","VR")}
    desired = {stem: interp(target,stem,t) for stem in ("U","UR","V","VR")}
    component = {}
    thresholds = {"U":1e-6,"UR":np.deg2rad(.01),"V":.02,"VR":.2}
    component_stems = ("U","UR","V","VR")
    for stem in component_stems:
        component[stem] = {str(i+1):metric(t,actual[stem][:,i]-desired[stem][:,i],thresholds[stem]) for i in range(3)}
    orientation = np.rad2deg((Rotation.from_rotvec(actual["UR"]).inv()*Rotation.from_rotvec(desired["UR"])).magnitude())
    crossings = {}
    for level in (.01,.05,.10,.50,1.00):
        idx=np.flatnonzero(orientation>level)
        if len(idx):
            i=int(idx[0]); crossings[str(level)]={"time_ms":float(t[i]*1e3),"target_UR":desired["UR"][i].tolist(),
                "actual_UR":actual["UR"][i].tolist(),"target_VR":desired["VR"][i].tolist(),"actual_VR":actual["VR"][i].tolist(),
                "amplitude_R":desired["UR"][i].tolist(),"step_time_s":float(t[i]),"total_time_s":float(t[i])}

    # Current componentwise 5-us rotvec interpolation versus orientation SLERP.
    native_t,_=vectors(target,"UR"); native_ur=interp(target,"UR",native_t)
    knots=np.r_[np.arange(0,STOP,5e-6),STOP]
    kur=interp(target,"UR",knots)
    sample=t[(t<=STOP)]
    linear=np.column_stack([np.interp(sample,knots,kur[:,i]) for i in range(3)])
    slerp=Slerp(knots,Rotation.from_rotvec(kur))(sample)
    interpolation_error=np.rad2deg((Rotation.from_rotvec(linear).inv()*slerp).magnitude())

    # Target VR convention: log-map angular velocity from successive orientations.
    mt=(native_t[:-1]+native_t[1:])/2; dt=np.diff(native_t); R=Rotation.from_rotvec(native_ur)
    spatial=(R[1:]*R[:-1].inv()).as_rotvec()/dt[:,None]
    body=(R[:-1].inv()*R[1:]).as_rotvec()/dt[:,None]
    native_vr=interp(target,"VR",mt)
    valid=(mt>1e-4)&(mt<STOP)
    def omega_stats(w):
        e=np.linalg.norm(w[valid]-native_vr[valid],axis=1)
        return {"rms_rad_s":float(np.sqrt(np.mean(e*e))),"p99_rad_s":float(np.percentile(e,99)),"max_rad_s":float(e.max())}
    convention={"spatial_global":omega_stats(spatial),"body_material":omega_stats(body)}
    convention["selected"]="GLOBAL_SPATIAL"

    # Figures for the invalid replay and target consistency.
    for name,stems,path in (("U/UR",("U","UR"),"A14P5_INVALID_REPLAY_U_UR_ERROR.png"),
                            ("V/VR",("V","VR"),"A14P5_INVALID_REPLAY_V_VR_ERROR.png")):
        fig,axs=plt.subplots(2,1,figsize=(11,7),sharex=True,layout="constrained")
        for ax,stem in zip(axs,stems):
            for i in range(3): ax.plot(t*1e3,actual[stem][:,i]-desired[stem][:,i],label=f"{stem}{i+1}")
            ax.set_ylabel(f"{stem} error");ax.grid(alpha=.2);ax.legend(ncol=3)
        axs[-1].set_xlabel("time (ms)");fig.suptitle(f"Invalid replay component errors: {name}")
        fig.savefig(ROOT/path,dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,4.5),layout="constrained");ax.plot(t*1e3,orientation)
    for level in (.01,.05,.1,.5,1):ax.axhline(level,ls="--",lw=.7)
    ax.set(xlabel="time (ms)",ylabel="geodesic orientation error (deg)");ax.grid(alpha=.2)
    fig.savefig(ROOT/"A14P5_ROTATION_MATRIX_ERROR_HISTORY.png",dpi=160);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(11,7),sharex=True,layout="constrained")
    for label,w in (("spatial/global",spatial),("body/material",body)):
        axs[0].plot(mt[valid]*1e3,np.linalg.norm(w[valid]-native_vr[valid],axis=1),label=label)
    axs[0].set_ylabel("|omega_from_R - VR| (rad/s)");axs[0].legend();axs[0].grid(alpha=.2)
    axs[1].plot(sample*1e3,interpolation_error);axs[1].set(xlabel="time (ms)",ylabel="rotvec-linear vs SLERP (deg)");axs[1].grid(alpha=.2)
    fig.savefig(ROOT/"A14P5_TARGET_UR_VR_CONSISTENCY.png",dpi=160);plt.close(fig)

    # Inspect and record exact submitted boundary semantics.
    submitted=(invalid_case/(INVALID+".inp")).read_text(encoding="latin1")
    cards=[]
    for dof,name in zip(range(1,7),("REPLAY_U1","REPLAY_U2","REPLAY_U3","REPLAY_R1","REPLAY_R2","REPLAY_R3")):
        pat=rf"\*Amplitude, name={name}, definition=TABULAR, time=STEP TIME\n(.*?)\*Boundary, amplitude={name}, type=DISPLACEMENT\nRP_ROBOT, {dof}, {dof}, 1\."
        m=re.search(pat,submitted,re.S)
        if not m:raise RuntimeError(f"Missing submitted card {name}")
        first=m.group(1).splitlines()[0];last=m.group(1).splitlines()[-1]
        cards.append({"dof":dof,"amplitude":name,"reference_magnitude":1.0,"type":"DISPLACEMENT","op":"default MOD",
                      "time_basis":"STEP TIME","initial_pair":first,"final_pair":last,
                      "abaqus_received_value":"1.0 x tabular amplitude"})

    # Exact robot-only deck. V and VR are prescribed from native target histories.
    source=(invalid_case/(INVALID+".inp")).read_text(encoding="latin1")
    robot_part=block(source,r"(?ms)^\*Part, name=Robot_SOLID\s*$.*?^\*End Part\s*$")
    robot_instance=block(source,r"(?ms)^\*Instance, name=Robot_SOLID-1, part=Robot_SOLID\s*$.*?^\*End Instance\s*$")
    robot_material=block(source,r"(?ms)^\*Material, name=MAT_ROBOT_RIGID\s*$.*?(?=^\*\*\s*$)").rstrip()
    ident=json.loads((invalid_case/"case_identity.json").read_text())
    rp0=np.asarray(ident["initial_center_aba_mm"])
    amp=[]
    target_t=target["V1"][:,0]; use=(target_t<=STOP); target_t=target_t[use]
    for stem,prefix,offset in (("V","TARGET_V",1),("VR","TARGET_VR",4)):
        vals=np.column_stack([target[f"{stem}{i}"][:,1][use] for i in (1,2,3)])
        for axis in range(3):
            amp.append(f"*Amplitude, name={prefix}{axis+1}, definition=TABULAR, time=STEP TIME")
            amp.extend(f"{x:.16g}, {y:.16g}" for x,y in zip(target_t,vals[:,axis]))
            amp.append(f"*Boundary, amplitude={prefix}{axis+1}, type=VELOCITY")
            amp.append(f"RP_ROBOT, {offset+axis}, {offset+axis}, 1.")
    deck=(f"*Heading\n{JOB}: robot-only prescribed global V/VR harness\n**\n{robot_part}\n**\n*Assembly, name=Assembly\n"
          f"{robot_instance}\n*Elset, elset=ROBOT_SOLID_CEL_ALL, instance=Robot_SOLID-1, generate\n1, 7302, 1\n"
          f"*Node\n2, {rp0[0]:.15g}, {rp0[1]:.15g}, {rp0[2]:.15g}\n*Nset, nset=RP_ROBOT\n2,\n"
          f"*Rigid Body, ref node=RP_ROBOT, elset=ROBOT_SOLID_CEL_ALL\n*End Assembly\n**\n{robot_material}\n**\n"
          f"*Step, name=Step_Kinematic, nlgeom=YES\n*Dynamic, Explicit, DIRECT USER CONTROL\n1e-6, {STOP}\n"
          +"\n".join(amp)+"\n*Output, history, time interval=5e-6\n*Node Output, nset=RP_ROBOT\nU, UR, V, VR\n*End Step\n")
    CASE.mkdir(parents=True); inp=CASE/(JOB+".inp");inp.write_text(deck,encoding="latin1")
    audit={"invalid_classification":"FIXED_TRAJECTORY_REPLAY_INVALID","component_errors":component,
           "orientation_threshold_crossings":crossings,"orientation_max_deg":float(orientation.max()),
           "componentwise_rotvec_interpolation_vs_SLERP_max_deg":float(interpolation_error.max()),
           "target_UR_VR_consistency":convention,"submitted_boundary_cards":cards,
           "initial_target_U":interp(target,"U",np.array([0.]))[0].tolist(),"initial_target_UR":interp(target,"UR",np.array([0.]))[0].tolist(),
           "harness_method":"global translational V plus global/spatial angular VR; no CEL, wall, contact, fluid, magnetic load or external load",
           "harness_input_sha256":hashlib.sha256(inp.read_bytes()).hexdigest().upper()}
    (ROOT/"A14P5_PRESCRIBED_MOTION_PREHARNESS_AUDIT.json").write_text(json.dumps(audit,indent=2)+"\n")
    print(json.dumps({"job":JOB,"orientation_max_deg":audit["orientation_max_deg"],
        "rotvec_linear_vs_slerp_max_deg":audit["componentwise_rotvec_interpolation_vs_SLERP_max_deg"],
        "VR_convention":convention,"input_sha256":audit["harness_input_sha256"]},indent=2))


if __name__=="__main__":main()
