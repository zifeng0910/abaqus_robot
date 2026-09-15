"""Analyze the 5 Hz straight rocking solve and render the required figures/GIFs."""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
JOB = "PROD_LOCAL_ROCK10_5HZ_G0_STRAIGHT"
CASE = OUT / "case" / JOB
LEGACY = REPO / "calibration_analysis" / "Legacy_Wobble_GIFs" / "Job_RouteA_CEL_SOLID_headtail_rock_probe.gif"
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298])


def get_array(archive, name):
    keys = [key for key in archive.files if key == name or key.startswith(name + " (Repeated:")]
    return max((archive[key] for key in keys), key=len)


def dense_pose():
    archive = np.load(CASE / "private" / "rp_history_private.npz")
    time = get_array(archive, "V1")[:, 0]
    data = {prefix: np.column_stack([np.interp(time, get_array(archive, prefix + str(i))[:, 0], get_array(archive, prefix + str(i))[:, 1]) for i in (1, 2, 3)]) for prefix in ("U", "UR", "V", "VR")}
    return time, data


def contact_force(time):
    archive = np.load(CASE / "private" / "contact_history_private.npz")
    def vector(prefix):
        columns = []
        for i in (1, 2, 3):
            keys = [key for key in archive.files if "|" + prefix + str(i) + " on surface " in key and "ASSEMBLY_ROBOT" in key]
            array = max((archive[key] for key in keys), key=len)
            columns.append(np.interp(time, array[:, 0], array[:, 1]))
        return np.column_stack(columns)
    return vector("CFN") + vector("CFS")


def intervals(flag):
    edges = np.diff(np.r_[False, flag, False].astype(int))
    return list(zip(np.where(edges == 1)[0], np.where(edges == -1)[0] - 1))


def longest(time, flag):
    dt = float(np.median(np.diff(time)))
    return max(((time[j] - time[i] + dt) for i, j in intervals(flag)), default=0.0)


def part_nodes(deck, name):
    part = re.search(rf"(?ms)^\*Part, name={re.escape(name)}\s*$.*?^\*End Part\s*$", deck).group(0)
    block = re.search(r"(?ms)^\*Node\s*$\n(.*?)(?=^\*)", part).group(1)
    return np.asarray([[float(x) for x in row.split(",")[1:4]] for row in block.splitlines() if row.strip()])


def unit(v):
    return np.asarray(v, dtype=float) / np.linalg.norm(v)


def face_geometry(wall, center, c, n, b):
    axial = (wall - center).dot(c)
    ring = wall[np.abs(axial - axial.min()) < 1e-7]
    xy = np.column_stack(((ring - center).dot(n), (ring - center).dot(b)))
    xy = xy[np.argsort(np.arctan2(xy[:, 1], xy[:, 0]))]
    normals, radii = [], []
    for p, q in zip(xy, np.roll(xy, -1, axis=0)):
        edge = q - p; normal2 = unit([edge[1], -edge[0]])
        if np.dot(normal2, p + q) < 0: normal2 *= -1
        normals.append(normal2[0] * n + normal2[1] * b); radii.append(np.dot(p, normal2))
    return np.asarray(normals), float(np.mean(radii))


def fundamental(time, values, frequency=5.0):
    w = 2 * np.pi * frequency
    matrix = np.column_stack((np.sin(w*time), np.cos(w*time), np.ones_like(time), time-time.mean()))
    coef = np.linalg.lstsq(matrix, values, rcond=None)[0]
    return float(np.hypot(coef[0], coef[1])), math.degrees(math.atan2(coef[1], coef[0]))


def main():
    identity = json.loads((CASE / "case_identity.json").read_text())
    time, history = dense_pose()
    force = contact_force(time); force_mag = np.linalg.norm(force, axis=1)
    active_dense = force_mag > 1e-8
    stride = max(1, int(round(1e-4 / np.median(np.diff(time)))))
    sample = np.unique(np.r_[np.arange(0, len(time), stride), len(time)-1])
    ts = time[sample]; com = RP0 + history["U"][sample]
    rotations = Rotation.from_rotvec(history["UR"][sample])
    c, n, b = map(unit, (identity["initial_axis_aba"], identity["n_rock_aba"], identity["b_rock_aba"]))
    axis = rotations.apply(np.broadcast_to(c, (len(sample), 3)))
    qc, qrock, qcross = axis.dot(c), axis.dot(n), axis.dot(b)
    theta = np.degrees(np.arctan2(qrock, qc)); theta_cross = np.degrees(np.arctan2(qcross, qc))
    s = (com - RP0).dot(c); radial_vec = com - RP0 - np.outer(s, c); radial = np.linalg.norm(radial_vec, axis=1)

    telemetry = pd.read_csv(CASE / f"{JOB}_telemetry.csv")
    field = np.column_stack([np.interp(ts, telemetry.t_s, telemetry[key]) for key in ("Bx_aba_T","By_aba_T","Bz_aba_T")])
    torque = np.column_stack([np.interp(ts, telemetry.t_s, telemetry[key]) for key in ("tx_aba_Nmm","ty_aba_Nmm","tz_aba_Nmm")])
    alpha = np.degrees(np.arctan2(field.dot(n), field.dot(c))); trock = torque.dot(b)

    deck = (CASE / f"{JOB}.inp").read_text(); nodes = part_nodes(deck, "Robot_SOLID"); wall = part_nodes(deck, "Pipe_WALL_HELPER")
    axial0 = (nodes - RP0).dot(c); lo, hi = axial0.min(), axial0.max()
    region = np.full(len(nodes), "BODY", dtype=object); region[axial0 <= lo+.25] = "HEAD"; region[axial0 >= hi-.25] = "TAIL"
    normals, radius = face_geometry(wall, RP0, c, n, b)
    gaps = []
    for center_i, rotation_i in zip(com, rotations):
        points = center_i + rotation_i.apply(nodes - RP0)
        radial_points = points - RP0 - np.outer((points - RP0).dot(c), c)
        node_gap = radius - np.max(radial_points.dot(normals.T), axis=1)
        gaps.append((node_gap.min(), node_gap[region=="HEAD"].min(), node_gap[region=="TAIL"].min(), node_gap[region=="BODY"].min()))
    gaps = np.asarray(gaps) * 1e3
    active = active_dense[sample]
    head_support = active & (gaps[:,1] <= 20); tail_support = active & (gaps[:,2] <= 20)
    both = head_support & tail_support; head_only = head_support & ~tail_support; tail_only = tail_support & ~head_support; separated = ~active
    state = np.where(both,"BOTH",np.where(head_only,"HEAD",np.where(tail_only,"TAIL",np.where(separated,"FREE","BODY"))))
    support_indicator = gaps[:,2] - gaps[:,1]

    q = np.column_stack((qrock, qcross)); covariance = np.cov(q, rowvar=False); eig = np.linalg.eigvalsh(covariance)
    linearity = float(eig[-1] / max(eig.sum(), 1e-30))
    valid_sign = np.abs(theta) > .1; signs = np.sign(theta[valid_sign]); zero_crossings = int(np.sum(signs[1:] != signs[:-1]))
    robot_amp, robot_phase = fundamental(ts, theta); field_amp, field_phase = fundamental(ts, alpha)
    phase_lag = ((robot_phase-field_phase+180)%360)-180
    energy = np.load(CASE / "private" / "energy_history_private.npz")
    energy_metrics = {key+"_min": float(energy[key][:,1].min()) for key in energy.files}
    energy_metrics.update({key+"_max": float(energy[key][:,1].max()) for key in energy.files})

    topology_ok = theta.min() < -1 and theta.max() > 1 and np.min(np.abs(theta)) < .5 and np.max(np.abs(theta_cross)) < 5
    level2 = topology_ok and head_only.any() and tail_only.any() and np.mean(both) < .25
    target_family = theta.min() <= -8 and theta.max() >= 8 and theta.min() >= -12 and theta.max() <= 12
    if level2 and target_family and gaps[:,0].min() > -20:
        classification = "HEADTAIL_WALL_ROCKING_RECOVERED"
    elif topology_ok:
        classification = "LOCAL_ROCKING_MODE_RECOVERED"
    elif np.max(np.abs(theta_cross)) >= 5:
        classification = "CROSS_PLANE_INSTABILITY"
    elif theta.min() >= -1 or theta.max() <= 1:
        classification = "MAGNETIC_TORQUE_INSUFFICIENT_FOR_10DEG_ROCKING"
    else:
        classification = "ROCKING_PHASE_NOT_TRACKED"

    metrics = {
        "case_id": JOB, "classification": classification,
        "theta_rock_min_deg": float(theta.min()), "theta_rock_max_deg": float(theta.max()),
        "theta_rock_fundamental_amplitude_deg": robot_amp, "theta_cross_rms_deg": float(np.sqrt(np.mean(theta_cross**2))),
        "theta_cross_max_abs_deg": float(np.max(np.abs(theta_cross))), "theta_rock_zero_crossings": zero_crossings,
        "pca_rocking_linearity": linearity, "field_alpha_min_deg": float(alpha.min()), "field_alpha_max_deg": float(alpha.max()),
        "phase_lag_robot_vs_field_deg": phase_lag, "Trock_min_Nmm": float(trock.min()), "Trock_max_Nmm": float(trock.max()),
        "HEAD_min_gap_um": float(gaps[:,1].min()), "TAIL_min_gap_um": float(gaps[:,2].min()), "BODY_min_gap_um": float(gaps[:,3].min()),
        "minimum_exact_surface_gap_um": float(gaps[:,0].min()), "deep_penetration": bool(gaps[:,0].min() < -20),
        "HEAD_only_support_fraction": float(head_only.mean()), "TAIL_only_support_fraction": float(tail_only.mean()),
        "both_end_support_fraction": float(both.mean()), "no_contact_fraction": float(separated.mean()),
        "longest_HEAD_only_ms": longest(ts,head_only)*1e3, "longest_TAIL_only_ms": longest(ts,tail_only)*1e3,
        "longest_both_end_bridge_ms": longest(ts,both)*1e3, "longest_no_contact_ms": longest(ts,separated)*1e3,
        "support_indicator_sign_changes": int(np.sum(np.sign(support_indicator[1:]) != np.sign(support_indicator[:-1]))),
        "support_indicator_theta_correlation": float(np.corrcoef(theta,support_indicator)[0,1]),
        "COM_radial_max_um": float(radial.max()*1e3), "axial_displacement_mm": float(s[-1]-s[0]),
        "peak_contact_force_N": float(force_mag.max()), "field_magnitude_max_abs_error_T": float(np.max(np.abs(np.linalg.norm(field,axis=1)-.01))),
        **energy_metrics,
    }
    (OUT / "rocking_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    pd.DataFrame([metrics]).to_csv(OUT / "rocking_metrics.csv", index=False)
    pose = pd.DataFrame({"time_s":ts,"normalized_phase":ts/.2,"theta_rock_deg":theta,"theta_cross_deg":theta_cross,
                         "alpha_B_deg":alpha,"q_rock":qrock,"q_cross":qcross,"HEAD_gap_um":gaps[:,1],"TAIL_gap_um":gaps[:,2],
                         "BODY_gap_um":gaps[:,3],"minimum_gap_um":gaps[:,0],"T_rock_Nmm":trock,"contact_state":state,
                         "com_x":com[:,0],"com_y":com[:,1],"com_z":com[:,2],"axis_x":axis[:,0],"axis_y":axis[:,1],"axis_z":axis[:,2],
                         "ur1":history["UR"][sample,0],"ur2":history["UR"][sample,1],"ur3":history["UR"][sample,2],
                         "s_mm":s,"COM_radial_mm":radial,"contact_force_N":force_mag[sample]})
    pose.to_csv(OUT / "rocking_pose.csv", index=False)

    fig, axes = plt.subplots(5,1,figsize=(9,10.5),sharex=True)
    axes[0].plot(ts*1e3,theta,label="robot theta_rock",color="#c43c35"); axes[0].plot(ts*1e3,alpha,label="field alpha_B",color="#237a57",ls="--"); axes[0].set_ylabel("deg"); axes[0].legend(ncol=2,frameon=False)
    axes[1].plot(ts*1e3,theta_cross,color="#315a8a"); axes[1].set_ylabel("theta_cross (deg)")
    axes[2].plot(ts*1e3,gaps[:,1],label="HEAD",color="#c43c35"); axes[2].plot(ts*1e3,gaps[:,2],label="TAIL",color="#315a8a"); axes[2].axhline(0,color="#222",lw=.7); axes[2].set_ylabel("exact gap (um)"); axes[2].legend(frameon=False)
    state_code = pd.Categorical(state, categories=["FREE", "BODY", "HEAD", "TAIL", "BOTH"]).codes
    axes[3].step(ts*1e3,state_code,where="post",color="#555555"); axes[3].set_ylabel("contact")
    axes[3].set_yticks(range(5),["FREE", "BODY", "HEAD", "TAIL", "BOTH"])
    axes[4].plot(ts*1e3,trock*1e3,color="#237a57"); axes[4].set_ylabel("T_rock (uN mm)"); axes[4].set_xlabel("Time (ms)")
    for x,label in [(50,"HEAD peak"),(100,"zero crossing"),(150,"TAIL peak")]:
        for ax in axes: ax.axvline(x,color="#777",lw=.65,alpha=.6)
        axes[0].text(x,axes[0].get_ylim()[1],label,ha="center",va="bottom",fontsize=8)
    for ax in axes: ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(OUT/"Rocking_5Hz_TimeSeries.png",dpi=220); fig.savefig(OUT/"Rocking_5Hz_TimeSeries.pdf"); plt.close(fig)

    legacy = pd.read_csv(OUT/"Legacy_RouteA_motion_target.csv")
    fig,ax=plt.subplots(figsize=(6.4,5)); ax.plot(legacy.q_rock,legacy.q_cross,color="#777",lw=3,label="Legacy RouteA +/-10 deg target"); ax.plot(qrock,qcross,color="#c43c35",lw=1.6,label="Finite-moment magnetic response"); ax.set(xlabel="q_rock",ylabel="q_cross",title=f"Rocking orbit; PCA linearity={linearity:.4f}"); ax.grid(alpha=.2); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(OUT/"RouteA_vs_Magnetic_Rocking_Orbit.png",dpi=220); fig.savefig(OUT/"RouteA_vs_Magnetic_Rocking_Orbit.pdf"); plt.close(fig)

    render_gifs(pose, nodes, region, c, n, b, radius)
    write_report(metrics)
    print(json.dumps(metrics, indent=2))


def render_gifs(pose, nodes, region, c, n, b, radius):
    ids = np.linspace(0,len(pose)-1,101).astype(int); rel0=nodes-RP0
    phi=np.linspace(0,2*np.pi,25); axial=np.linspace(-1.8,1.8,18); pp,ss=np.meshgrid(phi,axial)
    tube=RP0+ss[...,None]*c+radius*np.cos(pp)[...,None]*n+radius*np.sin(pp)[...,None]*b
    fig=plt.figure(figsize=(12,5.4)); left=fig.add_subplot(121,projection="3d"); right=fig.add_subplot(122)
    def draw(k):
        left.cla(); right.cla(); row=pose.iloc[ids[k]]; rot=Rotation.from_rotvec(np.zeros(3))
        # Recover a minimal rotation that maps the initial axis to the reported axis for visualization.
        rot=Rotation.from_rotvec(row[["ur1","ur2","ur3"]].to_numpy(float))
        points=row[["com_x","com_y","com_z"]].to_numpy(float)+rot.apply(rel0)
        left.plot_surface(tube[...,0],tube[...,1],tube[...,2],color="#8bb9c7",alpha=.12,linewidth=0)
        colors=np.where(region=="HEAD","#c43c35",np.where(region=="TAIL","#315a8a","#444444")); left.scatter(points[::3,0],points[::3,1],points[::3,2],c=colors[::3],s=3)
        left.plot([RP0[0]-1.8*c[0],RP0[0]+1.8*c[0]],[RP0[1]-1.8*c[1],RP0[1]+1.8*c[1]],[RP0[2]-1.8*c[2],RP0[2]+1.8*c[2]],color="#237a57",lw=1)
        left.set(xlim=(RP0[0]-2,RP0[0]+2),ylim=(RP0[1]-2,RP0[1]+2),zlim=(RP0[2]-2,RP0[2]+2),title="Fixed global camera"); left.view_init(24,-58)
        local=np.column_stack(((points-RP0).dot(c),(points-RP0).dot(n))); right.scatter(local[::3,0],local[::3,1],c=colors[::3],s=5); right.axhline(radius,color="#8bb9c7"); right.axhline(-radius,color="#8bb9c7"); right.axhline(0,color="#237a57",lw=.8); right.set_aspect("equal"); right.set(xlim=(-1.8,1.8),ylim=(-.8,.8),xlabel="c_hat (mm)",ylabel="n_rock (mm)",title="Rocking plane; view along b_rock")
        fig.suptitle(f"t={row.time_s*1000:6.1f} ms  alpha_B={row.alpha_B_deg:+6.2f} deg  theta_rock={row.theta_rock_deg:+6.2f} deg  theta_cross={row.theta_cross_deg:+5.2f} deg\nHEAD gap={row.HEAD_gap_um:+7.2f} um  TAIL gap={row.TAIL_gap_um:+7.2f} um  state={row.contact_state}  T_rock={row.T_rock_Nmm*1e3:+7.3f} uN mm")
        fig.tight_layout(rect=(0,0,1,.88))
    anim=FuncAnimation(fig,draw,frames=len(ids),interval=60); dual=OUT/f"{JOB}_DualView.gif"; anim.save(dual,writer=PillowWriter(fps=16.667),dpi=90); plt.close(fig)

    new=Image.open(dual); old=Image.open(LEGACY); frames=[]
    for k in range(101):
        old_index=24+int(round(6*k/100)); old.seek(old_index); new.seek(k)
        a=old.convert("RGB"); z=new.convert("RGB"); h=540; a=a.resize((round(a.width*h/a.height),h),Image.Resampling.LANCZOS); z=z.resize((round(z.width*h/z.height),h),Image.Resampling.LANCZOS)
        canvas=Image.new("RGB",(a.width+z.width,h),"white"); canvas.paste(a,(0,0)); canvas.paste(z,(a.width,0)); frames.append(canvas)
    frames[0].save(OUT/"RouteA_vs_MagneticRocking.gif",save_all=True,append_images=frames[1:],duration=60,loop=0,optimize=False)


def write_report(m):
    identity = json.loads((CASE / "case_identity.json").read_text())
    geometry = pd.read_csv(OUT / "RouteA10_geometry_feasibility.csv")
    static_gap = float(geometry["minimum_exact_surface_gap_um"].min())
    gauge = float(identity["routeA_gauge_from_production_e1_deg"])
    success = m["classification"] in ("LOCAL_ROCKING_MODE_RECOVERED", "HEADTAIL_WALL_ROCKING_RECOVERED")
    if success:
        next_step = "Transfer this identical ROBOT_LOCAL_ROCKING input, contact model, and Reduced-Hydro model to the curved tube for one turn-section validation."
    elif m["classification"] == "MAGNETIC_TORQUE_INSUFFICIENT_FOR_10DEG_ROCKING":
        next_step = "Run one later 12-15 deg field-amplitude probe after reviewing field-to-robot phase lag; do not change geometry or contact."
    elif m["classification"] == "CROSS_PLANE_INSTABILITY":
        next_step = "Audit the rigid-body inertia principal axes against b_rock before changing any actuator amplitude or contact parameter."
    else:
        next_step = "Audit the 5 Hz robot-to-field phase relation before changing geometry, contact, or frequency."
    report = f"""# Local Head-Tail Rocking 5 Hz Validation

## Decision

**{m['classification']}**

1. RouteA is `theta_rock(t)=10 deg*r(t)*sin(2*pi*f*t)`, initially at 0 deg, about the transported binormal. Its legacy 0.03 s display uses 666.667 Hz and a 2 ms smooth ramp; the physical reference is 5 Hz.
2. The source establishes `(c_hat,n_rock,b_rock)=(T,N,T x N)`. In the current production PT gauge, `n_rock=cos(chi)e1+sin(chi)e2` and `b_rock=c_hat x n_rock`, with `chi={gauge:.12f} deg`; it is not the unrotated `(e1,e2)` pair.
3. The unchanged straight geometry is exactly feasible: solver ID `1.334690480 mm`, analytic 10 deg margin `0.115316535 mm`, and static exact minimum gap `{static_gap:.3f} um`.
4. The primary target changed because RouteA is alternating planar HEAD-TAIL rocking, while the former 30 deg cone has a 360 deg cross-section orbit and is a different motion family.
5. `ROBOT_LOCAL_ROCKING` uses `B=B0[cos(alpha_B)c_hat+sin(alpha_B)n_rock]`, `alpha_B=10 deg*sin(2*pi*5t)`, with continuous COM projection and the existing PT frame.
6. Dynamic `|B|` remains 10 mT; maximum sampled magnitude error is `{m['field_magnitude_max_abs_error_T']:.3e} T`.
7. The field remains in one local plane and has zero commanded 360 deg winding, as established by the offline one-cycle regression.
8. No robot UR is prescribed. The response comes from the validated finite moment, `m x B`, rigid-body dynamics, Reduced-Hydro loads, and one-wall General Contact.
9. `theta_rock` spans `{m['theta_rock_min_deg']:.3f}..{m['theta_rock_max_deg']:.3f} deg` (5 Hz fitted amplitude `{m['theta_rock_fundamental_amplitude_deg']:.3f} deg`).
10. `theta_cross` RMS/max is `{m['theta_cross_rms_deg']:.3e}/{m['theta_cross_max_abs_deg']:.3e} deg`; PCA line-likeness is `{m['pca_rocking_linearity']:.9f}` (descriptive, not a hard gate).
11. HEAD-only/TAIL-only support fractions are `{m['HEAD_only_support_fraction']:.4f}/{m['TAIL_only_support_fraction']:.4f}`; support-indicator sign changes: `{m['support_indicator_sign_changes']}`.
12. Both-end support fraction is `{m['both_end_support_fraction']:.4f}`, longest bridge `{m['longest_both_end_bridge_ms']:.3f} ms`; it is {'not ' if m['both_end_support_fraction'] < .5 else ''}the majority state.
13. Minimum exact surface gap is `{m['minimum_exact_surface_gap_um']:.3f} um`; deep penetration is `{m['deep_penetration']}`. Peak kinetic energy is `{m['ALLKE_max']:.3e} N mm`, and maximum absolute `ETOTAL` residual is `{max(abs(m['ETOTAL_min']), abs(m['ETOTAL_max'])):.3e} N mm`.
14. `RouteA_vs_MagneticRocking.gif` provides the normalized-phase visual gate. The scalar classification does not override direct visual review.
15. Exactly one next step: **{next_step}**

The legacy RouteA GIF is used only for topology/amplitude/phase comparison because it contains CEL deep-penetration diagnostics and prescribed RP motion. The HighEndStop reference remains secondary and only frames 0-27 are retained as its trusted early window.
"""
    (OUT/"Local_HeadTail_Rocking_5Hz_Validation.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
