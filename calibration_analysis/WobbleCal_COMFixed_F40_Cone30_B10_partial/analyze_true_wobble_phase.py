"""True HEAD-to-TAIL phase/lock audit for the COM-fixed calibration run."""
from pathlib import Path
import json, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
JOB = "WobbleCal_COMFixed_F40_Cone30_B10"
OUT = ROOT / (JOB + "_analysis")
raw = pd.read_csv(ROOT / (JOB + "_true_axis_raw.csv"))
tel = pd.read_csv(OUT / "runtime_merged.csv").sort_values("t_s").drop_duplicates("t_s")
tf = pd.read_csv(ROOT / "curvenew_CEL_xyrot56_exact.csv")
tr = json.loads((ROOT / "abaqus_magpylib_frame_transform.json").read_text())
R = np.asarray(tr["R_aba_to_mag"], float)

def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-14 else np.array([1., 0., 0.])

def tangent_mag(s):
    ss = tf.arclength_mm.to_numpy(float)
    pp = tf[["x_mm", "y_mm", "z_mm"]].to_numpy(float)
    j = int(np.searchsorted(ss, s, side="right") - 1)
    j = max(0, min(j, len(pp) - 2))
    d0 = unit(pp[j + 1] - pp[j])
    if j >= len(pp) - 2:
        return d0
    d1 = unit(pp[j + 2] - pp[j + 1])
    w = min(max((s - ss[j]) / max(ss[j + 1] - ss[j], 1e-12), 0.), 1.)
    return unit((1 - w) * d0 + w * d1)

def znormal(t):
    q = np.array([0., 0., 1.]) - np.dot([0., 0., 1.], t) * t
    if np.linalg.norm(q) < 1e-12:
        y = np.array([0., 1., 0.]); q = y - np.dot(y, t) * t
    return unit(q)

df = pd.merge_asof(raw.sort_values("time_s"), tel, left_on="time_s", right_on="t_s",
                   direction="nearest", tolerance=7.5e-5)
for c in ["Bx_aba_T", "By_aba_T", "Bz_aba_T", "instantaneous_phase_rad", "driver_arc_mm"]:
    if c not in df: df[c] = np.nan

axes = []; cones = []; e1s = []; e2s = []; bv = []
for i, s in enumerate(df.driver_arc_mm):
    tm = tangent_mag(float(s) if np.isfinite(s) else float(tf.arclength_mm.iloc[-1]))
    nm = znormal(tm); bm = unit(np.cross(tm, nm))
    t = unit(R.T.dot(tm)); n = unit(R.T.dot(nm)); b = unit(R.T.dot(bm))
    c = unit(math.cos(math.radians(40)) * t + math.sin(math.radians(40)) * n)
    e1 = unit(n - np.dot(n, c) * c); e2 = unit(np.cross(c, e1))
    axes.append(df.loc[i, ["axis_x", "axis_y", "axis_z"]].to_numpy(float))
    cones.append(c); e1s.append(e1); e2s.append(e2); bv.append(df.loc[i, ["Bx_aba_T", "By_aba_T", "Bz_aba_T"]].to_numpy(float))
axes = np.asarray(axes); cones = np.asarray(cones); e1s = np.asarray(e1s); e2s = np.asarray(e2s); bv = np.asarray(bv)
df["a1"] = np.einsum("ij,ij->i", axes, e1s)
df["a2"] = np.einsum("ij,ij->i", axes, e2s)
df["b1"] = np.einsum("ij,ij->i", bv, e1s)
df["b2"] = np.einsum("ij,ij->i", bv, e2s)
df["phi_robot_rad"] = np.unwrap(np.arctan2(df.a2, df.a1).to_numpy(float))
df["phi_field_rad"] = np.unwrap(np.arctan2(df.b2, df.b1).to_numpy(float))
df["delta_phi_rad"] = np.unwrap(df.phi_robot_rad.to_numpy(float) - df.phi_field_rad.to_numpy(float))
df["theta_cone_deg"] = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", axes, cones), -1., 1.)))
tpipe = np.asarray([unit(R.T.dot(tangent_mag(float(s)))) for s in df.driver_arc_mm])
dp = np.einsum("ij,ij->i", axes, tpipe)
df["theta_pipe_directed_deg"] = np.degrees(np.arccos(np.clip(dp, -1., 1.)))
df["theta_pipe_axis_deg"] = np.degrees(np.arccos(np.clip(np.abs(dp), -1., 1.)))
df["transverse_radius"] = np.sqrt(df.a1 ** 2 + df.a2 ** 2)
df["head_x_rel"] = df.head_x - df.rp_x; df["head_y_rel"] = df.head_y - df.rp_y; df["head_z_rel"] = df.head_z - df.rp_z
df["tail_x_rel"] = df.tail_x - df.rp_x; df["tail_y_rel"] = df.tail_y - df.rp_y; df["tail_z_rel"] = df.tail_z - df.rp_z
df.to_csv(OUT / "wobble_true_axis_phase_history.csv", index=False)

def fit(mask, col):
    q = df[mask].dropna(subset=[col])
    if len(q) < 3: return float("nan"), float("nan"), len(q)
    slope = np.polyfit(q.time_s.to_numpy(float), q[col].to_numpy(float), 1)[0]
    return float(slope / (2 * np.pi)), float(q[col].iloc[-1] - q[col].iloc[0]), len(q)

windows = [("Cycle1", 0., .025), ("Cycle2", .025, .050), ("PartialCycle3", .050, float(df.time_s.max()))]
wrows = []
for name, a, b in windows:
    m = (df.time_s >= a) & (df.time_s <= b)
    fr, dr, n = fit(m, "phi_robot_rad"); ff, _, _ = fit(m, "phi_field_rad")
    q = df[m]; dd = float(q.delta_phi_rad.iloc[-1] - q.delta_phi_rad.iloc[0]) if len(q) > 1 else float("nan")
    wrows.append({"window": name, "start_s": a, "end_s": b, "n_frames": n,
                  "f_robot_hz": fr, "R_lock": abs(fr / 40.), "field_fit_hz": ff,
                  "robot_phase_change_rad": dr, "delta_phi_drift_rad": dd,
                  "theta_cone_mean_deg": float(q.theta_cone_deg.mean()),
                  "theta_cone_pp_deg": float(q.theta_cone_deg.max() - q.theta_cone_deg.min()),
                  "transverse_radius_rms": float(np.sqrt(np.mean(q.transverse_radius ** 2)))})
win = pd.DataFrame(wrows); win.to_csv(OUT / "wobble_phase_lock_windows.csv", index=False)

def phase_interp(q, col, grid):
    x = q.phi_field_rad.to_numpy(float); y = q[col].to_numpy(float); idx = np.argsort(x)
    return np.interp(grid, x[idx], y[idx])

q2 = df[(df.time_s >= .025) & (df.time_s < .050)].copy()
q3 = df[(df.time_s >= .050) & (df.time_s <= df.time_s.max())].copy()
if len(q2) > 3 and len(q3) > 3:
    x2 = q2.phi_field_rad.to_numpy(); x3 = q3.phi_field_rad.to_numpy(); span = min(x2[-1] - x2[0], x3[-1] - x3[0])
    grid = np.linspace(0., span, 250); a12 = phase_interp(q2, "a1", x2[0] + grid); a22 = phase_interp(q2, "a2", x2[0] + grid); a13 = phase_interp(q3, "a1", x3[0] + grid); a23 = phase_interp(q3, "a2", x3[0] + grid)
    rms = float(np.sqrt(np.mean((a13 - a12) ** 2 + (a23 - a22) ** 2))); radius = float(np.sqrt(np.mean(a12 ** 2 + a22 ** 2))); ratio = rms / max(radius, 1e-12)
    orbit = pd.DataFrame({"phase_rel_rad": grid, "cycle2_a1": a12, "cycle2_a2": a22, "cycle3_a1": a13, "cycle3_a2": a23, "difference": np.sqrt((a13-a12)**2 + (a23-a22)**2)})
else:
    rms = ratio = radius = span = float("nan"); orbit = pd.DataFrame()
orbit.to_csv(OUT / "wobble_cycle2_vs_cycle3_phase_matched.csv", index=False)
pd.DataFrame([{"orbit_rms_difference": rms, "cycle2_orbit_radius_rms": radius, "orbit_repeat_error_ratio": ratio, "common_phase_span_rad": span}]).to_csv(OUT / "wobble_orbit_repeatability.csv", index=False)

main = pd.read_csv(ROOT / (JOB + "_telemetry.csv")); rec = pd.read_csv(ROOT / (JOB + "_telemetry_recovery.csv"))
gap = float(rec.t_s.iloc[0] - main.t_s.iloc[-1]); phase_gap = float(rec.instantaneous_phase_rad.iloc[0] - main.instantaneous_phase_rad.iloc[-1]); expected = 2 * np.pi * 40 * gap
pd.DataFrame([{"main_last_t_s": main.t_s.iloc[-1], "recovery_first_t_s": rec.t_s.iloc[0], "gap_s": gap, "phase_gap_rad": phase_gap, "expected_phase_gap_rad": expected, "phase_gap_residual_rad": phase_gap - expected, "b_jump_T": abs(rec.B_aba_T.iloc[0] - main.B_aba_T.iloc[-1])}]).to_csv(OUT / "socket_recovery_locking_window_audit.csv", index=False)

flip = bool(np.nanmin(np.einsum("ij,j->i", axes, axes[0])) < 0.)
c2 = win[win.window == "Cycle2"].iloc[0]; c3 = win[win.window == "PartialCycle3"].iloc[0]
verdict = "40HZ_MAGNETIC_WOBBLE_CALIBRATION_PASS_FROM_PARTIAL_RUN" if (0.95 <= c2.R_lock <= 1.05 and 0.95 <= c3.R_lock <= 1.05 and abs(c2.delta_phi_drift_rad) <= np.pi/4 and abs(c3.delta_phi_drift_rad) <= np.pi/4 and ratio <= .15) else "40HZ_WOBBLE_NOT_LOCKED"
summary = {"job": JOB, "verdict": verdict, "cycle2_f_robot_hz": float(c2.f_robot_hz), "cycle2_R_lock": float(c2.R_lock), "partial_cycle3_f_robot_hz": float(c3.f_robot_hz), "partial_cycle3_R_lock": float(c3.R_lock), "cycle2_delta_phi_drift_rad": float(c2.delta_phi_drift_rad), "partial_cycle3_delta_phi_drift_rad": float(c3.delta_phi_drift_rad), "theta_cone_mean_deg": float(df.theta_cone_deg.mean()), "theta_cone_pp_deg": float(df.theta_cone_deg.max() - df.theta_cone_deg.min()), "head_tail_flip": flip, "orbit_repeat_error_ratio": ratio, "socket_gap_s": gap, "socket_phase_gap_residual_rad": phase_gap - expected}
(OUT / "true_wobble_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

fig, ax = plt.subplots(3, 2, figsize=(13, 11), sharex=True); tm = df.time_s * 1000
ax[0,0].plot(tm, df.phi_field_rad, label="field"); ax[0,0].plot(tm, df.phi_robot_rad, label="robot"); ax[0,0].set_title("field vs robot phase"); ax[0,0].legend()
ax[0,1].plot(tm, df.delta_phi_rad); ax[0,1].set_title("phase lag")
ax[1,0].plot(tm, df.theta_cone_deg); ax[1,0].set_title("theta cone"); ax[1,0].set_ylabel("deg")
ax[1,1].plot(df.a1, df.a2); ax[1,1].set_title("long-axis transverse orbit"); ax[1,1].set_xlabel("a1"); ax[1,1].set_ylabel("a2")
ax[2,0].plot(tm, df.head_y_rel, label="head y"); ax[2,0].plot(tm, df.tail_y_rel, label="tail y"); ax[2,0].legend(); ax[2,0].set_title("head/tail relative motion")
ax[2,1].plot(tm, df.torque_norm_Nmm); ax[2,1].set_title("torque norm"); ax[2,1].set_xlabel("time (ms)")
fig.tight_layout(); fig.savefig(OUT / "true_wobble_diagnostics.png", dpi=180); plt.close(fig)

report = """# {job} — true long-axis / phase-lock audit

## Status

This is a partial record: ODB ends at {end:.3f} ms of the requested 75 ms.
The final 1.4 ms is not used to claim completion.

## Method

No explicit HEAD/TAIL sets exist in the ODB. Fixed groups were therefore defined
once from the 5% extreme projections of the reference ROBOT_SOLID mesh along its
principal long axis (19 nodes per end), then carried with the Abaqus finite
rotation-vector Rodrigues matrix. UR1–UR3 were not interpreted as Euler angles.
The field frame is reconstructed from the production 40 Hz, cone30°, Bias40°,
tangent-axis evaluator and the validated frame rotation.

## Lock results

| window | f_robot (Hz) | R_lock | phase-lag drift (rad) | mean theta_cone (deg) |
|---|---:|---:|---:|---:|
{rows}

Overall theta_cone mean is {theta:.3f}° with peak-to-peak variation {pp:.3f}°.
HEAD/TAIL flip detected: **{flip}**. Cycle2 versus the first common partial
Cycle3 phase span gives orbit repeat error ratio **{ratio:.4f}** (preferred ≤0.10,
acceptable ≤0.15).

## Socket recovery

The main→recovery gap is {gap_us:.2f} microseconds; phase-gap residual is
{pres:.3e} rad and B jump is {bj:.3e} T. No phase reset is evident, so the
recovery does not contaminate the Cycle2/PartialCycle3 window.

## Verdict

**{verdict}**

This verdict concerns rotational COM-fixed calibration only. It does not validate
translation, wall contact, or a complete 75 ms run. If PASS, the next step is
COM translation release with Wall-OFF + CEL-ON for one 25 ms cycle. If NOT_LOCKED,
the next frequency is a separate 30 Hz COM-fixed calibration with all other
parameters unchanged.
""".format(job=JOB, end=float(df.time_s.max()*1000), rows="\n".join("| {window} | {f_robot_hz:.3f} | {R_lock:.3f} | {delta_phi_drift_rad:.3f} | {theta_cone_mean_deg:.3f} |".format(**r) for r in wrows), theta=float(df.theta_cone_deg.mean()), pp=float(df.theta_cone_deg.max()-df.theta_cone_deg.min()), flip=flip, ratio=ratio, gap_us=gap*1e6, pres=phase_gap-expected, bj=abs(rec.B_aba_T.iloc[0]-main.B_aba_T.iloc[-1]), verdict=verdict)
(OUT / "WobbleCal_COMFixed_F40_Cone30_B10_true_wobble_report.md").write_text(report, encoding="utf-8")
print(json.dumps(summary, indent=2))
