"""Analyze the three initial v_mean calibration candidates and choose one refinement G."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
sys.path.insert(0, str(REPO / "calibration_analysis" / "FastPrecomputedForwardFlowScreen" / "scripts"))
import analyze_render_forward_flow as base

JOBS = ["FAST_VMEAN5_G2_UFLOW0P5_2CYCLES", "FAST_VMEAN5_G3_UFLOW0P5_2CYCLES", "FAST_VMEAN5_G4_UFLOW0P5_2CYCLES"]
FORCE_THRESHOLD = 1e-10


def total_contact_history(case, time):
    archive = np.load(case / "private" / "contact_history_private.npz")
    vectors = []
    for prefix in ("CFN", "CFS"):
        columns = []
        for axis in (1, 2, 3):
            keys = [key for key in archive.files if "|" + prefix + str(axis) + " on surface " in key and "ASSEMBLY_ROBOT" in key]
            if keys:
                arr = max((archive[key] for key in keys), key=len)
                columns.append(np.interp(time, arr[:, 0], arr[:, 1]))
            else:
                columns.append(np.zeros_like(time))
        vectors.append(np.column_stack(columns))
    archive.close()
    return vectors[0] + vectors[1]


def candidate(job):
    case = OUT / "case" / job
    identity = json.loads((case / "case_identity.json").read_text())
    c, n, b = (base.unit(identity[key]) for key in ("canonical_plus_s_axis_aba", "n_routeA_aba", "b_routeA_aba"))
    rp0 = np.asarray(identity["initial_center_aba_mm"], float)
    center = rp0 - identity["radial_offset_n_mm"] * n
    deck = (case / (job + ".inp")).read_text(encoding="latin1")
    labels, nodes = base.part_nodes(deck, "Robot_SOLID")
    _, wall = base.part_nodes(deck, "Pipe_WALL_HELPER")
    normals, radius = base.wall_planes(wall, center, c, n, b)
    s0 = (nodes - rp0).dot(c)
    region = np.full(len(nodes), "BODY", dtype=object)
    region[s0 <= s0.min() + .25] = "TAIL"; region[s0 >= s0.max() - .25] = "HEAD"

    rp = np.load(case / "private" / "rp_history_private.npz")
    time = rp["V1"][:, 0]
    U, UR, V = (base.vector(rp, prefix, time) for prefix in ("U", "UR", "V"))
    rp.close()
    axial, vs = U.dot(c), V.dot(c)
    body_axis = Rotation.from_rotvec(UR).apply(np.broadcast_to(c, UR.shape))
    rocking = np.degrees(np.arctan2(body_axis.dot(n), body_axis.dot(c)))
    centered = rocking - np.mean(rocking)
    frequencies = np.fft.rfftfreq(len(time), np.median(np.diff(time)))
    spectrum = np.abs(np.fft.rfft(centered))
    band = (frequencies >= 60) & (frequencies <= 240)
    dominant = float(frequencies[band][np.argmax(spectrum[band])])

    field = np.load(case / "private" / "robot_contact_fields_private.npz")
    ctime = field["time"].astype(float)
    if not np.array_equal(field["node_labels"], labels):
        raise RuntimeError("node order mismatch: " + job)
    normal, shear = field["CNORMF"].astype(float), field["CSHEARF"].astype(float)
    field.close()
    Ui = np.column_stack([np.interp(ctime, time, U[:, i]) for i in range(3)])
    URi = np.column_stack([np.interp(ctime, time, UR[:, i]) for i in range(3)])
    positions = [rp0 + Ui[k] + Rotation.from_rotvec(URi[k]).apply(nodes - rp0) for k in range(len(ctime))]
    gaps, active, counts, impact, separation, same_end_recontact = {}, {}, {}, {}, {}, {}
    for name in ("HEAD", "TAIL"):
        mask = region == name
        values = []
        for pts in positions:
            rel = pts[mask] - center
            transverse = rel - np.outer(rel.dot(c), c)
            values.append(radius - transverse.dot(normals.T).max(axis=1).max())
        gaps[name] = np.asarray(values)
        force = np.linalg.norm((normal[:, mask, :] + shear[:, mask, :]).sum(axis=1), axis=1)
        active[name] = (force > FORCE_THRESHOLD) & (gaps[name] <= .015)
        spans = base.runs(active[name]); counts[name] = len(spans)
        gv = np.gradient(gaps[name], ctime)
        impact_values, separation_values, ratios = [], [], []
        for index, (i, j) in enumerate(spans):
            vin = float(min(0, np.min(gv[max(0, i-2):min(len(gv), i+2)])))
            end = spans[index+1][0] if index+1 < len(spans) else min(len(gv)-1, j+10)
            vout = float(max(0, np.max(gv[min(j+1, len(gv)-1):end+1])))
            impact_values.append(vin); separation_values.append(vout)
            if vin < 0: ratios.append(vout / abs(vin))
        impact[name] = min(impact_values, default=0.0)
        separation[name] = max(separation_values, default=0.0)
        half_counts = []
        for half in range(4):
            lo, hi = half / 240.0, (half+1) / 240.0
            half_counts.append(sum(ctime[i] < hi and ctime[j] >= lo for i, j in spans))
        same_end_recontact[name] = max(half_counts, default=0) > 1

    total_force = total_contact_history(case, time)
    peak_force = float(np.linalg.norm(total_force, axis=1).max())
    both = active["HEAD"] & active["TAIL"]
    both_longest = max((ctime[j] - ctime[i] for i, j in base.runs(both)), default=0.0)
    min_gap = min(float(gaps["HEAD"].min()), float(gaps["TAIL"].min()))
    log = pd.read_csv(case / "magnetic_hydro_increment.csv")
    phase = np.unwrap(np.radians(log.phase_deg.to_numpy()))
    phase_rate = np.gradient(phase, log.time_s.to_numpy())
    phase_error = float(np.max(np.abs(phase_rate[10:-10] / (2*np.pi*120) - 1)))
    duration = float(time[-1] - time[0]); delta = float(axial[-1] - axial[0]); mean = delta / duration
    catastrophic_force = peak_force > 0.5
    tumble = bool(np.any(body_axis.dot(c) < 0))
    gross_penetration = min_gap < -0.05
    persistent_bridge = both_longest > .001
    safety_passed = not (catastrophic_force or tumble or gross_penetration or persistent_bridge)
    return {
        "job": job, "G_mT": identity["gradient_mT"], "delta_s_mm": delta, "v_mean_mm_s": mean,
        "v_final_mm_s": float(vs[-1]), "v_max_mm_s": float(vs.max()),
        "fraction_v_s_positive": float(np.mean(vs > 0)),
        "rocking_min_deg": float(rocking.min()), "rocking_max_deg": float(rocking.max()),
        "dominant_rocking_frequency_Hz": dominant, "phase_rate_max_relative_error": phase_error,
        "TAIL_contact_episodes": counts["TAIL"], "HEAD_contact_episodes": counts["HEAD"],
        "peak_contact_resultant_N": peak_force,
        "TAIL_impact_wall_normal_velocity_mm_s": impact["TAIL"],
        "HEAD_impact_wall_normal_velocity_mm_s": impact["HEAD"],
        "TAIL_post_contact_separation_velocity_mm_s": separation["TAIL"],
        "HEAD_post_contact_separation_velocity_mm_s": separation["HEAD"],
        "TAIL_same_end_recontact": bool(same_end_recontact["TAIL"]), "HEAD_same_end_recontact": bool(same_end_recontact["HEAD"]),
        "obvious_rebound": bool(same_end_recontact["TAIL"] or same_end_recontact["HEAD"]),
        "minimum_wall_gap_mm": min_gap, "gross_penetration": bool(gross_penetration),
        "both_wall_bridge_longest_ms": both_longest * 1e3, "tumble": tumble,
        "contact_force_catastrophic": bool(catastrophic_force), "contact_safety_passed": bool(safety_passed),
        "background_flow_mm_s": 0.5, "v_mean_over_U_flow": mean / 0.5,
        "runtime_s": identity["wallclock_s"], "socket_calls": identity["socket_calls"],
    }


def main():
    rows = [candidate(job) for job in JOBS]
    safe = [row for row in rows if row["contact_safety_passed"]]
    if len(safe) < 2:
        raise RuntimeError("fewer than two contact-safe candidates; no interpolation")
    coefficients = np.polyfit([row["G_mT"] for row in safe], [row["v_mean_mm_s"] for row in safe], 1)
    predicted = float((5.0 - coefficients[1]) / coefficients[0])
    predicted = round(predicted, 3)
    plan = {"fit_slope_mm_s_per_mT": float(coefficients[0]), "fit_intercept_mm_s": float(coefficients[1]),
            "predicted_G_mT_for_5mm_s": predicted, "target_v_mean_mm_s": 5.0,
            "initial_candidate_count": 3, "refinement_dynamics_remaining": 1}
    (OUT / "initial_candidate_metrics.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="ascii")
    pd.DataFrame(rows).to_csv(OUT / "initial_candidate_metrics.csv", index=False)
    (OUT / "refinement_plan.json").write_text(json.dumps(plan, indent=2) + "\n", encoding="ascii")
    print(json.dumps({"candidates": rows, "refinement_plan": plan}, indent=2))


if __name__ == "__main__":
    main()
