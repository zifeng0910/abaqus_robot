"""Offline force/backtrack audit and open-loop first-order G sensitivity for F100."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
JOB = "TRUECEL_B0P11_A14P5_F100_FAST"
CASE = OUT / "case" / JOB
G0 = 2.0
CANDIDATES = np.arange(2.05, 2.401, 0.05)


def unit(x):
    x = np.asarray(x, dtype=float)
    return x / np.linalg.norm(x)


def vector(npz, stem):
    t = npz[stem + "1"][:, 0].astype(float)
    assert all(np.array_equal(t, npz[stem + str(i)][:, 0]) for i in (2, 3))
    return t, np.column_stack([npz[stem + str(i)][:, 1] for i in (1, 2, 3)])


def contact_vector(npz, surface):
    keys = [next(k for k in npz.files if f"|CFT{i} " in k and k.endswith(surface)) for i in (1, 2, 3)]
    t = npz[keys[0]][:, 0].astype(float)
    assert all(np.array_equal(t, npz[k][:, 0]) for k in keys)
    return t, np.column_stack([npz[k][:, 1] for k in keys])


def negative_runs(t, v):
    mask = v < 0
    starts = np.flatnonzero(np.diff(np.r_[False, mask].astype(int)) == 1)
    ends = np.flatnonzero(np.diff(np.r_[mask, False].astype(int)) == -1)
    return [{"start_s": float(t[a]), "end_s": float(t[b]),
             "duration_s": float(t[b] - t[a]), "minimum_v_s_mm_s": float(v[a:b+1].min())}
            for a, b in zip(starts, ends)]


def metrics(t, s, v):
    back = np.maximum.accumulate(s) - s
    return {
        "max_backtrack_mm": float(back.max()),
        "minimum_v_s_mm_s": float(v.min()),
        "cycle2_end_v_s_mm_s": float(np.interp(0.02, t, v)),
        "cycle2_delta_s_mm": float(np.interp(0.02, t, s) - np.interp(0.01, t, s)),
    }


def main():
    identity = json.loads((CASE / "case_identity.json").read_text(encoding="utf-8"))
    c = unit(identity["canonical_plus_s_axis_aba"])
    rp = np.load(CASE / "private" / "rp_history_private.npz")
    t, u = vector(rp, "U")
    vt, velocity = vector(rp, "V")
    assert np.array_equal(t, vt)
    rp.close()
    s = u @ c
    v_s = velocity @ c

    magnetic = np.genfromtxt(CASE / "magnetic_increment_f100.csv", delimiter=",", names=True)
    mt = magnetic["time_s"].astype(float)
    fmag = np.column_stack([magnetic[f"Fmag{i}_N"] for i in (1, 2, 3)])
    unique = np.r_[True, np.diff(mt) > 0]
    mt, fmag = mt[unique], fmag[unique]
    fmag_s = fmag @ c

    contact = np.load(CASE / "private" / "contact_history_private.npz")
    robot_surface = "on surface ASSEMBLY_ROBOT_SOLID-1_ROBOT_SOLID_SURF"
    pipe_surface = "on surface ASSEMBLY_PIPE_WALL_HELPER-1_PIPE_WALL_HELPER_SURF"
    ct, whole = contact_vector(contact, robot_surface)
    pt, pipe = contact_vector(contact, pipe_surface)
    contact.close()
    whole_s = whole @ c
    pipe_reaction_on_robot_proxy_s = -(pipe @ c)

    accumulated = np.maximum.accumulate(s)
    back = accumulated - s
    deficit_i = int(np.argmax(back))
    peak_i = int(np.argmax(s[:deficit_i + 1]))
    runs = negative_runs(t, v_s)
    recovery = np.flatnonzero((t > t[deficit_i]) & (s >= s[peak_i]))
    interval = {
        "maximum_forward_position_time_s": float(t[peak_i]),
        "maximum_forward_position_mm_relative": float(s[peak_i]),
        "backtracking_begin_time_s": float(t[peak_i]),
        "maximum_positional_deficit_time_s": float(t[deficit_i]),
        "maximum_positional_deficit_mm": float(back[deficit_i]),
        "forward_recovery_time_s": float(t[recovery[0]]) if len(recovery) else None,
        "negative_velocity_intervals": runs,
    }

    # Exact fixed-state scaling from the actual VUAMP implementation:
    # grad = G*grad_unit; F_i = ramp*sum_j(m_j*grad_ji).
    # The generated table uses grad_unit = d(s_eff)*(c outer c), hence
    # F_mag,s = ramp*G*d(s_eff)*(m dot c), exactly linear in G.
    mass_tonne = float(identity["robot_mass_mg"]) * 1e-9
    fmag_rp = np.interp(t, mt, fmag_s)
    rows = []
    replay = {}
    for g in CANDIDATES:
        delta_a = (g / G0 - 1.0) * fmag_rp / mass_tonne
        delta_v = cumulative_trapezoid(delta_a, t, initial=0.0)
        delta_s = cumulative_trapezoid(delta_v, t, initial=0.0)
        sg, vg = s + delta_s, v_s + delta_v
        row = {"G_mT": round(float(g), 2), **metrics(t, sg, vg)}
        rows.append(row)
        replay[round(float(g), 2)] = (sg, vg)
    threshold = next((r["G_mT"] for r in rows if r["max_backtrack_mm"] < 0.005), None)
    selected = None if threshold is None or threshold > 2.40 else min(2.40, round(threshold + 0.05, 2))

    screen_csv = OUT / "TRUECEL_F100_OPEN_LOOP_G_SCREEN.csv"
    with screen_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)

    timeline = np.linspace(0.014, 0.02, 1201)
    force_csv = OUT / "TRUECEL_F100_G_FORCE_BACKTRACK_HISTORY.csv"
    with force_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s", "s_relative_mm", "v_s_mm_s", "delta_s_backtrack_mm",
                         "F_mag_s_N", "whole_general_contact_force_s_N",
                         "pipe_surface_reaction_on_robot_proxy_s_N"])
        writer.writerows(zip(
            timeline, np.interp(timeline, t, s), np.interp(timeline, t, v_s),
            np.interp(timeline, t, back), np.interp(timeline, mt, fmag_s),
            np.interp(timeline, ct, whole_s), np.interp(timeline, pt, pipe_reaction_on_robot_proxy_s),
        ))

    result = {
        "classification": "OPEN_LOOP_FIRST_ORDER_G_ESTIMATE",
        "source_case": JOB,
        "exact_G_dependency": {
            "fortran": "grad=G*grad_unit; F_i=ramp*sum_j(moment_j*grad_ji)",
            "table_tensor": "grad_unit=d(s_eff)*(c outer c)",
            "axial_expression": "F_mag,s=ramp*G*d(s_eff)*(moment dot c)",
            "fixed_state_scaling": "F_mag,s(G)=G/G0*F_mag,s(G0), exact",
            "G_definition": "Gaussian-profile amplitude parameter in tesla; not a uniform mT/mm gradient",
            "B0_effect": "B0 scales magnetic torque only in this implementation; it is unchanged",
        },
        "backtrack_interval": interval,
        "baseline_metrics": metrics(t, s, v_s),
        "force_window_14_20ms": {
            "F_mag_s_min_N": float(fmag_s[(mt >= 0.014) & (mt <= 0.02)].min()),
            "F_mag_s_max_N": float(fmag_s[(mt >= 0.014) & (mt <= 0.02)].max()),
            "F_mag_s_impulse_Ns": float(np.trapz(fmag_s[(mt >= 0.014) & (mt <= 0.02)], mt[(mt >= 0.014) & (mt <= 0.02)])),
            "whole_contact_available": True,
            "direct_wall_force_available": False,
            "wall_proxy_caveat": "The F100 ODB has whole-surface CFT only. Pipe-surface CFT also contains pipe-fluid contact, so its opposite sign is exported only as a non-pair-isolated wall-reaction proxy, not claimed as direct robot-wall force.",
        },
        "screen": rows,
        "threshold_G_mT": threshold,
        "selected_G_mT": selected,
        "selection_rule": "lowest G with predicted MAX_BACKTRACK<0.005 mm plus one 0.05-mT robustness increment",
        "model_limitations": "Open-loop replay freezes all non-magnetic loads and the baseline trajectory-dependent magnetic field state; it is a screening estimate, not dynamics.",
    }
    (OUT / "TRUECEL_F100_OPEN_LOOP_G_ESTIMATE.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    table = "\n".join(
        f"| {r['G_mT']:.2f} | {r['max_backtrack_mm']:.6f} | {r['minimum_v_s_mm_s']:+.3f} | {r['cycle2_end_v_s_mm_s']:+.3f} | {r['cycle2_delta_s_mm']:+.6f} |"
        for r in rows)
    report = f"""# F100 axial-bias offline screen

Classification: **OPEN_LOOP_FIRST_ORDER_G_ESTIMATE**

The actual VUAMP source computes `grad = G*grad_unit` and `F_i = ramp*sum_j(moment_j*grad_ji)`. The table tensor is `grad_unit = d(s_eff)*(c outer c)`, therefore at fixed state:

`F_mag,s = ramp*G*d(s_eff)*(moment dot c)`

This is exactly linear in G. G is the amplitude of the Gaussian gradient profile in tesla, not a spatially uniform mT/mm gradient. B0 enters the torque term and remains unchanged.

The measured baseline reaches its last maximum at {interval['maximum_forward_position_time_s']*1e3:.6f} ms. Backtracking then continues to 20 ms, where the maximum deficit is {interval['maximum_positional_deficit_mm']:.6f} mm. No forward recovery occurs before the run ends. The negative-v interval is {runs[0]['start_s']*1e3:.6f}–{runs[0]['end_s']*1e3:.6f} ms.

| G (mT) | predicted MAX_BACKTRACK (mm) | predicted min v_s | predicted Cycle-2 end v_s | predicted Cycle-2 delta_s (mm) |
|---:|---:|---:|---:|---:|
{table}

The lowest predicted strong-pass value is G={threshold:.2f} mT. Applying the allowed one-step robustness margin selects **G={selected:.2f} mT** as the only dynamics candidate.

Important limitation: non-magnetic loads and the baseline trajectory are frozen. The pipe-surface force is exported only as a wall-reaction proxy because the existing ODB does not contain pair-isolated robot-wall CFT; it also includes pipe-fluid contact.
"""
    (OUT / "TRUECEL_F100_OPEN_LOOP_G_REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps({"threshold_G_mT": threshold, "selected_G_mT": selected,
                      "backtrack_interval": interval, "screen": rows}, indent=2))


if __name__ == "__main__":
    main()
