"""Post-process the sole TAIL-gap local mesh-resolution diagnostic."""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapezoid
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
COARSE = "TRUECEL_A14P5_FORCE_FLUX_DIAG_2CYCLES"
FINE = "TRUECEL_A14P5_TAILGAP_REFINE_2CYCLES"
T = 1.0 / 120.0
BASELINE_J_FLUID = -6.777448767817853e-7

sys.path.insert(0, str(REPO / "calibration_analysis" / "FastPrecomputedForwardFlowScreen" / "scripts"))
import analyze_render_forward_flow as base


def unit(values):
    values = np.asarray(values, dtype=float)
    return values / np.linalg.norm(values)


def vector(npz, prefix, time):
    return np.column_stack([
        np.interp(time, npz[f"{prefix}{axis}"][:, 0], npz[f"{prefix}{axis}"][:, 1])
        for axis in (1, 2, 3)
    ])


def integrate_window(time, values, start, end):
    inner = (time > start) & (time < end)
    sample_t = np.r_[start, time[inner], end]
    sample_v = np.r_[np.interp(start, time, values), values[inner], np.interp(end, time, values)]
    return float(trapezoid(sample_v, sample_t))


def duration_runs(time, active):
    return [float(time[j] - time[i]) for i, j in base.runs(active)]


def weighted_mean(values, weights):
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    return float(np.average(values[valid], weights=weights[valid])) if np.any(valid) else float("nan")


def weighted_percentile(values, weights, percentile):
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(valid):
        return float("nan")
    values, weights = values[valid], weights[valid]
    order = np.argsort(values)
    target = percentile / 100.0 * np.sum(weights)
    return float(values[order][np.searchsorted(np.cumsum(weights[order]), target, side="left")])


def hex_volumes(nodes, connectivity):
    xyz = nodes[connectivity]
    signs = np.asarray(((-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
                        (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)), dtype=float)
    gauss = 1.0 / np.sqrt(3.0)
    result = np.zeros(len(connectivity), dtype=float)
    for xi in (-gauss, gauss):
        for eta in (-gauss, gauss):
            for zeta in (-gauss, gauss):
                dndxi = 0.125 * signs[:, 0] * (1 + signs[:, 1]*eta) * (1 + signs[:, 2]*zeta)
                dndeta = 0.125 * signs[:, 1] * (1 + signs[:, 0]*xi) * (1 + signs[:, 2]*zeta)
                dndzeta = 0.125 * signs[:, 2] * (1 + signs[:, 0]*xi) * (1 + signs[:, 1]*eta)
                jac = np.stack((np.einsum("eni,n->ei", xyz, dndxi),
                                np.einsum("eni,n->ei", xyz, dndeta),
                                np.einsum("eni,n->ei", xyz, dndzeta)), axis=1)
                result += np.linalg.det(jac)
    if np.min(result) <= 0:
        raise RuntimeError(f"Nonpositive computed hex volume: {np.min(result)}")
    return result


def history_vector(npz, pairwise):
    values, keys = [], []
    for axis in (1, 2, 3):
        matches = [key for key in npz.files if f"|CFT{axis} " in key
                   and "on surface ASSEMBLY_ROBOT_SOLID-1_ROBOT_SOLID_SURF" in key
                   and (("/ASSEMBLY_PIPE_WALL_HELPER" in key) == pairwise)]
        if len(matches) != 1:
            raise RuntimeError(f"CFT{axis} pairwise={pairwise}: {matches}")
        keys.append(matches[0]); values.append(npz[matches[0]].astype(float))
    time = values[0][:, 0]
    force = np.column_stack([np.interp(time, row[:, 0], row[:, 1]) for row in values])
    return time, force, keys


def contact_metrics(case, identity, c, n, pipe, rp_time, displacement, rotation):
    deck = (case / f"{case.name}.inp").read_text(encoding="latin1")
    _, robot = base.part_nodes(deck, "Robot_SOLID")
    center0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    material_s = np.dot(robot - center0, c)
    regions = {
        "TAIL": material_s <= material_s.min() + 0.25,
        "HEAD": material_s >= material_s.max() - 0.25,
    }
    fields = np.load(case / "private" / "robot_contact_fields_private.npz")
    time = fields["time"].astype(float)
    force = fields["CNORMF"].astype(float) + fields["CSHEARF"].astype(float)
    u = np.column_stack([np.interp(time, rp_time, displacement[:, k]) for k in range(3)])
    ur = np.column_stack([np.interp(time, rp_time, rotation[:, k]) for k in range(3)])
    radius = float(identity["lumen_radius_mm"])
    gaps = {name: np.empty(len(time)) for name in regions}
    normal = {name: np.zeros(len(time)) for name in regions}
    for k in range(len(time)):
        points = center0 + u[k] + Rotation.from_rotvec(ur[k]).apply(robot - center0)
        relative = points - pipe
        qn = relative @ n
        qb = relative - np.outer(relative @ c, c) - np.outer(qn, n)
        radial = np.sqrt(qn*qn + np.einsum("ij,ij->i", qb, qb))
        outward = (relative - np.outer(relative @ c, c)) / np.maximum(radial[:, None], 1e-12)
        gap = radius - radial
        for name, mask in regions.items():
            gaps[name][k] = float(np.min(gap[mask]))
            selected = mask & (gap <= 0.02)
            if np.any(selected):
                normal[name][k] = float(np.sum(np.abs(np.einsum(
                    "ij,ij->i", force[k, selected], outward[selected]
                ))))
    active = {name: (gaps[name] <= 0.015) & (normal[name] > 1e-9) for name in regions}
    episodes = {name: base.runs(active[name]) for name in active}
    event = {}
    for name, opposite in (("TAIL", "HEAD"), ("HEAD", "TAIL")):
        if not episodes[name]:
            event[name] = {"episodes": 0, "first_R_v": None, "episode_windows_ms": []}
            continue
        i, j = episodes[name][0]
        rate = np.gradient(gaps[name], time)
        inward = max(0.0, -float(rate[max(0, i-1)]))
        next_opposite = np.flatnonzero(active[opposite] & (np.arange(len(time)) > j))
        stop = int(next_opposite[0]) if len(next_opposite) else min(len(time)-1, j+40)
        outward = max(0.0, float(np.max(rate[j:stop+1])))
        event[name] = {
            "episodes": len(episodes[name]),
            "first_R_v": None if inward <= 1e-12 else outward/inward,
            "episode_windows_ms": [[float(time[a]*1e3), float(time[z]*1e3)] for a,z in episodes[name]],
            "first_contact_ms": float(time[i]*1e3),
            "first_dwell_ms": float((time[j]-time[i])*1e3),
        }
    both = active["TAIL"] & active["HEAD"]
    event["both_bridge_ms"] = max(duration_runs(time, both), default=0.0)*1e3
    event["max_penetration_mm"] = max(0.0, -min(np.min(gaps["TAIL"]), np.min(gaps["HEAD"])))
    return event, time, active


def load_case(job):
    case = OUT / "case" / job
    identity = json.loads((case / "case_identity.json").read_text())
    c = unit(identity["canonical_plus_s_axis_aba"])
    n = unit(identity["n_routeA_aba"])
    b = unit(identity["b_routeA_aba"])
    center0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    pipe = center0 - float(identity["initial_axial_shift_mm"])*c - float(identity["radial_offset_n_mm"])*n
    rp = np.load(case / "private" / "rp_history_private.npz")
    time = rp["U1"][:, 0].astype(float)
    u, ur, velocity = (vector(rp, prefix, time) for prefix in ("U", "UR", "V"))
    axial = u @ c
    vs = velocity @ c
    axis = Rotation.from_rotvec(ur).apply(np.broadcast_to(c, ur.shape))
    rocking = np.degrees(np.arctan2(axis@n, axis@c))
    contact, contact_time, active = contact_metrics(case, identity, c, n, pipe, time, u, ur)
    return {"job":job, "case":case, "identity":identity, "c":c, "n":n, "b":b,
            "center0":center0, "pipe":pipe, "time":time, "u":u, "ur":ur,
            "velocity":velocity, "axial":axial, "vs":vs, "axis":axis,
            "rocking":rocking, "contact":contact, "contact_time":contact_time, "active":active}


def force_closure(data):
    histories = np.load(data["case"] / "private" / "force_flux_contact_history_private.npz")
    ht, whole, whole_keys = history_vector(histories, False)
    wt, wall, wall_keys = history_vector(histories, True)
    wall_i = np.column_stack([np.interp(ht, wt, wall[:, k]) for k in range(3)])
    fwhole = whole @ data["c"]
    fwall = wall_i @ data["c"]
    ffluid = fwhole - fwall
    magnetic_path = data["case"] / "magnetic_increment.csv"
    if not magnetic_path.exists() and data["job"] == COARSE:
        magnetic_path = OUT / "case" / "TRUECEL_A14P5_DIAG_2CYCLES" / "magnetic_increment.csv"
    magnetic = np.loadtxt(magnetic_path, delimiter=",", skiprows=1)
    mt = np.r_[0.0, magnetic[:, 0]]
    mf = np.vstack((np.zeros(3), magnetic[:, 3:6])) @ data["c"]
    start, end = 0.008, 0.0125
    jmag = integrate_window(mt, mf, start, end)
    jwall = integrate_window(ht, fwall, start, end)
    jfluid = integrate_window(ht, ffluid, start, end)
    mass = float(data["identity"]["robot_mass_mg"]) * 1e-9
    dp = mass * (np.interp(end, data["time"], data["vs"]) - np.interp(start, data["time"], data["vs"]))
    error = jmag + jwall + jfluid - dp
    scale = max(abs(dp), abs(jmag)+abs(jwall)+abs(jfluid), 1e-30)
    return {
        "history_time": ht, "fwhole": fwhole, "fwall": fwall, "ffluid": ffluid,
        "magnetic_time": mt, "fmag": mf,
        "critical_8p0_12p5_ms": {
            "J_mag_s_Ns_direct": jmag,
            "J_wall_s_Ns_direct": jwall,
            "J_fluid_s_Ns_inferred_whole_minus_wall": jfluid,
            "Delta_p_robot_s_Ns": dp,
            "closure_error_Ns": error,
            "relative_closure_error": error/scale,
        },
        "history_keys": {"whole":whole_keys, "wall_pair":wall_keys},
    }


def field_metrics(data):
    field = np.load(data["case"] / "private" / "truecel_field_private.npz")
    time = field["time"].astype(float)
    nodes = field["fluid_node_coordinates_mm"].astype(float)
    conn = field["fluid_connectivity_index"].astype(int)
    centers = field["fluid_element_centroids_mm"].astype(float)
    volumes = hex_volumes(nodes, conn)
    evf = np.clip(field["fluid_evf"].astype(float), 0.0, 1.0)
    pressure = field["fluid_pressure"].astype(float)
    velocity = np.nanmean(field["fluid_velocity_mm_s"][:, conn, :], axis=2)
    vs = velocity @ data["c"]
    speed = np.linalg.norm(velocity, axis=2)
    relative = centers - data["pipe"]
    s = relative @ data["c"]
    qn = relative @ data["n"]
    qb = relative @ data["b"]
    radius = np.sqrt(qn*qn + qb*qb)
    envelope = data["identity"].get("tail_swept_envelope", {
        "buffered_min_snb_mm":[-4.6868,-.75,-.75], "buffered_max_snb_mm":[-3.6511,.75,.70]
    })
    s0, s1 = envelope["buffered_min_snb_mm"][0], envelope["buffered_max_snb_mm"][0]
    corridor = (s >= s0) & (s <= s1) & (radius >= float(data["identity"]["lumen_radius_mm"])-0.25)
    rows = []
    strongest = {"value_mg_mm_s": 0.0}
    global_speed_p95_max = 0.0
    rho = float(data["identity"]["fluid_density_tonne_mm3"]) * 1e9
    for k, t in enumerate(time):
        weights = evf[k, corridor] * volumes[corridor]
        global_weights = evf[k] * volumes
        global_speed_p95_max = max(
            global_speed_p95_max, weighted_percentile(speed[k], global_weights, 95)
        )
        momentum_cell = rho * evf[k, corridor] * volumes[corridor] * vs[k, corridor]
        if 0.008 <= t <= 0.0125 and len(momentum_cell) and np.nanmin(momentum_cell) < strongest["value_mg_mm_s"]:
            local_index = int(np.nanargmin(momentum_cell))
            global_index = np.flatnonzero(corridor)[local_index]
            strongest = {
                "value_mg_mm_s": float(momentum_cell[local_index]),
                "time_ms": float(t*1e3),
                "element_label": int(field["fluid_element_labels"][global_index]),
                "s_n_b_mm": [float(s[global_index]), float(qn[global_index]), float(qb[global_index])],
            }
        rows.append({
            "time_s": float(t),
            "pressure_mean_N_mm2": weighted_mean(pressure[k,corridor], weights),
            "pressure_p95_abs_N_mm2": weighted_percentile(np.abs(pressure[k,corridor]), weights, 95),
            "pressure_max_abs_N_mm2": float(np.nanmax(np.abs(pressure[k,corridor][weights>0]))) if np.any(weights>0) else float("nan"),
            "velocity_s_mean_mm_s": weighted_mean(vs[k,corridor], weights),
            "velocity_s_p95_abs_mm_s": weighted_percentile(np.abs(vs[k,corridor]), weights, 95),
            "speed_p95_mm_s": weighted_percentile(speed[k,corridor], weights, 95),
            "axial_momentum_mg_mm_s": float(np.nansum(momentum_cell)),
            "EVF_volume_mm3": float(np.nansum(weights)),
        })
    inventory = evf @ volumes
    return {"time":time, "centers":centers, "volumes":volumes, "evf":evf, "pressure":pressure,
            "velocity":velocity, "vs":vs, "speed":speed, "s":s, "qn":qn, "qb":qb,
            "corridor":corridor, "rows":rows, "inventory":inventory,
            "retention":float(inventory[-1]/inventory[0]),
            "global_speed_p95_max_mm_s":float(global_speed_p95_max),
            "strongest_negative_momentum":strongest}


def cycle_and_recoil(data):
    cycles = []
    for number in (1,2):
        start, end = (number-1)*T, number*T
        s0, s1 = np.interp((start,end), data["time"], data["axial"])
        cycles.append({"cycle":number, "delta_s_mm":float(s1-s0),
                       "mean_v_s_mm_s":float((s1-s0)/T),
                       "end_v_s_mm_s":float(np.interp(end,data["time"],data["vs"]))})
    mask = data["time"] >= T
    negative = mask & (data["vs"] < 0)
    runs = base.runs(negative)
    return {"cycles":cycles, "minimum_v_s_mm_s":float(np.min(data["vs"][mask])),
            "negative_duration_ms":float(sum(duration_runs(data["time"],negative))*1e3),
            "negative_recoil_onset_ms":None if not runs else float(data["time"][runs[0][0]]*1e3),
            "cycle_end_v_s_mm_s":float(np.interp(2*T,data["time"],data["vs"]))}


def energy_metrics(data):
    energy = np.load(data["case"] / "private" / "energy_history_private.npz")
    result = {}
    for key in ("ETOTAL","ALLPW","ALLFD","ALLVD","ALLAE"):
        series = energy[key].astype(float)
        result[key] = {
            "initial_Nmm":float(series[0,1]), "final_Nmm":float(series[-1,1]),
            "max_abs_drift_Nmm":float(np.max(np.abs(series[:,1]-series[0,1]))),
            "cycle1_increment_Nmm":float(np.interp(T,series[:,0],series[:,1])-np.interp(0,series[:,0],series[:,1])),
            "cycle2_increment_Nmm":float(np.interp(2*T,series[:,0],series[:,1])-np.interp(T,series[:,0],series[:,1])),
        }
    sta = (data["case"] / f"{data['job']}.sta").read_text(encoding="latin1")
    rows = re.findall(r"(?m)^\s*\d+\s+[0-9.E+-]+\s+[0-9.E+-]+\s+\S+\s+([0-9.E+-]+)\s+\d+",sta)
    dt = np.asarray([float(value) for value in rows])
    result["stable_timestep_s"] = {"min":float(dt.min()),"max":float(dt.max()),"initial":float(dt[0])}
    return result


def pressure_imbalance(data, fields, closure):
    mask = (closure["history_time"] >= .008) & (closure["history_time"] <= .0125)
    peak_index = np.flatnonzero(mask)[np.argmin(closure["ffluid"][mask])]
    peak_time = float(closure["history_time"][peak_index])
    k = int(np.argmin(np.abs(fields["time"]-peak_time)))
    com = data["center0"] + np.asarray([np.interp(fields["time"][k],data["time"],data["u"][:,j]) for j in range(3)])
    ur = np.asarray([np.interp(fields["time"][k],data["time"],data["ur"][:,j]) for j in range(3)])
    axis = Rotation.from_rotvec(ur).apply(data["c"])
    tail = com - 1.3*axis
    rel = fields["centers"] - tail
    ds = rel @ axis
    tail_radial = tail - data["pipe"] - np.dot(tail-data["pipe"],data["c"])*data["c"]
    tail_radial /= max(np.linalg.norm(tail_radial),1e-12)
    side = rel @ tail_radial
    near = (np.abs(ds)<=0.5) & (np.linalg.norm(rel-np.outer(ds,axis),axis=1)<=0.8)
    w = fields["evf"][k]*fields["volumes"]
    groups = {"TAIL_forward":near&(ds>0.05), "TAIL_rear":near&(ds<-.05),
              "near_wall_gap_side":near&(side>0), "opposite_side":near&(side<0)}
    means = {name:weighted_mean(fields["pressure"][k,mask],w[mask]) for name,mask in groups.items()}
    axial_imbalance = means["TAIL_forward"]-means["TAIL_rear"]
    return {"peak_negative_fluid_force_time_ms":peak_time*1e3,
            "matched_field_time_ms":float(fields["time"][k]*1e3), "pressure_means_N_mm2":means,
            "forward_minus_rear_pressure_N_mm2":axial_imbalance,
            "implied_robot_axial_pressure_load":("canonical_-s" if axial_imbalance>0 else "canonical_+s")}


def main():
    coarse = load_case(COARSE); fine = load_case(FINE)
    coarse_closure = force_closure(coarse); fine_closure = force_closure(fine)
    coarse_fields = field_metrics(coarse); fine_fields = field_metrics(fine)
    coarse_recoil = cycle_and_recoil(coarse); fine_recoil = cycle_and_recoil(fine)
    coarse_energy = energy_metrics(coarse); fine_energy = energy_metrics(fine)
    imbalance = pressure_imbalance(fine, fine_fields, fine_closure)

    jfine = fine_closure["critical_8p0_12p5_ms"]["J_fluid_s_Ns_inferred_whole_minus_wall"]
    ratio = jfine / BASELINE_J_FLUID
    magnitude_change = abs(abs(jfine)-abs(BASELINE_J_FLUID))/abs(BASELINE_J_FLUID)
    common_t = fine_fields["time"]
    coarse_momentum = np.interp(common_t, coarse_fields["time"], [row["axial_momentum_mg_mm_s"] for row in coarse_fields["rows"]])
    fine_momentum = np.asarray([row["axial_momentum_mg_mm_s"] for row in fine_fields["rows"]])
    critical = (common_t>=.008)&(common_t<=.0125)
    correlation = float(np.corrcoef(coarse_momentum[critical],fine_momentum[critical])[0,1])
    recoil_improved = (fine_recoil["minimum_v_s_mm_s"] > 0.8*coarse_recoil["minimum_v_s_mm_s"]
                       and fine_recoil["negative_duration_ms"] < .8*coarse_recoil["negative_duration_ms"])
    pressure_drop = (max(row["pressure_p95_abs_N_mm2"] for row in fine_fields["rows"])
                     < .7*max(row["pressure_p95_abs_N_mm2"] for row in coarse_fields["rows"]))
    if magnitude_change <= .20 and correlation >= .70:
        classification = "TAIL_SQUEEZE_IMPULSE_MESH_ROBUST"
    elif abs(jfine) < .70*abs(BASELINE_J_FLUID) and recoil_improved and pressure_drop:
        classification = "TAIL_GAP_UNDERRESOLUTION_DOMINANT"
    else:
        classification = "TAIL_GAP_MESH_SENSITIVE"

    def status(new, old, lower_is_better=True, tolerance=.10):
        if abs(old)<1e-30: return "UNCHANGED"
        change=(new-old)/abs(old)
        if abs(change)<=tolerance: return "UNCHANGED"
        improved = change<0 if lower_is_better else change>0
        return "IMPROVED" if improved else "WORSE"
    recoil_status = status(abs(fine_recoil["minimum_v_s_mm_s"]),abs(coarse_recoil["minimum_v_s_mm_s"]))
    coarse_episodes=coarse["contact"]["TAIL"]["episodes"]+coarse["contact"]["HEAD"]["episodes"]
    fine_episodes=fine["contact"]["TAIL"]["episodes"]+fine["contact"]["HEAD"]["episodes"]
    contact_status = status(float(fine_episodes),float(coarse_episodes))
    numerics_status = status(fine_energy["ETOTAL"]["max_abs_drift_Nmm"],coarse_energy["ETOTAL"]["max_abs_drift_Nmm"])

    numerical_gates = {
        "EVF_retention_ge_0p98": bool(fine_fields["retention"] >= 0.98),
        "EVF_retention_preferred_ge_0p99": bool(fine_fields["retention"] >= 0.99),
        "max_abs_ETOTAL_drift_Nmm": fine_energy["ETOTAL"]["max_abs_drift_Nmm"],
        "Mach_p95": fine_fields["global_speed_p95_max_mm_s"] / float(fine["identity"]["fluid_EOS_c0_mm_s"]),
        "stable_timestep_s": fine_energy["stable_timestep_s"],
        "max_penetration_mm": fine["contact"]["max_penetration_mm"],
        "BOTH_bridge_ms": fine["contact"]["both_bridge_ms"],
        "tumble": bool(np.max(np.abs(fine["rocking"])) >= 45.0),
    }

    setup=json.loads((OUT/f"{FINE}_Setup_Audit.json").read_text())
    result={
        "case":FINE, "physical_parent":COARSE,
        "mesh":setup,
        "impulse_closure":{"coarse":coarse_closure["critical_8p0_12p5_ms"],
                           "refined":fine_closure["critical_8p0_12p5_ms"],
                           "R_J":ratio,"absolute_magnitude_change_fraction":magnitude_change,
                           "method":"robot-fluid force inferred as whole-robot General Contact minus direct pair-isolated robot-wall"},
        "local_tail_metrics":{"coarse":coarse_fields["rows"],"refined":fine_fields["rows"],
                              "refined_strongest_negative_axial_momentum":fine_fields["strongest_negative_momentum"],
                              "critical_momentum_history_correlation":correlation,
                              "pressure_imbalance_at_peak":imbalance},
        "axial_recoil":{"coarse":coarse_recoil,"refined":fine_recoil},
        "contact":{"coarse":coarse["contact"],"refined":fine["contact"]},
        "energy":{"coarse":coarse_energy,"refined":fine_energy},
        "numerical_gates":numerical_gates,
        "fluid_mass":{"coarse_initial_mg":coarse_fields["inventory"][0]*float(coarse["identity"]["fluid_density_tonne_mm3"])*1e9,
                      "coarse_final_mg":coarse_fields["inventory"][-1]*float(coarse["identity"]["fluid_density_tonne_mm3"])*1e9,
                      "coarse_EVF_retention":coarse_fields["retention"],
                      "refined_initial_mg":fine_fields["inventory"][0]*float(fine["identity"]["fluid_density_tonne_mm3"])*1e9,
                      "refined_final_mg":fine_fields["inventory"][-1]*float(fine["identity"]["fluid_density_tonne_mm3"])*1e9,
                      "refined_EVF_retention":fine_fields["retention"]},
        "runtime":{"coarse_wallclock_s":float(coarse["identity"]["wallclock_s"]),
                   "refined_wallclock_s":float(fine["identity"]["wallclock_s"]),
                   "ratio":float(fine["identity"]["wallclock_s"])/float(coarse["identity"]["wallclock_s"])},
        "classification":classification,"RECOIL":recoil_status,"CONTACT":contact_status,"NUMERICS":numerics_status,
        "no_physics_tuning":True,"dynamics_run_count":int(fine["identity"]["dynamics_run_count"]),"third_cycle":False,
    }
    (OUT/f"{FINE}_Mesh_Convergence.json").write_text(json.dumps(result,indent=2)+"\n",encoding="ascii")
    with (OUT/f"{FINE}_TAIL_LOCAL_HISTORY.csv").open("w",newline="",encoding="ascii") as handle:
        writer=csv.DictWriter(handle,fieldnames=fine_fields["rows"][0].keys());writer.writeheader();writer.writerows(fine_fields["rows"])

    lines=[f"# {FINE}: local TAIL-gap mesh diagnostic","",f"Primary classification: **{classification}**","",
           f"RECOIL: **{recoil_status}**  ",f"CONTACT: **{contact_status}**  ",f"NUMERICS: **{numerics_status}**","",
           "## Critical impulse closure","",f"- baseline J_fluid_s: `{BASELINE_J_FLUID:+.7e} N s`",
           f"- refined J_fluid_s: `{jfine:+.7e} N s`",f"- R_J: `{ratio:.6f}`",
           f"- refined closure error: `{fine_closure['critical_8p0_12p5_ms']['relative_closure_error']:.3%}`",
           "- robot-fluid force remains inferred (whole-robot contact minus direct robot-wall), not native pair-isolated output.","",
           "## Mesh and initialization","",f"- local refinement type: `{fine['identity']['mesh_refinement_type']}`",
           f"- minimum cell size: `{setup['minimum_cell_size_mm']:.6f} mm`; target refined cells: `{setup['target_corridor_refined_elements']}`",
           f"- old/new elements: `{setup['old_element_count']} / {setup['new_element_count']}`; max local aspect ratio `{setup['maximum_local_aspect_ratio']:.3f}`",
           f"- refined EVF retention: `{fine_fields['retention']:.6f}`; Mach p95 `{numerical_gates['Mach_p95']:.6g}`",
           f"- ETOTAL drift: `{numerical_gates['max_abs_ETOTAL_drift_Nmm']:.6g} N mm`; stable dt min `{numerical_gates['stable_timestep_s']['min']:.6g} s`",
           f"- penetration/BOTH/tumble: `{numerical_gates['max_penetration_mm']:.6g} mm` / `{numerical_gates['BOTH_bridge_ms']:.6g} ms` / `{numerical_gates['tumble']}`","",
           "## Pressure mechanism","",f"- peak inferred negative fluid force at `{imbalance['peak_negative_fluid_force_time_ms']:.6f} ms`",
           f"- forward-minus-rear pressure: `{imbalance['forward_minus_rear_pressure_N_mm2']:+.6e} N/mm2`",
           f"- implied pressure-load direction: `{imbalance['implied_robot_axial_pressure_load']}`","",
           "Exactly one two-cycle dynamics was run. No magnetic, contact, EOS, boundary, material, or timestep-scale parameter was changed; no third cycle was run."]
    (OUT/f"{FINE}_Mesh_Convergence.md").write_text("\n".join(lines)+"\n",encoding="ascii")
    print(json.dumps({k:result[k] for k in ("classification","RECOIL","CONTACT","NUMERICS")},indent=2))


if __name__=="__main__":
    main()
