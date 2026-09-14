"""Post-process the unique L2300 run with frozen-CAD gap authority."""
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import maximum_filter1d, minimum_filter1d
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT = REPO.parent
AUDIT = REPO / "calibration_analysis" / "ReducedHydro_hidden_impact_audit"
PROBE = AUDIT / "normal_contact_damping_probe"
sys.path[:0] = [str(AUDIT), str(PROBE)]
from audit_stage_a import mesh_properties
from exact_gap_audit import Wall, tests as gap_tests
from contact_probe_common import contact_force, dense_rp
from audit_frozen_cad_static import (EXPECTED_SHA256, STEP, clusters,
                                     certified_near_wall_query, exact_query,
                                     rotation_x_to_axis, surface_points)

JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L2300_D0815_WallSupported_0083"
DT = 1e-7
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
A0 = np.array([0.9647382600216, -0.1188742372140, 0.2348382536499])
A0 /= np.linalg.norm(A0)
CAD_COM = np.array([1.2060186937156343, 0.0, 0.0])
CAD_LENGTH = 2.3
HEAD_POLES = np.asarray(((0.0, 0.0), (0.0, .09664525), (.005, .18497527),
                         (.17271875, .4075), (.25, .4075), (.26149477, .4075)))
PROFILE = pd.read_csv(STEP.with_name("Robot_L2300_D0815_WallWobble_profile.csv"))


def unit(values):
    values = np.asarray(values, float)
    return values / np.maximum(np.linalg.norm(values, axis=-1, keepdims=True), 1e-30)


def interp_vector(frame, t, columns, time_column="t_s"):
    source_t = frame[time_column].to_numpy(float)
    assert np.all(np.diff(source_t) > 0)
    return np.column_stack([np.interp(t, source_t, frame[c]) for c in columns])


def canonical_projection(points, curve):
    xyz = curve[["x_mm", "y_mm", "z_mm"]].to_numpy(float)
    arc = curve.arclength_mm.to_numpy(float)
    edge = np.diff(xyz, axis=0)
    edge2 = np.einsum("ij,ij->i", edge, edge)
    tangents = unit(edge)
    rows = []
    previous = None
    for point in points:
        candidates = (np.arange(len(edge)) if previous is None else
                      np.arange(max(0, previous - 10), min(len(edge), previous + 11)))
        origin = xyz[candidates]
        fraction = np.einsum("ij,ij->i", point-origin, edge[candidates]) / edge2[candidates]
        fraction = np.clip(fraction, 0.0, 1.0)
        projected = origin + fraction[:, None] * edge[candidates]
        distance2 = np.einsum("ij,ij->i", point-projected, point-projected)
        best = int(np.argmin(distance2))
        if distance2[best] > 1.0:
            candidates = np.arange(len(edge)); origin = xyz[:-1]
            fraction = np.clip(np.einsum("ij,ij->i", point-origin, edge)/edge2, 0.0, 1.0)
            projected = origin + fraction[:, None] * edge
            distance2 = np.einsum("ij,ij->i", point-projected, point-projected)
            best = int(np.argmin(distance2))
        segment = int(candidates[best])
        s = arc[segment] + fraction[best] * (arc[segment+1]-arc[segment])
        rows.append((s, *projected[best], *tangents[segment], segment))
        previous = segment
    return pd.DataFrame(rows, columns=["raw_s_mm", "proj_x", "proj_y", "proj_z",
                                       "raw_t_x", "raw_t_y", "raw_t_z", "segment"])


def local_frames(tangent):
    e1 = np.empty_like(tangent)
    for i, direction in enumerate(tangent):
        candidate = np.array([0., 0., 1.]) - direction[2]*direction
        if np.linalg.norm(candidate) < 1e-8:
            candidate = np.array([0., 1., 0.]) - direction[1]*direction
        candidate /= np.linalg.norm(candidate)
        if i and np.dot(candidate, e1[i-1]) < 0:
            candidate *= -1
        e1[i] = candidate
    return e1, unit(np.cross(tangent, e1))


def robust_slope(t, values):
    x = t-t.mean()
    slope, intercept = np.polyfit(x, values, 1)
    residual = values-(slope*x+intercept)
    mad = np.median(np.abs(residual-np.median(residual)))
    keep = np.ones(len(x), bool) if mad < 1e-14 else np.abs(residual-np.median(residual)) <= 4*1.4826*mad
    if keep.sum() >= max(3, int(.7*len(x))):
        slope = np.polyfit(x[keep], values[keep], 1)[0]
    return float(slope)


def intervals(t, flag):
    edges = np.diff(np.r_[False, flag, False].astype(int))
    starts = np.where(edges == 1)[0]; ends = np.where(edges == -1)[0]-1
    rows = [(int(a), int(b), float(t[a]), float(t[b]), float((b-a+1)*np.median(np.diff(t))))
            for a, b in zip(starts, ends)]
    return rows


def phase_plateau(t, phase, field_phase, width=.0005, tolerance_deg=10.0):
    count = int(round(width/np.median(np.diff(t))))+1
    spread = maximum_filter1d(phase, count, mode="nearest")-minimum_filter1d(phase, count, mode="nearest")
    field_moves = np.abs(np.gradient(field_phase, t)) > math.radians(10.0)/width
    flag = (np.abs(spread) < math.radians(tolerance_deg)) & field_moves
    return flag, np.degrees(spread), intervals(t, flag)


def bezier(parameters):
    degree = len(HEAD_POLES)-1
    basis = np.column_stack([math.comb(degree, i)*(1-parameters)**(degree-i)*parameters**i
                             for i in range(degree+1)])
    return basis@HEAD_POLES


def refined_local_points(local, region, coarse_gaps, region_names=("HEAD", "BODY", "TAIL")):
    hp = PROFILE[PROFILE.region == "HEAD"]
    points = []; labels = []
    for name in region_names:
        ids = np.flatnonzero(region == name)
        seed = local[ids[int(np.argmin(coarse_gaps[ids]))]]
        phi0 = math.atan2(seed[2], seed[1])
        phi = phi0 + np.radians(np.arange(-2.0, 2.0001, .5))
        if name == "HEAD":
            k = int(np.argmin((hp.x_mm.to_numpy()-seed[0])**2 +
                              (hp.radius_mm.to_numpy()-math.hypot(seed[1], seed[2]))**2))
            curve = bezier(np.clip(np.arange(float(hp.parameter.iloc[k])-.016,
                                               float(hp.parameter.iloc[k])+.016001, .002), 0, 1))
            candidate = np.column_stack((np.repeat(curve[:, 0], len(phi)),
                np.repeat(curve[:, 1], len(phi))*np.tile(np.cos(phi), len(curve)),
                np.repeat(curve[:, 1], len(phi))*np.tile(np.sin(phi), len(curve))))
        elif name == "BODY":
            x = np.clip(np.arange(seed[0]-.032, seed[0]+.032001, .002), HEAD_POLES[-1, 0], CAD_LENGTH)
            candidate = np.column_stack((np.repeat(x, len(phi)), .4075*np.tile(np.cos(phi), len(x)),
                                         .4075*np.tile(np.sin(phi), len(x))))
        else:
            radius = np.clip(np.arange(math.hypot(seed[1], seed[2])-.022,
                                       math.hypot(seed[1], seed[2])+.022001, .002), 0, .4075)
            candidate = np.column_stack((np.full(len(radius)*len(phi), CAD_LENGTH),
                np.repeat(radius, len(phi))*np.tile(np.cos(phi), len(radius)),
                np.repeat(radius, len(phi))*np.tile(np.sin(phi), len(radius))))
        points.append(candidate); labels.extend([name]*len(candidate))
    return np.vstack(points), np.asarray(labels)


def pose_points(local, displacement, dynamic_rotation, initial_rotation):
    initial = RP0 + (local-CAD_COM)@initial_rotation.T
    return RP0 + displacement + dynamic_rotation.apply(initial-RP0)


def sector_values(gaps, triangles, regions, wall, e1, e2, threshold_um):
    result = {}
    selected = gaps*1e3 <= threshold_um
    all_angles = np.arctan2(wall.normals[triangles[selected]]@e2,
                            wall.normals[triangles[selected]]@e1) if selected.any() else np.array([])
    centers = clusters(all_angles)
    separation = max((math.degrees(math.acos(np.clip(math.cos(a-b), -1, 1)))
                      for i, a in enumerate(centers) for b in centers[i+1:]), default=0.0)
    result["cluster_count"] = len(centers)
    result["separation_deg"] = separation
    result["opposing"] = bool(len(centers) >= 2 and separation >= 120.0)
    for name in ("HEAD", "TAIL"):
        mask = selected & (regions == name)
        if mask.any():
            angles = np.arctan2(wall.normals[triangles[mask]]@e2, wall.normals[triangles[mask]]@e1)
            result[name.lower()+"_sector_deg"] = float(np.degrees(np.angle(np.mean(np.exp(1j*angles)))) % 360)
        else:
            result[name.lower()+"_sector_deg"] = np.nan
    return result


def exact_cad_trajectory(t, data, rotations, e1, e2, active):
    local, region, cad_com, _ = surface_points()
    assert np.linalg.norm(cad_com-CAD_COM) < 1e-10
    wall = Wall(); initial_rotation = rotation_x_to_axis(A0)
    # The global locator retains every fourth azimuthal sample (12 deg, at
    # most 2.24 um circumferential sag on the 0.4075 mm radius). Regional
    # minima are then refined on 0.25 deg x 1 um local parameter grids.
    locator_parts = []
    for name in ("HEAD", "BODY", "TAIL"):
        ids = np.flatnonzero(region == name)
        key = local[ids, 0] if name != "TAIL" else np.linalg.norm(local[ids, 1:], axis=1)
        levels = np.unique(np.round(key, 10))
        level_stride = 10 if name != "TAIL" else 3
        selected_levels = np.r_[levels[::level_stride], levels[-1]]
        for level in np.unique(selected_levels):
            ring = ids[np.isclose(key, level, atol=1e-9)]
            locator_parts.append(ring[::4])
    locator = np.unique(np.concatenate(locator_parts))
    locator_local = local[locator]; locator_region = region[locator]
    coarse = np.arange(0, len(t), 20, dtype=int)
    if coarse[-1] != len(t)-1: coarse = np.r_[coarse, len(t)-1]
    rows = []; seeds = {}
    for number, i in enumerate(coarse):
        points = pose_points(locator_local, data["U"][i], rotations[i], initial_rotation)
        gaps, triangles, _ = exact_query(wall, points)
        refined, refined_region = refined_local_points(locator_local, locator_region, gaps)
        refined_global = pose_points(refined, data["U"][i], rotations[i], initial_rotation)
        rgaps, rtriangles, _ = exact_query(wall, refined_global)
        all_gaps = np.r_[gaps, rgaps]; all_triangles = np.r_[triangles, rtriangles]
        all_region = np.r_[locator_region, refined_region]
        row = {"time_s": t[i], "increment": int(i), "sampling": "coarse_2us"}
        for name in ("HEAD", "BODY", "TAIL"):
            mask = all_region == name
            row[name+"_gap_um"] = float(all_gaps[mask].min()*1e3)
        row["min_gap_um"] = min(row[x+"_gap_um"] for x in ("HEAD", "BODY", "TAIL"))
        for threshold in (20.0, 10.0, 5.0, 0.0):
            sv = sector_values(all_gaps, all_triangles, all_region, wall, e1[i], e2[i], threshold)
            tag = str(int(threshold))
            row["clusters_"+tag+"um"] = sv["cluster_count"]
            row["sector_separation_"+tag+"um_deg"] = sv["separation_deg"]
            row["opposing_bridge_"+tag+"um"] = sv["opposing"]
            if threshold == 20.0:
                row["HEAD_sector_deg"] = sv["head_sector_deg"]
                row["TAIL_sector_deg"] = sv["tail_sector_deg"]
        rows.append(row)
        seeds[int(i)] = gaps
        if number % 500 == 0:
            print("EXACT_CAD_COARSE=%d/%d" % (number, len(coarse)), flush=True)

    # Resolve CAD minima at every solver-active increment using the nearest 1 us seed.
    dense_rows = []
    dense_indices = np.flatnonzero(active)
    for number, i in enumerate(dense_indices):
        nearest = int(coarse[np.argmin(np.abs(coarse-i))])
        seed_gaps = seeds[nearest]
        nearest_region = min(("HEAD", "BODY", "TAIL"),
                             key=lambda name: seed_gaps[locator_region == name].min())
        refined, refined_region = refined_local_points(locator_local, locator_region, seed_gaps,
                                                       region_names=(nearest_region,))
        points = pose_points(refined, data["U"][i], rotations[i], initial_rotation)
        gaps, triangles, _ = exact_query(wall, points)
        row = {"time_s": t[i], "increment": int(i), "sampling": "contact_0p1us"}
        for name in ("HEAD", "BODY", "TAIL"):
            mask = refined_region == name
            row[name+"_gap_um"] = float(gaps[mask].min()*1e3) if mask.any() else np.nan
        row["min_gap_um"] = float(np.nanmin([row[x+"_gap_um"] for x in ("HEAD", "BODY", "TAIL")]))
        dense_rows.append(row)
        if number % 5000 == 0:
            print("EXACT_CAD_CONTACT=%d/%d" % (number, len(dense_indices)), flush=True)
    frame = pd.concat((pd.DataFrame(rows), pd.DataFrame(dense_rows)), ignore_index=True)
    frame = frame.sort_values(["increment", "sampling"]).drop_duplicates("increment", keep="first")
    frame.to_csv(HERE/"L2300_exact_gap.csv", index=False)
    return frame, pd.DataFrame(rows), wall


def event_contact_kinematics(index, start, end, region_name, local, region, wall, initial_rotation,
                             t, data, rotations, com, vcom):
    region_local = local[region == region_name][::4]
    global_points = pose_points(region_local, data["U"][index], rotations[index], initial_rotation)
    gaps, _, _ = exact_query(wall, global_points)
    refined, _ = refined_local_points(region_local, np.repeat(region_name, len(region_local)), gaps,
                                      region_names=(region_name,))
    refined_global = pose_points(refined, data["U"][index], rotations[index], initial_rotation)
    refined_gaps, triangles, _ = exact_query(wall, refined_global)
    k = int(np.argmin(refined_gaps)); local_point = refined[k]; triangle = int(triangles[k])
    initial_point = RP0+(local_point-CAD_COM)@initial_rotation.T
    normal = wall.normals[triangle]
    def velocity(indices):
        points = RP0+data["U"][indices]+rotations[indices].apply(np.broadcast_to(initial_point-RP0, (len(indices),3)))
        return vcom[indices]+np.cross(data["VR"][indices], points-com[indices])
    pre_indices = np.arange(max(0, start-5), start)
    post_indices = np.arange(end+1, min(len(t), end+6))
    vn_in = float(np.mean(velocity(pre_indices)@normal)) if len(pre_indices) else np.nan
    vn_out = float(np.mean(velocity(post_indices)@normal)) if len(post_indices) else np.nan
    point = refined_global[k]
    return vn_in, vn_out, triangle, point, float(refined_gaps[k]*1e3)


def main():
    assert gap_tests().startswith("PASS")
    assert hashlib.sha256(STEP.read_bytes()).hexdigest() == EXPECTED_SHA256
    t, data = dense_rp(HERE/"candidate_private")
    assert len(t) == 83331 and np.allclose(np.diff(t), DT, atol=2e-12, rtol=0)
    meshes, rp, _ = mesh_properties((HERE/(JOB+".inp")).read_text())
    mesh = meshes["Robot_SOLID"]
    assert np.linalg.norm(np.asarray(rp)-RP0) < 1e-9
    rotations = Rotation.from_rotvec(data["UR"])
    initial_rotation = rotation_x_to_axis(A0)
    head0 = RP0 + (np.array([0., 0., 0.])-CAD_COM)@initial_rotation.T
    tail0 = RP0 + (np.array([CAD_LENGTH, 0., 0.])-CAD_COM)@initial_rotation.T
    head = RP0+data["U"]+rotations.apply(np.broadcast_to(head0-RP0, data["U"].shape))
    tail = RP0+data["U"]+rotations.apply(np.broadcast_to(tail0-RP0, data["U"].shape))
    axis = unit(tail-head)
    offset = rotations.apply(np.broadcast_to(mesh["com"]-RP0, data["U"].shape))
    com = RP0+data["U"]+offset
    vcom = data["V"]+np.cross(data["VR"], offset)

    curve = pd.read_csv(ROOT/"CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_TRUE_centerline_odb.csv")
    projection = canonical_projection(com, curve)
    raw_tangent = projection[["raw_t_x", "raw_t_y", "raw_t_z"]].to_numpy(float)
    forward_sign = 1.0 if np.dot(axis[0], raw_tangent[0]) >= 0 else -1.0
    tangent = forward_sign*raw_tangent
    canonical_s = forward_sign*projection.raw_s_mm.to_numpy(float)
    e1, e2 = local_frames(tangent)
    dot = np.einsum("ij,ij->i", axis, tangent)
    directed = np.degrees(np.arccos(np.clip(dot, -1, 1)))
    folded = np.degrees(np.arccos(np.clip(np.abs(dot), 0, 1)))
    transverse_axis = axis-dot[:, None]*tangent
    robot_phase = np.unwrap(np.arctan2(np.einsum("ij,ij->i", transverse_axis, e2),
                                       np.einsum("ij,ij->i", transverse_axis, e1)))

    telemetry = pd.read_csv(ROOT/(JOB+"_telemetry.csv"))
    B = interp_vector(telemetry, t, ["Bx_aba_T", "By_aba_T", "Bz_aba_T"])
    fmag = interp_vector(telemetry, t, ["fx_aba_N", "fy_aba_N", "fz_aba_N"])
    tmag = interp_vector(telemetry, t, ["tx_aba_Nmm", "ty_aba_Nmm", "tz_aba_Nmm"])
    B_phase = np.unwrap(np.arctan2(np.einsum("ij,ij->i", B, e2), np.einsum("ij,ij->i", B, e1)))
    cmd_phase = np.interp(t, telemetry.t_s, telemetry.instantaneous_phase_rad)
    plateau, phase_range, plateau_events = phase_plateau(t, robot_phase, B_phase)

    hydro = pd.read_csv(ROOT/(JOB+"_hydro_increment.csv"), skiprows=2,
        names=["time_s", "v1", "v2", "v3", "w1", "w2", "w3", "Fh1", "Fh2", "Fh3", "Th1", "Th2", "Th3"])
    fhyd = interp_vector(hydro, t, ["Fh1", "Fh2", "Fh3"], "time_s")
    thyd = interp_vector(hydro, t, ["Th1", "Th2", "Th3"], "time_s")
    fn, fs = contact_force(HERE/"candidate_private", t); fwall = fn+fs
    contact_active = np.linalg.norm(fwall, axis=1) > 1e-8
    contact_intervals = intervals(t, contact_active)

    if "--reuse-exact" in sys.argv and (HERE/"L2300_exact_gap.csv").exists():
        exact_all = pd.read_csv(HERE/"L2300_exact_gap.csv")
        exact_coarse = exact_all[exact_all.sampling.str.startswith("coarse")].copy()
        wall = Wall()
    else:
        exact_all, exact_coarse, wall = exact_cad_trajectory(t, data, rotations, e1, e2, contact_active)
    exact_by_increment = exact_all.set_index("increment")

    cad_local, cad_region, _, _ = surface_points()
    event_initial_rotation = rotation_x_to_axis(A0)
    event_rows = []
    for number, (start, end, start_s, end_s, duration_s) in enumerate(contact_intervals, 1):
        scope = exact_by_increment.loc[exact_by_increment.index.intersection(np.arange(start, end+1))]
        if scope.empty:
            nearest = exact_all.iloc[np.argmin(np.abs(exact_all.increment.to_numpy()-start)):]
            scope = nearest.iloc[:1]
        minimum = scope.loc[scope.min_gap_um.idxmin()]
        impulse_scope = slice(max(0, start-1), min(len(t), end+2))
        impulse = np.trapz(fwall[impulse_scope], t[impulse_scope], axis=0)
        region_names = ("HEAD", "BODY", "TAIL")
        region = region_names[int(np.nanargmin([minimum[name+"_gap_um"] for name in region_names]))]
        min_index = int(minimum.name)
        vn_in, vn_out, wall_triangle, point, refined_min_gap = event_contact_kinematics(
            min_index, start, end, region, cad_local, cad_region, wall, event_initial_rotation,
            t, data, rotations, com, vcom)
        short_isolated = duration_s <= 10e-6
        separable = bool(short_isolated and vn_in < 0 and vn_out > 0)
        event_rows.append({"event_number": number, "start_s": start_s, "end_s": end_s,
            "duration_us": duration_s*1e6, "peak_force_N": float(np.linalg.norm(fwall[start:end+1], axis=1).max()),
            "impulse_x_Ns": impulse[0], "impulse_y_Ns": impulse[1], "impulse_z_Ns": impulse[2],
            "impulse_norm_Ns": float(np.linalg.norm(impulse)), "contact_location": region,
            "min_exact_gap_um": min(float(minimum.min_gap_um), refined_min_gap),
            "wall_triangle": wall_triangle, "contact_x_mm": point[0], "contact_y_mm": point[1], "contact_z_mm": point[2],
            "separable_rebound": separable,
            "pre_normal_velocity_mm_s": vn_in, "post_normal_velocity_mm_s": vn_out,
            "effective_restitution": vn_out/-vn_in if separable else np.nan,
            "normal_basis": "actual wall triangle at frozen-CAD material-point minimum",
            "separability_reason": ("short isolated event with closing/opening sign change" if separable else
                                     "no clean short-event closing/opening sign change")})
    event_catalog = pd.DataFrame(event_rows)
    event_catalog.to_csv(HERE/"L2300_contact_events.csv", index=False)

    axis_frame = pd.DataFrame({"time_s": t, "head_x_mm": head[:,0], "head_y_mm": head[:,1], "head_z_mm": head[:,2],
        "tail_x_mm": tail[:,0], "tail_y_mm": tail[:,1], "tail_z_mm": tail[:,2],
        "axis_x": axis[:,0], "axis_y": axis[:,1], "axis_z": axis[:,2],
        "tangent_x": tangent[:,0], "tangent_y": tangent[:,1], "tangent_z": tangent[:,2],
        "e1_x": e1[:,0], "e1_y": e1[:,1], "e1_z": e1[:,2], "e2_x": e2[:,0], "e2_y": e2[:,1], "e2_z": e2[:,2]})
    axis_frame.to_csv(HERE/"L2300_true_axis.csv", index=False)
    pd.DataFrame({"time_s": t, "directed_tilt_deg": directed, "folded_tilt_deg": folded,
                  "axis_dot_tangent": dot}).to_csv(HERE/"L2300_directed_tilt.csv", index=False)
    phase = pd.DataFrame({"time_s": t, "robot_phase_rad": robot_phase,
        "robot_phase_advance_deg": np.degrees(robot_phase-robot_phase[0]), "B_local_phase_rad": B_phase,
        "command_phase_rad": cmd_phase, "robot_minus_B_rad": np.unwrap(robot_phase-B_phase),
        "phase_window_range_deg": phase_range, "phase_plateau": plateau.astype(int)})
    phase.to_csv(HERE/"L2300_phase.csv", index=False)

    coarse_i = exact_coarse.increment.to_numpy(int)
    support = exact_coarse[["time_s", "increment", "HEAD_gap_um", "BODY_gap_um", "TAIL_gap_um", "min_gap_um"]].copy()
    for threshold in (15, 20, 25):
        support[f"HEAD_support_{threshold}um"] = (support.HEAD_gap_um <= threshold).astype(int)
        support[f"TAIL_support_{threshold}um"] = (support.TAIL_gap_um <= threshold).astype(int)
        support[f"either_end_support_{threshold}um"] = ((support.HEAD_gap_um <= threshold)|(support.TAIL_gap_um <= threshold)).astype(int)
        support[f"both_end_support_{threshold}um"] = ((support.HEAD_gap_um <= threshold)&(support.TAIL_gap_um <= threshold)).astype(int)
    support.to_csv(HERE/"L2300_wall_support_timeline.csv", index=False)
    support_summary = []
    for threshold in (15, 20, 25):
        support_summary.append({"threshold_um": threshold,
            "HEAD_supported_fraction": float((support.HEAD_gap_um <= threshold).mean()),
            "TAIL_supported_fraction": float((support.TAIL_gap_um <= threshold).mean()),
            "either_end_supported_fraction": float(((support.HEAD_gap_um <= threshold)|(support.TAIL_gap_um <= threshold)).mean()),
            "both_end_supported_fraction": float(((support.HEAD_gap_um <= threshold)&(support.TAIL_gap_um <= threshold)).mean())})
    pd.DataFrame(support_summary).to_csv(HERE/"L2300_wall_support_15_20_25um.csv", index=False)
    sectors = exact_coarse[["time_s", "increment", "HEAD_sector_deg", "TAIL_sector_deg"]].copy()
    sectors["robot_phase_deg"] = np.degrees(robot_phase[coarse_i])
    sectors["HEAD_supported_20um"] = (exact_coarse.HEAD_gap_um <= 20).astype(int)
    sectors["TAIL_supported_20um"] = (exact_coarse.TAIL_gap_um <= 20).astype(int)
    sectors.to_csv(HERE/"L2300_wall_sector_timeline.csv", index=False)

    bridge = exact_coarse[["time_s", "increment"]].copy()
    for threshold in (20, 10, 5, 0):
        bridge[f"opposing_bridge_{threshold}um"] = exact_coarse[f"opposing_bridge_{threshold}um"].astype(int)
        bridge[f"sector_separation_{threshold}um_deg"] = exact_coarse[f"sector_separation_{threshold}um_deg"]
    bridge["phase_plateau"] = plateau[coarse_i].astype(int)
    bridge["directed_tilt_deg"] = directed[coarse_i]
    bridge.to_csv(HERE/"L2300_bridge_timeline.csv", index=False)

    omega = data["VR"]
    tmag_axis = np.einsum("ij,ij->i", tmag, axis)
    hydro_power = (np.einsum("ij,ij->i", fhyd, vcom)+np.einsum("ij,ij->i", thyd, omega))*1e-3
    mag_power = (np.einsum("ij,ij->i", fmag, vcom)+np.einsum("ij,ij->i", tmag, omega))*1e-3
    torque = pd.DataFrame({"time_s": t, "Tmag_x_Nmm": tmag[:,0], "Tmag_y_Nmm": tmag[:,1], "Tmag_z_Nmm": tmag[:,2],
        "Tmag_norm_Nmm": np.linalg.norm(tmag, axis=1), "Tmag_axis_Nmm": tmag_axis,
        "Tmag_tilt_Nmm": np.linalg.norm(tmag-tmag_axis[:,None]*axis, axis=1),
        "Thydro_x_Nmm": thyd[:,0], "Thydro_y_Nmm": thyd[:,1], "Thydro_z_Nmm": thyd[:,2],
        "Thydro_norm_Nmm": np.linalg.norm(thyd, axis=1), "hydro_power_W": hydro_power})
    torque.to_csv(HERE/"L2300_torque.csv", index=False)

    translation = pd.DataFrame({"time_s": t, "canonical_s_mm": canonical_s,
        "delta_s_mm": canonical_s-canonical_s[0],
        "projection_residual_mm": np.linalg.norm(com-projection[["proj_x","proj_y","proj_z"]].to_numpy(), axis=1),
        "Vt_mm_s": np.einsum("ij,ij->i", vcom, tangent)})
    translation.to_csv(HERE/"L2300_translation.csv", index=False)

    rotmat = rotations.as_matrix(); inertia_global = np.einsum("nij,jk,nlk->nil", rotmat, mesh["I"], rotmat)
    ktrans = .5*mesh["mass"]*np.einsum("ij,ij->i", vcom, vcom)*1e-3
    krot = .5*np.einsum("ni,nij,nj->n", omega, inertia_global, omega)*1e-3
    wmag = cumulative_trapezoid(mag_power, t, initial=0.0)
    whydro = cumulative_trapezoid(hydro_power, t, initial=0.0)
    energy = pd.DataFrame({"time_s": t, "Ktrans_J": ktrans, "Krot_J": krot,
        "Ktotal_rigid_J": ktrans+krot, "Wmag_J": wmag, "Whydro_J": whydro,
        "Wcontact_inferred_J": ktrans+krot-ktrans[0]-krot[0]-wmag-whydro})
    assembly = np.load(HERE/"candidate_private"/"assembly_history_private.npz")
    for key in assembly.files:
        energy["Abaqus_"+key+"_J"] = np.interp(t, assembly[key][:,0], assembly[key][:,1])*1e-3
    energy.to_csv(HERE/"L2300_energy.csv", index=False)

    first_contact = contact_intervals[0][2]
    first_index = contact_intervals[0][0]
    late = t >= .002
    bridge20_dense = np.interp(t, exact_coarse.time_s, exact_coarse.opposing_bridge_20um.astype(float)) >= .5
    bridge_events = intervals(t, bridge20_dense)
    max_bridge = max((row[4] for row in bridge_events), default=0.0)
    max_plateau = max((row[4] for row in plateau_events), default=0.0)
    support_values = pd.DataFrame(support_summary).set_index("threshold_um")
    head_sector = np.unwrap(np.radians(sectors.HEAD_sector_deg.interpolate(limit_direction="both")))
    tail_sector = np.unwrap(np.radians(sectors.TAIL_sector_deg.interpolate(limit_direction="both")))
    sector_movement = max(float(np.degrees(np.ptp(head_sector))), float(np.degrees(np.ptp(tail_sector))))
    post_phase_advance = float(np.degrees(robot_phase[-1]-robot_phase[first_index]))
    late_phase_rate = robust_slope(t[late], robot_phase[late])/(2*math.pi)
    late_field_phase_rate = robust_slope(t[late], B_phase[late])/(2*math.pi)
    post_direction = np.sign(post_phase_advance) if abs(post_phase_advance) > 1e-9 else np.sign(late_field_phase_rate)
    crossed = bool(np.any(dot < 0))
    persistent_jam = bool(max_bridge > .0005 and max_plateau > .0005 and np.ptp(directed[late]) < 10)
    repeated_support = bool(support_values.loc[20, "either_end_supported_fraction"] > .05)
    sector_moves = bool(sector_movement > 30)
    late_phase_reverse = bool(np.sign(late_phase_rate) != post_direction)
    phase_continues = bool(abs(post_phase_advance) > 20 and max_plateau <= .0005 and
                           not late_phase_reverse and late_phase_rate*late_field_phase_rate > 0)
    geometry_success = (not crossed and repeated_support and sector_moves and not persistent_jam and phase_continues)
    ds_total = float(translation.delta_s_mm.iloc[-1])
    ds_2ms = float(translation.delta_s_mm.iloc[-1]-np.interp(.002, t, translation.delta_s_mm))
    late_vt = float(np.median(translation.Vt_mm_s[late]))
    if crossed:
        classification = "L2300_STILL_TOO_SHORT_TUMBLE"
    elif persistent_jam:
        classification = "L2300_TOO_LONG_OR_GEOMETRIC_JAM"
    elif geometry_success and ds_2ms > 0 and late_vt > 0:
        classification = "L2300_WALL_SUPPORTED_FORWARD_WOBBLE"
    else:
        classification = "L2300_WOBBLE_GEOMETRY_SUCCESS_PROPULSION_NOT_RECOVERED"

    runtime = {"job": JOB, "datacheck_run": True, "dynamic_run": True, "run_completed_successfully": True,
        "completed_s": float(t[-1]), "increments": len(t)-1, "dt_s": DT, "field_frames": 335,
        "first_contact_s": first_contact, "contact_event_count": len(contact_intervals),
        "longest_contact_us": float(event_catalog.duration_us.max()), "min_exact_gap_um": float(exact_all.min_gap_um.min()),
        "strict_support_fraction": float(support_values.loc[15,"either_end_supported_fraction"]),
        "nominal_support_fraction": float(support_values.loc[20,"either_end_supported_fraction"]),
        "loose_support_fraction": float(support_values.loc[25,"either_end_supported_fraction"]),
        "HEAD_nominal_support_fraction": float(support_values.loc[20,"HEAD_supported_fraction"]),
        "TAIL_nominal_support_fraction": float(support_values.loc[20,"TAIL_supported_fraction"]),
        "both_end_nominal_support_fraction": float(support_values.loc[20,"both_end_supported_fraction"]),
        "wall_sector_movement_deg": sector_movement, "max_directed_tilt_deg": float(directed.max()),
        "final_directed_tilt_deg": float(directed[-1]), "axis_crossed_90deg": crossed,
        "polarity_reversal": bool(np.any(np.signbit(dot) != np.signbit(dot[0]))),
        "longest_opposing_bridge_ms": max_bridge*1e3, "persistent_jam": persistent_jam,
        "longest_phase_plateau_ms": max_plateau*1e3, "phase_plateau_over_0p5ms": bool(max_plateau > .0005),
        "postimpact_phase_advance_deg": post_phase_advance, "late_phase_rate_Hz": late_phase_rate,
        "late_field_phase_rate_Hz": late_field_phase_rate, "late_phase_reverse": late_phase_reverse,
        "Tmag_RMS_Nmm": float(np.sqrt(np.mean(np.linalg.norm(tmag,axis=1)**2))),
        "Tmag_collapsed": bool(np.sqrt(np.mean(np.linalg.norm(tmag[late],axis=1)**2)) < .5*np.sqrt(np.mean(np.linalg.norm(tmag[~late],axis=1)**2))),
        "hydro_power_max_W": float(hydro_power.max()), "hydro_power_positive_fraction": float(np.mean(hydro_power > 1e-15)),
        "Abaqus_ETOTAL_min_J": float(energy.Abaqus_ETOTAL_J.min()), "Abaqus_ETOTAL_max_J": float(energy.Abaqus_ETOTAL_J.max()),
        "Abaqus_ETOTAL_end_J": float(energy.Abaqus_ETOTAL_J.iloc[-1]), "canonical_delta_s_mm": ds_total,
        "delta_s_2ms_to_end_mm": ds_2ms, "late_median_Vt_mm_s": late_vt,
        "geometry_success": geometry_success, "classification": classification,
        "gap_authority": "frozen CAD parameterization reconstructed from dense RP motion",
        "contact_authority": "Abaqus robot-side CFN+CFS history every increment",
        "mesh_refinement_performed": False, "mesh_elements": 29141}
    pd.DataFrame([runtime]).to_csv(HERE/"L2300_runtime_identity.csv", index=False)
    (HERE/"L2300_dynamic_summary.json").write_text(json.dumps(runtime, indent=2), encoding="utf-8")
    np.savez_compressed(HERE/"candidate_private"/"L2300_visualization_private.npz", t=t, U=data["U"], UR=data["UR"],
        axis=axis, tangent=tangent, e1=e1, e2=e2, B=B, fwall=fwall, contact_active=contact_active,
        robot_phase=robot_phase, B_phase=B_phase, directed_tilt=directed, canonical_delta_s=translation.delta_s_mm.to_numpy())
    print(json.dumps(runtime, indent=2))


if __name__ == "__main__":
    main()
