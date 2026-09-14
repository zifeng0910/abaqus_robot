"""Exact-wall static audit of the frozen L2300 FreeCAD geometry; no Abaqus."""

from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT = REPO.parent
CAD = REPO / "cad" / "freecad_parametric_robot" / "variants" / "L2300_D0815_wallwobble"
AUDIT = REPO / "calibration_analysis" / "ReducedHydro_hidden_impact_audit"
sys.path.insert(0, str(AUDIT))
from exact_gap_audit import Wall, closest

STEP = CAD / "Robot_L2300_D0815_WallWobble.step"
EXPECTED_SHA256 = "461a36dabd0cc1f94390b3e3dc24b2e059b8a74a98d3080ca67694c3a5c0d0e2"
CENTERLINE = ROOT / "CEL_HighEnd83Geom_Z90_XYp2m6_Bias40_Lead5_PolMinus_D055_Forward_Probe006_R014_TRUE_centerline_odb.csv"
RP = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
A0 = np.array([0.9647382600216, -0.1188742372140, 0.2348382536499])
A0 /= np.linalg.norm(A0)
TILTS = (15, 20, 25, 28, 30, 31, 32, 33, 35, 38)
AZIMUTHS = tuple(range(0, 360, 10))
THRESHOLDS_UM = (50.0, 20.0, 10.0, 5.0, 0.0)
HEAD_POLES = np.asarray(((0.0, 0.0), (0.0, .09664525), (.005, .18497527),
                         (.17271875, .4075), (.25, .4075), (.26149477, .4075)))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def unit(vector):
    return vector/np.linalg.norm(vector)


def rotation_x_to_axis(axis):
    x = np.array([1.0, 0.0, 0.0])
    axis = unit(axis)
    v = np.cross(x, axis)
    c = float(np.dot(x, axis))
    if c < -1.0+1e-12:
        return np.diag([1.0, -1.0, -1.0])
    vx = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return np.eye(3)+vx+vx@vx/(1.0+c)


def local_frame():
    curve = pd.read_csv(CENTERLINE)[["x_mm", "y_mm", "z_mm"]].to_numpy(float)
    i = int(np.argmin(np.linalg.norm(curve-RP, axis=1)))
    tangent_station = unit(curve[min(i+2, len(curve)-1)]-curve[max(i-2, 0)])
    forward = -tangent_station if np.dot(A0, tangent_station) < 0.0 else tangent_station
    reference = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(reference, forward)) > .90:
        reference = np.array([0.0, 1.0, 0.0])
    e1 = unit(reference-np.dot(reference, forward)*forward)
    e2 = unit(np.cross(forward, e1))
    return {"station": i, "centerline_point_mm": curve[i], "tangent_station": tangent_station,
            "forward_tangent": forward, "e1": e1, "e2": e2,
            "a0_dot_forward": float(np.dot(A0, forward))}


def surface_points():
    geometry = json.loads((CAD/"Robot_L2300_D0815_WallWobble_geometry.json").read_text())
    props = json.loads((CAD/"Robot_L2300_D0815_WallWobble_mass_properties.json").read_text())
    profile = pd.read_csv(CAD/"Robot_L2300_D0815_WallWobble_profile.csv")
    head = profile[profile.region == "HEAD"].iloc[::20]
    if head.iloc[-1].parameter < 1.0:
        head = pd.concat((head, profile[profile.region == "HEAD"].iloc[[-1]]))
    phi = np.radians(np.arange(0.0, 360.0, 3.0))
    rows = []
    for row in head.itertuples():
        if row.radius_mm < 1e-12:
            rows.append((row.x_mm, 0.0, 0.0, "HEAD"))
        else:
            rows.extend((row.x_mm, row.radius_mm*math.cos(p), row.radius_mm*math.sin(p), "HEAD") for p in phi)
    xh = geometry["head"]["axial_extent_mm"]
    length = geometry["L_total_mm"]
    radius = geometry["R_body_mm"]
    for x in np.linspace(xh, length, 102)[1:-1]:
        rows.extend((x, radius*math.cos(p), radius*math.sin(p), "BODY") for p in phi)
    rows.append((length, 0.0, 0.0, "TAIL"))
    for r in np.linspace(0.0, radius, 32)[1:]:
        rows.extend((length, r*math.cos(p), r*math.sin(p), "TAIL") for p in phi)
    xyz = np.asarray([row[:3] for row in rows], float)
    region = np.asarray([row[3] for row in rows])
    return xyz, region, np.asarray(props["center_of_mass"]["value"], float), geometry


def clusters(angles):
    if not len(angles):
        return []
    values = np.sort(np.mod(angles, 2*math.pi))
    gaps = np.diff(np.r_[values, values[0]+2*math.pi])
    cut = int(np.argmax(gaps))
    ordered = np.r_[values[cut+1:], values[:cut+1]+2*math.pi]
    groups = np.split(ordered, np.where(np.diff(ordered) > math.radians(45.0))[0]+1)
    return [float(np.mod(np.angle(np.mean(np.exp(1j*group))), 2*math.pi)) for group in groups if len(group)]


def bezier(t):
    degree = len(HEAD_POLES)-1
    basis = np.column_stack([math.comb(degree, i)*(1-t)**(degree-i)*t**i
                             for i in range(degree+1)])
    return basis@HEAD_POLES


def sector_metrics(gaps_mm, wall_indices, wall, e1, e2, threshold_um):
    selected = gaps_mm*1e3 <= threshold_um
    if not selected.any():
        return 0, 0.0, False
    normals = wall.normals[wall_indices[selected]]
    angles = np.arctan2(normals@e2, normals@e1)
    centers = clusters(angles)
    separation = max((math.degrees(math.acos(np.clip(math.cos(a-b), -1.0, 1.0)))
                      for i, a in enumerate(centers) for b in centers[i+1:]), default=0.0)
    return len(centers), separation, bool(len(centers) >= 2 and separation >= 120.0)


def exact_query(wall, points):
    """Certified nearest-triangle query with adaptive candidate expansion."""
    points = np.asarray(points, float)
    total = len(wall.tri)
    signed = np.empty(len(points))
    triangle = np.empty(len(points), int)
    nearest = np.empty((len(points), 3))
    remaining = np.arange(len(points))
    k = min(8, total)
    while len(remaining):
        center_distance, indices = wall.tree.query(points[remaining], k=k)
        if k == 1:
            center_distance, indices = center_distance[:, None], indices[:, None]
        tri = wall.tri[indices]
        q, distance = closest(points[remaining, None, :], tri[:, :, 0], tri[:, :, 1], tri[:, :, 2])
        row = np.arange(len(remaining))
        best = np.argmin(distance, axis=1)
        dd = distance[row, best]
        jj = indices[row, best]
        qb = q[row, best]
        resolved = (np.ones(len(remaining), bool) if k == total
                    else center_distance[:, -1]-wall.radius > dd+1e-10)
        output = remaining[resolved]
        triangle[output] = jj[resolved]
        nearest[output] = qb[resolved]
        signed[output] = dd[resolved]*np.sign(np.sum(
            (points[output]-qb[resolved])*wall.normals[jj[resolved]], axis=1))
        remaining = remaining[~resolved]
        k = min(2*k, total)
    return signed, triangle, nearest


def certified_near_wall_query(wall, points, region):
    """Exact threshold counts and lower-bound-certified regional minima."""
    center_distance, nearest_center = wall.tree.query(points, k=1)
    lower_bound = center_distance-wall.radius
    plane_sign_proxy = np.sum(
        (points-wall.centers[nearest_center])*wall.normals[nearest_center], axis=1)
    candidate = (lower_bound <= .050000001) | (plane_sign_proxy <= 0.0)
    for name in ("HEAD", "BODY", "TAIL"):
        ids = np.flatnonzero(region == name)
        candidate[ids[np.argsort(lower_bound[ids])[:128]]] = True

    gaps = np.full(len(points), np.inf)
    triangles = nearest_center.copy()
    while True:
        ids = np.flatnonzero(candidate)
        exact_gap, exact_triangle, _ = exact_query(wall, points[ids])
        gaps[ids] = exact_gap
        triangles[ids] = exact_triangle
        expanded = candidate.copy()
        for name in ("HEAD", "BODY", "TAIL"):
            rid = np.flatnonzero(region == name)
            best_distance = float(np.min(np.abs(gaps[rid][np.isfinite(gaps[rid])])))
            expanded[rid[lower_bound[rid] <= best_distance+1e-10]] = True
        if np.array_equal(expanded, candidate):
            break
        candidate = expanded
    return gaps, triangles, int(candidate.sum())


def refine_region_minima(local, region, coarse_gaps, rotation, wall):
    """Locally refine each CAD parameterization around its coarse minimum."""
    refined_points = []
    refined_regions = []
    full_profile = pd.read_csv(CAD/"Robot_L2300_D0815_WallWobble_profile.csv")
    hp = full_profile[full_profile.region == "HEAD"]
    for name in ("HEAD", "BODY", "TAIL"):
        ids = np.flatnonzero(region == name)
        point = local[ids[int(np.argmin(coarse_gaps[ids]))]]
        phi0 = math.atan2(point[2], point[1])
        phi = phi0+np.radians(np.arange(-2.0, 2.0001, .25))
        if name == "HEAD":
            nearest = int(np.argmin((hp.x_mm.to_numpy()-point[0])**2
                                    +(hp.radius_mm.to_numpy()-math.hypot(point[1], point[2]))**2))
            t0 = float(hp.parameter.iloc[nearest])
            curve = bezier(np.clip(np.arange(t0-.015, t0+.015001, .001), 0.0, 1.0))
            points = np.column_stack((
                np.repeat(curve[:, 0], len(phi)),
                np.repeat(curve[:, 1], len(phi))*np.tile(np.cos(phi), len(curve)),
                np.repeat(curve[:, 1], len(phi))*np.tile(np.sin(phi), len(curve))))
        elif name == "BODY":
            x = np.clip(np.arange(point[0]-.030, point[0]+.030001, .001),
                        HEAD_POLES[-1, 0], 2.300)
            points = np.column_stack((np.repeat(x, len(phi)),
                .4075*np.tile(np.cos(phi), len(x)), .4075*np.tile(np.sin(phi), len(x))))
        else:
            radius0 = math.hypot(point[1], point[2])
            radius = np.clip(np.arange(radius0-.020, radius0+.020001, .001), 0.0, .4075)
            points = np.column_stack((np.full(len(radius)*len(phi), 2.300),
                np.repeat(radius, len(phi))*np.tile(np.cos(phi), len(radius)),
                np.repeat(radius, len(phi))*np.tile(np.sin(phi), len(radius))))
        refined_points.append(points)
        refined_regions.extend([name]*len(points))
    local_refined = np.vstack(refined_points)
    global_refined = RP+(local_refined-np.array([1.2060186937156343, 0.0, 0.0]))@rotation.T
    gaps, triangles, _ = exact_query(wall, global_refined)
    return gaps, triangles, np.asarray(refined_regions)


def pose_row(local, region, cad_com, wall, frame, tilt, azimuth):
    radial = math.cos(math.radians(azimuth))*frame["e1"]+math.sin(math.radians(azimuth))*frame["e2"]
    axis = math.cos(math.radians(tilt))*frame["forward_tangent"]+math.sin(math.radians(tilt))*radial
    rotation = rotation_x_to_axis(axis)
    points = RP+(local-cad_com)@rotation.T
    gaps, wall_indices, queried_count = certified_near_wall_query(wall, points, region)
    refined_gaps, refined_wall_indices, refined_region = refine_region_minima(
        local, region, gaps, rotation, wall)
    all_gaps = np.r_[gaps, refined_gaps]
    all_wall_indices = np.r_[wall_indices, refined_wall_indices]
    all_region = np.r_[region, refined_region]
    row = {
        "tilt_deg": tilt,
        "azimuth_deg": azimuth,
        "axis_x": axis[0], "axis_y": axis[1], "axis_z": axis[2],
        "min_gap_um": float(all_gaps.min()*1e3),
        "HEAD_min_gap_um": float(all_gaps[all_region == "HEAD"].min()*1e3),
        "BODY_min_gap_um": float(all_gaps[all_region == "BODY"].min()*1e3),
        "TAIL_min_gap_um": float(all_gaps[all_region == "TAIL"].min()*1e3),
        "exactly_queried_surface_point_count": queried_count,
    }
    for threshold in THRESHOLDS_UM:
        tag = str(int(threshold))
        selected = gaps*1e3 <= threshold
        count, separation, opposing = sector_metrics(
            all_gaps, all_wall_indices, wall, frame["e1"], frame["e2"], threshold)
        row["point_count_le_%sum" % tag] = int(selected.sum())
        row["sector_clusters_%sum" % tag] = count
        row["sector_separation_%sum_deg" % tag] = separation
        row["opposing_bridge_%sum" % tag] = opposing
    row["penetration_count"] = int((gaps < 0.0).sum())
    end_support = min(row["HEAD_min_gap_um"], row["TAIL_min_gap_um"]) <= 20.0
    row["WALL_SUPPORT"] = bool(end_support and not row["opposing_bridge_20um"])
    return row


def first_discrete(group, column, predicate):
    hit = group[predicate(group[column])]
    return float(hit.tilt_deg.iloc[0]) if len(hit) else float("nan")


def main():
    actual_sha = sha256(STEP)
    if actual_sha != EXPECTED_SHA256:
        raise RuntimeError("Frozen STEP SHA256 mismatch: %s" % actual_sha)
    local, region, cad_com, geometry = surface_points()
    frame = local_frame()
    wall = Wall()
    rows = []
    for tilt in TILTS:
        for azimuth in AZIMUTHS:
            rows.append(pose_row(local, region, cad_com, wall, frame, tilt, azimuth))
        print("STATIC_TILT_COMPLETE=%s" % tilt, flush=True)
    result = pd.DataFrame(rows)
    result.to_csv(HERE/"L2300_static_clearance_map.csv", index=False)

    ceiling = []
    for azimuth, group in result.groupby("azimuth_deg"):
        group = group.sort_values("tilt_deg")
        ceiling.append({
            "azimuth_deg": azimuth,
            "first_end_gap_le_50um_tilt_deg": first_discrete(group, "HEAD_min_gap_um",
                lambda _: np.minimum(group.HEAD_min_gap_um, group.TAIL_min_gap_um) <= 50.0),
            "first_wall_support_tilt_deg": first_discrete(group, "WALL_SUPPORT", lambda x: x.astype(bool)),
            "first_contact_tilt_deg": first_discrete(group, "min_gap_um", lambda x: x <= 0.0),
            "first_opposing_20um_tilt_deg": first_discrete(group, "opposing_bridge_20um", lambda x: x.astype(bool)),
            "first_opposing_10um_tilt_deg": first_discrete(group, "opposing_bridge_10um", lambda x: x.astype(bool)),
            "first_opposing_5um_tilt_deg": first_discrete(group, "opposing_bridge_5um", lambda x: x.astype(bool)),
            "first_opposing_0um_tilt_deg": first_discrete(group, "opposing_bridge_0um", lambda x: x.astype(bool)),
        })
    ceiling = pd.DataFrame(ceiling)
    ceiling.to_csv(HERE/"L2300_static_tilt_ceiling.csv", index=False)
    fractions = result.groupby("tilt_deg").agg(
        wall_support_fraction=("WALL_SUPPORT", "mean"),
        opposing_20um_fraction=("opposing_bridge_20um", "mean"),
        opposing_0um_fraction=("opposing_bridge_0um", "mean"),
        penetrating_pose_fraction=("penetration_count", lambda x: float(np.mean(x > 0))),
        gap_min_um=("min_gap_um", "min"), gap_median_um=("min_gap_um", "median"),
        gap_max_um=("min_gap_um", "max"),
    ).reset_index()
    fractions.to_csv(HERE/"L2300_static_wall_support_fraction_vs_tilt.csv", index=False)

    low = fractions[fractions.tilt_deg <= 20]
    mid = fractions[fractions.tilt_deg.between(28, 33)]
    high = fractions[fractions.tilt_deg >= 35]
    gates = {
        "frozen_STEP_SHA256_matches": True,
        "low_tilt_not_widely_hard_clamped": bool(low.opposing_0um_fraction.max() <= .10),
        "20deg_primarily_clear": bool(fractions.loc[fractions.tilt_deg == 20, "wall_support_fraction"].iloc[0] <= .35),
        "28_33deg_wall_support_present": bool(mid.wall_support_fraction.max() >= .20),
        "35_38deg_wall_limiting_present": bool((high.wall_support_fraction+high.opposing_20um_fraction).max() >= .50),
    }
    if not gates["low_tilt_not_widely_hard_clamped"]:
        classification = "L2300_STATIC_TOO_LONG"
    elif not gates["35_38deg_wall_limiting_present"]:
        classification = "L2300_STATIC_TOO_SHORT"
    elif all(gates.values()):
        classification = "L2300_STATIC_WALL_SUPPORT_WINDOW_PRESENT"
    else:
        classification = "L2300_STATIC_GEOMETRY_INVALID"
    summary = {
        "classification": classification,
        "source_STEP": str(STEP.relative_to(REPO)).replace("\\", "/"),
        "source_STEP_SHA256": actual_sha,
        "surface_representation": "frozen degree-5 Bezier/BRep profile surface; 3-degree azimuth sampling",
        "surface_point_count": int(len(local)),
        "surface_point_count_by_region": {name: int((region == name).sum()) for name in ("HEAD", "BODY", "TAIL")},
        "wall": "SmoothWall114 actual triangles",
        "wall_triangle_count": int(len(wall.tri)),
        "RP_mm": RP.tolist(),
        "local_frame": {key: (value.tolist() if isinstance(value, np.ndarray) else value) for key, value in frame.items()},
        "tilts_deg": list(TILTS), "azimuths_deg": list(AZIMUTHS),
        "counts_are_sampled_surface_point_counts": True,
        "near_wall_does_not_mean_solver_contact": True,
        "gates": gates,
    }
    (HERE/"L2300_static_geometry_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    identity = {
        "source_STEP": summary["source_STEP"], "SHA256": actual_sha,
        "SHA256_matches_manifest": True, "L_total_mm": geometry["L_total_mm"],
        "D_body_mm": geometry["D_body_mm"], "volume_mm3": geometry["master_shape_gate"]["volume_mm3"],
        "classification": classification,
    }
    pd.DataFrame([identity]).to_csv(HERE/"L2300_frozen_step_identity.csv", index=False)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
