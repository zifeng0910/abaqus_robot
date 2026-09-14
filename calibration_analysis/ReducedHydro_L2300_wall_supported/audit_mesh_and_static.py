"""Audit the sole L2300 mesh, build its deck, and run the exact-wall static pose gate."""
from pathlib import Path
import csv
import json
import math
import re
import sys

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
from scipy.special import comb


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT = REPO.parent
OLD = REPO / "calibration_analysis" / "ReducedHydro_geometry_L1800_validation"
FREECAD = REPO / "calibration_analysis" / "ReducedHydro_FreeCAD_L1800_validation"
AUDIT = REPO / "calibration_analysis" / "ReducedHydro_hidden_impact_audit"
sys.path[:0] = [str(FREECAD), str(AUDIT)]
from build_and_audit_freecad import block_for_mesh, exterior, rotation_x_to_axis, tetra_properties
from exact_gap_audit import Wall
from cross_geometry_trajectory_replay import clusters

STEM = "Robot_parametric_L2p300_D0p815_WallWobble"
JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L2300_D0815_WallSupported_0083"
OLD_JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083"
RP = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
DENSITY = 7.80906654321e-9
L_TOTAL = 2.300
R_BODY = 0.4075
TILTS = (20, 25, 28, 30, 31, 32, 33, 35)
AZIMUTHS = tuple(range(0, 360, 10))


def bezier(t, poles):
    degree = len(poles) - 1
    b = np.column_stack([comb(degree, i) * (1-t)**(degree-i) * t**i
                         for i in range(degree + 1)])
    db = degree * np.column_stack([
        comb(degree-1, i) * (1-t)**(degree-1-i) * t**i
        for i in range(degree)
    ])
    return b @ poles, db @ np.diff(poles, axis=0)


def exact_normal(point, poles, head_extent):
    x, y, z = point; rr = math.hypot(y, z)
    if x >= L_TOTAL - 1e-7:
        return np.array([1., 0., 0.]), "tail"
    if x >= head_extent - 1e-8:
        return np.array([0., y, z]) / max(rr, 1e-30), "body"
    tt = np.linspace(0, 1, 10001)
    curve, derivative = bezier(tt, poles)
    i = int(np.argmin(abs(curve[:, 0] - x)))
    dx, dr = derivative[i]
    n = np.array([-dr, dx*y/max(rr, 1e-30), dx*z/max(rr, 1e-30)])
    return n / np.linalg.norm(n), "head"


def rotation_x_to_vector(axis):
    x = np.array([1., 0., 0.]); axis = axis / np.linalg.norm(axis)
    v = np.cross(x, axis); c = np.dot(x, axis)
    vx = np.array([[0., -v[2], v[1]], [v[2], 0., -v[0]], [-v[1], v[0], 0.]])
    return np.eye(3) + vx + vx @ vx / (1. + c)


def bridge_metrics(gaps, triangles, wall, e1, e2, threshold_um):
    selected = gaps * 1e3 <= threshold_um
    angles = (np.arctan2(wall.normals[triangles[selected]] @ e2,
                         wall.normals[triangles[selected]] @ e1)
              if selected.any() else np.array([]))
    centers = clusters(angles)
    separation = max((np.degrees(np.arccos(np.clip(np.cos(a-b), -1, 1)))
                      for i, a in enumerate(centers) for b in centers[i+1:]), default=0.0)
    return len(centers), separation, bool(len(centers) >= 2 and separation >= 120.0)


def pose_row(local_surface, local_x, wall, tangent, e1, e2, tilt, azimuth):
    radial = math.cos(math.radians(azimuth))*e1 + math.sin(math.radians(azimuth))*e2
    axis = -math.cos(math.radians(tilt))*tangent + math.sin(math.radians(tilt))*radial
    rot = rotation_x_to_vector(axis)
    points = RP + local_surface @ rot.T
    gaps, triangles, _ = wall.query(points)
    row = {"tilt_deg": tilt, "azimuth_deg": azimuth,
           "min_gap_um": float(gaps.min()*1e3),
           "head_min_gap_um": float(gaps[local_x <= .30].min()*1e3),
           "tail_min_gap_um": float(gaps[local_x >= L_TOTAL-.05].min()*1e3)}
    for threshold in (50.0, 20.0, 0.0):
        count, separation, bridge = bridge_metrics(gaps, triangles, wall, e1, e2, threshold)
        tag = str(int(threshold)); row[f"sector_count_{tag}um"] = count
        row[f"sector_separation_{tag}um_deg"] = separation
        row[f"opposing_bridge_{tag}um"] = bridge
    return row


def first_crossing(theta, values, threshold=0.0):
    values = np.asarray(values)
    hit = np.flatnonzero(values <= threshold)
    if not len(hit): return float("nan")
    i = int(hit[0])
    if i == 0: return float(theta[0])
    return float(np.interp(threshold, [values[i], values[i-1]], [theta[i], theta[i-1]]))


def main():
    raw = json.loads((HERE / "L2300_mesh060_raw.json").read_text())
    geometry = json.loads((HERE / f"{STEM}_geometry.json").read_text())
    cad = json.loads((HERE / f"{STEM}_mass_properties.json").read_text())
    fit = json.loads((HERE / "L2300_head_fit.json").read_text())
    poles = np.asarray(fit["poles_mm"], float); head_extent = fit["head_axial_extent_mm"]
    ndf = pd.read_csv(HERE / "L2300_mesh060_nodes_local.csv")
    edf = pd.read_csv(HERE / "L2300_mesh060_elements.csv")
    labels = ndf.node.to_numpy(int); local = ndf[["x_mm","y_mm","z_mm"]].to_numpy(float)
    lookup = {label:i for i,label in enumerate(labels)}
    element_labels = edf.element.to_numpy(int)
    elements = np.array([[lookup[int(v)] for v in row]
        for row in edf[["n1","n2","n3","n4"]].to_numpy()], int)
    rotation = rotation_x_to_axis()
    cad_com = np.asarray(cad["center_of_mass"]["value"], float)
    xyz = RP + (local-cad_com) @ rotation.T
    volume, mass, com, inertia, _ = tetra_properties(xyz, elements)
    principal = np.sort(np.linalg.eigvalsh(inertia))
    cad_principal = np.sort(np.asarray(cad["principal_mass_moments"]["value_tonne_mm2"], float))
    faces = exterior(elements)
    normal_rows=[]
    for ei, side, face in faces:
        p = local[list(face)]; center = p.mean(axis=0)
        normal = np.cross(p[1]-p[0], p[2]-p[0]); normal /= np.linalg.norm(normal)
        exact, region = exact_normal(center, poles, head_extent)
        angle = math.degrees(math.acos(np.clip(abs(np.dot(normal, exact)), -1, 1)))
        normal_rows.append({"element":int(element_labels[ei]), "side":side,
                            "region":region, "center_x_mm":float(center[0]),
                            "angle_error_deg":angle})
    normals = pd.DataFrame(normal_rows)
    region_stats = {name:{"triangle_count":int(len(g)),
        "median_deg":float(g.angle_error_deg.median()),
        "P95_deg":float(g.angle_error_deg.quantile(.95)),
        "max_deg":float(g.angle_error_deg.max())} for name,g in normals.groupby("region")}
    surface_index = np.unique(np.asarray([f for _,_,f in faces]).ravel())
    axial = local[surface_index,0]; radial = np.linalg.norm(local[surface_index,1:],axis=1)
    length=float(axial.max()-axial.min()); diameter=float(2*radial.max())
    mass_error=abs(mass-float(cad["mass"]["tonne"]))/float(cad["mass"]["tonne"])
    inertia_error=np.abs(principal-cad_principal)/cad_principal
    com_error=float(np.linalg.norm(com-RP)); p95=float(normals.angle_error_deg.quantile(.95))
    gates = {"element_count_expected":25000 <= len(elements) <= 35000,
             "overall_normal_P95_lt_2deg":p95 < 2.0,
             "head_normal_P95_lt_2deg":region_stats["head"]["P95_deg"] < 2.0,
             "mass_error_lt_0p5pct":mass_error < .005,
             "COM_error_lt_5um":com_error < .005,
             "inertia_error_lt_1p5pct":float(inertia_error.max()) < .015,
             "L_error_lt_0p005mm":abs(length-L_TOTAL) < .005,
             "D_error_lt_0p005mm":abs(diameter-.815) < .005}
    mesh_result={"job":JOB,"mesh_size_mm":.060,"node_count":len(labels),"element_count":len(elements),
        "exterior_node_count":len(surface_index),"exterior_triangle_count":len(faces),
        "length_mm":length,"diameter_mm":diameter,"volume_mm3":volume,"mass_mg":mass*1e9,
        "mass_relative_error":mass_error,"mesh_COM_mm":com.tolist(),"COM_error_um":com_error*1e3,
        "principal_inertia_tonne_mm2":principal.tolist(),"inertia_relative_error":inertia_error.tolist(),
        "surface_normal_P95_deg":p95,"surface_normal_by_region":region_stats,
        "gates":{k:bool(v) for k,v in gates.items()},
        "element_count_expected_is_preferred_not_hard":True,
        "all_mesh_gates_pass":bool(all(v for k,v in gates.items() if k != "element_count_expected"))}
    normals.to_csv(HERE/"L2300_mesh060_surface_normal_audit.csv",index=False)
    pd.DataFrame([mesh_result]).to_csv(HERE/"L2300_mesh_quality.csv",index=False)
    pd.DataFrame([{"quantity":"mass_mg","CAD":cad["mass"]["mg"],"mesh":mass*1e9,"relative_error":mass_error}]
        + [{"quantity":f"principal_inertia_{i+1}_tonne_mm2","CAD":cad_principal[i],"mesh":principal[i],"relative_error":inertia_error[i]} for i in range(3)]
        + [{"quantity":"COM_error_um","CAD":0.0,"mesh":com_error*1e3,"relative_error":np.nan}]
        ).to_csv(HERE/"L2300_mass_inertia.csv",index=False)
    (HERE/"L2300_mesh_gate.json").write_text(json.dumps(mesh_result,indent=2)+"\n")
    if not mesh_result["all_mesh_gates_pass"]:
        print(json.dumps(mesh_result,indent=2)); raise RuntimeError("L2300 mesh gate failed; STOP")

    surface_local = local[surface_index] - cad_com
    local_x = local[surface_index,0]
    frame = pd.read_csv(OLD/"L1800_8p333_true_axis.csv").iloc[0]
    tangent=frame[["tangent_x","tangent_y","tangent_z"]].to_numpy(float)
    e1=frame[["e1_x","e1_y","e1_z"]].to_numpy(float); e2=frame[["e2_x","e2_y","e2_z"]].to_numpy(float)
    wall=Wall(); rows=[]
    for tilt in TILTS:
        for azimuth in AZIMUTHS:
            rows.append(pose_row(surface_local,local_x,wall,tangent,e1,e2,tilt,azimuth))
    static=pd.DataFrame(rows); static.to_csv(HERE/"L2300_static_tilt_azimuth_clearance.csv",index=False)

    ceiling=[]
    fine_theta=np.arange(0.,60.0001,.25)
    for azimuth in AZIMUTHS:
        rr=[pose_row(surface_local,local_x,wall,tangent,e1,e2,t,azimuth) for t in fine_theta]
        gaps=np.array([x["min_gap_um"] for x in rr]); bridge20=np.array([x["opposing_bridge_20um"] for x in rr])
        bridge0=np.array([x["opposing_bridge_0um"] for x in rr])
        ceiling.append({"azimuth_deg":azimuth,
            "theta_near50_first_deg":first_crossing(fine_theta,gaps,50.),
            "theta_near20_first_deg":first_crossing(fine_theta,gaps,20.),
            "theta_contact_first_deg":first_crossing(fine_theta,gaps,0.),
            "theta_opposing_bridge_deg":float(fine_theta[np.flatnonzero(bridge20)[0]]) if bridge20.any() else np.nan,
            "theta_opposing_hard_bridge_deg":float(fine_theta[np.flatnonzero(bridge0)[0]]) if bridge0.any() else np.nan})
    ceiling=pd.DataFrame(ceiling); ceiling.to_csv(HERE/"L2300_geometric_tilt_ceiling.csv",index=False)
    low=static[static.tilt_deg.isin([20,25])]; mid=static[static.tilt_deg.isin([28,30,31,32])]
    high=static[static.tilt_deg.isin([33,35])]
    static_gates={
        "20_25_no_severe_opposing_penetration":bool(not low.opposing_bridge_0um.any()),
        "28_32_some_azimuth_nearwall_0_50um":bool(((mid.min_gap_um>=0)&(mid.min_gap_um<=50)).any()),
        "below_28_not_all_azimuth_opposing_hard":bool(not low.groupby("tilt_deg").opposing_bridge_0um.all().any()),
        "33_35_contact_limited":bool(high.groupby("tilt_deg").apply(lambda g:(g.min_gap_um<=0).any()).all()),
        "opposing_hard_bridge_not_below_25":bool(np.nanmin(ceiling.theta_opposing_hard_bridge_deg) >= 25.0),
    }
    static_result={"wall":"SmoothWall114 exact triangles","COM_mm":RP.tolist(),
        "tilts_deg":list(TILTS),"azimuths_deg":list(AZIMUTHS),
        "threshold_definitions":{"near_wall_um":[0,50],"opposing_bridge":"two wall-normal sectors separated >=120 deg"},
        "theta_contact_first_deg":{"min":float(ceiling.theta_contact_first_deg.min()),"median":float(ceiling.theta_contact_first_deg.median()),"max":float(ceiling.theta_contact_first_deg.max())},
        "theta_opposing_bridge_deg":{"min":float(ceiling.theta_opposing_bridge_deg.min()),"median":float(ceiling.theta_opposing_bridge_deg.median()),"max":float(ceiling.theta_opposing_bridge_deg.max())},
        "gates":static_gates,"all_static_gates_pass":bool(all(static_gates.values()))}
    (HERE/"L2300_static_gate.json").write_text(json.dumps(static_result,indent=2)+"\n")
    if not static_result["all_static_gates_pass"]:
        print(json.dumps(static_result,indent=2)); raise RuntimeError("L2300 static gate failed; STOP before dynamic")

    block=block_for_mesh(labels,xyz,element_labels,elements,faces)
    deck=(OLD/f"{OLD_JOB}.inp").read_text(); start=deck.index("** ROBOT GEOMETRY CANDIDATE:")
    end=deck.index("*End Part",start)+len("*End Part")
    deck=deck[:start]+block+deck[end:]; deck=deck.replace(OLD_JOB,JOB)
    deck=re.sub(r"(\*Instance, name=Robot_SOLID-1, part=Robot_SOLID\n)[^*]+?(\*End Instance)",r"\1\2",deck,count=1)
    deck=re.sub(r"(\*Elset, elset=ROBOT_SOLID_CEL_ALL, instance=Robot_SOLID-1, generate\n)\d+,\s*\d+,\s*1",
                rf"\g<1>{int(element_labels.min())}, {int(element_labels.max())}, 1",deck,count=1)
    (HERE/f"{JOB}.inp").write_text(deck); (ROOT/f"{JOB}.inp").write_text(deck)
    with (HERE/"L2300_mesh060_surface_triangles.csv").open("w",newline="") as stream:
        w=csv.writer(stream);w.writerow(["element","side","n1","n2","n3"])
        for ei,side,face in faces:w.writerow([int(element_labels[ei]),side]+[int(labels[i]) for i in face])
    pd.DataFrame({"node":labels[surface_index],"x_mm":xyz[surface_index,0],"y_mm":xyz[surface_index,1],"z_mm":xyz[surface_index,2]}).to_csv(HERE/"L2300_mesh060_surface_nodes_global.csv",index=False)
    old_id=json.loads((OLD/"L1800_input_identity.json").read_text())
    old_volume=old_id["candidate_mass_mg"]/7.80906654321
    moment_density=old_id["candidate_magnetic_moment_Am2"]/old_volume
    magnetic={"constant_magnetization_rule":True,"reference_volume_mm3":old_volume,
        "reference_moment_Am2":old_id["candidate_magnetic_moment_Am2"],"moment_density_Am2_per_mm3":moment_density,
        "L2300_CAD_volume_mm3":geometry["volume_mm3"],"L2300_magnetic_moment_Am2":moment_density*geometry["volume_mm3"]}
    pd.DataFrame([magnetic]).to_csv(HERE/"L2300_magnetic_moment.csv",index=False)
    (HERE/"L2300_dynamic_eligibility.json").write_text(json.dumps({"job":JOB,"mesh_gate":True,"static_gate":True,
        "eligible_for_unique_dynamic_run":True,"magnetic":magnetic},indent=2)+"\n")
    print(json.dumps({"mesh":mesh_result,"static":static_result,"magnetic":magnetic},indent=2))


if __name__ == "__main__":
    main()
