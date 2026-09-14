"""Build the exact-CAD input deck and enforce all zero-solve geometry gates."""
from pathlib import Path
import csv
import hashlib
import json
import math
import re

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT = REPO.parent
CAD = REPO / 'cad' / 'freecad_parametric_robot'
OLD = REPO / 'calibration_analysis' / 'ReducedHydro_geometry_L1800_validation'
OLD_JOB = 'Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083'
JOB = 'Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_WallOn_Free_0083'
AXIS = np.array([0.9647382600216, -0.1188742372140, 0.2348382536499], float)
AXIS /= np.linalg.norm(AXIS)
RP = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
DENSITY = 7.80906654321e-9
CAD_COM = np.array([0.966978181129355, 0., 0.])
CAD_VOLUME = 0.866162685064448
CAD_MASS = 6.763922044913722e-9
CAD_INERTIA = np.array([5.429470082050118e-10, 1.8566445332010883e-9, 1.8566445332010866e-9])
OLD_MASS_MG = 6.161278579
OLD_MOMENT = 0.000719639311013
R_HEAD = 0.4142916666666666
HEAD_EXTENT = 0.4142359935270392
R_BODY = 0.4075


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rotation_x_to_axis():
    x = np.array([1., 0., 0.]); v = np.cross(x, AXIS); c = np.dot(x, AXIS)
    vx = np.array([[0., -v[2], v[1]], [v[2], 0., -v[0]], [-v[1], v[0], 0.]])
    return np.eye(3) + vx + vx @ vx / (1. + c)


def tetra_properties(nodes, elements):
    verts = nodes[elements]
    signed = np.linalg.det(verts[:, 1:] - verts[:, :1]) / 6.
    vol = np.abs(signed)
    centers = verts.mean(axis=1)
    volume = vol.sum(); mass = volume * DENSITY
    com = np.average(centers, axis=0, weights=vol)
    second = np.zeros((3, 3))
    for vv, v in zip(verts, vol):
        s = vv.sum(axis=0)
        second += v / 20. * (np.outer(s, s) + vv.T @ vv)
    second_com = DENSITY * second - mass * np.outer(com, com)
    inertia = np.eye(3) * np.trace(second_com) - second_com
    return volume, mass, com, inertia, signed


def exterior(elements):
    definitions = [(1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)]
    found = {}
    for ei, tet in enumerate(elements):
        for side, ids in enumerate(definitions, 1):
            face = tuple(int(tet[i]) for i in ids)
            key = tuple(sorted(face))
            if key in found: found[key] = None
            else: found[key] = (ei, side, face)
    return [v for v in found.values() if v is not None]


def exact_normal_local(point):
    x, y, z = point; radial = np.array([0., y, z]); rr = np.linalg.norm(radial)
    if x >= 1.8 - 1e-7:
        return np.array([1., 0., 0.]), 'tail'
    if x >= HEAD_EXTENT - 1e-8:
        return radial / max(rr, 1e-30), 'body'
    xc = HEAD_EXTENT
    radial_center = R_BODY - R_HEAD
    # Normal of the revolved generating circle, including its radial offset.
    n = np.array([x-xc, (rr-radial_center)*y/max(rr,1e-30),
                  (rr-radial_center)*z/max(rr,1e-30)])
    return n / np.linalg.norm(n), 'head'


def wall_gap(points):
    sys_path = str(REPO / 'calibration_analysis' / 'ReducedHydro_hidden_impact_audit')
    import sys
    if sys_path not in sys.path: sys.path.insert(0, sys_path)
    from exact_gap_audit import Wall
    wall = Wall(); gaps, tri, closest = wall.query(points)
    i = int(np.argmin(gaps))
    return float(gaps[i]), i, int(tri[i]), closest[i]


def block_for_mesh(labels, xyz, element_labels, elements, faces):
    side_sets = {i: [] for i in range(1, 5)}
    for ei, side, _ in faces: side_sets[side].append(int(element_labels[ei]))
    lines = ['** ROBOT SOURCE: authoritative FreeCAD STEP; no coordinate scaling', '*Part, name=Robot_SOLID', '*Node']
    lines += [f'{int(label)}, {p[0]:.12g}, {p[1]:.12g}, {p[2]:.12g}' for label, p in zip(labels, xyz)]
    lines.append('*Element, type=C3D4')
    lines += [f'{int(el)}, ' + ', '.join(str(int(labels[i])) for i in tet) for el, tet in zip(element_labels, elements)]
    for side in range(1, 5):
        lines.append(f'*Elset, elset=ROBOT_EXTERIOR_S{side}')
        vals = side_sets[side]
        lines += [', '.join(str(x) for x in vals[i:i+16]) for i in range(0, len(vals), 16)]
    lines += [f'*Nset, nset=ROBOT_CEL_BODY, generate', f'{int(labels.min())}, {int(labels.max())}, 1',
              '*Elset, elset=ROBOT_CEL_BODY, generate', f'{int(element_labels.min())}, {int(element_labels.max())}, 1',
              '*Elset, elset=ROBOT_SOLID_ALL, generate', f'{int(element_labels.min())}, {int(element_labels.max())}, 1',
              '*Surface, type=ELEMENT, name=ROBOT_SOLID_SURF']
    lines += [f'ROBOT_EXTERIOR_S{i}, S{i}' for i in range(1, 5)]
    lines += ['** Section: SEC_ROBOT_RIGID', '*Solid Section, elset=ROBOT_CEL_BODY, material=MAT_ROBOT_RIGID', ',', '*End Part']
    return '\n'.join(lines)


def main():
    raw = json.loads((HERE / 'abaqus_step_import_raw.json').read_text())
    ndf = pd.read_csv(HERE / 'abaqus_step_nodes_local.csv')
    edf = pd.read_csv(HERE / 'abaqus_step_elements.csv')
    labels = ndf.node.to_numpy(int); local = ndf[['x_mm','y_mm','z_mm']].to_numpy(float)
    lookup = {label: i for i, label in enumerate(labels)}
    element_labels = edf.element.to_numpy(int)
    elements = np.array([[lookup[int(x)] for x in row] for row in edf[['n1','n2','n3','n4']].to_numpy()], int)
    R = rotation_x_to_axis(); xyz = RP + (local - CAD_COM) @ R.T
    volume, mass, com, inertia, signed = tetra_properties(xyz, elements)
    eig = np.sort(np.linalg.eigvalsh(inertia)); target = np.sort(CAD_INERTIA)
    faces = exterior(elements)
    normal_rows=[]
    for ei, side, face in faces:
        p = local[list(face)]; center = p.mean(axis=0)
        n = np.cross(p[1]-p[0], p[2]-p[0]); n /= np.linalg.norm(n)
        exact, region = exact_normal_local(center)
        angle = math.degrees(math.acos(np.clip(abs(np.dot(n, exact)), -1., 1.)))
        normal_rows.append({'element':int(element_labels[ei]),'side':side,'region':region,'angle_error_deg':angle})
    normals = pd.DataFrame(normal_rows); normals.to_csv(HERE/'freecad_robot_surface_normal_audit.csv',index=False)
    surface_idx = np.unique(np.array([f for _,_,f in faces]).ravel())
    surf = xyz[surface_idx]
    axial = (surf-RP)@AXIS
    radial = np.linalg.norm((surf-RP)-np.outer(axial,AXIS),axis=1)
    length=float(axial.max()-axial.min()); diameter=float(2*radial.max())
    dot=float(np.dot(R@np.array([1.,0.,0.]),AXIS))
    gap, near_i, wall_tri, closest = wall_gap(surf)
    old_volume=OLD_MASS_MG/7.80906654321
    moment_density=OLD_MOMENT/old_volume; moment=moment_density*CAD_VOLUME
    mass_err=abs(mass-CAD_MASS)/CAD_MASS; inertia_err=np.abs(eig-target)/target
    com_err=np.linalg.norm(com-RP)
    gates = {k: bool(v) for k, v in {
             'step_identity':abs(raw['bbox_extents_mm'][0]-1.8)<1e-6 and max(abs(raw['bbox_extents_mm'][i]-.815) for i in (1,2))<1e-6,
             'head_tail_polarity':dot>0.999999,'com_placement':com_err<.005,
             'surface_faceting':float(normals.angle_error_deg.quantile(.95))<2.5,
             'mesh_dimensions':abs(length-1.8)<=.005 and abs(diameter-.815)<=.005,
             'mesh_mass':mass_err<.005,'mesh_inertia':float(inertia_err.max())<.02,
             'magnetic_moment':abs(moment/OLD_MOMENT-CAD_VOLUME/old_volume)<1e-12,
             'initial_exact_gap':gap>0}.items()}
    old_identity=json.loads((OLD/'L1800_input_identity.json').read_text())
    old_principal=np.sort(np.linalg.eigvalsh(np.asarray(old_identity['candidate_inertia_tonne_mm2'])))
    identity={'job':JOB,'source_step':str(CAD/'Robot_parametric_L1p800_D0p815.step'),'source_step_sha256':sha(CAD/'Robot_parametric_L1p800_D0p815.step'),
      'abaqus_import':raw,'rotation_matrix':R.tolist(),'head_to_tail_dot_a0':dot,'RP_mm':RP.tolist(),'mesh_COM_mm':com.tolist(),'RP_mesh_COM_error_um':com_err*1e3,
      'mesh_length_mm':length,'mesh_diameter_mm':diameter,'mesh_volume_mm3':volume,'mesh_mass_mg':mass*1e9,'mass_relative_error':mass_err,
      'mesh_inertia_principal_tonne_mm2':eig.tolist(),'CAD_inertia_principal_tonne_mm2':target.tolist(),'inertia_relative_error':inertia_err.tolist(),
      'surface_triangle_count':len(faces),'surface_normal_median_deg':float(normals.angle_error_deg.median()),'surface_normal_P95_deg':float(normals.angle_error_deg.quantile(.95)),'surface_normal_max_deg':float(normals.angle_error_deg.max()),
      'initial_exact_gap_um':gap*1e3,'nearest_robot_node':int(labels[surface_idx[near_i]]),'nearest_wall_triangle':wall_tri,'nearest_wall_point_mm':closest.tolist(),
      'old_volume_mm3':old_volume,'old_moment_Am2':OLD_MOMENT,'moment_density_Am2_per_mm3':moment_density,'new_CAD_volume_mm3':CAD_VOLUME,'new_moment_Am2':moment,
      'moment_ratio':moment/OLD_MOMENT,'old_mass_mg':OLD_MASS_MG,
      'old_transverse_inertia_tonne_mm2':float(old_principal[-1]),
      'new_old_transverse_inertia_ratio':float(target[-1]/old_principal[-1]),
      'new_old_moment_per_Iperp_ratio':float((moment/target[-1])/(OLD_MOMENT/old_principal[-1])),
      'gates':gates}
    identity['all_pre_datacheck_gates_pass']=all(gates.values())
    (HERE/'freecad_preflight_identity.json').write_text(json.dumps(identity,indent=2))
    pd.DataFrame([{'metric':k,'value':v} for k,v in identity.items() if np.isscalar(v)]).to_csv(HERE/'freecad_cad_import_identity.csv',index=False)
    pd.DataFrame([{'metric':'length_mm','CAD':1.8,'Abaqus_mesh':length},{'metric':'diameter_mm','CAD':.815,'Abaqus_mesh':diameter},{'metric':'volume_mm3','CAD':CAD_VOLUME,'Abaqus_mesh':volume}]).to_csv(HERE/'freecad_vs_abaqus_geometry_audit.csv',index=False)
    pd.DataFrame([{'metric':'mass_mg','CAD':CAD_MASS*1e9,'Abaqus_mesh':mass*1e9,'relative_error':mass_err}]+[{'metric':f'inertia_{i}','CAD':target[i],'Abaqus_mesh':eig[i],'relative_error':inertia_err[i]} for i in range(3)]).to_csv(HERE/'freecad_vs_abaqus_mass_inertia.csv',index=False)
    pd.DataFrame([{'RP_x':RP[0],'RP_y':RP[1],'RP_z':RP[2],'mesh_COM_x':com[0],'mesh_COM_y':com[1],'mesh_COM_z':com[2],'error_um':com_err*1e3,'axis_dot':dot}]).to_csv(HERE/'freecad_robot_placement_audit.csv',index=False)
    pd.DataFrame([{'old_volume_mm3':old_volume,'old_moment_Am2':OLD_MOMENT,'moment_density_Am2_per_mm3':moment_density,'new_volume_mm3':CAD_VOLUME,'new_moment_Am2':moment,'ratio':moment/OLD_MOMENT}]).to_csv(HERE/'freecad_magnetic_moment_scaling.csv',index=False)
    pd.DataFrame([{'gap_um':gap*1e3,'robot_node':int(labels[surface_idx[near_i]]),'wall_triangle':wall_tri}]).to_csv(HERE/'freecad_initial_gap_audit.csv',index=False)
    if not identity['all_pre_datacheck_gates_pass']:
        raise RuntimeError('Preflight gate failed: '+json.dumps(gates))
    old_text=(OLD/f'{OLD_JOB}.inp').read_text()
    block=block_for_mesh(labels,xyz,element_labels,elements,faces)
    start=old_text.index('** ROBOT GEOMETRY CANDIDATE:')
    end=old_text.index('*End Part',start)+len('*End Part')
    deck=old_text[:start]+block+old_text[end:]
    deck=deck.replace(OLD_JOB,JOB)
    deck=re.sub(r'(\*Instance, name=Robot_SOLID-1, part=Robot_SOLID\n)[^*]+?(\*End Instance)',r'\1\2',deck, count=1)
    deck=re.sub(r'(\*Elset, elset=ROBOT_SOLID_CEL_ALL, instance=Robot_SOLID-1, generate\n)\d+,\s*\d+,\s*1',rf'\g<1>{int(element_labels.min())}, {int(element_labels.max())}, 1',deck,count=1)
    (HERE/f'{JOB}.inp').write_text(deck)
    (ROOT/f'{JOB}.inp').write_text(deck)
    pd.DataFrame([{'node':int(labels[i]),'x_mm':xyz[i,0],'y_mm':xyz[i,1],'z_mm':xyz[i,2]} for i in surface_idx]).to_csv(HERE/'freecad_robot_surface_nodes_global.csv',index=False)
    with (HERE/'freecad_robot_surface_triangles_exact.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['element','side','n1','n2','n3'])
        for ei,side,face in faces:w.writerow([int(element_labels[ei]),side]+[int(labels[i]) for i in face])
    print(json.dumps(identity,indent=2))


if __name__=='__main__': main()
