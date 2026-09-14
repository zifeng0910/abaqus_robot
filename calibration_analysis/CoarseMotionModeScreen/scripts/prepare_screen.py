"""Build all self-contained coarse-screen decks and per-case bridge sources."""

from pathlib import Path
import csv
import json
import re
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SCREEN = HERE.parent
REPO = SCREEN.parents[1]
ROOT = REPO.parent
CAD = REPO / "cad/freecad_parametric_robot/screening/coarse_motion_mode"
MESH = SCREEN / "cases/_meshes"
REFERENCE = REPO / "calibration_analysis/ReducedHydro_L2300_frozenCAD_validation"
TEMPLATE_JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L2300_D0815_WallSupported_0083"
TEMPLATE = REFERENCE / (TEMPLATE_JOB + ".inp")
BRIDGE = REPO / "calibration_analysis/Wobble30Hz_reduced_hydrodynamics/vuforc_socket_bridge_reduced_hydro.f"
RP = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
A0 = np.array([0.9647382600216, -0.1188742372140, 0.2348382536499])
DENSITY = 7.80906654321e-9
MOMENT_PER_VOLUME = 0.0010400401426980412 / 1.1402680223988795
sys.path.insert(0, str(HERE))
from case_manifest import CASES, ALIASES


def rotation_x_to_axis(axis):
    x = np.array([1., 0., 0.]); axis = axis / np.linalg.norm(axis)
    v = np.cross(x, axis); c = np.dot(x, axis)
    k = np.array([[0.,-v[2],v[1]],[v[2],0.,-v[0]],[-v[1],v[0],0.]])
    return np.eye(3) + k + k@k/(1+c)


def exterior(elements):
    definitions=((1,2,3),(0,3,2),(0,1,3),(0,2,1)); found={}
    for ei,tet in enumerate(elements):
        for side,ids in enumerate(definitions,1):
            face=tuple(int(tet[i]) for i in ids); key=tuple(sorted(face))
            found[key]=None if key in found else (ei,side,face)
    return [v for v in found.values() if v is not None]


def mesh_data(length):
    code="%04d" % round(length*1000); stem="Robot_SCREENING_L%s_D0815" % code
    nodes=pd.read_csv(MESH/(stem+"_nodes.csv")); elems=pd.read_csv(MESH/(stem+"_elements.csv"))
    props=json.loads((CAD/(stem+"_geometry.json")).read_text())
    labels=nodes.node.to_numpy(int); local=nodes[["x_mm","y_mm","z_mm"]].to_numpy(float)
    lookup={label:i for i,label in enumerate(labels)}
    element_labels=elems.element.to_numpy(int)
    elements=np.asarray([[lookup[int(v)] for v in row] for row in elems[["n1","n2","n3","n4"]].to_numpy()],int)
    rotation=rotation_x_to_axis(A0); cad_com=np.asarray(props["center_of_mass_mm"])
    xyz=RP+(local-cad_com)@rotation.T
    verts=xyz[elements]; volume=np.abs(np.linalg.det(verts[:,1:]-verts[:,:1])/6).sum()
    faces=exterior(elements)
    return labels,xyz,element_labels,elements,faces,props,volume


def mesh_block(labels,xyz,element_labels,elements,faces):
    side_sets={i:[] for i in range(1,5)}
    for ei,side,_ in faces: side_sets[side].append(int(element_labels[ei]))
    lines=["** SCREENING ONLY: exact STEP import; no node scaling", "*Part, name=Robot_SOLID", "*Node"]
    lines += ["%d, %.12g, %.12g, %.12g" % (label,p[0],p[1],p[2]) for label,p in zip(labels,xyz)]
    lines.append("*Element, type=C3D4")
    lines += ["%d, %s" % (label,", ".join(str(int(labels[i])) for i in tet)) for label,tet in zip(element_labels,elements)]
    for side in range(1,5):
        lines.append("*Elset, elset=ROBOT_EXTERIOR_S%d" % side)
        values=side_sets[side]
        lines += [", ".join(str(v) for v in values[i:i+16]) for i in range(0,len(values),16)]
    lines += ["*Nset, nset=ROBOT_CEL_BODY, generate", "%d, %d, 1"%(labels.min(),labels.max()),
              "*Elset, elset=ROBOT_CEL_BODY, generate", "%d, %d, 1"%(element_labels.min(),element_labels.max()),
              "*Elset, elset=ROBOT_SOLID_ALL, generate", "%d, %d, 1"%(element_labels.min(),element_labels.max()),
              "*Surface, type=ELEMENT, name=ROBOT_SOLID_SURF",
              "ROBOT_EXTERIOR_S1, S1","ROBOT_EXTERIOR_S2, S2","ROBOT_EXTERIOR_S3, S3","ROBOT_EXTERIOR_S4, S4",
              "*Solid Section, elset=ROBOT_CEL_BODY, material=MAT_ROBOT_RIGID",",","*End Part"]
    return "\n".join(lines)


def build_deck(item, mesh):
    labels,xyz,elabels,elements,faces,props,volume=mesh
    deck=TEMPLATE.read_text()
    start=deck.index("** ROBOT SOURCE:"); end=deck.index("*End Part",start)+len("*End Part")
    deck=deck[:start]+mesh_block(labels,xyz,elabels,elements,faces)+deck[end:]
    deck=deck.replace(TEMPLATE_JOB,item["job_name"])
    deck=re.sub(r"(\*Elset, elset=ROBOT_SOLID_CEL_ALL, instance=Robot_SOLID-1, generate\n)\d+,\s*\d+,\s*1",
                r"\g<1>%d, %d, 1"%(elabels.min(),elabels.max()),deck,count=1)
    deck=re.sub(r"(\*Friction\n)[^\n]+",r"\g<1>%.8g,"%item["mu"],deck,count=1)
    deck=re.sub(r"(\*Contact Damping, definition=CRITICAL DAMPING FRACTION, tangent fraction=0\.0\n)[^\n]+",
                r"\g<1>%.8g"%item["zeta"],deck,count=1)
    deck=re.sub(r"(\*Dynamic, Explicit, DIRECT\n)1\.0e-7,\s*[^\n]+",
                r"\g<1>1.0e-7, %.9g"%item["duration_s"],deck,count=1)
    deck=re.sub(r"\*Output, field, time interval=[^,\n]+, time marks=NO",
                "*Output, field, time interval=1.0e-4, time marks=NO",deck,count=1)
    header=("** SCREENING ONLY - NOT AN AUTHORITATIVE DESIGN\n"
            "** CASE %s; global 0.10 mm C3D4 rigid mesh; fixed direct dt=1e-7 s\n"%item["case_id"])
    return header+deck, props, volume, len(labels), len(elements), len(faces)


def main():
    (SCREEN/"cases").mkdir(parents=True,exist_ok=True); (SCREEN/"metrics").mkdir(exist_ok=True)
    meshes={length:mesh_data(length) for length in sorted(set(c["length_mm"] for c in CASES))}
    bridge_template=BRIDGE.read_text(); public=[]
    for item in CASES:
        folder=SCREEN/"cases"/item["case_id"]; (folder/"private").mkdir(parents=True,exist_ok=True)
        deck,props,mesh_volume,nodes,elements,faces=build_deck(item,meshes[item["length_mm"]])
        (folder/(item["job_name"]+".inp")).write_text(deck)
        bridge=bridge_template.replace("DATA cpar /4.0D-9/", "DATA cpar /%.8g/"%(4e-9*item["cparallel_scale"]))
        bridge=bridge.replace("kwob /3.0D-9/", "kwob /%.8g/"%(3e-9*item["kwobble_scale"]))
        bridge=bridge.replace("FILE='J:\\\\abaqusfangzhen\\\\reduced_hydro_runtime_load.csv'", "FILE='hydro_increment.csv'")
        (folder/"vuforc_screen.f").write_text(bridge)
        identity=dict(item)
        identity.update(mesh_size_mm=.100,node_count=nodes,element_count=elements,exterior_triangles=faces,
                        mesh_volume_mm3=mesh_volume,mesh_mass_mg=mesh_volume*DENSITY*1e9,
                        cad_volume_mm3=props["volume_mm3"],magnetic_moment_Am2=props["volume_mm3"]*MOMENT_PER_VOLUME,
                        field_interval_s=1e-4,direct_dt_s=1e-7,status="PREPARED")
        (folder/"case_identity.json").write_text(json.dumps(identity,indent=2)+"\n")
        public.append(identity)
    pd.DataFrame(public).to_csv(SCREEN/"metrics/case_manifest.csv",index=False)
    with (SCREEN/"metrics/case_aliases.csv").open("w",newline="") as handle:
        writer=csv.writer(handle); writer.writerow(("requested_case_id","executed_as")); writer.writerows(ALIASES.items())
    print("PREPARED %d unique cases"%len(public))


if __name__=="__main__": main()

