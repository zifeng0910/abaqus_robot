"""Prepare Case B or C using the independently validated global V + spatial VR driver."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TARGET="TRUECEL_A14P5_DIAG_2CYCLES"
COARSE_SOURCE="TRUECEL_A14P5_FORCE_FLUX_DIAG_2CYCLES"
INVALID_REFINED="TRUECEL_A14P5_COARSETRAJ_CEL_HREFINE_REPLAY"
CASE_B="TRUECEL_A14P5_COARSETRAJ_CEL_COARSE_REPLAY_VVR"
CASE_C="TRUECEL_A14P5_COARSETRAJ_CEL_HREFINE_REPLAY_VVR"
STOP=.0108

def digest_bytes(data):return hashlib.sha256(data).hexdigest().upper()
def digest(path):return digest_bytes(path.read_bytes())

def replace_once(text,old,new):
    if text.count(old)!=1:raise RuntimeError(f"Expected one occurrence, got {text.count(old)}: {old[:80]}")
    return text.replace(old,new,1)

def driver_block():
    z=np.load(ROOT/"case"/TARGET/"private"/"rp_history_private.npz")
    t=z["V1"][:,0];use=t<=STOP;t=t[use]
    out=["** VALIDATED GLOBAL V + SPATIAL/GLOBAL VR REPLAY DRIVER"]
    for stem,prefix,offset in (("V","REPLAY_V",1),("VR","REPLAY_VR",4)):
        for axis in range(3):
            values=z[f"{stem}{axis+1}"][:,1][use]
            out.append(f"*Amplitude, name={prefix}{axis+1}, definition=TABULAR, time=STEP TIME")
            out.extend(f"{x:.16g}, {y:.16g}" for x,y in zip(t,values))
            out.append(f"*Boundary, amplitude={prefix}{axis+1}, type=VELOCITY")
            out.append(f"RP_ROBOT, {offset+axis}, {offset+axis}, 1.")
    return "\n".join(out)+"\n"

def robot_part(deck):return re.search(r"(?ms)^\*Part, name=Robot_SOLID\s*$.*?^\*End Part\s*$",deck).group(0)
def wall_parts(deck):
    return [m.group(0) for m in re.finditer(r"(?ms)^\*Part, name=(?:Pipe_SOLID|Pipe_WALL_HELPER|Pipe_TERMINAL_STOP|Pipe_TERMINAL_STOP_HIGH)\s*$.*?^\*End Part\s*$",deck)]

def main():
    ap=argparse.ArgumentParser();ap.add_argument("case",choices=("B","C"));args=ap.parse_args()
    validation=json.loads((ROOT/"A14P5_PRESCRIBED_MOTION_REPLAY_DEBUG.json").read_text())
    if validation["replay_classification"]!="KINEMATIC_REPLAY_VALIDATED":raise RuntimeError("V/VR driver not validated")
    job=CASE_B if args.case=="B" else CASE_C
    case=ROOT/"case"/job
    if case.exists():raise RuntimeError(f"Case exists; refusing overwrite: {case}")
    if args.case=="C":
        gate=json.loads((ROOT/(CASE_B+"_Kinematic_Gate.json")).read_text())
        if not gate["passed"]:raise RuntimeError("Case B kinematic gate did not pass")
    source_job=COARSE_SOURCE if args.case=="B" else INVALID_REFINED
    source_case=ROOT/"case"/source_job
    deck=(source_case/(source_job+".inp")).read_text(encoding="latin1")
    deck=deck.replace(source_job,job)
    if args.case=="B":
        deck=replace_once(deck,", 0.016666666667\n*Bulk Viscosity",", 0.0108\n*Bulk Viscosity")
        marker="** OUTPUT REQUESTS\n"
        deck=replace_once(deck,marker,driver_block()+marker)
    else:
        start="** PRESCRIBED ORIGINAL COARSE RP TRAJECTORY: U AND TOTAL ROTATION VECTOR UR\n"
        end="** OUTPUT REQUESTS\n"
        if deck.count(start)!=1 or deck.count(end)!=1:raise RuntimeError("Invalid refined driver markers")
        deck=deck.split(start,1)[0]+driver_block()+end+deck.split(end,1)[1]
    if "REPLAY_R1" in deck or "REPLAY_U1" in deck:raise RuntimeError("Legacy U/UR driver remains")
    for required in ("REPLAY_V1","REPLAY_VR1","SECOND SURFACE=Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF","*Contact Inclusions"):
        if required not in deck:raise RuntimeError(f"Missing {required}")
    case.mkdir(parents=True);inp=case/(job+".inp");inp.write_text(deck,encoding="latin1")
    for name in ("magnetic_field_gradient_table_A14P5.dat","magnetic_field_gradient_table_A14P5.json"):
        shutil.copy2(source_case/name,case/name)
    vuamp=(source_case/"vuamp_precomputed_truecel.f90").read_text(encoding="ascii")
    vuamp=re.sub(
        r"open\(unit=77,file='[^']*magnetic_increment\.csv'",
        lambda _: f"open(unit=77,file='{case}\\magnetic_increment.csv'",
        vuamp,
        count=1,
    )
    (case/"vuamp_precomputed_truecel.f90").write_text(vuamp,encoding="ascii")
    source_identity=json.loads((source_case/"case_identity.json").read_text())
    elements=17600 if args.case=="B" else 86640
    dims=[44,20,20] if args.case=="B" else [57,40,38]
    identity=dict(source_identity);identity.update(case_id=job,status="PREPARED",dynamics_run_count=0,
        diagnostic_type="VALIDATED_GLOBAL_V_SPATIAL_VR_PRESCRIBED_CEL_REPLAY",duration_s=STOP,
        Eulerian_element_count=elements,Eulerian_dimensions=dims,input_sha256=digest(inp),
        driver_sha256=digest_bytes(driver_block().encode("latin1")),validated_driver_commit="a5d466b")
    (case/"case_identity.json").write_text(json.dumps(identity,indent=2)+"\n",encoding="ascii")
    target_identity=json.loads((ROOT/"case"/TARGET/"case_identity.json").read_text())
    initial={"target_U":[0.,0.,0.],"target_UR":[0.,0.,0.],"replay_U":[0.,0.,0.],"replay_UR":[0.,0.,0.],"position_error_mm":0.,"orientation_error_deg":0.}
    gate={"job":job,"case_label":args.case,"CEL_elements":elements,"CEL_dimensions":dims,
          "driver_sha256":identity["driver_sha256"],"driver_validated":True,"initial_pose_gate":initial,
          "robot_part_byte_identical_to_coarse_original":robot_part(deck)==robot_part((ROOT/"case"/COARSE_SOURCE/(COARSE_SOURCE+".inp")).read_text(encoding="latin1")),
          "wall_parts_byte_identical_to_coarse_original":wall_parts(deck)==wall_parts((ROOT/"case"/COARSE_SOURCE/(COARSE_SOURCE+".inp")).read_text(encoding="latin1")),
          "magnetic_table_sha256":digest(case/"magnetic_field_gradient_table_A14P5.dat"),
          "physics":{"B0_mT":target_identity["B0_mT"],"G_mT":target_identity["gradient_mT"],"frequency_Hz":target_identity["frequency_Hz"],"mu":target_identity["mu"],"zeta":target_identity["zeta"]},
          "input_sha256":digest(inp)}
    if not gate["robot_part_byte_identical_to_coarse_original"] or not gate["wall_parts_byte_identical_to_coarse_original"]:raise RuntimeError("Robot/wall identity gate failed")
    (ROOT/(job+"_PreRun_Gate.json")).write_text(json.dumps(gate,indent=2)+"\n",encoding="ascii")
    print(json.dumps(gate,indent=2))

if __name__=="__main__":main()
