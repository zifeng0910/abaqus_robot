"""Prepare the single B0=11 mT, A_main=14.5 deg coarse TRUE-CEL full-transit screen."""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
PARENT_JOB = "TRUECEL_NO_RECOIL_B0P11_2CYCLES"
JOB = "TRUECEL_B0P11_A14P5_FAST_FULL_TRANSIT"
TABLE_NAME = "magnetic_field_gradient_table_B0P11_A14P5.dat"
MAX_DURATION_S = 2.0
FLUID_FIELD_INTERVAL_S = 1.0 / 120.0

BACKEND = REPO / "calibration_analysis" / "PrecomputedStraightMagneticBackend"
sys.path.insert(0, str(BACKEND / "scripts"))
import table_model


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def part_nodes(deck: str, name: str) -> np.ndarray:
    part = re.search(
        rf"^\*Part, name={re.escape(name)}\s*$([\s\S]*?)^\*End Part\s*$",
        deck,
        re.MULTILINE | re.IGNORECASE,
    )
    if not part:
        raise RuntimeError(f"Part not found: {name}")
    node_block = re.search(
        r"^\*Node\s*$([\s\S]*?)(?=^\*)",
        part.group(1),
        re.MULTILINE | re.IGNORECASE,
    )
    if not node_block:
        raise RuntimeError(f"Node block not found: {name}")
    rows = []
    for line in node_block.group(1).splitlines():
        fields = [item.strip() for item in line.split(",")]
        if len(fields) >= 4:
            rows.append([float(value) for value in fields[1:4]])
    return np.asarray(rows, dtype=float)


def replace_once(text: str, old: str, new: str, label: str) -> str:
    if text.count(old) != 1:
        raise RuntimeError(f"Expected exactly one {label}; found {text.count(old)}")
    return text.replace(old, new, 1)


def add_online_monitor(sub: str, target_displacement_mm: float) -> str:
    declarations = """  logical :: terminate_requested=.false., event_open=.false., monitor_open=.false.\n  character(len=64) :: terminate_reason='NONE'\n  real(8) :: neg_start=-1.0d0, loss_start=-1.0d0, collapse_start=-1.0d0\n  real(8) :: max_disp=0.0d0, max_dt_seen=0.0d0, next_monitor_sample=0.0d0\n  integer, parameter :: monitor_hist_n=512\n  real(8) :: monitor_t(monitor_hist_n)=0.0d0, monitor_s(monitor_hist_n)=0.0d0\n  integer :: monitor_hist_count=0, monitor_hist_pos=0\n"""
    sub = replace_once(
        sub,
        "  logical :: table_loaded=.false., log_open=.false.\n",
        "  logical :: table_loaded=.false., log_open=.false.\n" + declarations,
        "module monitor declarations",
    )
    helpers = f"""
  subroutine request_termination(reason,t,s_abs,disp,vs,dt)
    character(len=*),intent(in) :: reason
    real(8),intent(in) :: t,s_abs,disp,vs,dt
    if (terminate_requested) return
    terminate_requested=.true.
    terminate_reason=reason
    open(unit=78,file=trim('J:\\abaqusfangzhen\\abaqus_robot\\calibration_analysis\\' // &
      'TrueCELLongForwardTransit\\case\\{JOB}\\fast_screen_event.txt'),status='replace')
    write(78,'(A)') trim(reason)
    write(78,'(A,ES18.10)') 'time_s=',t
    write(78,'(A,ES18.10)') 's_abs_mm=',s_abs
    write(78,'(A,ES18.10)') 'displacement_mm=',disp
    write(78,'(A,ES18.10)') 'v_s_mm_s=',vs
    write(78,'(A,ES18.10)') 'stable_dt_s=',dt
    close(78)
    write(6,*) 'FAST_SCREEN_TERMINATION ',trim(reason),' time=',t,' s=',s_abs,' vs=',vs
  end subroutine request_termination

  subroutine monitor_motion(t,dt,disp,vs)
    real(8),intent(in) :: t,dt,disp,vs
    real(8) :: s_abs,old_t,old_s,best_t
    integer :: i,old_idx
    s_abs=-3.2d0+disp
    max_disp=max(max_disp,disp)
    max_dt_seen=max(max_dt_seen,dt)
    if (vs < -2.0d0) then
      if (neg_start < 0.0d0) neg_start=t
      if (t-neg_start >= 2.5d-4) then
        call request_termination('SUSTAINED_RECOIL_FAIL_VELOCITY',t,s_abs,disp,vs,dt)
      endif
    else
      neg_start=-1.0d0
    endif
    if (max_disp-disp > 2.0d-2) then
      if (loss_start < 0.0d0) loss_start=t
      if (t-loss_start >= 5.0d-4) then
        call request_termination('SUSTAINED_RECOIL_FAIL_POSITION',t,s_abs,disp,vs,dt)
      endif
    else
      loss_start=-1.0d0
    endif
    if (t > 1.0d-3 .and. dt < 0.25d0*max_dt_seen) then
      if (collapse_start < 0.0d0) collapse_start=t
      if (t-collapse_start >= 2.5d-4) then
        call request_termination('NUMERICALLY_INVALID_STABLE_DT_COLLAPSE',t,s_abs,disp,vs,dt)
      endif
    else
      collapse_start=-1.0d0
    endif
    if (disp >= {target_displacement_mm:.16g}d0) then
      call request_termination('FULL_TRANSIT_FORWARD_PASS',t,s_abs,disp,vs,dt)
    endif
    if (t+1.0d-15 >= next_monitor_sample) then
      monitor_hist_pos=mod(monitor_hist_pos,monitor_hist_n)+1
      monitor_t(monitor_hist_pos)=t
      monitor_s(monitor_hist_pos)=disp
      monitor_hist_count=min(monitor_hist_count+1,monitor_hist_n)
      old_idx=0; best_t=-1.0d99
      do i=1,monitor_hist_count
        old_t=monitor_t(i)
        if (old_t <= t-2.0d-3 .and. old_t > best_t) then
          best_t=old_t; old_idx=i
        endif
      enddo
      if (t >= 1.0d0/120.0d0 .and. old_idx > 0) then
        old_s=monitor_s(old_idx)
        if (disp-old_s < 5.0d-3 .and. disp < {target_displacement_mm:.16g}d0) then
          call request_termination('FORWARD_PROGRESS_STALL',t,s_abs,disp,vs,dt)
        endif
      endif
      if (.not.monitor_open) then
        open(unit=79,file=trim('J:\\abaqusfangzhen\\abaqus_robot\\calibration_analysis\\' // &
          'TrueCELLongForwardTransit\\case\\{JOB}\\fast_screen_online.csv'),status='replace')
        write(79,'(A)') 'time_s,s_abs_mm,displacement_mm,v_s_mm_s,max_displacement_mm,stable_dt_s,termination_reason'
        monitor_open=.true.
      endif
      write(79,'(6(ES18.10,:,","),A)') t,s_abs,disp,vs,max_disp,dt,trim(terminate_reason)
      next_monitor_sample=t+1.0d-5
    endif
  end subroutine monitor_motion
"""
    sub = replace_once(
        sub,
        "  end subroutine open_increment_log\n",
        "  end subroutine open_increment_log\n" + helpers,
        "monitor helper insertion point",
    )
    sub = replace_once(
        sub,
        "  real(8) :: loads(6),last_time,t,u(3),ur(3)\n",
        "  real(8) :: loads(6),last_time,t,u(3),ur(3),v(3),vs\n",
        "VUAMP variables",
    )
    sub = replace_once(
        sub,
        "    ur(3)=VGETSENSORVALUE('RP_UR3',jSensorLookUpTable,sensorValues)\n",
        "    ur(3)=VGETSENSORVALUE('RP_UR3',jSensorLookUpTable,sensorValues)\n"
        "    v(1)=VGETSENSORVALUE('RP_V1',jSensorLookUpTable,sensorValues)\n"
        "    v(2)=VGETSENSORVALUE('RP_V2',jSensorLookUpTable,sensorValues)\n"
        "    v(3)=VGETSENSORVALUE('RP_V3',jSensorLookUpTable,sensorValues)\n",
        "RP velocity sensor reads",
    )
    sub = replace_once(
        sub,
        "    dotu=dot_product(u,c)\n",
        "    dotu=dot_product(u,c)\n    vs=dot_product(v,c)\n    call monitor_motion(t,dble(dt),dotu,vs)\n",
        "online monitor call",
    )
    sub = replace_once(
        sub,
        "  AmpValueNew=0.0d0\n",
        "  if (terminate_requested) lFlagsDefine(iConcludeStep)=1\n  AmpValueNew=0.0d0\n",
        "graceful conclude flag",
    )
    sub = replace_once(
        sub,
        "      flush(77); close(77); log_open=.false.\n",
        "      flush(77); close(77); log_open=.false.\n"
        "    endif\n"
        "    if (monitor_open) then\n"
        "      flush(79); close(79); monitor_open=.false.\n",
        "monitor log close",
    )
    return sub


def main() -> None:
    parent = OUT / "case" / PARENT_JOB
    case = OUT / "case" / JOB
    if case.exists():
        raise RuntimeError("FAST candidate directory already exists; refusing another preparation")

    identity = json.loads((parent / "case_identity.json").read_text())
    expected = {
        "B0_mT": 11.0,
        "gradient_mT": 2.0,
        "frequency_Hz": 120.0,
        "rocking_main_amplitude_deg": 14.8,
        "rocking_cross_amplitude_deg": 2.5,
        "fluid_EOS_c0_mm_s": 100000.0,
        "Eulerian_dimensions": [44, 20, 20],
        "explicit_stable_time_scale_factor": 0.4,
    }
    for key, value in expected.items():
        if identity.get(key) != value:
            raise RuntimeError(f"Authoritative parent mismatch: {key}={identity.get(key)!r}, expected {value!r}")

    parent_inp = parent / f"{PARENT_JOB}.inp"
    parent_deck = parent_inp.read_text(encoding="latin1")
    c = np.asarray(identity["canonical_plus_s_axis_aba"], dtype=float)
    rp0 = np.asarray(identity["initial_center_aba_mm"], dtype=float)
    robot = part_nodes(parent_deck, "Robot_SOLID")
    head_extent = float(np.max((robot - rp0).dot(c)))
    s_start = float(identity["initial_center_s_relative_to_window_mm"])
    s_domain_finish = float(identity["axial_CEL_extent_s_mm"][1])
    end_clearance = float(identity["Eulerian_spacing_mm"][0])
    s_finish = s_domain_finish - head_extent - end_clearance
    travel = s_finish - s_start
    if not (7.0 < travel < 8.0):
        raise RuntimeError(f"Unexpected full-transit travel derived from geometry: {travel:.9g} mm")

    thresholds = json.loads((REPO / "calibration_analysis" / "RefinedDualEndTrueCELGate" / "selected_refined_geometry.json").read_text())
    selected_geometry = thresholds["best_near_candidate_not_selected"]
    head_touch = float(selected_geometry["HEAD_touch_deg"])
    tail_touch = float(selected_geometry["TAIL_touch_deg"])
    if not (14.5 > abs(head_touch) and 14.5 > abs(tail_touch)):
        raise RuntimeError("A_main=14.5 deg is not above both authoritative geometric touch thresholds")

    cfg = table_model.config()
    verified_case = OUT / "case" / "TRUECEL_A14P5_STEADY_2CYCLES"
    verified_table = verified_case / "magnetic_field_gradient_table_A14P5.dat"
    verified_meta_path = verified_case / "magnetic_field_gradient_table_A14P5.json"
    verified_meta = json.loads(verified_meta_path.read_text())
    if verified_meta["rocking_main_amplitude_deg"] != 14.5:
        raise RuntimeError("Existing verified table is not A_main=14.5 deg")
    if verified_meta["table_sha256"] != sha(verified_table):
        raise RuntimeError("Existing verified A14.5 table hash mismatch")
    if verified_meta["authoritative_regeneration_max_abs_error"] > 1e-12:
        raise RuntimeError("Existing A14.5 table failed its authoritative regeneration gate")
    case.mkdir(parents=True)
    table_path = case / TABLE_NAME
    shutil.copy2(verified_table, table_path)
    table_meta = dict(verified_meta)
    table_meta.update({
        "reused_verified_source": str(verified_table),
        "only_physics_change": "rocking_main_amplitude_deg 14.8 -> 14.5",
    })
    (case / f"{TABLE_NAME}.json").write_text(json.dumps(table_meta, indent=2) + "\n", encoding="ascii")

    deck = parent_deck.replace(PARENT_JOB, JOB)
    deck, n_duration = re.subn(
        r"(\*Dynamic, Explicit, SCALE FACTOR=0\.4\s*\r?\n)\s*,\s*0\.016666666667",
        rf"\g<1>, {MAX_DURATION_S:.12f}",
        deck,
        count=1,
    )
    if n_duration != 1:
        raise RuntimeError("Could not replace parent step duration")
    old_field = re.search(
        r"\*Output, field, time interval=0\.00025, time marks=NO[\s\S]*?(?=\*Output, history, frequency=1)",
        deck,
    )
    if not old_field:
        raise RuntimeError("Could not locate parent field-output block")
    sparse_field = f"""*Output, field, time interval={FLUID_FIELD_INTERVAL_S:.12g}, time marks=NO
*Node Output, nset=Robot_SOLID-1.ROBOT_CEL_BODY
U, V
*Node Output, nset=RP_ROBOT
U, UR, V, VR
*Contact Output, surface=ROBOT_SOLID-1.ROBOT_SOLID_SURF
CDISP
*Node Output, nset=Fluid_EULERIAN-1.FLUID_CEL_NODES
V
*Element Output, elset=FLUID_CEL_ALL
EVF, PRESS
"""
    deck = deck[: old_field.start()] + sparse_field + deck[old_field.end() :]
    deck = (
        "** REDUCED_SOUND_SPEED_COARSE_CEL_FAST_SCREEN\n"
        "** ONE PHYSICS CHANGE FROM AUTHORITATIVE PARENT: A_main 14.8 -> 14.5 deg\n"
        f"** FULL TRANSIT s_start={s_start:.12f} s_finish={s_finish:.12f} travel={travel:.12f} mm\n"
        + deck
    )
    inp = case / f"{JOB}.inp"
    inp.write_text(deck, encoding="latin1")

    parent_sub = (parent / "vuamp_precomputed_truecel.f90").read_text(encoding="ascii")
    old_table = str(BACKEND / cfg["table"]["path"])
    table_assignment = (
        "trim('J:\\abaqusfangzhen\\abaqus_robot\\calibration_analysis\\' // &\n"
        f"      'TrueCELLongForwardTransit\\case\\{JOB}\\' // &\n"
        f"      '{TABLE_NAME}')"
    )
    log_assignment = (
        "trim('J:\\abaqusfangzhen\\abaqus_robot\\calibration_analysis\\' // &\n"
        f"      'TrueCELLongForwardTransit\\case\\{JOB}\\' // &\n"
        "      'magnetic_increment.csv')"
    )
    sub = replace_once(parent_sub, "'" + old_table + "'", table_assignment, "magnetic table path")
    sub = replace_once(sub, "'" + str(parent / "magnetic_increment.csv") + "'", log_assignment, "increment log path")
    sub = add_online_monitor(sub, travel)
    sub_path = case / "vuamp_precomputed_truecel.f90"
    sub_path.write_text(sub, encoding="ascii")

    new_identity = dict(identity)
    for key in ("wallclock_s", "cpus", "socket_calls", "failure"):
        new_identity.pop(key, None)
    new_identity.update({
        "case_id": JOB,
        "status": "PREPARED",
        "classification": "REDUCED_SOUND_SPEED_COARSE_CEL_FAST_SCREEN",
        "physical_parent": PARENT_JOB,
        "rocking_main_amplitude_deg": 14.5,
        "single_physics_change": "A_main 14.8 -> 14.5 deg",
        "duration_s": MAX_DURATION_S,
        "dynamics_run_count": 0,
        "magnetic_table_file": TABLE_NAME,
        "magnetic_table_sha256": sha(table_path),
        "input_sha256": sha(inp),
        "fortran_sha256": sha(sub_path),
        "s_start_mm": s_start,
        "s_finish_mm": s_finish,
        "travel_distance_mm": travel,
        "head_extent_from_RP_mm": head_extent,
        "end_clearance_mm": end_clearance,
        "end_clearance_basis": "one actual coarse axial CEL element",
        "geometric_touch_thresholds_deg": {"HEAD": head_touch, "TAIL": tail_touch},
        "online_fail_fast": {
            "sustained_recoil_velocity": "v_s < -2 mm/s continuously for >=0.25 ms",
            "unrecovered_position_loss": ">0.02 mm for >=0.5 ms",
            "stall": "rolling 2 ms gain <0.005 mm after first 120-Hz cycle",
            "stable_dt_collapse": "dt <25% of prior max continuously for >=0.25 ms",
            "termination": "VUAMP iConcludeStep; partial ODB preserved",
        },
        "field_output": {
            "combined_sparse_interval_s": FLUID_FIELD_INTERVAL_S,
            "motion_gif_source": "high-frequency RP history; field frames are not required",
        },
        "frozen_for_candidate": [
            "B0=11mT", "G=2mT", "f=120Hz", "A_cross=2.5deg",
            "TRUE-CEL 44x20x20 mesh/domain", "c0=100000mm/s",
            "water density/viscosity/EOS", "robot/wall geometry",
            "one-wall General Contact; mu=0.03; zeta=1.0",
            "Explicit scale factor=0.4; no mass scaling", "initial pose and magnetic ramp",
        ],
    })
    (case / "case_identity.json").write_text(json.dumps(new_identity, indent=2) + "\n", encoding="ascii")
    setup = {
        "case": JOB,
        "classification": new_identity["classification"],
        "parent": PARENT_JOB,
        "single_physics_change": new_identity["single_physics_change"],
        "B0_mT": 11.0,
        "G_mT": 2.0,
        "frequency_Hz": 120.0,
        "A_main_deg": 14.5,
        "A_cross_deg": 2.5,
        "c0_mm_s": 100000.0,
        "CEL_mesh": [44, 20, 20],
        "contact": "parent unchanged",
        "explicit_controls": "SCALE FACTOR=0.4; no mass scaling",
        "s_start_mm": s_start,
        "s_finish_mm": s_finish,
        "travel_distance_mm": travel,
        "max_duration_s": MAX_DURATION_S,
        "head_touch_deg": head_touch,
        "tail_touch_deg": tail_touch,
        "input_sha256": sha(inp),
        "fortran_sha256": sha(sub_path),
        "magnetic_table_sha256": sha(table_path),
        "dynamics_candidates_authorized": 1,
    }
    (OUT / f"{JOB}_Setup_Audit.json").write_text(json.dumps(setup, indent=2) + "\n", encoding="ascii")
    print(json.dumps(setup, indent=2))


if __name__ == "__main__":
    main()
