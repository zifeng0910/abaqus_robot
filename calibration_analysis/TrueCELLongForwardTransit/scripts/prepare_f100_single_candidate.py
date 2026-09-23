"""Prepare the one authorized fresh-start 100-Hz TRUE-CEL candidate."""
from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
CASE_ROOT = OUT / "case"
SOURCE_JOB = "TRUECEL_B0P11_A14P5_FAST_FULL_TRANSIT"
JOB = "TRUECEL_B0P11_A14P5_F100_FAST"
SOURCE = CASE_ROOT / SOURCE_JOB
DEST = CASE_ROOT / JOB
PERIOD = 0.01
STAGE1_END = 2 * PERIOD


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def replace_once(text: str, old: str, new: str, name: str) -> str:
    if text.count(old) != 1:
        raise RuntimeError(f"Expected one {name}; found {text.count(old)}")
    return text.replace(old, new, 1)


MONITOR = r"""
  subroutine monitor_motion(t,dt,disp,vs)
    real(8),intent(in) :: t,dt,disp,vs
    real(8) :: frac,bd,bv,delta,back
    integer :: unit
    character(len=48) :: reason
    reason='NONE'
    monitor_max_dt=max(monitor_max_dt,dt)
    monitor_max_s=max(monitor_max_s,disp)
    back=monitor_max_s-disp
    if (back > 2.0d-2) then
      if (monitor_loss_start < 0.0d0) monitor_loss_start=t
    else if (back <= 1.0d-2) then
      monitor_loss_start=-1.0d0
    endif
    if (monitor_loss_start >= 0.0d0 .and. t-monitor_loss_start >= 5.0d-3) &
      reason='F100_RECOIL_FAIL_POSITION'
    if (t > 1.0d-3 .and. dt < 0.25d0*monitor_max_dt) then
      if (monitor_collapse_start < 0.0d0) monitor_collapse_start=t
      if (t-monitor_collapse_start >= 2.5d-4) reason='NUMERICALLY_INVALID_STABLE_DT_COLLAPSE'
    else
      monitor_collapse_start=-1.0d0
    endif
    if (t >= monitor_next_boundary .and. monitor_last_t < monitor_next_boundary) then
      frac=(monitor_next_boundary-monitor_last_t)/max(t-monitor_last_t,1.0d-30)
      bd=monitor_last_s+frac*(disp-monitor_last_s)
      bv=monitor_last_v+frac*(vs-monitor_last_v)
      delta=bd-monitor_cycle_start_s
      if (delta <= -1.0d-2) reason='F100_RECOIL_FAIL_CYCLE'
      open(unit=80,file=trim('J:\abaqusfangzhen\abaqus_robot\calibration_analysis\' // &
        'TrueCELLongForwardTransit\case\TRUECEL_B0P11_A14P5_F100_FAST\online_cycle_gate.csv'), &
        status='unknown',position='append')
      if (monitor_cycle==1) write(80,'(A)') 'cycle,end_time_s,delta_s_mm,end_v_s_mm_s'
      write(80,'(I0,A,ES18.10,A,ES18.10,A,ES18.10)') monitor_cycle,',',monitor_next_boundary,',',delta,',',bv
      close(80)
      monitor_cycle=monitor_cycle+1
      monitor_cycle_start_t=monitor_next_boundary
      monitor_cycle_start_s=bd
      monitor_next_boundary=monitor_next_boundary+1.0d-2
    endif
    if (t+1.0d-15 >= monitor_next_sample) then
      if (.not.monitor_open) then
        open(unit=79,file=trim('J:\abaqusfangzhen\abaqus_robot\calibration_analysis\' // &
          'TrueCELLongForwardTransit\case\TRUECEL_B0P11_A14P5_F100_FAST\online_motion.csv'),status='replace')
        write(79,'(A)') 'time_s,displacement_mm,v_s_mm_s,max_backtrack_mm,stable_dt_s'
        monitor_open=.true.
      endif
      write(79,'(4(ES18.10,:,","),ES18.10)') t,disp,vs,back,dt
      flush(79)
      monitor_next_sample=t+1.0d-5
    endif
    if (reason/='NONE' .and. .not.terminate_requested) then
      terminate_requested=.true.; terminate_reason=reason
      open(unit=78,file=trim('J:\abaqusfangzhen\abaqus_robot\calibration_analysis\' // &
        'TrueCELLongForwardTransit\case\TRUECEL_B0P11_A14P5_F100_FAST\f100_event.txt'),status='replace')
      write(78,'(A)') trim(reason)
      write(78,'(A,ES18.10)') 'time_s=',t
      write(78,'(A,ES18.10)') 'displacement_mm=',disp
      write(78,'(A,ES18.10)') 'v_s_mm_s=',vs
      write(78,'(A,ES18.10)') 'max_backtrack_mm=',back
      close(78)
    endif
    monitor_last_t=t; monitor_last_s=disp; monitor_last_v=vs
  end subroutine monitor_motion
"""


def main() -> None:
    if DEST.exists():
        raise RuntimeError(f"Candidate already exists; refusing to overwrite: {DEST}")
    source_id = json.loads((SOURCE / "case_identity.json").read_text(encoding="utf-8"))
    required = {
        "B0_mT": 11.0, "gradient_mT": 2.0,
        "frequency_Hz": 120.0, "rocking_main_amplitude_deg": 14.5,
        "rocking_cross_amplitude_deg": 2.5, "fluid_EOS_c0_mm_s": 100000.0,
        "Eulerian_dimensions": [44, 20, 20],
        "explicit_stable_time_scale_factor": 0.4,
    }
    for key, value in required.items():
        if source_id.get(key) != value:
            raise RuntimeError(f"Authoritative F120 source mismatch for {key}: {source_id.get(key)!r}")
    src_deck_path = SOURCE / f"{SOURCE_JOB}.inp"
    deck = src_deck_path.read_text(encoding="latin1")
    deck = deck.replace(SOURCE_JOB, JOB)
    deck, count = re.subn(
        r"(\*Dynamic, Explicit, SCALE FACTOR=0\.4\s*\r?\n)\s*,\s*[0-9.eE+-]+",
        rf"\g<1>, {STAGE1_END:.12f}", deck, count=1, flags=re.IGNORECASE,
    )
    if count != 1:
        raise RuntimeError("Could not set the fresh-start Stage-1 duration to exactly 20 ms")
    deck = "** ONE AUTHORIZED PHYSICS CHANGE FROM F120: frequency 120 -> 100 Hz\n" + deck

    source_f90 = (SOURCE / "vuamp_precomputed_truecel.f90").read_text(encoding="ascii")
    source_f90 = source_f90.replace(SOURCE_JOB, JOB)
    source_f90 = replace_once(source_f90, "43200.0d0*t", "36000.0d0*t", "100-Hz phase law")
    # Keep the force/torque lookup, interpolation, ramp, and scaling byte-for-byte unchanged.
    decl = """  logical :: monitor_initialized=.false.
  real(8) :: monitor_last_t=-1.0d0,monitor_last_s=0.0d0,monitor_last_v=0.0d0
  real(8) :: monitor_max_s=0.0d0,monitor_cycle_start_t=0.0d0,monitor_cycle_start_s=0.0d0
  real(8) :: monitor_next_boundary=1.0d-2,monitor_next_sample=0.0d0
  real(8) :: monitor_loss_start=-1.0d0,monitor_collapse_start=-1.0d0,monitor_max_dt=0.0d0
  integer :: monitor_cycle=1
"""
    source_f90 = replace_once(source_f90,
        "  logical :: terminate_requested=.false., event_open=.false., monitor_open=.false.\n",
        "  logical :: terminate_requested=.false., event_open=.false., monitor_open=.false.\n" + decl,
        "monitor state declarations")
    monitor_pat = re.compile(r"  subroutine monitor_motion\(t,dt,disp,vs\)[\s\S]*?  end subroutine monitor_motion\n")
    source_f90, count = monitor_pat.subn(lambda _m: MONITOR, source_f90, count=1)
    if count != 1:
        raise RuntimeError("Could not replace exactly one outdated 120-Hz online monitor")
    source_f90 = replace_once(source_f90,
        "  real(8) :: loads(6),last_time,t,u(3),ur(3),v(3),vs\n",
        "  real(8) :: loads(6),last_time,t,u(3),ur(3),v(3),vs\n",
        "RP sensor declaration")
    # Initialize Stage-1 gate state on the first master VUAMP call.
    source_f90 = replace_once(source_f90,
        "  if (is_master==1 .and. abs(t-last_time)>1.0d-15) then\n",
        "  if (is_master==1 .and. abs(t-last_time)>1.0d-15) then\n"
        "    if (.not.monitor_initialized) then\n"
        "      monitor_initialized=.true.; monitor_last_t=t; monitor_max_s=0.0d0\n"
        "      monitor_cycle_start_s=0.0d0; monitor_next_boundary=1.0d-2\n"
        "    endif\n",
        "monitor initialization")
    source_f90 = replace_once(source_f90,
        "    call monitor_motion(t,dble(dt),dotu,vs)\n",
        "    call monitor_motion(t,dble(dt),dotu,vs)\n",
        "monitor call")
    source_f90 = source_f90.replace("magnetic_increment.csv", "magnetic_increment_f100.csv")
    if "SUSTAINED_RECOIL_FAIL_VELOCITY" in source_f90 or "1.0d0/120.0d0" in source_f90:
        raise RuntimeError("An outdated velocity fail-fast or 120-Hz monitor remains")

    DEST.mkdir(parents=True)
    shutil.copy2(SOURCE / "magnetic_field_gradient_table_B0P11_A14P5.dat", DEST)
    meta_src = SOURCE / "magnetic_field_gradient_table_B0P11_A14P5.dat.json"
    shutil.copy2(meta_src, DEST / meta_src.name)
    inp = DEST / f"{JOB}.inp"
    inp.write_text(deck, encoding="latin1")
    f90 = DEST / "vuamp_precomputed_truecel.f90"
    f90.write_text(source_f90, encoding="ascii")
    identity = dict(source_id)
    identity.update({
        "case_id": JOB, "status": "PREPARED", "classification": "F100_TWO_CYCLE_GATE",
        "physical_parent": SOURCE_JOB, "frequency_Hz": 100.0,
        "duration_s": STAGE1_END, "dynamics_run_count": 0,
        "single_physics_change": "frequency 120 -> 100 Hz",
        "frequency_change_only": True,
        "magnetic_table_file": "magnetic_field_gradient_table_B0P11_A14P5.dat",
        "magnetic_table_sha256": sha(DEST / "magnetic_field_gradient_table_B0P11_A14P5.dat"),
        "input_sha256": sha(inp), "fortran_sha256": sha(f90),
        "stage1_duration_s": STAGE1_END, "period_s": PERIOD,
        "stage2_authorized_only_if_stage1_passes": True,
        "frozen_for_candidate": [
            "B0=11mT", "G=2mT", "A_main=14.5deg", "A_cross=2.5deg",
            "TRUE-CEL 44x20x20 mesh/domain", "c0=100000mm/s",
            "water density/viscosity/EOS", "robot/wall geometry/contact",
            "Explicit scale factor=0.4; no mass scaling", "clean initial state",
        ],
        "online_fail_fast": {
            "position": "backtrack >0.020 mm and remains >0.010 mm below maximum for >=5 ms",
            "cycle": "completed cycle delta_s <= -0.010 mm",
            "numerical": "stable dt <25% of prior max for >=0.25 ms",
            "termination": "VUAMP iConcludeStep; partial ODB preserved",
        },
        "field_output": {"inherited_sparse_interval_s": 0.008333333333333333,
                         "RP_history_interval_s": 1e-5},
    })
    (DEST / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n", encoding="ascii")
    print(json.dumps({"job": JOB, "fresh_start": True, "frequency_Hz": 100,
        "cycle_period_s": PERIOD, "stage1_end_s": STAGE1_END,
        "physics_change_only": "frequency 120 -> 100 Hz",
        "input_sha256": sha(inp), "fortran_sha256": sha(f90)}, indent=2))


if __name__ == "__main__":
    main()
