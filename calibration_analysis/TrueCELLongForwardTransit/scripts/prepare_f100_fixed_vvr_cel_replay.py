"""Prepare one full-CEL force replay driven by native no-fluid wall-contact V/VR."""
from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PARENT = "TRUECEL_B0P11_G2P20_A14P5_F100_CLEAN50"
SOURCE = "F100_G2P20_NOFLUID_REALWALL50"
JOB = "F100_G2P20_REALWALL_FIXEDVVR_FULLCEL40"
TABLE = "magnetic_field_gradient_table_B0P11_A14P5.dat"
STOP = .04


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def replace_once(deck: str, old: str, new: str) -> str:
    if deck.count(old) != 1:
        raise RuntimeError(f"Expected one occurrence of {old[:80]!r}; got {deck.count(old)}")
    return deck.replace(old, new, 1)


def main() -> None:
    parent = ROOT / "case" / PARENT
    source = ROOT / "case" / SOURCE
    target = ROOT / "case" / JOB
    if target.exists():
        raise RuntimeError(f"Refusing to overwrite {target}")
    pi = json.loads((parent / "case_identity.json").read_text(encoding="utf-8-sig"))
    si = json.loads((source / "case_identity.json").read_text(encoding="utf-8-sig"))
    for key in ("B0_mT", "gradient_mT", "frequency_Hz", "rocking_main_amplitude_deg",
                "rocking_cross_amplitude_deg", "initial_center_aba_mm", "canonical_plus_s_axis_aba",
                "n_routeA_aba", "b_routeA_aba", "robot_moment_Am2", "lumen_radius_mm"):
        if pi[key] != si[key]:
            raise RuntimeError(f"Parent/source identity differs: {key}")
    if sha(parent / TABLE) != sha(source / TABLE) or sha(parent / TABLE) != pi["magnetic_table_sha256"]:
        raise RuntimeError("Magnetic table identity failed")
    if si["status"] != "SOLVED" or si["dynamics_run_count"] != 1:
        raise RuntimeError("Source trajectory is not the completed single clean run")
    source_history = source / "private" / "rp_history_private.npz"
    with np.load(source_history) as z:
        t = z["V1"][:, 0].astype(float)
        use = t <= STOP + 1e-10
        t = t[use]
        if abs(t[0]) > 1e-12 or abs(t[-1]-STOP) > 1e-8 or np.any(np.diff(t) <= 0):
            raise RuntimeError("Native V/VR history does not span 0-40 ms monotonically")
        for stem in ("V", "VR"):
            for axis in (1, 2, 3):
                a = z[f"{stem}{axis}"]
                if a.shape[0] < len(t) or not np.allclose(a[:len(t), 0], t, atol=1e-9):
                    raise RuntimeError(f"Native time mismatch for {stem}{axis}")
    parent_deck = (parent / f"{PARENT}.inp").read_text(encoding="latin1")
    deck = parent_deck.replace(PARENT, JOB)
    deck = replace_once(deck, ", 0.050000000000\n*Bulk Viscosity", ", 0.040000000000\n*Bulk Viscosity")
    deck = replace_once(deck, "** OUTPUT REQUESTS\n",
                        "** VALIDATED GLOBAL TRANSLATIONAL V + GLOBAL/SPATIAL VR REPLAY\n"
                        "*Include, input=driver_vvr_amplitudes.inp\n** OUTPUT REQUESTS\n")
    old_history = ("*Contact Output, surface=ROBOT_SOLID-1.ROBOT_SOLID_SURF\nCFN, CFS, CFT\n"
                   "*Contact Output, surface=Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF\nCFN, CFS, CFT\n")
    new_history = ("*Contact Output, surface=ROBOT_SOLID-1.ROBOT_SOLID_SURF\nCFN, CFS, CFT\n"
                   "*Contact Output, surface=ROBOT_SOLID-1.ROBOT_SOLID_SURF, "
                   "SECOND SURFACE=Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF\nCFN, CFS, CFT\n"
                   "*Contact Output, surface=Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF\nCFN, CFS, CFT\n")
    deck = replace_once(deck, old_history, new_history)
    deck = replace_once(deck, "*Output, history, frequency=1\n",
                        "*Output, field, time interval=2.5e-5, time marks=NO\n"
                        "*Contact Output, surface=ROBOT_SOLID-1.ROBOT_SOLID_SURF\nCDISP\n"
                        "*Output, history, frequency=1\n")
    if "*Restart, read" in deck or deck.count("*Step, name=Step_Drive") != 1:
        raise RuntimeError("Replay must have one fresh Explicit step")
    if any(token in deck for token in ("REPLAY_U1", "REPLAY_R1", "REPLAY_UR1")):
        raise RuntimeError("Invalid displacement rotation driver")
    target.mkdir(parents=True)
    inp = target / f"{JOB}.inp"
    inp.write_text(deck, encoding="latin1")
    driver = target / "driver_vvr_amplitudes.inp"
    with np.load(source_history) as z, driver.open("w", encoding="ascii", newline="\n") as handle:
        handle.write("** Native source RP history, validated V+spatial VR method; 0-40 ms\n")
        for stem, start in (("V", 1), ("VR", 4)):
            for axis in (1, 2, 3):
                name = f"REPLAY_{stem}{axis}"
                handle.write(f"*Amplitude, name={name}, definition=TABULAR, time=STEP TIME\n")
                values = z[f"{stem}{axis}"][:len(t), 1]
                for ti, value in zip(t, values):
                    handle.write(f"{ti:.16g}, {value:.16g}\n")
                dof = start + axis - 1
                handle.write(f"*Boundary, amplitude={name}, type=VELOCITY\n"
                             f"RP_ROBOT, {dof}, {dof}, 1.\n")
    original_sub = (parent / "vuamp_precomputed_truecel.f90").read_bytes()
    if original_sub.count(PARENT.encode()) != 6:
        raise RuntimeError("Unexpected VUAMP path count")
    f90 = target / "vuamp_precomputed_truecel.f90"
    f90.write_bytes(original_sub.replace(PARENT.encode(), JOB.encode()))
    if f90.read_bytes().replace(JOB.encode(), PARENT.encode()) != original_sub:
        raise RuntimeError("VUAMP changed beyond case paths")
    for name in (TABLE, TABLE + ".json"):
        shutil.copy2(parent / name, target / name)
    stripped = deck.replace("*Include, input=driver_vvr_amplitudes.inp\n", "").replace(
        ", 0.040000000000\n*Bulk Viscosity", ", 0.050000000000\n*Bulk Viscosity")
    stripped = stripped.replace(new_history, old_history).replace(
        "*Output, field, time interval=2.5e-5, time marks=NO\n"
        "*Contact Output, surface=ROBOT_SOLID-1.ROBOT_SOLID_SURF\nCDISP\n", "")
    stripped = stripped.replace("** VALIDATED GLOBAL TRANSLATIONAL V + GLOBAL/SPATIAL VR REPLAY\n", "")
    if stripped != parent_deck.replace(PARENT, JOB):
        raise RuntimeError("Parent physics changed beyond prescribed motion/duration/output")
    audit = {
        "job": JOB, "source_trajectory": SOURCE, "parent_CEL": PARENT,
        "source_RP_history_sha256": sha(source_history), "parent_input_sha256": sha(parent / f"{PARENT}.inp"),
        "input_sha256": sha(inp), "driver_sha256": sha(driver), "driver_points_per_component": len(t),
        "driver_method": "validated prescribed global V(t) plus global/spatial VR(t); unit velocity BC references",
        "driver_time_s": [float(t[0]), float(t[-1])], "single_step_duration_s": STOP,
        "magnetic_table_sha256": sha(target / TABLE), "fortran_sha256": sha(f90),
        "physics_parent_equivalent_except_prescribed_motion": True,
        "contact_outputs": "whole robot, pair-specific robot-wall, RP reaction; robot-CEL inferred as whole minus wall",
        "restart_read": False, "rotational_damping_added": False,
        "validated_thresholds": {"max_position_mm": .001, "max_orientation_deg": .02,
                                 "p99_abs_delta_v_s_mm_s": .05, "p99_abs_delta_omega_rock_rad_s": .5},
    }
    (target / "setup_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    identity = dict(pi)
    identity.update(case_id=JOB, source_trajectory=SOURCE, status="PREPARED",
                    classification="PENDING", dynamics_run_count=0, duration_s=STOP,
                    input_sha256=sha(inp), driver_sha256=sha(driver),
                    restart_read=False, prescribed_motion="GLOBAL_V_PLUS_SPATIAL_VR")
    (target / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
