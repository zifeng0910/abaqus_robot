"""Prepare the four fixed-geometry fast straight dynamic screening cases."""
from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
SOURCE_JOB = "PROD_LOCAL_ROCK10_5HZ_G0_STRAIGHT"
SOURCE = REPO / "calibration_analysis" / "LocalHeadTailRocking" / "case" / SOURCE_JOB
SERVER = OUT / "production" / "magpylib_socket_server_fast.py"

CASES = (
    ("FAST_PLANAR_15HZ_A12_FORWARD", 15.0, 12.0, 0.0, 0.100, "ROBOT_LOCAL_ROCKING"),
    ("FAST_PLANAR_20HZ_A12_FORWARD", 20.0, 12.0, 0.0, 0.075, "ROBOT_LOCAL_ROCKING"),
    ("FAST_PLANAR_20HZ_A14P343_FORWARD", 20.0, 14.343111711438091, 0.0, 0.075, "ROBOT_LOCAL_ROCKING"),
    ("FAST_ELLIPTIC_20HZ_A14P343_X2P5_FORWARD", 20.0, 14.343111711438091, 2.5, 0.075,
     "ROBOT_LOCAL_ELLIPTIC_ROCKING"),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    source_identity = json.loads((SOURCE / "case_identity.json").read_text())
    source_deck_path = SOURCE / f"{SOURCE_JOB}.inp"
    source_deck = source_deck_path.read_text()
    if source_identity["status"] != "SOLVED":
        raise RuntimeError("ROCK10 source is not solved")
    if "** Reduced-Hydro: CEL fluid part, material, initialization and CEL contact were removed." not in source_deck:
        raise RuntimeError("Source is not the verified non-CEL ReducedHydro deck")

    cases_root = OUT / "cases"
    cases_root.mkdir(parents=True, exist_ok=True)
    for job, frequency, main_amp, cross_amp, duration, mode in CASES:
        case = cases_root / job
        case.mkdir(parents=True, exist_ok=True)
        deck = source_deck.replace(SOURCE_JOB, job)
        deck, count = re.subn(
            r"(\*Dynamic, Explicit, DIRECT\s*\n)1\.0e-7,\s*0\.200000000",
            rf"\g<1>1.0e-7, {duration:.9f}", deck, count=1)
        if count != 1:
            raise RuntimeError(f"{job}: duration replacement count={count}")
        header = (
            "** FAST STRAIGHT DYNAMIC SCREEN\n"
            f"** {job}; mode={mode}; f={frequency:.15g}Hz; A_main={main_amp:.15g}deg; "
            f"A_cross={cross_amp:.15g}deg; G=6mT; L=45mm; duration={duration:.9f}s\n"
            "** NON-CEL; SAME RIGID ROBOT/TUBE, ONE-WALL GENERAL CONTACT, MU=0.03, "
            "ZETA=0.50, REDUCEDHYDRO, INITIAL POSE, ROUTEA GAUGE, AND DT=1E-7\n"
        )
        deck_path = case / f"{job}.inp"
        deck_path.write_text(header + deck, encoding="ascii")
        shutil.copyfile(SOURCE / "vuforc_production_local.f", case / "vuforc_production_local.f")
        shutil.copyfile(SOURCE / "straight_control_centerline.dxf", case / "straight_control_centerline.dxf")

        identity = {
            "case_id": job,
            "status": "PREPARED",
            "source_job": SOURCE_JOB,
            "field_frame_mode": mode,
            "duration_s": duration,
            "direct_dt_s": 1e-7,
            "B0_mT": 10.0,
            "frequency_Hz": frequency,
            "rocking_main_amplitude_deg": main_amp,
            "rocking_cross_amplitude_deg": cross_amp,
            "gradient_mT": 6.0,
            "gradient_length_mm": 45.0,
            "gradient_profile": "legacy",
            "forward_definition": "canonical increasing centerline arclength +s",
            "mu": 0.03,
            "zeta": 0.50,
            "robot_length_mm": source_identity["robot_length_mm"],
            "robot_diameter_mm": source_identity["robot_diameter_mm"],
            "robot_moment_Am2": source_identity["robot_moment_Am2"],
            "robot_mass_mg": source_identity["robot_mass_mg"],
            "initial_center_aba_mm": source_identity["initial_center_aba_mm"],
            "initial_axis_aba": source_identity["initial_axis_aba"],
            "n_routeA_aba": source_identity["n_rock_aba"],
            "b_routeA_aba": source_identity["b_rock_aba"],
            "routeA_gauge_from_production_e1_deg": source_identity["routeA_gauge_from_production_e1_deg"],
            "source_input_sha256": sha256(source_deck_path),
            "input_sha256": sha256(deck_path),
            "fortran_sha256": sha256(case / "vuforc_production_local.f"),
            "centerline_sha256": sha256(case / "straight_control_centerline.dxf"),
            "fast_server_sha256": sha256(SERVER),
        }
        (case / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n", encoding="ascii")

        invariants = (
            "*Rigid Body, ref node=RP_ROBOT, elset=ROBOT_SOLID_CEL_ALL",
            "ROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF",
            "0.03,",
            "\n0.5\n",
            "HYDRO_FX",
        )
        if not all(token in deck for token in invariants):
            raise RuntimeError(f"{job}: frozen invariant missing")
        if "*Eulerian" in deck or "*Eulerian Section" in deck:
            raise RuntimeError(f"{job}: unexpected true CEL definition")
        print(f"PREPARED {job}")


if __name__ == "__main__":
    main()
