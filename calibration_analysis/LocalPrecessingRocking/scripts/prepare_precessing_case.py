"""Prepare the single authorized precessing-rocking dynamic case."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
SOURCE_JOB = "PROD_LOCAL_ROCK_TOUCH_5HZ_G0_STRAIGHT"
JOB = "PROD_LOCAL_PRECESSROCK_10HZ_PREC2P5_G0_STRAIGHT"
SOURCE = REPO / "calibration_analysis" / "LocalHeadTailRocking" / "case" / SOURCE_JOB
CASE = OUT / "case" / JOB
SERVER = REPO / "calibration_analysis" / "ProductionLocalFrameValidation" / "production" / "magpylib_socket_server.py"
AMPLITUDE_DEG = 14.343111711438091
ROUTE_A_DEG = -61.37284757596327
PSI0_DEG = -83.87284757596327


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    source_identity = json.loads((SOURCE / "case_identity.json").read_text())
    if source_identity["status"] != "SOLVED":
        raise RuntimeError("Authoritative LIGHTCONTACT source case is not SOLVED")
    if float(source_identity["rocking_field_amplitude_deg"]) != AMPLITUDE_DEG:
        raise RuntimeError("LIGHTCONTACT amplitude identity drift")
    if CASE.exists() and any(CASE.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty case directory: {CASE}")
    CASE.mkdir(parents=True, exist_ok=True)

    source_deck = SOURCE / f"{SOURCE_JOB}.inp"
    deck = source_deck.read_text().replace(SOURCE_JOB, JOB)
    deck = deck.replace(
        f"** ROBOT_LOCAL_ROCKING; B0=10mT; alpha_B0={AMPLITUDE_DEG:.12f}deg; f=5Hz; phase=0; G=0",
        f"** ROBOT_LOCAL_PRECESSING_ROCKING; B0=10mT; alpha_B0={AMPLITUDE_DEG:.12f}deg; f_rock=10Hz; f_prec=2.5Hz; psi0={PSI0_DEG:.12f}deg; G=0")
    deck_path = CASE / f"{JOB}.inp"
    deck_path.write_text(deck)
    for name in ("vuforc_production_local.f", "straight_control_centerline.dxf"):
        shutil.copyfile(SOURCE / name, CASE / name)

    identity = dict(source_identity)
    identity.update({
        "case_id": JOB, "status": "PREPARED", "source_job": SOURCE_JOB,
        "field_frame_mode": "ROBOT_LOCAL_PRECESSING_ROCKING",
        "frequency_Hz": 10.0, "rocking_frequency_Hz": 10.0,
        "precession_frequency_Hz": 2.5,
        "rocking_field_amplitude_deg": AMPLITUDE_DEG,
        "precession_phase_deg": PSI0_DEG,
        "routeA_gauge_from_production_e1_deg": ROUTE_A_DEG,
        "first_positive_rocking_peak_s": 0.025,
        "precession_azimuth_at_first_positive_peak_deg": ROUTE_A_DEG,
        "source_input_sha256": sha256(source_deck),
        "input_sha256": sha256(deck_path),
        "fortran_sha256": sha256(CASE / "vuforc_production_local.f"),
        "centerline_sha256": sha256(CASE / "straight_control_centerline.dxf"),
        "production_server_sha256": sha256(SERVER),
        "implementation_commit": "e7632dec103ac8ad47ec7dac5b736e1dabe61a50",
    })
    identity.pop("wallclock_s", None)
    (CASE / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n")

    exact = {
        "duration_s": .2, "direct_dt_s": 1e-7, "B0_mT": 10.0,
        "gradient_mT": 0.0, "mu": .03, "zeta": .5,
        "robot_length_mm": 2.4, "robot_diameter_mm": .815,
        "robot_moment_Am2": .0010876227522174417,
    }
    for key, expected in exact.items():
        if float(identity[key]) != expected:
            raise RuntimeError(f"Invariant drift: {key}={identity[key]!r}, expected {expected!r}")
    if sha256(CASE / "vuforc_production_local.f") != sha256(SOURCE / "vuforc_production_local.f"):
        raise RuntimeError("Fortran bridge drift")
    if sha256(CASE / "straight_control_centerline.dxf") != sha256(SOURCE / "straight_control_centerline.dxf"):
        raise RuntimeError("Centerline drift")
    if "ROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF" not in deck:
        raise RuntimeError("One-wall General Contact topology missing")
    print(json.dumps(identity, indent=2))


if __name__ == "__main__":
    main()
