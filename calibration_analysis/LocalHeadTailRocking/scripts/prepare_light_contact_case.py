"""Prepare the single amplitude-derived light-contact dynamic case."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
SOURCE_JOB = "PROD_LOCAL_ROCK10_5HZ_G0_STRAIGHT"
JOB = "PROD_LOCAL_ROCK_TOUCH_5HZ_G0_STRAIGHT"
SOURCE = OUT / "case" / SOURCE_JOB
CASE = OUT / "case" / JOB
SERVER = OUT.parents[1] / "calibration_analysis" / "ProductionLocalFrameValidation" / "production" / "magpylib_socket_server.py"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    selection_path = OUT / "Selected_LightContact_Amplitude.json"
    selection = json.loads(selection_path.read_text())
    amplitude = float(selection["selected_A_B_required_deg"])
    if selection["classification"] != "LIGHT_CONTACT_AMPLITUDE_SELECTED" or amplitude > 20.0:
        raise RuntimeError("Light-contact amplitude gate did not pass")
    source_identity = json.loads((SOURCE / "case_identity.json").read_text())
    if source_identity["status"] != "SOLVED":
        raise RuntimeError("Authoritative ROCK10 source case is not SOLVED")
    if CASE.exists() and any(CASE.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty case directory: {CASE}")
    CASE.mkdir(parents=True, exist_ok=True)

    source_deck = SOURCE / f"{SOURCE_JOB}.inp"
    deck = source_deck.read_text().replace(SOURCE_JOB, JOB)
    deck = deck.replace("alpha_B0=10deg", f"alpha_B0={amplitude:.12f}deg")
    deck_path = CASE / f"{JOB}.inp"
    deck_path.write_text(deck)
    for name in ("vuforc_production_local.f", "straight_control_centerline.dxf"):
        shutil.copyfile(SOURCE / name, CASE / name)

    identity = dict(source_identity)
    identity.update({
        "case_id": JOB,
        "status": "PREPARED",
        "source_job": SOURCE_JOB,
        "rocking_field_amplitude_deg": amplitude,
        "light_contact_design_manifest_sha256": sha256(selection_path),
        "source_input_sha256": sha256(source_deck),
        "input_sha256": sha256(deck_path),
        "fortran_sha256": sha256(CASE / "vuforc_production_local.f"),
        "centerline_sha256": sha256(CASE / "straight_control_centerline.dxf"),
        "production_server_sha256": sha256(SERVER),
    })
    identity.pop("wallclock_s", None)
    (CASE / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n")

    if identity["B0_mT"] != 10.0 or identity["frequency_Hz"] != 5.0 or identity["gradient_mT"] != 0.0:
        raise RuntimeError("Magnetic invariant drift")
    if identity["mu"] != 0.03 or identity["zeta"] != 0.5 or identity["direct_dt_s"] != 1e-7 or identity["duration_s"] != 0.2:
        raise RuntimeError("Contact or numerical invariant drift")
    if identity["robot_length_mm"] != 2.4 or identity["robot_diameter_mm"] != 0.815:
        raise RuntimeError("Geometry invariant drift")
    if "ROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF" not in deck:
        raise RuntimeError("One-wall General Contact topology missing")
    print(json.dumps(identity, indent=2))


if __name__ == "__main__":
    main()
