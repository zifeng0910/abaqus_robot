"""Prepare the one authorized 100 ms production local-frame Abaqus case."""
import hashlib
import json
import re
import shutil
from pathlib import Path


HERE = Path(__file__).resolve().parent
VALIDATION = HERE.parent
SOURCE = VALIDATION.parent / "MagneticDriveFrameAudit" / "cases" / "FRAME_LOCAL30"
SOURCE_JOB = "MDA_FRAME_LOCAL30"
JOB = "PROD_LOCAL30_G0_3CYCLE"
CASE = VALIDATION / "case" / JOB


def main():
    CASE.mkdir(parents=True, exist_ok=True)
    source_deck = SOURCE / f"{SOURCE_JOB}.inp"
    source_fortran = SOURCE / "vuforc_audit.f"
    deck = source_deck.read_text()
    deck = deck.replace(SOURCE_JOB, JOB)
    deck, count = re.subn(r"(\*Dynamic, Explicit, DIRECT\n)1\.0e-7,\s*0\.033333333",
                          r"\g<1>1.0e-7, 0.100000000", deck, count=1)
    if count != 1:
        raise RuntimeError("Dynamic duration replacement failed")
    header = ("** PRODUCTION ROBOT-LOCAL FRAME VALIDATION; EXACTLY THREE 30 HZ CYCLES\n"
              "** ROBOT_LOCAL_TANGENT; B0=10mT; cone=30deg; Bias=0; G=0; ramp=1ms\n")
    deck_path = CASE / f"{JOB}.inp"
    deck_path.write_text(header + deck)
    shutil.copyfile(source_fortran, CASE / "vuforc_production_local.f")
    config = json.loads((VALIDATION / "config" / f"{JOB}.json").read_text())
    config.update({
        "job_name": JOB, "status": "PREPARED", "source_case": "FRAME_LOCAL30",
        "source_input_sha256": hashlib.sha256(source_deck.read_bytes()).hexdigest(),
        "input_sha256": hashlib.sha256(deck_path.read_bytes()).hexdigest(),
        "fortran_sha256": hashlib.sha256(source_fortran.read_bytes()).hexdigest(),
        "production_server_sha256": hashlib.sha256((VALIDATION / "production" / "magpylib_socket_server.py").read_bytes()).hexdigest(),
    })
    (CASE / "case_identity.json").write_text(json.dumps(config, indent=2) + "\n")
    if len(re.findall(r"(?im)^\*node\s*$", deck)) < 1:
        raise RuntimeError("Missing nodes")
    node_count = 1588
    element_count = 7302
    if config["node_count"] != node_count or config["element_count"] != element_count:
        raise RuntimeError("Configured coarse mesh identity mismatch")
    if "1.0e-7, 0.100000000" not in deck or "0.033333333" in deck:
        raise RuntimeError("Three-cycle duration identity failed")
    print(json.dumps(config, indent=2))


if __name__ == "__main__":
    main()
