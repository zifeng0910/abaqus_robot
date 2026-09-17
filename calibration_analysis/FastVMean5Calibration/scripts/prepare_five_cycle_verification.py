"""Freeze the selected G=1.184 mT configuration for five-cycle verification."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
SOURCE_JOB = "FAST_VMEAN5_G1P184_UFLOW0P5_2CYCLES"
SOURCE = OUT / "case" / SOURCE_JOB
JOB = "FAST_PRECOMP_F120_VMEAN5_SELECTED"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest().upper()


def main():
    source_identity = json.loads((SOURCE / "case_identity.json").read_text())
    if source_identity["status"] != "SOLVED" or abs(source_identity["gradient_mT"] - 1.184) > 1e-12:
        raise RuntimeError("selected source is not the solved G=1.184 refinement")
    deck = (SOURCE / (SOURCE_JOB + ".inp")).read_text(encoding="latin1")
    deck = deck.replace(SOURCE_JOB, JOB)
    deck = deck.replace("TWO 120 HZ CYCLES", "FIVE 120 HZ CYCLES", 1)
    deck = deck.replace("2.0e-7, 0.016666666667", "2.0e-7, 0.041666666667", 1)
    deck = deck.replace("*Output, field, time interval=5.0e-5", "*Output, field, time interval=1.0e-4", 1)
    sub = (SOURCE / "vuamp_vmean5.f90").read_text(encoding="ascii")
    sub = sub.replace("FastVMean5Calibration\\case\\" + SOURCE_JOB,
                      "FastVMean5Calibration\\case\\" + JOB, 1)
    if "2.0e-7, 0.041666666667" not in deck or "*Output, field, time interval=1.0e-4" not in deck:
        raise RuntimeError("five-cycle deck replacement failed")
    case = OUT / "case" / JOB
    case.mkdir(parents=True, exist_ok=True)
    inp = case / (JOB + ".inp")
    f90 = case / "vuamp_vmean5.f90"
    inp.write_text(deck, encoding="latin1")
    f90.write_text(sub, encoding="ascii")
    identity = dict(source_identity)
    for key in ("wallclock_s", "cpus", "socket_calls", "input_sha256", "fortran_sha256"):
        identity.pop(key, None)
    identity.update({
        "case_id": JOB, "source_case": SOURCE_JOB, "status": "PREPARED",
        "duration_s": 0.041666666667, "commanded_cycles": 5.0, "dynamics_run_count": 0,
        "calibration_role": "SELECTED_FIVE_CYCLE_VERIFICATION",
        "selection_basis": "G=1.184 mT yielded 5.068059 mm/s over two cycles and passed the stated hard contact gate",
        "field_output_interval_s": 1.0e-4, "history_output_interval_s": 5.0e-5,
        "input_sha256": sha(inp), "fortran_sha256": sha(f90),
    })
    (case / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n", encoding="ascii")
    manifest = {"job": JOB, "selected_G_mT": 1.184, "duration_s": 0.041666666667,
                "input_sha256": sha(inp), "fortran_sha256": sha(f90)}
    (OUT / "selected_case_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="ascii")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
