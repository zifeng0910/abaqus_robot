"""Prepare one short G0.30 mT precomputed-magnetic FAST_SURROGATE screen."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
PARENT_JOB = "FAST_PRECOMP_F120_DUALEND_NOREBOUND"
JOB = "FAST_PRECOMP_F120_G0P30_FLOW_1CYCLE"
PARENT = REPO / "calibration_analysis" / "FastPrecomputedNoReboundScreen" / "case" / PARENT_JOB
CASE = OUT / "case" / JOB
TABLE_DIR = REPO / "calibration_analysis" / "PrecomputedStraightMagneticBackend"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def main() -> None:
    source_inp = PARENT / (PARENT_JOB + ".inp")
    deck = source_inp.read_text(encoding="latin1").replace(PARENT_JOB, JOB)
    replacements = {
        "** AUTHORITATIVE CASE: " + JOB: "** AUTHORITATIVE CASE: " + JOB + "\n** SOLE PHYSICS CHANGE: MAGNETIC GRADIENT G=0.15 -> 0.30 mT",
        "2.0e-7, 0.016666666667": "2.0e-7, 0.008333333333",
        "*Output, field, time interval=1.0e-5, time marks=NO": "*Output, field, time interval=5.0e-5, time marks=NO",
        "** NON-CEL FAST_SURROGATE_FLUID; PRECOMPUTED_TABLE; TWO 120 HZ CYCLES": "** NON-CEL FAST_SURROGATE_FLUID; PRECOMPUTED_TABLE; ONE 120 HZ CYCLE",
    }
    for old, new in replacements.items():
        if deck.count(old) != 1:
            raise RuntimeError("expected one deck token: {}".format(old))
        deck = deck.replace(old, new, 1)

    source_fortran = REPO / "calibration_analysis" / "FastPrecomputedNoReboundScreen" / "vuamp_precomputed_fast_surrogate.f90"
    fortran = source_fortran.read_text(encoding="ascii")
    if fortran.count("grad=0.00015d0*grad") != 1:
        raise RuntimeError("gradient scale token changed")
    fortran = fortran.replace("grad=0.00015d0*grad", "grad=0.00030d0*grad", 1)
    fortran = fortran.replace("FastPrecomputedNoReboundScreen\\case\\FAST_PRECOMP_F120_DUALEND_NOREBOUND",
                              "FastPrecomputedForwardFlowScreen\\case\\" + JOB, 1)

    CASE.mkdir(parents=True, exist_ok=True)
    inp = CASE / (JOB + ".inp")
    sub = CASE / "vuamp_precomputed_g0p30_flow.f90"
    inp.write_text(deck, encoding="latin1")
    sub.write_text(fortran, encoding="ascii")

    identity = json.loads((PARENT / "case_identity.json").read_text())
    for key in ("solved_input_sha256", "artifact_input_sha256", "solved_fortran_sha256", "maintained_fortran_sha256",
                "post_run_nonphysics_changes", "wallclock_s", "lookup_log_line"):
        identity.pop(key, None)
    identity.update({
        "case_id": JOB, "source_case": PARENT_JOB, "status": "PREPARED",
        "gradient_mT": 0.30, "duration_s": 0.008333333333, "commanded_cycles": 1.0,
        "field_output_interval_s": 5.0e-5, "dynamics_run_count": 0,
        "change_scope": ["gradient G 0.15 -> 0.30 mT", "duration 2 -> 1 cycle", "field output 10 -> 50 us"],
        "frozen_physics": ["geometry", "rigid robot", "General Contact law", "mu=0.03", "zeta=1.0",
                           "B0=10mT", "f=120Hz", "A_main=14.8deg", "A_cross=2.5deg", "initial pose",
                           "FAST_SURROGATE_FLUID coefficients", "U_flow=+10mm/s"],
        "input_sha256": sha(inp), "fortran_sha256": sha(sub), "socket_calls_expected": 0,
    })

    sys.path.insert(0, str(TABLE_DIR / "scripts"))
    import table_model
    cfg = table_model.config()
    cfg["gradient_G_T"] = 0.00030
    phase, spatial = table_model.read_table(TABLE_DIR / "magnetic_field_gradient_table.dat")
    sample_t = 0.001
    _, _, force, _, _, _ = table_model.table_loads(sample_t, cfg["initial_center_aba_mm"], [0, 0, 0], cfg, phase, spatial)
    c = np.asarray(cfg["canonical_plus_s_axis_aba"], float)
    projection = float(np.dot(force, c))
    identity["forward_sign_sanity"] = {"sample_time_s": sample_t, "F_dot_plus_s_N": projection, "passed": projection > 0}
    if projection <= 0:
        raise RuntimeError("forward gradient sign sanity failed")
    (CASE / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n", encoding="ascii")
    (OUT / "forward_sign_sanity.json").write_text(json.dumps(identity["forward_sign_sanity"], indent=2) + "\n", encoding="ascii")
    print(json.dumps({"job": JOB, "F_dot_plus_s_N": projection, "input_sha256": sha(inp), "fortran_sha256": sha(sub)}, indent=2))


if __name__ == "__main__":
    main()
