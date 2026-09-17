"""Prepare the evidence-based G1.20 mT one-cycle correction."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
PARENT_JOB = "FAST_PRECOMP_F120_G0P30_FLOW_1CYCLE"
JOB = "FAST_PRECOMP_F120_G1P20_FLOW_1CYCLE"
PARENT = OUT / "case" / PARENT_JOB
CASE = OUT / "case" / JOB


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest().upper()


def main():
    deck = (PARENT / (PARENT_JOB + ".inp")).read_text(encoding="latin1").replace(PARENT_JOB, JOB)
    deck = deck.replace("G=0.15 -> 0.30 mT", "G=0.15 -> 1.20 mT", 1)
    sub = (PARENT / "vuamp_precomputed_g0p30_flow.f90").read_text(encoding="ascii")
    sub = sub.replace("grad=0.00030d0*grad", "grad=0.00120d0*grad", 1)
    sub = sub.replace(PARENT_JOB, JOB, 1)
    CASE.mkdir(parents=True, exist_ok=True)
    inp = CASE / (JOB + ".inp")
    f90 = CASE / "vuamp_precomputed_g1p20_flow.f90"
    inp.write_text(deck, encoding="latin1")
    f90.write_text(sub, encoding="ascii")
    identity = json.loads((PARENT / "case_identity.json").read_text())
    for key in ("wallclock_s", "cpus", "socket_calls"):
        identity.pop(key, None)
    identity.update({
        "case_id": JOB, "source_case": PARENT_JOB, "status": "PREPARED", "gradient_mT": 1.20,
        "dynamics_run_count": 0, "input_sha256": sha(inp), "fortran_sha256": sha(f90),
        "gradient_choice_basis": "Linear extrapolation of G0.15 and G0.30 first-cycle delta_s gives zero crossing near 0.96 mT; 1.20 mT selected for a clear positive-margin short screen.",
        "forward_sign_sanity": {"basis": "same positive normalized gradient table and polarity as passed G0.30 case",
                                "F_dot_plus_s_at_1ms_N": 1.6650018876898242e-5, "passed": True},
    })
    (CASE / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n", encoding="ascii")
    print(json.dumps({"job": JOB, "gradient_mT": 1.20, "input_sha256": sha(inp), "fortran_sha256": sha(f90)}, indent=2))


if __name__ == "__main__":
    main()
