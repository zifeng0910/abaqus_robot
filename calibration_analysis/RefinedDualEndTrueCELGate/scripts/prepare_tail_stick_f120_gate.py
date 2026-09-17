"""Prepare one full-cycle F120 strict-CEL, critically damped wall-contact gate."""
from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
PARENT_JOB = "TAIL_NO_REBOUND_ZETA080_8MS_TRUECEL"
JOB = "TAIL_STICK_F120_ZETA100_8P333MS_TRUECEL"
PARENT = OUT / "case" / PARENT_JOB
CASE = OUT / "case" / JOB


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    source = (PARENT / (PARENT_JOB + ".inp")).read_text(encoding="latin1")
    deck = source.replace(PARENT_JOB, JOB)
    deck, n_damping = re.subn(
        r"(\*Contact Damping, definition=CRITICAL DAMPING FRACTION, tangent fraction=0\.0\s*\n)0\.8\s*\n",
        r"\g<1>1.0\n", deck, count=1, flags=re.I)
    deck, n_duration = re.subn(r"\*Dynamic, Explicit\s*\n, 0\.008000000000",
                                "*Dynamic, Explicit\n, 0.008333333333", deck, count=1)
    deck = deck.replace("f=60Hz; A_main=14.8deg", "f=120Hz; A_main=14.8deg", 1)
    deck = deck.replace("duration=0.008000000000s", "duration=0.008333333333s", 1)
    if n_damping != 1 or n_duration != 1:
        raise RuntimeError("damping/duration replacement was not unique")

    CASE.mkdir(parents=True, exist_ok=True)
    inp = CASE / (JOB + ".inp")
    inp.write_text(deck, encoding="latin1")
    for name in ("vuforc_magnetic_only.f", "straight_control_centerline.dxf"):
        shutil.copy2(PARENT / name, CASE / name)

    identity = json.loads((PARENT / "case_identity.json").read_text())
    for stale in ("failure", "simulated_time_s", "stop_reason", "solver_end_state", "wallclock_s"):
        identity.pop(stale, None)
    identity.update({
        "case_id": JOB, "source_case": PARENT_JOB, "status": "PREPARED",
        "frequency_Hz": 120.0, "duration_s": 0.008333333333,
        "zeta": 1.0, "dynamics_run_count": 0,
        "contact_force_intent": "very weak wall support; no adhesion or artificial velocity clamp",
        "physics_changes": ["frequency 60 -> 120 Hz", "wall contact critical damping 0.80 -> 1.00",
                            "duration 8.000 -> 8.333333 ms (one complete commanded cycle)"],
        "input_sha256": sha(inp), "fortran_sha256": sha(CASE / "vuforc_magnetic_only.f"),
    })
    (CASE / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n")
    audit = {
        "job": JOB, "source_job": PARENT_JOB,
        "frequency_Hz": 120.0, "period_ms": 1000.0 / 120.0,
        "duration_ms": 8.333333333, "commanded_cycles": 1.0,
        "zeta": 1.0, "mu": identity["mu"],
        "contact_model": "hard nonpenetration with critical normal damping; no adhesion/clamp",
        "frozen": ["geometry", "CEL mesh/material", "A_main", "A_cross", "B0", "G",
                   "flow", "mu", "phase", "RouteA gauge", "HEAD/TAIL semantics"],
        "input_sha256": sha(inp),
    }
    (OUT / "TAIL_STICK_F120_setup_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
