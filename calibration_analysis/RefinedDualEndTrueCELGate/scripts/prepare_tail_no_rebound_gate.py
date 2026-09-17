"""Prepare the one 8 ms zeta=0.80 strict-CEL contact damping gate."""
from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
PARENT_JOB = "REFINED_DUALEND_F60_A14P8_G0P15_TRUECEL_10MS"
JOB = "TAIL_NO_REBOUND_ZETA080_8MS_TRUECEL"
PARENT = OUT / "case" / PARENT_JOB
CASE = OUT / "case" / JOB


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    table = pd.read_csv(PARENT / (PARENT_JOB + "_TAIL_contact_diagnosis.csv"))
    threshold = 0.005
    active = table.TAIL_gap_mm.to_numpy() <= threshold
    starts = np.where(active & np.r_[True, ~active[:-1]])[0]
    ends = np.where(active & np.r_[~active[1:], True])[0]
    if len(starts) < 2:
        raise RuntimeError("existing result does not show required rebound/recontact diagnosis")
    first_start, first_end, second_start = int(starts[0]), int(ends[0]), int(starts[1])
    time = table.time_s.to_numpy(float)
    gap = table.TAIL_gap_mm.to_numpy(float)
    impact_v = (gap[first_start] - gap[first_start - 1]) / (time[first_start] - time[first_start - 1])
    rebound_v = (gap[first_end + 1] - gap[first_end]) / (time[first_end + 1] - time[first_end])
    diagnosis = {
        "case": PARENT_JOB, "gap_threshold_mm": threshold,
        "first_TAIL_touch_ms": float(time[first_start] * 1e3),
        "first_TAIL_release_ms": float(time[first_end + 1] * 1e3),
        "second_TAIL_touch_ms": float(time[second_start] * 1e3),
        "N_tail_contacts_same_halfcycle": int(len(starts)),
        "V_impact_mm_s": float(impact_v), "V_rebound_mm_s": float(rebound_v),
        "R_v": float(abs(rebound_v) / abs(impact_v)),
        "classification": "REBOUND_RECONTACT_FAILURE",
    }
    (OUT / "existing_10ms_rebound_diagnosis.json").write_text(json.dumps(diagnosis, indent=2) + "\n")

    source = (PARENT / (PARENT_JOB + ".inp")).read_text(encoding="latin1")
    deck = source.replace(PARENT_JOB, JOB)
    deck, count_damping = re.subn(
        r"(\*Contact Damping, definition=CRITICAL DAMPING FRACTION, tangent fraction=0\.0\s*\n)0\.5\s*\n",
        r"\g<1>0.8\n", deck, count=1, flags=re.I)
    deck, count_duration = re.subn(r"\*Dynamic, Explicit\s*\n, 0\.010000000000",
                                    "*Dynamic, Explicit\n, 0.008000000000", deck, count=1)
    deck = deck.replace("duration=0.010000000000s", "duration=0.008000000000s", 1)
    if count_damping != 1 or count_duration != 1:
        raise RuntimeError("contact damping/duration replacement was not unique")
    CASE.mkdir(parents=True, exist_ok=True)
    inp = CASE / (JOB + ".inp")
    inp.write_text(deck, encoding="latin1")
    shutil.copy2(PARENT / "vuforc_magnetic_only.f", CASE / "vuforc_magnetic_only.f")
    shutil.copy2(PARENT / "straight_control_centerline.dxf", CASE / "straight_control_centerline.dxf")
    identity = json.loads((PARENT / "case_identity.json").read_text())
    identity.update({"case_id": JOB, "source_case": PARENT_JOB, "status": "PREPARED",
                     "duration_s": 0.008, "zeta": 0.8, "dynamics_run_count": 0,
                     "only_physics_change": "General Contact critical damping fraction 0.50 -> 0.80",
                     "input_sha256": sha(inp), "fortran_sha256": sha(CASE / "vuforc_magnetic_only.f")})
    (CASE / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n")
    audit = {"job": JOB, "source_job": PARENT_JOB, "zeta_old": 0.5, "zeta_new": 0.8,
             "duration_old_ms": 10.0, "duration_new_ms": 8.0,
             "frozen": ["geometry", "CEL mesh", "fluid", "A_main", "A_cross", "frequency", "B0", "G",
                        "flow", "mu", "phase", "RouteA gauge", "HEAD/TAIL semantics"],
             "existing_rebound_diagnosis": diagnosis}
    (OUT / "TAIL_NO_REBOUND_8MS_setup_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
