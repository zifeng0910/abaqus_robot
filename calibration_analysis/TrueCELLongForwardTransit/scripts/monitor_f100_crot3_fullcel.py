"""Read-only status summary for the running clean FULL CEL damping job."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
JOB = "TRUECEL_B0P11_G2P20_A14P5_F100_CROT3_FULLCEL50"
CASE = ROOT / "case" / JOB


def read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        try:
            return list(csv.DictReader(handle))
        except (csv.Error, ValueError):
            return []


def main() -> None:
    identity = json.loads((CASE / "case_identity.json").read_text(encoding="utf-8-sig"))
    sta_path = CASE / f"{JOB}.sta"
    sta = sta_path.read_text(encoding="latin1") if sta_path.exists() else ""
    progress = []
    for line in sta.splitlines():
        match = re.match(r"\s*(\d+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)\s+\d\d:\d\d:\d\d\s+([0-9.E+-]+)", line)
        if match:
            progress.append({"increment": int(match.group(1)),
                             "step_time_s": float(match.group(2)),
                             "total_time_s": float(match.group(3)),
                             "stable_dt_s": float(match.group(4))})
    motion = read_csv(CASE / "online_motion.csv")
    cycles = read_csv(CASE / "online_cycle_gate.csv")
    online = None
    if motion:
        try:
            online = {k: float(v) for k, v in motion[-1].items()}
        except (TypeError, ValueError):
            online = None
    result = {"job": JOB, "identity_status": identity["status"],
              "solver_completed": "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in sta,
              "solver_time_s": progress[-1]["step_time_s"] if progress else None,
              "solver_increment": progress[-1]["increment"] if progress else None,
              "latest_stable_dt_s": progress[-1]["stable_dt_s"] if progress else None,
              "minimum_reported_stable_dt_s": min((p["stable_dt_s"] for p in progress), default=None),
              "latest_online": online, "completed_cycles": cycles,
              "restart_read": False,
              "fatal_keyword_in_sta": "***ERROR" in sta or "Abaqus Error" in sta,
              "run_log_exists": (CASE / f"{JOB}.run.log").exists(),
              "odb_exists": (CASE / f"{JOB}.odb").exists()}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
