"""Gate a 0-40 ms V+spatial-VR replay before any force interpretation."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "case" / "F100_G2P20_NOFLUID_REALWALL50"
JOB = "F100_G2P20_REALWALL_FIXEDVVR_FULLCEL40"
CASE = ROOT / "case" / JOB


def vector(z, stem, t):
    return np.column_stack([np.interp(t, z[f"{stem}{i}"][:, 0], z[f"{stem}{i}"][:, 1])
                            for i in (1, 2, 3)])


def main() -> None:
    source_identity = json.loads((SOURCE / "case_identity.json").read_text(encoding="utf-8-sig"))
    with np.load(SOURCE / "private" / "rp_history_private.npz") as source, \
         np.load(CASE / "private" / "rp_replay_private.npz") as replay:
        t = replay["U1"][:, 0].astype(float)
        if abs(t[0]) > 1e-10 or t[-1] < .04 - 1e-8 or np.any(np.diff(t) <= 0):
            raise RuntimeError("Replay does not cover a monotonic 0-40 ms trajectory")
        s = {stem: vector(source, stem, t) for stem in ("U", "UR", "V", "VR")}
        r = {stem: vector(replay, stem, t) for stem in ("U", "UR", "V", "VR")}
    c = np.asarray(source_identity["canonical_plus_s_axis_aba"], float)
    b = np.asarray(source_identity["b_routeA_aba"], float)
    pos = np.linalg.norm(r["U"]-s["U"], axis=1)
    ori = np.degrees((Rotation.from_rotvec(r["UR"]).inv() *
                      Rotation.from_rotvec(s["UR"])).magnitude())
    dvs = (r["V"]-s["V"]) @ c
    dwr = (r["VR"]-s["VR"]) @ b
    thresholds = {"max_position_mm": .001, "max_orientation_deg": .02,
                  "p99_abs_delta_v_s_mm_s": .05, "p99_abs_delta_omega_rock_rad_s": .5}
    results = {"max_position_mm": float(pos.max()),
               "max_orientation_deg": float(ori.max()),
               "p99_abs_delta_v_s_mm_s": float(np.percentile(np.abs(dvs), 99)),
               "p99_abs_delta_omega_rock_rad_s": float(np.percentile(np.abs(dwr), 99)),
               "max_abs_delta_v_s_mm_s": float(np.abs(dvs).max()),
               "max_abs_delta_omega_rock_rad_s": float(np.abs(dwr).max())}
    passed = all(results[key] < limit for key, limit in thresholds.items())
    summary = {"job": JOB, "source": SOURCE.name, "time_ms": [float(t[0]*1000), float(t[-1]*1000)],
               "method": "geodesic orientation error; native source global V and spatial/global VR",
               "thresholds": thresholds, "results": results,
               "initial_pose": {"source_U_mm": s["U"][0].tolist(), "replay_U_mm": r["U"][0].tolist(),
                                "source_UR_rad": s["UR"][0].tolist(), "replay_UR_rad": r["UR"][0].tolist()},
               "passed": passed, "classification": "KINEMATIC_REPLAY_VALID" if passed else "KINEMATIC_REPLAY_INVALID"}
    (ROOT / f"{JOB}_KINEMATIC_GATE.json").write_text(json.dumps(summary, indent=2) + "\n")
    csv_path = ROOT / f"{JOB}_REPLAY_VALIDATION.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("time_s", "position_error_mm", "orientation_error_deg",
                         "delta_v_s_mm_s", "delta_omega_rock_rad_s"))
        rows = np.unique(np.r_[np.arange(0, len(t), max(1, len(t)//4000)), len(t)-1])
        writer.writerows(zip(t[rows], pos[rows], ori[rows], dvs[rows], dwr[rows]))
    print(json.dumps(summary, indent=2))
    if not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
