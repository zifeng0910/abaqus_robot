"""Prepare the single authorized prescribed-wobble replay."""
from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SOURCE_JOB = "PROD_LOCAL30_G0_STRAIGHT_CTRL"
JOB = "STRAIGHT_PRESCRIBED_WOBBLE_CONTACT_AUDIT"
SOURCE = REPO / "calibration_analysis/StraightPipeControl/case" / SOURCE_JOB
CASE = HERE / "case" / JOB
DURATION = 0.016667
DT = 1e-7
RP0 = np.array([-7.468174204284, -3.676918015967, -9.550745259298])
A0 = np.array([0.9647382600216, -0.1188742372140, 0.2348382536499])
A0 /= np.linalg.norm(A0)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unit(a):
    return a / np.linalg.norm(a)


def align_rotation(source, target):
    cross = np.cross(source, target)
    sine = np.linalg.norm(cross)
    cosine = np.clip(np.dot(source, target), -1.0, 1.0)
    if sine < 1e-14:
        return Rotation.identity()
    return Rotation.from_rotvec(cross / sine * np.arctan2(sine, cosine))


def amplitude(name, time, values):
    lines = [f"*Amplitude, name={name}, time=TOTAL TIME"]
    pairs = [f"{t:.9g}, {v:.12g}" for t, v in zip(time, values)]
    for i in range(0, len(pairs), 4):
        lines.append(", ".join(pairs[i:i + 4]))
    return "\n".join(lines)


def main():
    CASE.mkdir(parents=True, exist_ok=True)
    identity = json.loads((SOURCE / "case_identity.json").read_text())
    tangent = unit(np.asarray(identity["straight_tangent_aba"], float))
    e1 = unit(np.asarray(identity["initial_e1_aba"], float))
    e2 = unit(np.asarray(identity["initial_e2_aba"], float))

    # Resolve the 30 Hz cone at 0.1 ms knots. A quintic ramp makes the imposed
    # rotation and its first derivative continuous at start and at 1 ms.
    time = np.unique(np.r_[np.arange(0.0, DURATION, 1e-4), DURATION])
    x = np.clip(time / 0.001, 0.0, 1.0)
    ramp = 10*x**3 - 15*x**4 + 6*x**5
    phase = np.deg2rad(248.0) + 2*np.pi*30.0*time
    cone_axis = (np.cos(np.deg2rad(30.0))*tangent[None, :]
                 + np.sin(np.deg2rad(30.0))*(np.cos(phase)[:, None]*e1[None, :]
                                             + np.sin(phase)[:, None]*e2[None, :]))
    rotations = []
    for scale, target in zip(ramp, cone_axis):
        full = align_rotation(A0, unit(target))
        rotations.append(Rotation.from_rotvec(scale * full.as_rotvec()).as_rotvec())
    rotvec = np.asarray(rotations)

    deck = (SOURCE / f"{SOURCE_JOB}.inp").read_text()
    deck = deck.replace(SOURCE_JOB, JOB)
    deck, count_duration = re.subn(
        r"(\*Dynamic, Explicit, DIRECT\s*\n)1\.0e-7,\s*0\.033333333",
        rf"\g<1>{DT:.1e}, {DURATION:.6f}", deck, count=1)
    socket_hydro = re.compile(
        r"(?ms)^\*Amplitude, name=SOCKET_FX, definition=USER\s*$.*?"
        r"^\*Cload, amplitude=HYDRO_MZ\s*$\nRP_ROBOT, 6, 1\.\s*$")
    replacement = "\n".join([
        "** PRESCRIBED REPLAY: no magnetic Socket loads; translations remain free.",
        "** Reduced-Hydro force and moment loads are retained.",
        "*Amplitude, name=HYDRO_FX, definition=USER",
        "*Amplitude, name=HYDRO_FY, definition=USER",
        "*Amplitude, name=HYDRO_FZ, definition=USER",
        "*Amplitude, name=HYDRO_MX, definition=USER",
        "*Amplitude, name=HYDRO_MY, definition=USER",
        "*Amplitude, name=HYDRO_MZ, definition=USER",
        "*Cload, amplitude=HYDRO_FX", "RP_ROBOT, 1, 1.",
        "*Cload, amplitude=HYDRO_FY", "RP_ROBOT, 2, 1.",
        "*Cload, amplitude=HYDRO_FZ", "RP_ROBOT, 3, 1.",
        "*Cload, amplitude=HYDRO_MX", "RP_ROBOT, 4, 1.",
        "*Cload, amplitude=HYDRO_MY", "RP_ROBOT, 5, 1.",
        "*Cload, amplitude=HYDRO_MZ", "RP_ROBOT, 6, 1.",
        amplitude("PRESCRIBED_UR1", time, rotvec[:, 0]),
        amplitude("PRESCRIBED_UR2", time, rotvec[:, 1]),
        amplitude("PRESCRIBED_UR3", time, rotvec[:, 2]),
        "*Boundary, type=DISPLACEMENT, amplitude=PRESCRIBED_UR1",
        "RP_ROBOT, 4, 4, 1.",
        "*Boundary, type=DISPLACEMENT, amplitude=PRESCRIBED_UR2",
        "RP_ROBOT, 5, 5, 1.",
        "*Boundary, type=DISPLACEMENT, amplitude=PRESCRIBED_UR3",
        "RP_ROBOT, 6, 6, 1.",
    ])
    deck, count_loads = socket_hydro.subn(replacement, deck, count=1)
    if count_duration != 1 or count_loads != 1:
        raise RuntimeError(f"deck replacement failed: duration={count_duration}, loads={count_loads}")
    deck = ("** SINGLE AUTHORIZED CURRENT-GEOMETRY PRESCRIBED-WOBBLE CONTACT AUDIT\n"
            "** ROTATION PRESCRIBED; TRANSLATION FREE; MAGPYLIB LOADS REMOVED\n" + deck)
    deck_path = CASE / f"{JOB}.inp"
    deck_path.write_text(deck)
    shutil.copyfile(HERE / "vuforc_prescribed_hydro.f", CASE / "vuforc_prescribed_hydro.f")
    np.savetxt(CASE / "prescribed_rotation_history.csv",
               np.column_stack((time, rotvec)), delimiter=",",
               header="time_s,UR1_rad,UR2_rad,UR3_rad", comments="")
    replay = {
        "case_id": JOB, "status": "PREPARED", "source_job": SOURCE_JOB,
        "source_commit": "5f0de0ee8c29ee9045f13eddee086e7c61bc19c9",
        "duration_s": DURATION, "direct_dt_s": DT, "translation": "FREE",
        "rotation": "PRESCRIBED_30HZ_30DEG_CONE_WITH_1MS_QUINTIC_RAMP",
        "phase_deg": 248.0, "frequency_Hz": 30.0, "cone_deg": 30.0,
        "magpylib_loads": False, "reduced_hydro": True,
        "mu": 0.03, "zeta": 0.50, "tangent_fraction": 0.0,
        "input_sha256": sha256(deck_path),
        "fortran_sha256": sha256(CASE / "vuforc_prescribed_hydro.f"),
    }
    (CASE / "case_identity.json").write_text(json.dumps(replay, indent=2) + "\n")
    print(json.dumps(replay, indent=2))


if __name__ == "__main__":
    main()
