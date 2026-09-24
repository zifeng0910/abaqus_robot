"""Prepare the single uninterrupted 0-30 ms F100/G2.20 validation run."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT = "TRUECEL_B0P11_G2P20_A14P5_F100_FAST"
JOB = "TRUECEL_B0P11_G2P20_A14P5_F100_CLEAN30"
SOURCE = ROOT / "case" / PARENT
DEST = ROOT / "case" / JOB


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def once(data: str, old: str, new: str) -> str:
    count = data.count(old)
    if count != 1:
        raise RuntimeError(f"Expected one {old!r}, found {count}")
    return data.replace(old, new, 1)


def main() -> None:
    if DEST.exists():
        raise RuntimeError(f"CLEAN30 already exists: {DEST}")
    ident = json.loads((SOURCE / "case_identity.json").read_text(encoding="utf-8"))
    required = {
        "status": "SOLVED", "duration_s": 0.02, "B0_mT": 11.0,
        "gradient_mT": 2.2, "frequency_Hz": 100.0,
        "rocking_main_amplitude_deg": 14.5,
        "rocking_cross_amplitude_deg": 2.5,
        "fluid_EOS_c0_mm_s": 100000.0,
        "Eulerian_dimensions": [44, 20, 20],
        "explicit_stable_time_scale_factor": 0.4,
        "cpus": 1,
    }
    for key, expected in required.items():
        if ident.get(key) != expected:
            raise RuntimeError(f"Parent {key}: {ident.get(key)!r} != {expected!r}")

    old_inp = SOURCE / f"{PARENT}.inp"
    old_sub = SOURCE / "vuamp_precomputed_truecel.f90"
    table = "magnetic_field_gradient_table_B0P11_A14P5.dat"
    if sha(SOURCE / table) != ident["magnetic_table_sha256"]:
        raise RuntimeError("Parent magnetic table hash changed")

    deck = old_inp.read_text(encoding="latin1")
    deck = once(deck, ", 0.020000000000", ", 0.030000000000")
    deck = deck.replace(PARENT, JOB)
    if deck.count("*Step, name=Step_Drive") != 1 or "*RESTART, READ" in deck.upper():
        raise RuntimeError("CLEAN30 is not a one-step fresh-start deck")
    if "*Fixed Mass Scaling" in deck or "*Variable Mass Scaling" in deck:
        raise RuntimeError("Unexpected mass scaling")

    sub = old_sub.read_text(encoding="latin1")
    sub = sub.replace(PARENT, JOB)
    sub = once(sub, "if (reason/='NONE' .and. .not.terminate_requested) then",
               "if (reason=='NUMERICALLY_INVALID_STABLE_DT_COLLAPSE' .and. &\n"
               "        .not.terminate_requested) then")
    if "t=dble(time(iTotalTime))" not in sub or \
            "phase_deg=modulo(36000.0d0*t,360.0d0)" not in sub or \
            "b=0.011d0*b; grad=0.002200d0*grad" not in sub:
        raise RuntimeError("Frozen magnetic implementation changed")

    DEST.mkdir(parents=True)
    inp = DEST / f"{JOB}.inp"
    f90 = DEST / "vuamp_precomputed_truecel.f90"
    inp.write_text(deck, encoding="latin1")
    f90.write_text(sub, encoding="latin1")
    shutil.copy2(SOURCE / table, DEST / table)
    shutil.copy2(SOURCE / (table + ".json"), DEST / (table + ".json"))
    new_ident = dict(ident)
    for key in ("wallclock_s", "failure"):
        new_ident.pop(key, None)
    new_ident.update({"case_id": JOB, "status": "PREPARED",
                      "classification": "CLEAN30_UNINTERRUPTED_C3_VALIDATION",
                      "duration_s": 0.03, "dynamics_run_count": 0,
                      "parent_clean20": PARENT,
                      "physical_early_stop": "DISABLED; numerical timestep collapse retained",
                      "input_sha256": sha(inp), "fortran_sha256": sha(f90)})
    (DEST / "case_identity.json").write_text(
        json.dumps(new_ident, indent=2) + "\n", encoding="ascii")
    audit = {"job": JOB, "parent": PARENT, "fresh_start": True,
             "single_step_duration_s": .03, "restart_read": False,
             "physics_change": "NONE",
             "output_change": "same output keywords; number interval=4 restart writes shift from 5/10/15/20 ms to 7.5/15/22.5/30 ms with step duration",
             "vuamp_change": "job paths and physical early-stop gate only",
             "parent_input_sha256": sha(old_inp), "parent_vuamp_sha256": sha(old_sub),
             "input_sha256": sha(inp), "vuamp_sha256": sha(f90),
             "magnetic_table_sha256": sha(DEST / table)}
    (ROOT / f"{JOB}_SETUP.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="ascii")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
