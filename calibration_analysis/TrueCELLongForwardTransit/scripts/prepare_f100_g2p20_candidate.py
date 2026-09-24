"""Prepare the one authorized fresh-start G=2.20 mT F100 dynamics candidate."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
SOURCE_JOB = "TRUECEL_B0P11_A14P5_F100_FAST"
JOB = "TRUECEL_B0P11_G2P20_A14P5_F100_FAST"
SOURCE = OUT / "case" / SOURCE_JOB
DEST = OUT / "case" / JOB
G_OLD = 2.0
G_NEW = 2.2


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def replace_once(text: str, old: str, new: str, label: str) -> str:
    if text.count(old) != 1:
        raise RuntimeError(f"Expected exactly one {label}; found {text.count(old)}")
    return text.replace(old, new, 1)


def main():
    if DEST.exists():
        raise RuntimeError(f"Candidate already exists; refusing to overwrite: {DEST}")
    identity = json.loads((SOURCE / "case_identity.json").read_text(encoding="utf-8"))
    expected = {
        "B0_mT": 11.0, "gradient_mT": G_OLD, "frequency_Hz": 100.0,
        "rocking_main_amplitude_deg": 14.5, "rocking_cross_amplitude_deg": 2.5,
        "fluid_EOS_c0_mm_s": 100000.0, "Eulerian_dimensions": [44, 20, 20],
        "explicit_stable_time_scale_factor": 0.4,
    }
    for key, value in expected.items():
        if identity.get(key) != value:
            raise RuntimeError(f"Authoritative F100 mismatch: {key}={identity.get(key)!r}, expected {value!r}")

    deck = (SOURCE / f"{SOURCE_JOB}.inp").read_text(encoding="latin1")
    deck = deck.replace(SOURCE_JOB, JOB)
    deck = replace_once(deck,
        "** ONE AUTHORIZED PHYSICS CHANGE FROM F120: frequency 120 -> 100 Hz",
        "** FROZEN F100 BASELINE; ONLY NEW PHYSICS CHANGE: G 2.0 -> 2.2 mT",
        "lineage comment")
    if ", 0.020000000000" not in deck:
        raise RuntimeError("F100 source is not the required 20-ms fresh-start deck")

    source = (SOURCE / "vuamp_precomputed_truecel.f90").read_text(encoding="ascii")
    source = source.replace(SOURCE_JOB, JOB)
    source = replace_once(source, "grad=0.002000d0*grad", "grad=0.002200d0*grad", "G multiplier")
    if "b=0.011d0*b" not in source or "phase_deg=modulo(36000.0d0*t,360.0d0)" not in source:
        raise RuntimeError("B0 or 100-Hz phase law is not frozen")
    source = source.replace("magnetic_increment_f100.csv", "magnetic_increment_g2p20_f100.csv")

    DEST.mkdir(parents=True)
    table = "magnetic_field_gradient_table_B0P11_A14P5.dat"
    shutil.copy2(SOURCE / table, DEST / table)
    shutil.copy2(SOURCE / (table + ".json"), DEST / (table + ".json"))
    inp = DEST / f"{JOB}.inp"
    f90 = DEST / "vuamp_precomputed_truecel.f90"
    inp.write_text(deck, encoding="latin1")
    f90.write_text(source, encoding="ascii")

    new_identity = dict(identity)
    for key in ("wallclock_s", "cpus", "socket_calls", "failure"):
        new_identity.pop(key, None)
    new_identity.update({
        "case_id": JOB, "status": "PREPARED", "classification": "G2P20_F100_TWO_CYCLE_GATE",
        "physical_parent": SOURCE_JOB, "gradient_mT": G_NEW, "gradient_frozen_mT": G_NEW,
        "duration_s": 0.02, "dynamics_run_count": 0,
        "single_physics_change": "G 2.0 -> 2.2 mT",
        "only_new_physics_change": "axial magnetic profile amplitude G 2.0 -> 2.2 mT",
        "selection_basis": "OPEN_LOOP_FIRST_ORDER_G_ESTIMATE: 2.15 mT first predicted <0.005 mm; one 0.05-mT robustness increment",
        "input_sha256": sha(inp), "fortran_sha256": sha(f90),
        "magnetic_table_sha256": sha(DEST / table),
        "frozen_for_candidate": [
            "B0=11mT", "f=100Hz", "A_main=14.5deg", "A_cross=2.5deg",
            "TRUE-CEL 44x20x20 mesh/domain", "c0=100000mm/s",
            "water density/viscosity/EOS", "robot/wall geometry/contact",
            "Explicit scale factor=0.4; no mass scaling", "clean initial state",
        ],
    })
    (DEST / "case_identity.json").write_text(json.dumps(new_identity, indent=2) + "\n", encoding="ascii")
    audit = {
        "case": JOB, "parent": SOURCE_JOB, "fresh_start": True,
        "single_physics_change": "G 2.0 -> 2.2 mT",
        "B0_mT": 11.0, "G_mT": G_NEW, "frequency_Hz": 100.0,
        "A_main_deg": 14.5, "A_cross_deg": 2.5,
        "duration_s": 0.02, "CEL_mesh": [44, 20, 20], "c0_mm_s": 100000.0,
        "explicit_controls": "SCALE FACTOR=0.4; no mass scaling",
        "input_sha256": sha(inp), "fortran_sha256": sha(f90),
        "magnetic_table_sha256": sha(DEST / table), "dynamics_candidates_authorized": 1,
    }
    (OUT / f"{JOB}_Setup_Audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
