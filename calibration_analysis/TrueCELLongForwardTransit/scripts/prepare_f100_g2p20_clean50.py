"""Prepare one uninterrupted five-cycle TRUE-CEL run with 5-ms writes."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = "TRUECEL_B0P11_G2P20_A14P5_F100_FAST"
SOURCE_SUB = "TRUECEL_B0P11_G2P20_A14P5_F100_CLEAN30"
JOB = "TRUECEL_B0P11_G2P20_A14P5_F100_CLEAN50"
SOURCE = ROOT / "case" / BASE
DEST = ROOT / "case" / JOB
TABLE = "magnetic_field_gradient_table_B0P11_A14P5.dat"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def replace_once(data: bytes, old: bytes, new: bytes) -> bytes:
    if data.count(old) != 1:
        raise RuntimeError(f"Expected exactly one {old!r}")
    return data.replace(old, new, 1)


def main() -> None:
    if DEST.exists():
        raise RuntimeError(f"CLEAN50 case already exists: {DEST}")
    identity = json.loads((SOURCE / "case_identity.json").read_text(encoding="utf-8"))
    expected = {"duration_s": .02, "B0_mT": 11.0, "gradient_mT": 2.2,
                "frequency_Hz": 100.0, "rocking_main_amplitude_deg": 14.5,
                "rocking_cross_amplitude_deg": 2.5,
                "fluid_EOS_c0_mm_s": 100000.0,
                "Eulerian_dimensions": [44, 20, 20],
                "explicit_stable_time_scale_factor": .4, "cpus": 1}
    for key, value in expected.items():
        if identity.get(key) != value:
            raise RuntimeError(f"Baseline {key} differs: {identity.get(key)!r}")
    if sha(SOURCE / TABLE) != identity["magnetic_table_sha256"]:
        raise RuntimeError("Magnetic table hash differs from authoritative baseline")

    original = (SOURCE / f"{BASE}.inp").read_bytes()
    deck = original.replace(BASE.encode(), JOB.encode())
    deck = replace_once(deck, b", 0.020000000000", b", 0.050000000000")
    deck = replace_once(deck,
                        b"*Restart, write, number interval=4, time marks=YES",
                        b"*Restart, write, number interval=10, time marks=YES")
    if deck.count(b"*Step, name=Step_Drive") != 1 or b"*RESTART, READ" in deck.upper():
        raise RuntimeError("Expected one fresh-start dynamics step")
    if b"*Fixed Mass Scaling" in deck or b"*Variable Mass Scaling" in deck:
        raise RuntimeError("Unexpected mass scaling")
    restored = deck.replace(JOB.encode(), BASE.encode())
    restored = replace_once(restored, b", 0.050000000000", b", 0.020000000000")
    restored = replace_once(restored,
                            b"*Restart, write, number interval=10, time marks=YES",
                            b"*Restart, write, number interval=4, time marks=YES")
    if restored != original:
        raise RuntimeError("Other submitted INP bytes changed")

    clean30_sub = (ROOT / "case" / SOURCE_SUB / "vuamp_precomputed_truecel.f90").read_bytes()
    if clean30_sub.count(SOURCE_SUB.encode()) != 6:
        raise RuntimeError("Unexpected VUAMP path count")
    sub = clean30_sub.replace(SOURCE_SUB.encode(), JOB.encode())
    if b"reason=='NUMERICALLY_INVALID_STABLE_DT_COLLAPSE'" not in sub or \
            b"t=dble(time(iTotalTime))" not in sub or \
            b"b=0.011d0*b; grad=0.002200d0*grad" not in sub:
        raise RuntimeError("Required frozen magnetic/numerical-stop logic missing")

    DEST.mkdir(parents=True)
    inp = DEST / f"{JOB}.inp"
    f90 = DEST / "vuamp_precomputed_truecel.f90"
    inp.write_bytes(deck)
    f90.write_bytes(sub)
    shutil.copy2(SOURCE / TABLE, DEST / TABLE)
    shutil.copy2(SOURCE / (TABLE + ".json"), DEST / (TABLE + ".json"))
    new_identity = dict(identity)
    for key in ("wallclock_s", "failure"):
        new_identity.pop(key, None)
    new_identity.update({
        "case_id": JOB, "status": "PREPARED", "dynamics_run_count": 0,
        "classification": "CLEAN50_PENDING", "duration_s": .05,
        "parent_clean20": BASE, "restart_read": False,
        "restart_write_intervals": 10, "restart_period_s": .005,
        "restart_write_times_ms": list(range(5, 51, 5)),
        "physical_early_stop": "DISABLED; numerical timestep collapse retained",
        "input_sha256": sha(inp), "fortran_sha256": sha(f90),
    })
    (DEST / "case_identity.json").write_text(
        json.dumps(new_identity, indent=2) + "\n", encoding="ascii")
    audit = {
        "job": JOB, "baseline": BASE, "fresh_start": True,
        "single_explicit_step_s": .05, "restart_read": False,
        "restart_write_strategy": "number interval=10, time marks=YES; every 5 ms; never read",
        "restart_write_times_ms": list(range(5, 51, 5)),
        "input_differences": ["job name in two comments", "step period .020 -> .050 s",
                              "restart number interval 4 -> 10"],
        "input_normalized_to_baseline_after_these_changes": True,
        "vuamp_source": SOURCE_SUB,
        "vuamp_change": "job paths only",
        "physical_early_stop": "disabled; numerical collapse stop retained",
        "baseline_input_sha256": sha(SOURCE / f"{BASE}.inp"),
        "input_sha256": sha(inp),
        "vuamp_source_sha256": sha(ROOT / "case" / SOURCE_SUB / "vuamp_precomputed_truecel.f90"),
        "vuamp_sha256": sha(f90),
        "magnetic_table_sha256": sha(DEST / TABLE),
    }
    (ROOT / f"{JOB}_SETUP.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="ascii")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
