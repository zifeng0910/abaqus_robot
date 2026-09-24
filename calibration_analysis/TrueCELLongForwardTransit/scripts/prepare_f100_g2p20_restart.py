"""Prepare the authorized 20->30 ms restart for the sole G=2.20 candidate."""
from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "case" / "TRUECEL_B0P11_G2P20_A14P5_F100_FAST"
OLD = "TRUECEL_B0P11_G2P20_A14P5_F100_FAST"
NEW = "TRUECEL_B0P11_G2P20_A14P5_F100_CYCLE3_RESTART"
END_S = 0.23549538141218107


def main() -> None:
    initial = (CASE / f"{OLD}.inp").read_text(encoding="latin1")
    step = initial[initial.index("*Step, name=Step_Drive"):]
    step = re.sub(r"\*Step, name=Step_Drive", "*Step, name=Step_Cycle3_Restart", step, count=1)
    step = step.replace(", 0.020000000000", ", 0.010000000000", 1)
    step = re.sub(r"(?:\*Amplitude, name=SOCKET_(?:F[XYZ]|M[XYZ]), definition=USER\n){6}", "", step, count=1)
    step = step.replace("*Restart, write, number interval=4, time marks=YES",
                        "*Restart, write, number interval=2, time marks=YES", 1)
    deck = (
        "*HEADING\n"
        "TRUE-CEL G2P20 F100 third-cycle restart; same saved physical state.\n"
        "*RESTART, READ, STEP=1, INTERVAL=4, END STEP\n"
        + step
    )
    (CASE / f"{NEW}.inp").write_text(deck, encoding="latin1")

    src = (CASE / "vuamp_precomputed_truecel.f90").read_text(encoding="latin1")
    src = src.replace("magnetic_increment_g2p20_f100.csv", "magnetic_increment_g2p20_f100_cycle3.csv")
    src = src.replace("online_cycle_gate.csv", "online_cycle_gate_cycle3.csv")
    src = src.replace("online_motion.csv", "online_motion_cycle3.csv")
    src = src.replace("f100_event.txt", "f100_event_cycle3.txt")
    old_init = (
        "monitor_initialized=.true.; monitor_last_t=t; monitor_max_s=0.0d0\n"
        "      monitor_cycle_start_s=0.0d0; monitor_next_boundary=1.0d-2"
    )
    new_init = (
        f"monitor_initialized=.true.; monitor_last_t=t; monitor_max_s={END_S:.17g}d0\n"
        f"      monitor_last_s={END_S:.17g}d0; monitor_cycle_start_s={END_S:.17g}d0\n"
        "      monitor_cycle=3; monitor_next_boundary=3.0d-2"
    )
    if old_init not in src:
        raise RuntimeError("restart monitor initialization token not found")
    src = src.replace(old_init, new_init, 1)
    (CASE / "vuamp_precomputed_truecel_cycle3.f90").write_text(src, encoding="latin1")

    required = [
        "*RESTART, READ, STEP=1, INTERVAL=4, END STEP",
        ", 0.010000000000",
        "grad=0.002200d0*grad",
        "phase_deg=modulo(36000.0d0*t,360.0d0)",
        "monitor_cycle=3; monitor_next_boundary=3.0d-2",
    ]
    combined = deck + src
    missing = [token for token in required if token not in combined]
    if missing:
        raise RuntimeError(f"restart verification failed: {missing}")
    print(f"PREPARED {NEW}; oldjob={OLD}; continuation=10 ms; G=2.20 unchanged")


if __name__ == "__main__":
    main()
