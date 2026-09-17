"""Derive isolated 0.5 ms LIVE/TABLE cases from the frozen non-CEL deck."""
from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path

from table_model import OUT, REPO, config


SOURCE = REPO / "calibration_analysis" / "F60LowGBackgroundFlowScreen" / "case" / "S4_HEADFORWARD_F60_G0P10_FLOW"
SOURCE_DECK = SOURCE / "S4_HEADFORWARD_F60_G0P10_FLOW.inp"
CASES = {
    "AB_LIVE_SOCKET_0P5MS": ("LIVE_MAGPYLIB_SOCKET", "vuamp_live_socket.f"),
    "AB_PRECOMPUTED_TABLE_0P5MS_CPU1": ("PRECOMPUTED_TABLE", "vuamp_precomputed_table.f90"),
    "BENCH_PRECOMPUTED_0P5MS_CPU4": ("PRECOMPUTED_TABLE", "vuamp_precomputed_table.f90"),
    "BENCH_PRECOMPUTED_0P5MS_CPU6": ("PRECOMPUTED_TABLE", "vuamp_precomputed_table.f90"),
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_deck(job, cfg):
    text = SOURCE_DECK.read_text(encoding="ascii")
    text = "** PRECOMPUTED MAGNETIC BACKEND 0.5 MS A/B\n** JOB={}\n".format(job) + text
    text = text.replace("2.0e-7, 0.037500000", "2.0e-7, 0.000500000")
    text = text.replace("*Restart, write, number interval=1, time marks=NO", "** Restart output disabled for 0.5 ms benchmark")
    text = text.replace("*Output, field, time interval=5.0e-4, time marks=NO", "*Output, field, time interval=5.0e-5, time marks=NO")
    start = text.find("*Amplitude, name=HYDRO_FX")
    end = text.find("** \n** OUTPUT REQUESTS", start)
    if start < 0 or end < 0:
        start = text.find("*Amplitude, name=HYDRO_FX")
        end = text.find("**\n** OUTPUT REQUESTS", start)
    if start < 0 or end < 0:
        raise RuntimeError("could not isolate ReducedHydro load block")
    text = text[:start] + "** FLUID_MODE=NONE; ReducedHydro amplitudes removed for magnetic A/B.\n" + text[end:]
    if "*Eulerian" in text:
        raise RuntimeError("source deck unexpectedly contains CEL")
    if "0.000500000" not in text or "HYDRO_FX, definition=USER" in text:
        raise RuntimeError("deck rewrite gate failed")
    return text


def main():
    cfg = config()
    cases_root = OUT / "cases"
    cases_root.mkdir(exist_ok=True)
    for job, (backend, user_source) in CASES.items():
        case = cases_root / job
        case.mkdir(exist_ok=True)
        deck = case / (job + ".inp")
        deck.write_text(make_deck(job, cfg), encoding="ascii", newline="\n")
        shutil.copy2(OUT / user_source, case / user_source)
        shutil.copy2(SOURCE / "straight_control_centerline.dxf", case / "straight_control_centerline.dxf")
        if backend == "PRECOMPUTED_TABLE":
            shutil.copy2(OUT / cfg["table"]["path"], case / cfg["table"]["path"])
        identity = {
            "case_id": job, "status": "PREPARED", "MAGNETIC_BACKEND": backend,
            "FLUID_MODE": "NONE", "duration_s": cfg["dynamic_ab_duration_s"],
            "direct_dt_s": cfg["direct_dt_s"], "frequency_Hz": cfg["frequency_Hz"],
            "B0_T": cfg["B0_T"], "gradient_G_T": cfg["gradient_G_T"],
            "source_case": cfg["source_case"], "source_input_sha256": digest(SOURCE_DECK),
            "input_sha256": digest(deck), "user_source": user_source,
            "user_source_sha256": digest(case / user_source), "dynamics_run_count": 0,
        }
        (case / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n", encoding="ascii")
        print("PREPARED", job, backend)


if __name__ == "__main__":
    main()
