"""Prepare gated one-cycle magnetic-frame diagnostic decks."""
import argparse
from pathlib import Path
import json
import re
import shutil


HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent
REPO = AUDIT.parents[1]
SOURCE = REPO / "calibration_analysis/CoarseMotionModeScreen/cases/GEO_240"
SOURCE_JOB = "SCR_GEO240_Wobble"
CASE_DEFINITIONS = [
    {"case_id": "FRAME_CURRENT", "job_name": "MDA_FRAME_CURRENT", "field_mode": "CURRENT_PRODUCTION", "cone_deg": 30.0, "cone_axis_bias_deg": 40.0},
    {"case_id": "FRAME_LOCAL30", "job_name": "MDA_FRAME_LOCAL30", "field_mode": "LOCAL_TANGENT_CONE", "cone_deg": 30.0, "cone_axis_bias_deg": 0.0},
    {"case_id": "FRAME_LOCAL40", "job_name": "MDA_FRAME_LOCAL40", "field_mode": "LOCAL_TANGENT_CONE", "cone_deg": 40.0, "cone_axis_bias_deg": 0.0},
    {"case_id": "FRAME_LOCAL30_RAMP", "job_name": "MDA_FRAME_LOCAL30_RAMP", "field_mode": "LOCAL_TANGENT_CONE", "cone_deg": 30.0, "cone_axis_bias_deg": 0.0, "ramp_time_s": 0.005},
]
CASES = {item["case_id"]: item for item in CASE_DEFINITIONS}
MANDATORY = ("FRAME_CURRENT", "FRAME_LOCAL30", "FRAME_LOCAL40")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-id", action="append", choices=tuple(CASES),
                        help="Prepare only the explicitly gated case; defaults to the three mandatory cases")
    args = parser.parse_args()
    source_deck = (SOURCE / (SOURCE_JOB + ".inp")).read_text()
    for case_id in args.case_id or MANDATORY:
        item = CASES[case_id]
        folder = AUDIT / "cases" / item["case_id"]
        (folder / "private").mkdir(parents=True, exist_ok=True)
        deck = source_deck.replace(SOURCE_JOB, item["job_name"])
        deck = re.sub(r"(\*Dynamic, Explicit, DIRECT\n)1\.0e-7,\s*[^\n]+", r"\g<1>1.0e-7, 0.033333333", deck, count=1)
        header = ("** MAGNETIC DRIVE FRAME AUDIT - SCREENING ONLY\n"
                  "** L2.4 D0.815 coarse rigid carrier; one 30 Hz cycle; G=0\n"
                  "** FIELD_MODE=%s; cone=%.1f deg\n" % (item["field_mode"], item["cone_deg"]))
        (folder / (item["job_name"] + ".inp")).write_text(header + deck)
        shutil.copyfile(SOURCE / "vuforc_screen.f", folder / "vuforc_audit.f")
        identity = dict(item, ramp_time_s=item.get("ramp_time_s", 0.001),
                        length_mm=2.4, diameter_mm=0.815, B0_mT=10.0, frequency_Hz=30.0,
                        gradient_mT=0.0, mu=0.03, zeta=0.50, duration_s=0.033333333,
                        direct_dt_s=1.0e-7, field_interval_s=1.0e-4, mesh_size_mm=0.1,
                        node_count=1588, element_count=7302, magnetic_moment_Am2=0.0010876227522174417,
                        cparallel_scale=1.0, kwobble_scale=1.0, status="PREPARED")
        (folder / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n")
        print("PREPARED", item["case_id"], item["job_name"])


if __name__ == "__main__":
    main()
