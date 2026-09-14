from pathlib import Path
import json
import pandas as pd

HERE=Path(__file__).resolve().parent; SCREEN=HERE.parent

def main():
    manifest=pd.read_csv(SCREEN/"metrics/case_manifest.csv"); metrics=pd.read_csv(SCREEN/"metrics/all_cases_metrics.csv")
    assert len(manifest)==32 and manifest.case_id.is_unique
    assert set(manifest.case_id)==set(metrics.case_id)
    assert (manifest.element_count.between(5000,19999)).all()
    assert (manifest.direct_dt_s==1e-7).all()
    assert all(abs(row.duration_s*row.frequency_Hz-.2)<1e-9 for row in manifest.itertuples())
    required=("robot_phase_advance_deg","delta_s_mm","rotation_translation_index","impact_fraction","sliding_fraction","stuck_fraction")
    assert all(c in metrics for c in required) and metrics[list(required)].notna().all().all()
    assert len(list((SCREEN/"gifs").glob("*_motion.gif")))==32
    assert len(list((SCREEN/"family_mosaics").glob("*.gif")))==9
    review=pd.read_csv(SCREEN/"human_review_template.csv",keep_default_na=False)
    assert len(review)==32 and (review.drop(columns="case_id")=="").all().all()
    report=(SCREEN/"Coarse_Motion_Mode_Screen_Report.md").read_text()
    assert report.rstrip().endswith("AWAITING HUMAN GIF REVIEW")
    print("PASS: 32 unique cases, metrics, GIFs, mosaics, blank review template, terminal report state")

if __name__=="__main__":main()

