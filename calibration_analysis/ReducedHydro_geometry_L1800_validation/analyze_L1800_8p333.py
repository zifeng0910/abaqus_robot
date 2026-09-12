"""Analyze the single L=1.800 mm geometry run with the validated RH pipeline."""
from pathlib import Path
import json
import sys

import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BASE = REPO / "calibration_analysis" / "ReducedHydro_zeta050_8p333_validation"
sys.path.insert(0, str(BASE))
import analyze_zeta050_8p333 as audit


JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083"
SOURCE_PREFIX = "zeta050_8p333"
OUTPUT_PREFIX = "L1800_8p333"


def main():
    audit.HERE = HERE
    audit.JOB = JOB
    runtime, _ = audit.analyze_new_case()

    identity = json.loads((HERE / "L1800_input_identity.json").read_text())
    runtime.update({
        "geometry_change": "PCA axial span only",
        "base_pca_axial_span_mm": identity["base_pca_axial_span_mm"],
        "candidate_pca_axial_span_mm": identity["candidate_pca_axial_span_mm"],
        "candidate_pca_transverse_diameter_mm": identity["candidate_pca_transverse_diameter_mm"],
        "candidate_mass_mg": identity["candidate_mass_mg"],
        "candidate_magnetic_moment_Am2": identity["candidate_magnetic_moment_Am2"],
        "baseline_mesh_identity_correction": identity["baseline_mesh_identity_correction"],
    })

    for source in HERE.glob(SOURCE_PREFIX + "_*"):
        source.replace(HERE / source.name.replace(SOURCE_PREFIX, OUTPUT_PREFIX, 1))

    summary_path = HERE / f"{OUTPUT_PREFIX}_summary.json"
    summary_path.write_text(json.dumps(runtime, indent=2))
    public = {key: value for key, value in runtime.items()
              if key not in ("head_node_ids", "tail_node_ids")}
    pd.DataFrame([public]).to_csv(HERE / f"{OUTPUT_PREFIX}_runtime_identity.csv", index=False)
    print(json.dumps(runtime, indent=2))


if __name__ == "__main__":
    main()
