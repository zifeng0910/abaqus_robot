"""Create field-only and old-versus-continuous temporal validation datasets."""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from continuous_segment_reference import ContinuousSegmentReference, unit
from regress_production_local_frame import DXF, TRANSFORM, TRAJECTORY, RP0, VENDORED, build, load_module


HERE = Path(__file__).resolve().parent
OUT = HERE.parent


def angle_step(vectors):
    left = vectors[:-1] / np.linalg.norm(vectors[:-1], axis=1)[:, None]
    right = vectors[1:] / np.linalg.norm(vectors[1:], axis=1)[:, None]
    return np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", left, right), -1.0, 1.0)))


def field_from_frame(reference, time_s, tangent, e1, e2):
    phase = reference.phase0 + reference.sense * 2.0 * math.pi * reference.frequency * time_s
    return reference.magnitude * (math.cos(reference.cone) * tangent +
        math.sin(reference.cone) * (math.cos(phase) * e1 + reference.sense * math.sin(phase) * e2))


def main():
    production = load_module("continuity_production", VENDORED)
    model, wrapped, _ = build(production, field_mode="ROBOT_LOCAL_TANGENT")
    reference = ContinuousSegmentReference(DXF, TRANSFORM, model.robot_axis_global,
                                            polarity=model.robot_polarity)
    one_cycle = []
    reference.reset()
    for time_s in np.linspace(0.0, 1.0 / 30.0, 181):
        oracle = reference.evaluate(time_s, RP0, np.zeros(3))
        production_result = wrapped.evaluate(time_s, RP0, np.zeros(3))
        b = oracle["B_aba_T"]
        one_cycle.append((time_s, oracle["s_mm"], *reference.vector_to_aba(oracle["center_mag_mm"] - np.zeros(3)),
                          *reference.vector_to_aba(oracle["tangent_mag"]),
                          *reference.vector_to_aba(oracle["e1_mag"]), *reference.vector_to_aba(oracle["e2_mag"]),
                          *b, *production_result["B_aba_vec_T"]))
    columns = ["time_s", "s_mm", "center_vx", "center_vy", "center_vz",
               "t_x", "t_y", "t_z", "e1_x", "e1_y", "e1_z", "e2_x", "e2_y", "e2_z",
               "Bx_T", "By_T", "Bz_T", "prod_Bx_T", "prod_By_T", "prod_Bz_T"]
    field_table = pd.DataFrame(one_cycle, columns=columns)
    center_mag = reference.evaluate(0.0, RP0, np.zeros(3))["center_mag_mm"]
    center_aba = reference.origin_aba + reference.R.T.dot(center_mag)
    field_table[["center_vx", "center_vy", "center_vz"]] = center_aba
    field_table.to_csv(OUT / "field_one_cycle.csv", index=False)

    history = np.load(TRAJECTORY)
    indices = np.unique(np.r_[np.arange(0, len(history["U1"]), 100), len(history["U1"]) - 1])
    time = history["U1"][indices, 0]
    position = RP0 + np.column_stack([history[key][indices, 1] for key in ("U1", "U2", "U3")])
    reference.reset()
    rows = []
    for time_s, position_aba in zip(time, position):
        point_mag = reference.position_to_mag(position_aba)
        new_s, _ = reference.project_mag(point_mag)
        new_t, new_e1, new_e2 = reference.frame_at(new_s)
        new_b = field_from_frame(reference, time_s, new_t, new_e1, new_e2)
        old_index = int(np.argmin(np.linalg.norm(reference.curve - point_mag, axis=1)))
        old_s = float(reference.arc[old_index])
        old_t, old_e1, old_e2 = reference.frame_at(old_s)
        old_b = field_from_frame(reference, time_s, old_t, old_e1, old_e2)
        rows.append((time_s, new_s, old_s, *new_t, *old_t, *new_b, *old_b))
    table = pd.DataFrame(rows, columns=["time_s", "new_s_mm", "old_s_mm",
        "new_tx", "new_ty", "new_tz", "old_tx", "old_ty", "old_tz",
        "new_Bx", "new_By", "new_Bz", "old_Bx", "old_By", "old_Bz"])
    new_b = table[["new_Bx", "new_By", "new_Bz"]].to_numpy()
    old_b = table[["old_Bx", "old_By", "old_Bz"]].to_numpy()
    new_t = table[["new_tx", "new_ty", "new_tz"]].to_numpy()
    old_t = table[["old_tx", "old_ty", "old_tz"]].to_numpy()
    table["new_B_step_deg"] = np.r_[0.0, angle_step(new_b)]
    table["old_B_step_deg"] = np.r_[0.0, angle_step(old_b)]
    table["new_t_step_deg"] = np.r_[0.0, angle_step(new_t)]
    table["old_t_step_deg"] = np.r_[0.0, angle_step(old_t)]
    table.to_csv(OUT / "temporal_continuity.csv", index=False)
    metrics = {
        "samples": int(len(table)), "sample_interval_us": float(np.median(np.diff(time)) * 1e6),
        "new_B_step_max_deg": float(table.new_B_step_deg.max()),
        "old_B_step_max_deg": float(table.old_B_step_deg.max()),
        "new_t_step_max_deg": float(table.new_t_step_deg.max()),
        "old_t_step_max_deg": float(table.old_t_step_deg.max()),
        "new_s_step_max_mm": float(np.max(np.abs(np.diff(table.new_s_mm)))),
        "old_s_step_max_mm": float(np.max(np.abs(np.diff(table.old_s_mm)))),
        "new_s_reversals": int(np.sum(np.diff(table.new_s_mm) < 0.0)),
        "initial_phase_component_error_T": float(np.max(np.abs(field_table.iloc[0][["Bx_T", "By_T", "Bz_T"]].to_numpy(float) - field_table.iloc[0][["prod_Bx_T", "prod_By_T", "prod_Bz_T"]].to_numpy(float)))),
    }
    (OUT / "temporal_continuity_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
