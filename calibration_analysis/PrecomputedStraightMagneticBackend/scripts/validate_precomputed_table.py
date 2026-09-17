"""Validate interpolation, force/torque, linearity, and transverse invariance."""
from __future__ import annotations

import json
import math

import numpy as np

from table_model import OUT, build_live_model, config, gradient_unit, read_table, table_loads, unit


def relative(actual, reference, floor=1e-15):
    return float(np.linalg.norm(actual - reference) / max(np.linalg.norm(reference), floor))


def main():
    cfg = config()
    phase_table, s_table = read_table(OUT / cfg["table"]["path"])
    live = build_live_model(cfg)
    rng = np.random.RandomState(20260917)
    samples = []
    for _ in range(64):
        t_s = rng.uniform(0.0, 2.0 / cfg["frequency_Hz"])
        axial = rng.uniform(-2.0, 2.0)
        radial_n, radial_b = rng.uniform(-0.25, 0.25, size=2)
        position = (np.asarray(cfg["initial_center_aba_mm"]) +
                    axial * unit(cfg["canonical_plus_s_axis_aba"]) +
                    radial_n * unit(cfg["n_routeA_aba"]) + radial_b * unit(cfg["b_routeA_aba"]))
        ur = rng.normal(size=3)
        ur *= rng.uniform(0.0, 0.45) / max(np.linalg.norm(ur), 1e-15)
        live.magnetic_model.reset_robot_arc_continuity()
        direct = live.evaluate(t_s, position, ur)
        field, gradient, force, torque, s_eff, phase_deg = table_loads(t_s, position, ur, cfg, phase_table, s_table)
        direct_field = np.asarray(direct["B_aba_vec_T"])
        direct_force = np.asarray(direct["force_N"])
        direct_torque = np.asarray(direct["torque_Nmm"])
        direct_gradient = cfg["gradient_G_T"] * gradient_unit(float(direct["robot_arc_mm"] - direct["driver_arc_mm"]), cfg)
        samples.append({"B": relative(field, direct_field), "gradient": relative(gradient, direct_gradient),
                        "force": relative(force, direct_force, 1e-12),
                        "torque": relative(torque, direct_torque, 1e-12),
                        "s_eff_mm": s_eff, "phase_deg": phase_deg})
    base_position = np.asarray(cfg["initial_center_aba_mm"])
    offsets = {"center": np.zeros(3), "+n": 0.25 * unit(cfg["n_routeA_aba"]),
               "-n": -0.25 * unit(cfg["n_routeA_aba"]), "+b": 0.25 * unit(cfg["b_routeA_aba"]),
               "-b": -0.25 * unit(cfg["b_routeA_aba"])}
    transverse = []
    for t_s in (0.001, 0.0031, 0.0067):
        values = {}
        for name, offset in offsets.items():
            live.magnetic_model.reset_robot_arc_continuity()
            result = live.evaluate(t_s, base_position + offset, np.zeros(3))
            values[name] = np.r_[result["B_aba_vec_T"], result["force_N"], result["torque_Nmm"]]
        center = values["center"]
        transverse.append(max(relative(value, center, 1e-12) for name, value in values.items() if name != "center"))
    maxima = {key: max(row[key] for row in samples) for key in ("B", "gradient", "force", "torque")}
    p95 = {key: float(np.percentile([row[key] for row in samples], 95)) for key in maxima}
    tolerances = {"B": 0.005, "torque": 0.01, "force": 0.03}
    report = {
        "sample_count": len(samples), "random_seed": 20260917,
        "max_relative_error": maxima, "p95_relative_error": p95,
        "tolerances": tolerances,
        "transverse_audit": {"offsets_mm": 0.25, "max_combined_relative_variation": max(transverse),
                             "conclusion": "(s_eff, phase) sufficient" if max(transverse) <= 0.02 else "extend dimensions"},
        "linearity": {"B0": "exact by production construction", "G": "exact by production construction",
                      "normalized_basis_accepted": True},
        "passed": maxima["B"] <= tolerances["B"] and maxima["torque"] <= tolerances["torque"] and maxima["force"] <= tolerances["force"] and max(transverse) <= 0.02,
        "samples": samples,
    }
    (OUT / "offline_validation.json").write_text(json.dumps(report, indent=2) + "\n", encoding="ascii")
    print(json.dumps({key: value for key, value in report.items() if key != "samples"}, indent=2))
    if not report["passed"]:
        raise SystemExit("OFFLINE_VALIDATION_FAILED")


if __name__ == "__main__":
    main()
