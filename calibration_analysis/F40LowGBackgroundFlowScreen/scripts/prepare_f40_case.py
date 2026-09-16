"""Prepare and gate the single F40 background-flow screening case."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
CASE = OUT / "case" / "S4_HEADFORWARD_F40_G0P05_FLOW"
JOB = "S4_HEADFORWARD_F40_G0P05_FLOW"
SOURCE_ID = REPO / "calibration_analysis" / "S4HeadForwardLowGScreen" / "cases" / "S4_HEADFORWARD_G0" / "case_identity.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unit(value):
    value = np.asarray(value, dtype=float)
    return value / np.linalg.norm(value)


def sign_gates(source):
    c = unit(source["canonical_plus_s_axis_aba"])
    a = unit(source["reduced_hydro_body_axis_aba"])
    flow = 10.0 * c
    cpar, cperp = 4.0e-9, 1.2e-8
    def hydro(v):
        rel = np.asarray(v, dtype=float) - flow
        parallel = float(np.dot(rel, a))
        transverse = rel - parallel * a
        return -cpar * parallel * a - cperp * transverse
    states = {
        "A_stationary": hydro([0.0, 0.0, 0.0]),
        "B_with_flow": hydro(flow),
        "C_faster_than_flow": hydro(15.0 * c),
    }
    hydro_s = {name: float(np.dot(force, c)) for name, force in states.items()}
    gates = {
        "A_stationary_pushes_plus_s": hydro_s["A_stationary"] > 0.0,
        "B_with_flow_is_zero": abs(hydro_s["B_with_flow"]) < 1.0e-20,
        "C_faster_drag_opposes_motion": hydro_s["C_faster_than_flow"] < 0.0,
    }
    driver_arc = 18.899960626
    robot_arc = 53.799881868759286
    delta = robot_arc - driver_arc
    d_b_ds = -(delta / 45.0**2) * math.exp(-0.5 * (delta / 45.0) ** 2) * 0.05e-3 * 1000.0
    gradient_force_s = -float(source["robot_moment_Am2"]) * d_b_ds
    gates["F_gradient_dot_c_positive"] = gradient_force_s > 0.0
    gates["all_passed"] = all(gates.values())
    return {
        "flow_axis": c.tolist(), "flow_velocity_mm_s": 10.0,
        "hydro_force_s_N": hydro_s, "gates": gates,
        "gradient_gate": {"driver_arc_mm": driver_arc, "robot_arc_mm": robot_arc,
                           "delta_s_mm": delta, "dB_ds_T_per_m": d_b_ds,
                           "F_gradient_dot_c_N": gradient_force_s,
                           "passed": gates["F_gradient_dot_c_positive"]},
    }


def main():
    source = json.loads(SOURCE_ID.read_text())
    deck = CASE / f"{JOB}.inp"
    fortran = CASE / "vuforc_background_flow.f"
    centerline = CASE / "straight_control_centerline.dxf"
    if not all(path.exists() for path in (deck, fortran, centerline)):
        raise RuntimeError("prepared case files are incomplete")
    deck_text = deck.read_text()
    if "*Eulerian" in deck_text or "*Eulerian Section" in deck_text:
        raise RuntimeError("unexpected CEL definition")
    if "2.0e-7, 0.037500000" not in deck_text:
        raise RuntimeError("F40 deck does not use the requested dt/duration")
    checks = sign_gates(source)
    if not checks["gates"]["all_passed"]:
        raise RuntimeError("sign sanity gate failed: {}".format(checks["gates"]))
    identity = {
        "case_id": JOB, "status": "PREPARED", "classification": "LOW_PRECISION_SCREENING",
        "source_case": source["case_id"], "source_commit": "8ccddd8",
        "field_frame_mode": "ROBOT_LOCAL_ELLIPTIC_ROCKING", "duration_s": 0.0375,
        "direct_dt_s": 2.0e-7, "frequency_Hz": 40.0,
        "rocking_main_amplitude_deg": source["rocking_main_amplitude_deg"],
        "rocking_cross_amplitude_deg": source["rocking_cross_amplitude_deg"],
        "B0_mT": source["B0_mT"], "gradient_mT": 0.05, "gradient_length_mm": 45.0,
        "gradient_profile": "legacy", "forward_definition": source["forward_definition"],
        "flow_model": "prescribed background flow + ReducedHydro relative-velocity drag",
        "flow_status": "DIAGNOSTIC_SCREENING_FLOW_ONLY",
        "U_flow_mm_s": 10.0, "flow_axis_aba": checks["flow_axis"],
        "flow_source": "No authoritative straight-tube flow value found in repository; prescribed by task.",
        "mu": source["mu"], "zeta": source["zeta"],
        "Cparallel_Ns_per_mm": 4.0e-9, "Cperp_Ns_per_mm": 1.2e-8,
        "Kspin_Nmm_s": 1.0e-9, "Kwobble_Nmm_s": 3.0e-9,
        "robot_length_mm": source["robot_length_mm"], "robot_diameter_mm": source["robot_diameter_mm"],
        "robot_moment_Am2": source["robot_moment_Am2"], "robot_mass_mg": source["robot_mass_mg"],
        "initial_center_aba_mm": source["initial_center_aba_mm"],
        "canonical_plus_s_axis_aba": source["canonical_plus_s_axis_aba"],
        "head_tail_axis_aba": source["head_tail_axis_aba"],
        "initial_magnetic_moment_axis_aba": source["initial_magnetic_moment_axis_aba"],
        "reduced_hydro_body_axis_aba": source["reduced_hydro_body_axis_aba"],
        "n_routeA_aba": source["n_routeA_aba"], "b_routeA_aba": source["b_routeA_aba"],
        "routeA_gauge_for_flipped_body_deg": source["routeA_gauge_for_flipped_body_deg"],
        "dt_gate": {"attempted_dt_s": 2.0e-7, "fallback_dt_s": 1.0e-7,
                    "accepted": None, "gate_status": "PENDING_DATACHECK_AND_DYNAMIC"},
        "sign_sanity_check": checks,
        "source_input_sha256": sha256(REPO / "calibration_analysis" / "S4HeadForwardLowGScreen" / "cases" / "S4_HEADFORWARD_G0" / "S4_HEADFORWARD_G0.inp"),
        "input_sha256": sha256(deck), "fortran_sha256": sha256(fortran),
        "centerline_sha256": sha256(centerline),
        "fast_server_sha256": sha256(REPO / "calibration_analysis" / "FastStraightDynamicScreen" / "production" / "magpylib_socket_server_fast.py"),
    }
    (CASE / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n", encoding="ascii")
    (OUT / "sign_sanity_check.json").write_text(json.dumps(checks, indent=2) + "\n", encoding="ascii")
    print(json.dumps(checks, indent=2))
    print("PREPARED {}".format(JOB))


if __name__ == "__main__":
    main()
