"""Prepare exactly one clean-start, waveform-only F100/G2.20 candidate."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from audit_f100_wall_exposure import ROOT, CASE, unit

NAME = "TRUECEL_B0P11_G2P20_F100_ASYMROCK_FAST"
DEST = ROOT / "case" / NAME
TABLE = "magnetic_field_gradient_table_B0P11_ASYMROCK.dat"
CENTER = -0.4822655
AMPLITUDE = 14.0255795
CROSS = 2.5


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def once(text, old, new):
    if text.count(old) != 1:
        raise RuntimeError(f"Expected one {old!r}, found {text.count(old)}")
    return text.replace(old, new, 1)


def main():
    if DEST.exists():
        raise RuntimeError(f"Candidate already exists: {DEST}")
    ident = json.loads((CASE / "case_identity.json").read_text())
    for key, expected in {"B0_mT": 11., "gradient_mT": 2.2, "frequency_Hz": 100.,
                          "rocking_main_amplitude_deg": 14.5, "rocking_cross_amplitude_deg": 2.5,
                          "fluid_EOS_c0_mm_s": 100000., "Eulerian_dimensions": [44, 20, 20],
                          "explicit_stable_time_scale_factor": .4}.items():
        if ident.get(key) != expected:
            raise RuntimeError(f"Frozen parent mismatch: {key}")
    exposure = ROOT / "F100_G2P20_WALL_EXPOSURE_NATIVE.csv"
    if not exposure.is_file():
        raise RuntimeError("Phase-0 geometry must be audited first")
    old = CASE / "magnetic_field_gradient_table_B0P11_A14P5.dat"
    with old.open(encoding="ascii") as file:
        header = file.readline()
        rows = file.readlines()
    ns, nph = map(int, header.split()[:2])
    phase = np.array([[float(v) for v in row.split()] for row in rows[:nph]])
    spatial = np.array([[float(v) for v in row.split()] for row in rows[nph:]])
    if len(spatial) != ns or nph != 361:
        raise RuntimeError("Unexpected magnetic table dimensions")
    c, n, rock = [unit(ident[key]) for key in ("canonical_plus_s_axis_aba", "n_routeA_aba", "b_routeA_aba")]
    radians = np.deg2rad(phase[:, 0])
    old_alpha = np.deg2rad(14.5) * np.sin(radians)
    old_direction = (-c[None, :] - np.tan(old_alpha)[:, None] * n[None, :]
                     + np.tan(np.deg2rad(CROSS) * np.cos(radians))[:, None] * rock[None, :])
    old_direction /= np.linalg.norm(old_direction, axis=1)[:, None]
    if np.max(np.abs(old_direction - phase[:, 1:4])) > 1e-12:
        raise RuntimeError("Analytic field expression does not match authoritative table")
    alpha = np.deg2rad(CENTER + AMPLITUDE * np.sin(radians))
    direction = (-c[None, :] - np.tan(alpha)[:, None] * n[None, :]
                 + np.tan(np.deg2rad(CROSS) * np.cos(radians))[:, None] * rock[None, :])
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    phase[:, 1:4] = direction
    deck = (CASE / f"{CASE.name}.inp").read_text(encoding="latin1").replace(CASE.name, NAME)
    deck = once(deck, ", 0.020000000000", ", 0.030000000000")
    sub = (CASE / "vuamp_precomputed_truecel.f90").read_text(encoding="latin1").replace(CASE.name, NAME)
    sub = once(sub, "magnetic_field_gradient_table_B0P11_A14P5.dat", TABLE)
    sub = sub.replace("magnetic_increment_g2p20_f100.csv", "magnetic_increment_asymrock_f100.csv")
    sub = sub.replace("f100_event.txt", "asymrock_event.txt")
    if "phase_deg=modulo(36000.0d0*t,360.0d0)" not in sub or "b=0.011d0*b; grad=0.002200d0*grad" not in sub:
        raise RuntimeError("B0/G/frequency implementation changed")
    DEST.mkdir(parents=True)
    inp = DEST / f"{NAME}.inp"
    f90 = DEST / "vuamp_precomputed_truecel.f90"
    inp.write_text(deck, encoding="latin1")
    f90.write_text(sub, encoding="latin1")
    tab = DEST / TABLE
    with tab.open("w", encoding="ascii", newline="\n") as file:
        file.write(header)
        np.savetxt(file, phase, fmt="%.17e")
        np.savetxt(file, spatial, fmt="%.17e")
    for key in ("wallclock_s", "cpus", "socket_calls", "failure", "change_scope",
                "frozen_physics", "frozen", "single_change", "frequency_change_only"):
        ident.pop(key, None)
    ident.update({"case_id": NAME, "status": "PREPARED", "classification": "ASYMROCK_THREE_CYCLE_GATE",
                  "physical_parent": CASE.name, "duration_s": .03, "dynamics_run_count": 0,
                  "rocking_main_amplitude_deg": AMPLITUDE, "rocking_main_center_deg": CENTER,
                  "rocking_cross_amplitude_deg": CROSS,
                  "single_physics_change": "symmetric to geometric-margin-balanced elliptic rocking field waveform",
                  "only_new_physics_change": "symmetric to geometrically balanced main rocking field-angle waveform",
                  "magnetic_table_file": TABLE, "stage1_duration_s": .03,
                  "geometric_touch_thresholds_deg": {"HEAD": -14.387845, "TAIL": 13.423314},
                  "frozen_for_candidate": ["B0=11mT", "G=2.20mT", "f=100Hz", "A_cross=2.5deg",
                      "TRUE-CEL 44x20x20 mesh/domain", "c0=100000mm/s", "water density/viscosity/EOS",
                      "robot/wall geometry/contact", "Explicit scale factor=0.4; no mass scaling", "clean initial state"],
                  "field_expression": "B_unit=normalize(-c-tan((center+amp*sin(2*pi*100*t))*deg)*n+tan(2.5*deg*cos(2*pi*100*t))*rock)",
                  "input_sha256": sha(inp), "fortran_sha256": sha(f90), "magnetic_table_sha256": sha(tab),
                  "fresh_t0": True})
    (DEST / "case_identity.json").write_text(json.dumps(ident, indent=2) + "\n", encoding="utf-8")
    (DEST / (TABLE + ".json")).write_text(json.dumps({"source_table_sha256": sha(old),
        "field_expression": ident["field_expression"], "center_deg": CENTER, "amplitude_deg": AMPLITUDE,
        "cross_deg": CROSS, "B0_mT": 11., "gradient_mT": 2.2, "frequency_Hz": 100.,
        "max_parent_formula_error": float(np.max(np.abs(old_direction - np.array([[float(v) for v in row.split()] for row in rows[:nph]])[:, 1:4]))),
        "table_sha256": sha(tab)}, indent=2) + "\n", encoding="ascii")
    print(json.dumps({"candidate": NAME, "table_sha256": sha(tab), "input_sha256": sha(inp),
                      "phase0_exposure_csv": str(exposure)}, indent=2))


if __name__ == "__main__":
    main()
