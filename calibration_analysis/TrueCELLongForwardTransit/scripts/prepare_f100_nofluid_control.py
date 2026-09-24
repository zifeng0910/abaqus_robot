"""Remove only fluid and contact physics from the frozen F100 CLEAN50 deck."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_NAME = "TRUECEL_B0P11_G2P20_A14P5_F100_CLEAN50"
JOB = "TRUECEL_B0P11_G2P20_F100_NOFLUID_CONTROL"
SOURCE = ROOT / "case" / SOURCE_NAME
DEST = ROOT / "case" / JOB
TABLE = "magnetic_field_gradient_table_B0P11_A14P5.dat"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def keyword(line: str) -> str:
    return line.strip().lower() if line.startswith("*") and not line.startswith("**") else ""


def remove_section(lines: list[str], start: str, end: str, end_inclusive: bool,
                   expected: int = 1) -> tuple[list[str], int]:
    out, count, i = [], 0, 0
    while i < len(lines):
        current = keyword(lines[i])
        matches = current == start if start == "*contact" else current.startswith(start)
        if matches:
            count += 1
            j = i + 1
            while j < len(lines) and not keyword(lines[j]).startswith(end):
                j += 1
            if j == len(lines):
                raise RuntimeError(f"Missing closing keyword {end!r} after {start!r}")
            i = j + int(end_inclusive)
        else:
            out.append(lines[i])
            i += 1
    if count != expected:
        raise RuntimeError(f"Expected {expected} {start!r} sections, found {count}")
    return out, count


def remove_keyword_blocks(lines: list[str], predicate, expected: int) -> tuple[list[str], int]:
    out, count, i = [], 0, 0
    while i < len(lines):
        k = keyword(lines[i])
        if k and predicate(k):
            count += 1
            i += 1
            while i < len(lines) and not keyword(lines[i]):
                i += 1
        else:
            out.append(lines[i])
            i += 1
    if count != expected:
        raise RuntimeError(f"Expected {expected} matching keyword blocks, found {count}")
    return out, count


def main() -> None:
    if DEST.exists() and any(DEST.iterdir()):
        prior_path = DEST / "case_identity.json"
        if not prior_path.exists():
            raise RuntimeError(f"Control case already exists: {DEST}")
        prior = json.loads(prior_path.read_text(encoding="utf-8"))
        prior_input = DEST / f"{JOB}.inp"
        if (prior.get("status") != "PREPARED" or
                prior.get("dynamics_run_count") != 0 or
                not prior_input.exists() or
                sha(prior_input) != prior.get("input_sha256")):
            raise RuntimeError("Refusing to replace a changed or submitted control case")
    identity = json.loads((SOURCE / "case_identity.json").read_text(encoding="utf-8"))
    for key, expected in {
        "B0_mT": 11.0, "gradient_mT": 2.2, "frequency_Hz": 100.0,
        "rocking_main_amplitude_deg": 14.5,
        "rocking_cross_amplitude_deg": 2.5, "duration_s": .05,
    }.items():
        if identity[key] != expected:
            raise RuntimeError(f"Frozen {key} changed: {identity[key]}")
    if sha(SOURCE / TABLE) != identity["magnetic_table_sha256"]:
        raise RuntimeError("Frozen magnetic table hash differs")

    original = (SOURCE / f"{SOURCE_NAME}.inp").read_text(encoding="latin1")
    lines = original.splitlines(keepends=True)
    removed = {}
    lines, removed["Eulerian_part"] = remove_section(
        lines, "*part, name=fluid_eulerian", "*end part", True)
    lines, removed["Eulerian_instance"] = remove_section(
        lines, "*instance, name=fluid_eulerian-1", "*end instance", True)
    lines, removed["fluid_material"] = remove_section(
        lines, "*material, name=mat_fluid_cel", "*material, name=mat_pipe_rigid", False)
    lines, removed["interaction_properties"] = remove_section(
        lines, "*surface interaction", "*boundary", False)
    lines, removed["contact_definition"] = remove_section(
        lines, "*contact", "*step, name=step_drive", False)
    filters = (
        ("Eulerian_elsets", lambda k: k.startswith("*elset, elset=fluid_cel"), 2),
        ("Eulerian_surface", lambda k: k.startswith("*surface, type=eulerian material"), 1),
        ("fluid_initial_conditions", lambda k: k.startswith("*initial conditions"), 2),
        ("contact_output", lambda k: k.startswith("*contact output"), 3),
        ("Eulerian_field_output", lambda k: (
            k.startswith("*node output, nset=fluid_eulerian-1") or
            k.startswith("*element output, elset=fluid_cel_all")), 2),
    )
    for label, predicate, expected in filters:
        lines, removed[label] = remove_keyword_blocks(lines, predicate, expected)
    deck = "** AUTHORITATIVE MAGNETIC-ONLY CONTROL; WALL GEOMETRY PASSIVE; NO CONTACT\n" + \
           "".join(lines).replace(SOURCE_NAME, JOB)
    old_step = "*Dynamic, Explicit, SCALE FACTOR=0.4\n, 0.050000000000"
    new_step = "*Dynamic, Explicit, DIRECT\n1.0e-7, 0.050000000000"
    if deck.count(old_step) != 1:
        raise RuntimeError("Expected one automatic-time-increment step")
    deck = deck.replace(old_step, new_step, 1)
    active = "\n".join(line.lower() for line in deck.splitlines()
                       if line.strip() and not line.startswith("**"))
    forbidden = ("fluid_eulerian", "fluid_cel", "mat_fluid", "prop_cel",
                 "*contact", "*eulerian", "*eos", "*viscosity", "*initial conditions")
    if any(token in active for token in forbidden):
        raise RuntimeError("Fluid/contact keyword remains active in control deck")
    if active.count("*step, name=step_drive") != 1 or "*restart, read" in active:
        raise RuntimeError("Control must have one fresh-start Explicit step")
    for part in ("Robot_SOLID", "Pipe_WALL_HELPER", "Pipe_TERMINAL_STOP",
                 "Pipe_TERMINAL_STOP_HIGH", "Pipe_SOLID"):
        start, end = f"*Part, name={part}", "*End Part"
        source_chunk = original[original.index(start):original.index(end, original.index(start)) + len(end)]
        if source_chunk not in deck:
            raise RuntimeError(f"Frozen geometry changed: {part}")
    for token in ("*Rigid Body, ref node=RP_ROBOT, elset=ROBOT_SOLID_CEL_ALL",
                  "*Dynamic, Explicit, DIRECT",
                  "1.0e-7, 0.050000000000"):
        if token not in deck:
            raise RuntimeError(f"Frozen robot/step setting missing: {token}")

    sub = (SOURCE / "vuamp_precomputed_truecel.f90").read_bytes()
    if sub.count(SOURCE_NAME.encode()) != 6:
        raise RuntimeError("Unexpected VUAMP job path count")
    new_sub = sub.replace(SOURCE_NAME.encode(), JOB.encode())
    if new_sub.replace(JOB.encode(), SOURCE_NAME.encode()) != sub:
        raise RuntimeError("VUAMP changed beyond case paths")
    DEST.mkdir(parents=True, exist_ok=True)
    inp = DEST / f"{JOB}.inp"
    f90 = DEST / "vuamp_precomputed_truecel.f90"
    inp.write_bytes(deck.encode("latin1"))
    f90.write_bytes(new_sub)
    shutil.copy2(SOURCE / TABLE, DEST / TABLE)
    shutil.copy2(SOURCE / f"{TABLE}.json", DEST / f"{TABLE}.json")
    new_identity = {key: identity[key] for key in (
        "B0_mT", "gradient_mT", "frequency_Hz", "rocking_main_amplitude_deg",
        "rocking_cross_amplitude_deg", "duration_s", "initial_center_aba_mm",
        "canonical_plus_s_axis_aba", "n_routeA_aba", "b_routeA_aba",
        "head_tail_axis_aba", "initial_magnetic_moment_axis_aba",
        "robot_moment_Am2", "robot_mass_mg", "robot_length_mm",
        "robot_diameter_mm", "radial_offset_n_mm", "s_start_mm",
        "lumen_radius_mm", "axial_CEL_extent_s_mm", "magnetic_table_sha256")}
    new_identity.update({
        "case_id": JOB, "source_case": SOURCE_NAME, "status": "PREPARED",
        "classification": "PENDING", "dynamics_run_count": 0,
        "fluid_mode": "NONE", "wall_geometry": "PRESENT_PASSIVE",
        "contact_mode": "NONE", "artificial_damping_added": False,
        "magnetic_backend": "SAME_PRECOMPUTED_MAGPYLIB_TABLE_AND_VUAMP",
        "time_step_control": "DIRECT 1.0e-7 s; required because all remaining bodies are rigid",
        "direct_dt_s": 1e-7,
        "restart_read": False, "restart_write_times_ms": list(range(5, 51, 5)),
        "input_sha256": sha(inp), "fortran_sha256": sha(f90),
        "wall_reference_extent_s_mm": identity["axial_CEL_extent_s_mm"],
    })
    (DEST / "case_identity.json").write_text(
        json.dumps(new_identity, indent=2) + "\n", encoding="utf-8")
    audit = {
        "case": JOB, "source": SOURCE_NAME,
        "source_input_sha256": sha(SOURCE / f"{SOURCE_NAME}.inp"),
        "control_input_sha256": sha(inp),
        "source_fortran_sha256": sha(SOURCE / "vuamp_precomputed_truecel.f90"),
        "control_fortran_sha256": sha(f90),
        "magnetic_table_sha256": sha(DEST / TABLE),
        "removed_blocks": removed,
        "wall_geometry_present": True, "wall_contact_present": False,
        "fluid_domain_present": False, "fluid_contact_present": False,
        "single_explicit_step_s": .05, "restart_read": False,
        "time_step_change": "Automatic stable increment -> DIRECT 1.0e-7 s; all-rigid Explicit model requires it",
        "restart_write_times_ms": list(range(5, 51, 5)),
    }
    (ROOT / f"{JOB}_SETUP.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
