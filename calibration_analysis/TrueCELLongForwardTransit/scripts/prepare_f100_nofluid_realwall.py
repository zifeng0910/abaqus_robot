"""Prepare the sole F100 no-fluid, real-wall isolation case."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

from audit_a14p5_rotational_closure import mesh_mass_properties

ROOT = Path(__file__).resolve().parents[1]
PARENT = "TRUECEL_B0P11_G2P20_A14P5_F100_CLEAN50"
BASE = "TRUECEL_B0P11_G2P20_F100_NOFLUID_CONTROL"
JOB = "F100_G2P20_NOFLUID_REALWALL50"
TABLE = "magnetic_field_gradient_table_B0P11_A14P5.dat"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def section(deck: str, begin: str, end: str) -> str:
    assert deck.count(begin) == 1 and deck.count(end) == 1, (begin, end)
    return deck.split(begin, 1)[1].split(end, 1)[0]


def main() -> None:
    source = ROOT / "case" / PARENT
    base = ROOT / "case" / BASE
    target = ROOT / "case" / JOB
    if target.exists():
        raise RuntimeError(f"Refusing to replace existing case: {target}")
    identity = json.loads((source / "case_identity.json").read_text())
    assert [identity[k] for k in ("B0_mT", "gradient_mT", "frequency_Hz",
            "rocking_main_amplitude_deg", "rocking_cross_amplitude_deg")] == [11., 2.2, 100., 14.5, 2.5]
    assert digest(source / TABLE) == identity["magnetic_table_sha256"]
    parent = (source / f"{PARENT}.inp").read_text(encoding="latin1")
    nofluid = (base / f"{BASE}.inp").read_text(encoding="latin1")
    mass, com, inertia, principal, _, elements, volume = mesh_mass_properties(parent)
    assert elements == 7302 and 10.0e-9 < mass < 10.1e-9
    assert "*Mass" not in parent and "*Rotary Inertia" not in parent

    law = "*Surface Interaction, name=PROP_CEL_HARD\n" + section(parent,
        "*Surface Interaction, name=PROP_CEL_HARD\n", "*Surface Interaction, name=PROP_CEL_FLUID_HARD")
    pair = "*Contact\n" + section(parent, "*Contact\n", "** FAST_SURROGATE, PRECOMPUTED_TABLE")
    pair = "\n".join(line for line in pair.splitlines()
                     if "FLUID_CEL_SURF" not in line) + "\n"
    assert pair.count("ROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF") == 2
    assert "FLUID" not in pair and "PROP_CEL_HARD" in pair
    deck = nofluid.replace(BASE, JOB)
    boundary = "*Boundary\nRP_PIPE, ENCASTRE\n"
    assert deck.count(boundary) == 1
    deck = deck.replace(boundary, law + boundary)
    step = "*Step, name=Step_Drive, nlgeom=YES\n"
    assert deck.count(step) == 1
    deck = deck.replace(step, pair + step)
    field = "*Node Output, nset=RP_ROBOT\nU, UR, V, VR\n"
    assert deck.count(field) == 1
    deck = deck.replace(field, field +
        "*Contact Output, surface=ROBOT_SOLID-1.ROBOT_SOLID_SURF\nCDISP, CNORMF, CSHEARF\n" +
        "*Contact Output, surface=Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF\nCDISP, CNORMF, CSHEARF\n")
    hist = "*Node Output, nset=RP_ROBOT\nCOORD, U, UR, V, VR, A, AR, CF, RF, RM\n"
    assert deck.count(hist) == 1
    deck = deck.replace(hist, hist +
        "*Contact Output, surface=ROBOT_SOLID-1.ROBOT_SOLID_SURF\nCFN, CFS, CFT\n" +
        "*Contact Output, surface=Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF\nCFN, CFS, CFT\n")
    # Preserve the parent's 25-us contact field cadence without changing dynamics.
    assert deck.count("*Output, field, time interval=0.00833333333333, time marks=NO") == 1
    deck = deck.replace("*Output, field, time interval=0.00833333333333, time marks=NO",
                        "*Output, field, time interval=2.5e-5, time marks=NO")
    for forbidden in ("FLUID_EULERIAN", "FLUID_CEL", "MAT_FLUID", "*EOS", "*Initial Conditions"):
        assert forbidden.lower() not in "\n".join(x.lower() for x in deck.splitlines() if not x.startswith("**")), forbidden
    assert deck.count("*Step, name=Step_Drive") == 1 and "*Restart, read" not in deck
    assert "*Dynamic, Explicit, DIRECT\n1.0e-7, 0.050000000000" in deck
    target.mkdir(parents=True)
    inp = target / f"{JOB}.inp"
    inp.write_text(deck, encoding="latin1")
    original_sub = (source / "vuamp_precomputed_truecel.f90").read_bytes()
    assert original_sub.count(PARENT.encode()) == 6
    f90 = target / "vuamp_precomputed_truecel.f90"
    f90.write_bytes(original_sub.replace(PARENT.encode(), JOB.encode()))
    assert f90.read_bytes().replace(JOB.encode(), PARENT.encode()) == original_sub
    for name in (TABLE, f"{TABLE}.json"):
        shutil.copy2(source / name, target / name)
    audit = {
        "case_id": JOB, "parent": PARENT, "parent_input_sha256": digest(source / f"{PARENT}.inp"),
        "input_sha256": digest(inp), "fortran_sha256": digest(f90),
        "magnetic_table_sha256": digest(target / TABLE),
        "robot_mesh_elements": elements, "robot_mesh_volume_mm3": volume,
        "robot_mass_mg": mass * 1e9, "robot_com_aba_mm": com.tolist(),
        "robot_inertia_about_com_Nmm_s2": inertia.tolist(),
        "I11_I22_I33_about_com_Nmm_s2": inertia.diagonal().tolist(),
        "principal_inertia_Nmm_s2": principal.tolist(),
        "mass_source": "parent Robot_SOLID C3D4 tetrahedral integration * MAT_ROBOT_RIGID density; Abaqus rigid body inherits this mass/inertia",
        "stale_parent_metadata_robot_mass_mg": identity["robot_mass_mg"],
        "contact_law": law.strip(), "contact_pair": pair.strip(),
        "direct_dt_s": 1e-7, "duration_s": .05,
        "restart_read": False, "restart_write": "inherited 10 intervals; never used for continuation",
        "fluid": False, "rotational_damping_added": False,
    }
    (target / "setup_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    new_identity = dict(identity)
    new_identity.update(case_id=JOB, source_case=PARENT, status="PREPARED",
                        classification="PENDING", dynamics_run_count=0,
                        fluid_mode="NONE", contact_mode="REAL_ROBOT_WALL_ONLY",
                        robot_mass_mg=mass * 1e9, input_sha256=digest(inp),
                        direct_dt_s=1e-7, restart_read=False)
    (target / "case_identity.json").write_text(json.dumps(new_identity, indent=2) + "\n")
    print(json.dumps({"case": JOB, "mass_mg": mass * 1e9,
                      "inertia_diagonal": inertia.diagonal().tolist(),
                      "input_sha256": digest(inp)}, indent=2))


if __name__ == "__main__":
    main()
