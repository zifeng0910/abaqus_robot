"""Prepare the clean 50-ms full-CEL case with only RP rotational damping added."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_NAME = "TRUECEL_B0P11_G2P20_A14P5_F100_CLEAN50"
JOB = "TRUECEL_B0P11_G2P20_A14P5_F100_CROT3_FULLCEL50"
SOURCE = ROOT / "case" / SOURCE_NAME
DEST = ROOT / "case" / JOB
TABLE = "magnetic_field_gradient_table_B0P11_A14P5.dat"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def replace_once(text: str, old: str, new: str) -> str:
    if text.count(old) != 1:
        raise RuntimeError(f"Expected one insertion anchor: {old[:70]}")
    return text.replace(old, new, 1)


def main() -> None:
    if DEST.exists():
        raise RuntimeError(f"Refusing to overwrite prepared or running case: {DEST}")
    parent = json.loads((SOURCE / "case_identity.json").read_text(encoding="utf-8"))
    expected = {"B0_mT": 11.0, "gradient_mT": 2.2, "frequency_Hz": 100.0,
                "rocking_main_amplitude_deg": 14.5,
                "rocking_cross_amplitude_deg": 2.5, "duration_s": .05}
    if any(parent.get(k) != v for k, v in expected.items()):
        raise RuntimeError("Frozen FULL CEL parent parameters differ")
    if sha(SOURCE / TABLE) != parent["magnetic_table_sha256"]:
        raise RuntimeError("Frozen magnetic table hash differs")
    screen = json.loads((ROOT / "TRUECEL_MAGNETIC_ONLY_DAMPING_OFFLINE_ESTIMATE.json").read_text())
    selected = next(x for x in screen["cases"] if x["factor"] == 3.0)
    c_rot = selected["c_rot_N_mm_s"]
    deck_parent = (SOURCE / f"{SOURCE_NAME}.inp").read_text(encoding="latin1")
    sub_parent = (SOURCE / "vuamp_precomputed_truecel.f90").read_text(encoding="latin1")
    if "*Dynamic, Explicit, SCALE FACTOR=0.4\n, 0.050000000000" not in deck_parent:
        raise RuntimeError("Expected clean uninterrupted Explicit step")
    for token in ("*Element, type=EC3D8R", "*Contact Inclusions",
                  "PROP_CEL_FLUID_HARD", "*Viscosity", "*Eos, type=USUP"):
        if token not in deck_parent:
            raise RuntimeError(f"Missing frozen CEL physics: {token}")
    deck = deck_parent.replace(SOURCE_NAME, JOB)
    if deck.replace(JOB, SOURCE_NAME) != deck_parent:
        raise RuntimeError("Deck changed beyond job label")
    sub = sub_parent.replace(SOURCE_NAME, JOB)
    decl = "  real(8) :: loads(6),last_time,t,u(3),ur(3),v(3),vs"
    sub = replace_once(sub, decl, decl + "\n  real(8) :: vr(3)\n" +
                       f"  real(8), parameter :: c_rot={f'{c_rot:.16e}'.replace('e', 'd')}")
    anchor = "    dotu=dot_product(u,c)"
    vr_lines = "\n".join(
        f"    vr({i})=VGETSENSORVALUE('RP_VR{i}',jSensorLookUpTable,sensorValues)"
        for i in range(1, 4))
    sub = replace_once(sub, anchor, vr_lines + "\n" + anchor)
    anchor = "    write(77,'(9(ES18.10,:,\",\"))') t,s_eff,phase_deg,loads"
    sub = replace_once(sub, anchor, anchor + "\n    loads(4:6)=loads(4:6)-c_rot*vr")
    if "loads(1:3)=loads(1:3)-" in sub:
        raise RuntimeError("Unexpected translation damping")
    DEST.mkdir(parents=True)
    inp = DEST / f"{JOB}.inp"
    f90 = DEST / "vuamp_precomputed_truecel.f90"
    inp.write_text(deck, encoding="latin1")
    f90.write_text(sub, encoding="latin1")
    (DEST / TABLE).write_bytes((SOURCE / TABLE).read_bytes())
    identity = dict(parent)
    identity.update({"case_id": JOB, "source_case": SOURCE_NAME, "status": "PREPARED",
                     "classification": "PENDING", "dynamics_run_count": 0,
                     "single_physics_change": "isotropic RP rotational damping torque -c_rot * VR",
                     "rotational_damping_c_N_mm_s": c_rot,
                     "rotational_damping_screen_factor": 3.0,
                     "translational_damping_added": False,
                     "fluid_and_contact_preserved": True,
                     "restart_read": False,
                     "input_sha256": sha(inp), "fortran_sha256": sha(f90),
                     "magnetic_table_sha256": sha(DEST / TABLE)})
    (DEST / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n")
    setup = {"case": JOB, "parent": SOURCE_NAME,
             "parent_input_sha256": sha(SOURCE / f"{SOURCE_NAME}.inp"),
             "candidate_input_sha256": sha(inp),
             "parent_fortran_sha256": sha(SOURCE / "vuamp_precomputed_truecel.f90"),
             "candidate_fortran_sha256": sha(f90),
             "magnetic_table_sha256": sha(DEST / TABLE),
             "c_rot_N_mm_s": c_rot,
             "source_screen": selected,
             "full_CEL_start_s": 0.0, "full_CEL_end_s": .05,
             "single_explicit_step": True,
             "restart_write_strategy": "number interval=10, time marks=YES; every 5 ms; never read",
             "fluid_contact_mesh_EOS_c0_and_magnetic_parameters_unchanged": True,
             "damping_scope": "RP rotational moment components 4-6 only"}
    (ROOT / f"{JOB}_SETUP.json").write_text(json.dumps(setup, indent=2) + "\n")
    print(json.dumps(setup, indent=2))


if __name__ == "__main__":
    main()
