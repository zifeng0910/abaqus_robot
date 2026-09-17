"""Prepare the sole non-CEL F120 precomputed-magnetic contact screen."""
from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
JOB = "FAST_PRECOMP_F120_DUALEND_NOREBOUND"
SOURCE_JOB = "TAIL_STICK_F120_ZETA100_8P333MS_TRUECEL"
SOURCE = REPO / "calibration_analysis" / "RefinedDualEndTrueCELGate" / "case" / SOURCE_JOB
CASE = OUT / "case" / JOB
TABLE = REPO / "calibration_analysis" / "PrecomputedStraightMagneticBackend" / "magnetic_field_gradient_table.dat"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def replace_once(text: str, pattern: str, replacement: str, label: str) -> str:
    result, count = re.subn(pattern, replacement, text, count=1, flags=re.I | re.M | re.S)
    if count != 1:
        raise RuntimeError("expected one {} block, found {}".format(label, count))
    return result


def main() -> None:
    source_path = SOURCE / (SOURCE_JOB + ".inp")
    source_identity = json.loads((SOURCE / "case_identity.json").read_text())
    deck = source_path.read_text(encoding="latin1").replace(SOURCE_JOB, JOB)
    deck = ("** AUTHORITATIVE CASE: FAST_PRECOMP_F120_DUALEND_NOREBOUND\n"
            "** NON-CEL FAST_SURROGATE_FLUID; PRECOMPUTED_TABLE; TWO 120 HZ CYCLES\n"
            "** HEAD FORWARD/RIGHT; TAIL REAR/LEFT; CANONICAL +s LEFT TO RIGHT\n"
            "** INHERITED PROVENANCE COMMENTS BELOW ARE NON-AUTHORITATIVE\n" + deck)

    deck = replace_once(deck, r"^\*Part, name=FLUID_EULERIAN\s*$.*?^\*End Part\s*$\n?", "", "Eulerian part")
    deck = replace_once(deck, r"^\*Instance, name=Fluid_EULERIAN-1.*?(?=^\*Elset, elset=PIPE_SOLID_CEL_ALL)", "", "Eulerian instance/sets")
    deck = replace_once(deck, r"^\*Surface, type=EULERIAN MATERIAL, name=FLUID_CEL_SURF\s*$\nFluid_EULERIAN-1_WATER\s*$\n?", "", "Eulerian surface")
    deck = replace_once(deck, r"^\*\* TRUE CEL:.*?(?=^\*\* MATERIALS)", "** FAST_SURROGATE_FLUID: no CEL mesh, EOS, advection, or fluid contact.\n", "CEL initialization")
    deck = replace_once(deck, r"^\*Material, name=MAT_FLUID_CEL\s*$.*?(?=^\*Material, name=MAT_PIPE_RIGID)", "", "fluid material")

    old_property = ("*Surface Interaction, name=PROP_CEL_HARD\n"
                    "*Friction\n0.03,\n"
                    "*Surface Behavior, pressure-overclosure=HARD\n"
                    "*Contact Damping, definition=CRITICAL DAMPING FRACTION, tangent fraction=0.0\n1.0")
    new_property = ("*Surface Interaction, name=PROP_CEL_HARD\n"
                    "*Friction\n0.03,\n"
                    "*Surface Behavior, pressure-overclosure=SCALE FACTOR\n"
                    "5.0, , 10.0, 0.01\n"
                    "*Contact Damping, definition=CRITICAL DAMPING FRACTION, tangent fraction=0.0\n1.0")
    if old_property not in deck:
        raise RuntimeError("source contact property changed")
    deck = deck.replace(old_property, new_property, 1)
    deck = replace_once(
        deck,
        r"^\*Contact Inclusions\s*$.*?^\*Contact Property Assignment\s*$.*?(?=^\*\* TRUE-CEL BASELINE)",
        "*Contact Inclusions\nROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF\n"
        "*Contact Property Assignment\nROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF, PROP_CEL_HARD\n",
        "contact pairs")
    deck = replace_once(deck, r"^\*\* TRUE-CEL BASELINE.*?(?=^\*Step, name=Step_Drive)",
                        "** FAST_SURROGATE, PRECOMPUTED_TABLE, fixed refined geometry.\n", "legacy CEL comments")
    deck = replace_once(deck, r"^\*Dynamic, Explicit\s*$\n, 0\.008333333333\s*$",
                        "*Dynamic, Explicit, DIRECT\n2.0e-7, 0.016666666667", "step duration")
    deck = replace_once(deck, r"^\*\* MAGPYLIB_SOCKET_LOADS.*?(?=^\*Amplitude, name=SOCKET_FX)",
                        "** PRECOMPUTED_TABLE magnetic loads plus FAST_SURROGATE_FLUID.\n"
                        "** No Python, TCP, Magpylib, CEL, prescribed motion, or velocity clamp during dynamics.\n",
                        "load comments")
    hydro = """*Amplitude, name=HYDRO_FX, definition=USER
*Amplitude, name=HYDRO_FY, definition=USER
*Amplitude, name=HYDRO_FZ, definition=USER
*Amplitude, name=HYDRO_MX, definition=USER
*Amplitude, name=HYDRO_MY, definition=USER
*Amplitude, name=HYDRO_MZ, definition=USER
*Cload, amplitude=HYDRO_FX
RP_ROBOT, 1, 1.
*Cload, amplitude=HYDRO_FY
RP_ROBOT, 2, 1.
*Cload, amplitude=HYDRO_FZ
RP_ROBOT, 3, 1.
*Cload, amplitude=HYDRO_MX
RP_ROBOT, 4, 1.
*Cload, amplitude=HYDRO_MY
RP_ROBOT, 5, 1.
*Cload, amplitude=HYDRO_MZ
RP_ROBOT, 6, 1.
** FAST_SURROGATE_FLUID uses V_rel=V_robot-10*c_hat mm/s.
"""
    deck = replace_once(deck, r"^\*\* REDUCEDHYDRO_EXTERNAL_LOADS = OFF;.*?$", hydro.rstrip(), "hydro insertion")
    deck = re.sub(r"^\*Element Output, elset=FLUID_CEL_ALL\s*$\nEVF\s*$\n?", "", deck, flags=re.I | re.M)
    deck = replace_once(deck, r"^\*Output, field, time interval=2\.5e-4, time marks=NO$",
                        "*Output, field, time interval=1.0e-5, time marks=NO", "field interval")
    deck = replace_once(deck, r"^\*Node Output\s*$\nU, V\s*$",
                        "*Node Output, nset=Robot_SOLID-1.ROBOT_CEL_BODY\nU, V", "robot-only field output")

    forbidden = ("FLUID_EULERIAN", "FLUID_CEL", "MAT_FLUID_CEL", "*Eulerian", "*Eos", "EVF")
    present = [token for token in forbidden if token.lower() in deck.lower()]
    if present:
        raise RuntimeError("CEL residue in FAST deck: {}".format(present))
    required = ("*Dynamic, Explicit, DIRECT\n2.0e-7, 0.016666666667",
                "pressure-overclosure=SCALE FACTOR", "5.0, , 10.0, 0.01",
                "HYDRO_FX", "RP_VR3", "ROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF")
    missing = [token for token in required if token not in deck]
    if missing:
        raise RuntimeError("missing FAST deck requirements: {}".format(missing))

    CASE.mkdir(parents=True, exist_ok=True)
    inp = CASE / (JOB + ".inp")
    inp.write_text(deck, encoding="latin1")
    shutil.copy2(SOURCE / "straight_control_centerline.dxf", CASE / "straight_control_centerline.dxf")
    shutil.copy2(TABLE, CASE / "magnetic_field_gradient_table.dat")
    fortran_source = OUT / "vuamp_precomputed_fast_surrogate.f90"
    shutil.copy2(fortran_source, CASE / fortran_source.name)

    identity = dict(source_identity)
    for stale in ("failure", "simulated_time_s", "available_commanded_cycles", "stop_reason", "solver_end_state", "wallclock_s"):
        identity.pop(stale, None)
    for cel_only in ("flow_classification", "time_step_control", "U_fluid_initial_mm_s", "fluid_density_tonne_mm3",
                     "fluid_viscosity_N_s_mm2", "fluid_EOS_c0_mm_s", "fluid_filled_volume_mm3", "fluid_mass_mg",
                     "fluid_element_count", "Eulerian_element_count", "Eulerian_transverse_bounds_mm",
                     "failed_technical_attempt_count", "failed_technical_attempt_reason", "only_physics_change", "physics_changes"):
        identity.pop(cel_only, None)
    identity.update({
        "case_id": JOB, "source_case": SOURCE_JOB, "status": "PREPARED",
        "classification": "FAST_SURROGATE_FLUID", "magnetic_backend": "PRECOMPUTED_TABLE",
        "precomputed_magnetic_baseline_commit": "74b8205a39092943a76bb499b11387be8d78dff1",
        "magnetic_table_sha256": sha(TABLE), "frequency_Hz": 120.0,
        "duration_s": 0.016666666667, "commanded_cycles": 2.0, "direct_dt_s": 2.0e-7,
        "U_flow_mm_s": 10.0, "Cparallel_Ns_per_mm": 4.0e-9,
        "Cperp_Ns_per_mm": 1.2e-8, "Kspin_Nmm_s": 1.0e-9, "Kwobble_Nmm_s": 3.0e-9,
        "contact_normal_behavior_old": "HARD penalty, zeta=1.0",
        "contact_normal_behavior_new": "SCALE FACTOR r=5%, geometric scale=10, initial scale=0.01; zeta=1.0",
        "contact_law_intent": "weak initial impact impulse with progressive nonpenetration; no adhesion/no-separation/clamp",
        "dynamics_run_count": 0, "socket_calls_expected": 0,
        "time_step_control": "Abaqus/Explicit DIRECT 2e-7 s",
        "flow_model": "FAST_SURROGATE_FLUID relative-velocity drag; not CEL/FSI",
        "contact_output": "direct robot-wall General Contact fields partitioned by actual HEAD/TAIL node regions",
        "source_input_sha256": sha(source_path), "input_sha256": sha(inp),
        "fortran_sha256": sha(CASE / fortran_source.name),
    })
    (CASE / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n", encoding="ascii")
    audit = {
        "old": {"pressure_overclosure": "HARD", "enforcement": "General Contact penalty", "mu": 0.03,
                "normal_critical_damping_fraction": 1.0, "tangent_fraction": 0.0},
        "new": {"pressure_overclosure": "SCALE FACTOR", "overclosure_factor_percent": 5.0,
                "geometric_stiffness_scale": 10.0, "initial_stiffness_scale": 0.01,
                "mu": 0.03, "normal_critical_damping_fraction": 1.0, "tangent_fraction": 0.0},
        "minimum_robot_tet_edge_mm": 0.0449331796,
        "first_stiffness_transition_overclosure_mm": 0.00224665898,
        "rationale": "Low initial penalty stiffness weakens impact impulse; geometric progression restores nonpenetration as overclosure grows.",
        "architecture": "General Contact only; robot-wall assignment only",
    }
    (OUT / "contact_law_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="ascii")
    (OUT / "contact_law_audit.md").write_text(
        "# General Contact audit\n\nOld: General Contact penalty, HARD pressure-overclosure, mu=0.03, normal zeta=1.0.\n\n"
        "New (single modification): `SCALE FACTOR` with `r=5%`, geometric stiffness multiplier `10`, and initial stiffness scale `0.01`; damping and friction are unchanged. "
        "For the measured 0.0449331796 mm minimum robot tetra edge, the first transition is about 0.00224666 mm. "
        "This targets a weak initial impulse while progressively recovering nonpenetration. It does not add adhesion, no-separation, a velocity clamp, Contact Pair, or SplitWall.\n",
        encoding="ascii")
    print(json.dumps({"job": JOB, "input": str(inp), "input_sha256": sha(inp), "table_sha256": sha(TABLE)}, indent=2))


if __name__ == "__main__":
    main()
