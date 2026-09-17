"""Prepare the single fitted refinement candidate for the v_mean screen."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
SOURCE_JOB = "FAST_PRECOMP_F120_G1P20_FLOW_1CYCLE"
SOURCE = REPO / "calibration_analysis" / "FastPrecomputedForwardFlowScreen" / "case" / SOURCE_JOB
JOB = "FAST_VMEAN5_G1P184_UFLOW0P5_2CYCLES"
G = 1.184
FLOW = (0.4881399801232148, -0.03091601237234885, 0.1037475782093262)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest().upper()


def main():
    plan = json.loads((OUT / "refinement_plan.json").read_text())
    if abs(plan["predicted_G_mT_for_5mm_s"] - G) > 1e-12:
        raise RuntimeError("refinement G does not match fitted plan")
    source_deck = (SOURCE / (SOURCE_JOB + ".inp")).read_text(encoding="latin1")
    source_sub = (SOURCE / "vuamp_precomputed_g1p20_flow.f90").read_text(encoding="ascii")
    source_identity = json.loads((SOURCE / "case_identity.json").read_text())
    case = OUT / "case" / JOB
    deck = source_deck.replace(SOURCE_JOB, JOB)
    deck = deck.replace("G=0.15 -> 1.20 mT", "G=0.15 -> 1.184 mT", 1)
    deck = deck.replace("ONE 120 HZ CYCLE", "TWO 120 HZ CYCLES", 1)
    deck = deck.replace("2.0e-7, 0.008333333333", "2.0e-7, 0.016666666667", 1)
    sub = source_sub.replace("grad=0.00120d0*grad", "grad=0.001184d0*grad", 1)
    sub = sub.replace(
        "data flow /9.762799602464296d0,-0.618320247446977d0,2.074951564186524d0/",
        "data flow /{:.16f}d0,{:.16f}d0,{:.16f}d0/".format(*FLOW), 1)
    sub = sub.replace(
        "FastPrecomputedForwardFlowScreen\\case\\" + SOURCE_JOB,
        "FastVMean5Calibration\\case\\" + JOB, 1)
    if "grad=0.001184d0*grad" not in sub or "data flow /0.4881399801232148d0" not in sub:
        raise RuntimeError("Fortran replacement failed")
    case.mkdir(parents=True, exist_ok=True)
    inp = case / (JOB + ".inp")
    f90 = case / "vuamp_vmean5.f90"
    inp.write_text(deck, encoding="latin1")
    f90.write_text(sub, encoding="ascii")
    identity = dict(source_identity)
    for key in ("wallclock_s", "cpus", "socket_calls", "gradient_choice_basis", "input_sha256", "fortran_sha256"):
        identity.pop(key, None)
    identity.update({
        "case_id": JOB, "source_case": SOURCE_JOB, "status": "PREPARED", "gradient_mT": G,
        "duration_s": 0.016666666667, "commanded_cycles": 2.0, "U_flow_mm_s": 0.5,
        "flow_vector_aba_mm_s": list(FLOW), "dynamics_run_count": 0,
        "calibration_role": "SINGLE_FITTED_REFINEMENT",
        "gradient_choice_basis": "linear fit to contact-safe G=2,3,4 mT two-cycle results",
        "input_sha256": sha(inp), "fortran_sha256": sha(f90),
        "frozen": ["B0=10mT", "f=120Hz", "A_main=14.8deg", "A_cross=2.5deg", "geometry",
                   "contact law", "initial pose", "hydrodynamic coefficients"],
        "forward_sign_sanity": {"passed": True, "basis": "positive normalized gradient basis and canonical +s polarity"},
    })
    (case / "case_identity.json").write_text(json.dumps(identity, indent=2) + "\n", encoding="ascii")
    manifest = {"job": JOB, "G_mT": G, "input_sha256": sha(inp), "fortran_sha256": sha(f90)}
    (OUT / "refinement_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="ascii")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
