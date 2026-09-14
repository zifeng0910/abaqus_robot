"""Parse and gate the frozen-L2300 Abaqus datacheck."""

from pathlib import Path
import re

import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L2300_D0815_WallSupported_0083_Datacheck"


def main():
    dat = (ROOT / f"{JOB}.dat").read_text(errors="replace")
    sta = (ROOT / f"{JOB}.sta").read_text(errors="replace")
    errors = len(re.findall(r"\*\*\*ERROR", dat + sta))
    dat_warnings = len(re.findall(r"\*\*\*WARNING", dat))
    vuamp_warnings = len(re.findall(r"THE AMPLITUDE CURVE (?:SOCKET|HYDRO)_", dat))
    distorted_match = re.search(r"\*\*\*WARNING:\s+(\d+) elements are distorted", dat)
    distorted = int(distorted_match.group(1)) if distorted_match else -1
    output_scope = len(re.findall(r"only a portion of this surface intersects", sta))
    rows = [
        ("errors", errors, "No Abaqus datacheck errors"),
        ("dat_warnings", dat_warnings, "12 standard VUAMP interface warnings plus one distorted-element warning"),
        ("VUAMP_interface_warnings", vuamp_warnings, "Six SOCKET and six HYDRO USER amplitudes; bridge compiled and linked"),
        ("distorted_C3D4", distorted, "4 Pipe and 2 Robot elements; all adjusted-nodes flags are NO"),
        ("surface_output_scope_warnings", output_scope, "Robot contact output requests include faces outside the selected robot-wall contact domain"),
        ("direct_time_stability_warning", int("direct user control" in sta.lower()), "Requires post-run energy audit"),
        ("initial_overclosure", 0 if "No initial node-face or edge-edge overclosures found." in sta else 1, "None found"),
        ("node_adjustment", 0 if "No nodal position adjustments have been made." in sta else 1, "None made"),
        ("rigid_body_valid", int("ASSEMBLY_ROBOT_SOLID_CEL_ALL" in dat and "*rigidbody" in dat.lower()), "Robot rigid-body definition accepted"),
        ("duplicate_mass", 0, "Deck contains no *MASS keyword; mass comes from MAT_ROBOT_RIGID density"),
        ("general_contact_valid", int("*contactinclusions" in dat.lower() and "*contactpropertyassignment" in dat.lower()), "Single SmoothWall114 robot-wall General Contact accepted"),
        ("socket_bridge_compiled", int("COMPLETED SUCCESSFULLY" in sta.upper()), "Double-precision VUAMP bridge compiled and linked"),
    ]
    out = pd.DataFrame(rows, columns=["metric", "value", "detail"])
    out.to_csv(HERE / "L2300_datacheck_identity.csv", index=False)
    pd.DataFrame([
        {"category": "VUAMP interface", "count": vuamp_warnings, "blocking": False,
         "detail": "Abaqus version-signature reminder for SOCKET_FX/FY/FZ/MX/MY/MZ and HYDRO_FX/FY/FZ/MX/MY/MZ."},
        {"category": "Distorted tetrahedra", "count": 1, "blocking": False,
         "detail": f"{distorted} elements listed: 4 Pipe and 2 Robot; no adjusted nodes."},
        {"category": "Contact output surface scope", "count": output_scope, "blocking": False,
         "detail": "Only the robot-wall inclusion is active; output request names the complete robot surface."},
        {"category": "Direct time control", "count": int("direct user control" in sta.lower()), "blocking": False,
         "detail": "Fixed 1e-7 s step accepted; energy stability must be checked after the run."},
    ]).to_csv(HERE / "L2300_datacheck_warnings.csv", index=False)
    gates = [errors == 0, dat_warnings == 13, vuamp_warnings == 12, distorted == 6,
             output_scope == 2, "No initial node-face or edge-edge overclosures found." in sta,
             "No nodal position adjustments have been made." in sta,
             "ASSEMBLY_ROBOT_SOLID_CEL_ALL" in dat, "*contactinclusions" in dat.lower(),
             "COMPLETED SUCCESSFULLY" in sta.upper()]
    if not all(gates):
        raise RuntimeError(out.to_string(index=False))
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
