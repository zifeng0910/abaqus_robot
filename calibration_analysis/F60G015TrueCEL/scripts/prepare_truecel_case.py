"""Prepare and hard-gate the single F60/G0.15 true-CEL screening case."""
from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE.parent
REPO = OUT.parents[1]
JOB = "S4_HEADFORWARD_F60_G0P15_TRUECEL"
CASE = OUT / "case" / JOB
PARENT_DIR = REPO / "calibration_analysis" / "F60LowGBackgroundFlowScreen"
PARENT_CASE = PARENT_DIR / "case" / "S4_HEADFORWARD_F60_G0P10_FLOW"
PARENT_INP = PARENT_CASE / "S4_HEADFORWARD_F60_G0P10_FLOW.inp"
PARENT_ID = PARENT_CASE / "case_identity.json"
LEGACY_CEL = REPO / "calibration_analysis" / "WobbleCal_COMFixed_F30_Cone30_B10" / "WobbleCal_COMFixed_F30_Cone30_B10.inp"
LEGACY_RUN_STA = Path(r"J:\abaqusfangzhen\WobbleCal_COMFixed_F30_Cone30_B10.sta")

RHO_TONNE_MM3 = 1.0007e-9
MU_N_S_MM2 = 7.1e-10
C0_MM_S = 1_480_000.0
CELL_MM = 0.15
S_MIN, S_MAX = -6.0, 6.0
Q_MIN, Q_MAX = -0.45, 0.45
LUMEN_RADIUS_MM = 0.4075
ROBOT_RADIUS_MM = 0.4075
ROBOT_LENGTH_MM = 2.4
U_FLOW_MM_S = 10.0
DURATION_S = 0.033333333333


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unit(values):
    values = np.asarray(values, dtype=float)
    return values / np.linalg.norm(values)


def chunks(values, size=16):
    values = list(values)
    return [", ".join(str(v) for v in values[i:i + size]) for i in range(0, len(values), size)]


def cell_interior_intersects_capsule(s, n, b):
    """Exact local AABB-to-capsule test; face-only touching is retained."""
    half_cylinder = 0.5 * ROBOT_LENGTH_MM - ROBOT_RADIUS_MM
    half_cell = 0.5 * CELL_MM
    ds = max(abs(s) - half_cell - half_cylinder, 0.0)
    dn = max(abs(n) - half_cell, 0.0)
    db = max(abs(b) - half_cell, 0.0)
    return ds * ds + dn * dn + db * db < (ROBOT_RADIUS_MM - 1.0e-9) ** 2


def make_fluid_part(rp0, c, n, b):
    ns = int(round((S_MAX - S_MIN) / CELL_MM))
    nq = int(round((Q_MAX - Q_MIN) / CELL_MM))
    if (ns, nq) != (80, 6):
        raise RuntimeError("unexpected CEL dimensions")

    def node_id(i, j, k):
        return 1 + i + (ns + 1) * (j + (nq + 1) * k)

    lines = [
        "** TRUE CEL LOCAL STRAIGHT DOMAIN: 80 x 6 x 6 EC3D8R cells; h=0.15 mm",
        "*Part, name=FLUID_EULERIAN",
        "*Node",
    ]
    for k in range(nq + 1):
        qb = Q_MIN + k * CELL_MM
        for j in range(nq + 1):
            qn = Q_MIN + j * CELL_MM
            for i in range(ns + 1):
                s = S_MIN + i * CELL_MM
                xyz = rp0 + s * c + qn * n + qb * b
                lines.append("{0}, {1:.12f}, {2:.12f}, {3:.12f}".format(node_id(i, j, k), *xyz))

    lines.append("*Element, type=EC3D8R")
    filled = []
    eid = 0
    for k in range(nq):
        qb = Q_MIN + (k + 0.5) * CELL_MM
        for j in range(nq):
            qn = Q_MIN + (j + 0.5) * CELL_MM
            for i in range(ns):
                s = S_MIN + (i + 0.5) * CELL_MM
                eid += 1
                conn = [node_id(i, j, k), node_id(i + 1, j, k),
                        node_id(i + 1, j + 1, k), node_id(i, j + 1, k),
                        node_id(i, j, k + 1), node_id(i + 1, j, k + 1),
                        node_id(i + 1, j + 1, k + 1), node_id(i, j + 1, k + 1)]
                lines.append("{}, {}".format(eid, ", ".join(str(value) for value in conn)))
                in_lumen = qn * qn + qb * qb <= LUMEN_RADIUS_MM ** 2
                if in_lumen and not cell_interior_intersects_capsule(s, qn, qb):
                    filled.append(eid)
    total_elements = eid
    total_nodes = (ns + 1) * (nq + 1) * (nq + 1)
    lines.extend([
        "*Elset, elset=FLUID_CEL_BODY, generate",
        "1, {}, 1".format(total_elements),
        "*Nset, nset=FLUID_CEL_NODES, generate",
        "1, {}, 1".format(total_nodes),
        "*Eulerian Section, elset=FLUID_CEL_BODY",
        "MAT_FLUID_CEL, WATER",
        "*End Part",
        "**",
    ])
    domain_volume = (S_MAX - S_MIN) * (Q_MAX - Q_MIN) ** 2
    filled_volume = len(filled) * CELL_MM ** 3
    return "\n".join(lines) + "\n", filled, total_elements, total_nodes, domain_volume, filled_volume


def replace_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise RuntimeError("{} expected once, found {}".format(label, count))
    return text.replace(old, new, 1)


def regex_once(text, pattern, replacement, label):
    text, count = re.subn(pattern, replacement, text, count=1, flags=re.I | re.S | re.M)
    if count != 1:
        raise RuntimeError("{} expected once, found {}".format(label, count))
    return text


def main():
    parent = json.loads(PARENT_ID.read_text(encoding="utf-8"))
    deck = PARENT_INP.read_text(encoding="latin1")
    legacy = LEGACY_CEL.read_text(encoding="latin1")
    if parent["status"] != "SOLVED" or parent["frequency_Hz"] != 60.0:
        raise RuntimeError("authoritative F60 parent is not solved")
    if "*Eulerian Section" in deck:
        raise RuntimeError("authoritative parent unexpectedly contains CEL")
    for marker in ("*Eulerian Section, elset=FLUID_CEL_BODY", "*Initial Conditions, type=VOLUME FRACTION",
                   "1.0007e-09", "1480000, 0., 0.", "7.1e-10"):
        if marker not in legacy:
            raise RuntimeError("legacy CEL source missing {}".format(marker))

    rp0 = np.asarray(parent["initial_center_aba_mm"], dtype=float)
    c = unit(parent["canonical_plus_s_axis_aba"])
    n = unit(parent["n_routeA_aba"])
    b = unit(parent["b_routeA_aba"])
    if np.dot(np.cross(c, n), b) < 0.999999:
        raise RuntimeError("RouteA basis is not right handed")
    fluid_part, filled, n_elem, n_node, domain_volume, filled_volume = make_fluid_part(rp0, c, n, b)
    fluid_mass_tonne = filled_volume * RHO_TONNE_MM3
    fluid_mass_mg = fluid_mass_tonne * 1.0e9
    if not filled or filled_volume <= 0.0 or fluid_mass_tonne <= 0.0:
        raise RuntimeError("EMPTY_EULERIAN_DOMAIN")

    deck = replace_once(deck, "*Part, name=Pipe_WALL_HELPER", fluid_part + "*Part, name=Pipe_WALL_HELPER", "fluid part insertion")
    assembly_insert = """*Instance, name=Fluid_EULERIAN-1, part=FLUID_EULERIAN
*End Instance
*Elset, elset=FLUID_CEL_ALL, instance=Fluid_EULERIAN-1, generate
1, {n_elem}, 1
*Elset, elset=FLUID_CEL_INIT_ALL, instance=Fluid_EULERIAN-1
{filled}
**
""".format(n_elem=n_elem, filled="\n".join(chunks(filled)))
    deck = replace_once(deck, "*Elset, elset=PIPE_SOLID_CEL_ALL", assembly_insert + "*Elset, elset=PIPE_SOLID_CEL_ALL", "fluid instance insertion")
    fluid_surface = """*Surface, type=EULERIAN MATERIAL, name=FLUID_CEL_SURF
Fluid_EULERIAN-1_WATER
*End Assembly"""
    deck = replace_once(deck, "*End Assembly", fluid_surface, "fluid material surface")
    flow = U_FLOW_MM_S * c
    initial_and_material = """** STRICT TRUE-CEL INITIAL MATERIAL AND UNIFORM INITIAL FLOW
*Initial Conditions, type=VOLUME FRACTION
FLUID_CEL_INIT_ALL, Fluid_EULERIAN-1_WATER, 1.0
*Initial Conditions, type=VELOCITY
Fluid_EULERIAN-1.FLUID_CEL_NODES, 1, {vx:.15g}
Fluid_EULERIAN-1.FLUID_CEL_NODES, 2, {vy:.15g}
Fluid_EULERIAN-1.FLUID_CEL_NODES, 3, {vz:.15g}
**
** MATERIALS
**
*Material, name=MAT_FLUID_CEL
*Density
1.0007e-09,
*Eos, type=USUP
1480000, 0., 0.
*Viscosity
7.1e-10,
""".format(vx=flow[0], vy=flow[1], vz=flow[2])
    deck = replace_once(deck, "** MATERIALS\n**", initial_and_material.rstrip("\n"), "fluid material insertion")

    deck = regex_once(
        deck,
        r"\*Contact Inclusions\s*\nROBOT_SOLID-1\.ROBOT_SOLID_SURF,\s*Pipe_WALL_HELPER-1\.PIPE_WALL_HELPER_SURF\s*\n",
        "*Contact Inclusions\n"
        "ROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF\n"
        "ROBOT_SOLID-1.ROBOT_SOLID_SURF, FLUID_CEL_SURF\n"
        "Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF, FLUID_CEL_SURF\n",
        "CEL contact inclusions",
    )
    deck = regex_once(
        deck,
        r"\*Contact Property Assignment\s*\nROBOT_SOLID-1\.ROBOT_SOLID_SURF,\s*Pipe_WALL_HELPER-1\.PIPE_WALL_HELPER_SURF,\s*PROP_CEL_HARD\s*\n",
        "*Contact Property Assignment\n"
        "ROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF, PROP_CEL_HARD\n"
        "ROBOT_SOLID-1.ROBOT_SOLID_SURF, FLUID_CEL_SURF, PROP_CEL_HARD\n"
        "Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF, FLUID_CEL_SURF, PROP_CEL_HARD\n",
        "CEL contact assignments",
    )
    deck = regex_once(deck, r"\*Dynamic, Explicit, DIRECT\s*\n2\.0e-7,\s*0\.037500000", "*Dynamic, Explicit\n, 0.033333333333", "automatic CEL time step")
    deck = regex_once(
        deck,
        r"\*Amplitude, name=HYDRO_FX, definition=USER.*?\*Cload, amplitude=HYDRO_MZ\s*\nRP_ROBOT, 6, 1\.\s*\n",
        "** REDUCEDHYDRO_EXTERNAL_LOADS = OFF; CEL is the only hydrodynamic load source.\n",
        "ReducedHydro removal",
    )
    deck = deck.replace("** F60 LOW-G BACKGROUND-FLOW DYNAMIC SCREEN", "** F60 G0.15 STRICT TRUE-CEL DYNAMIC SCREEN", 1)
    deck = deck.replace("** S4_HEADFORWARD_F60_G0P10_FLOW; f=60Hz; A_main=14.343111711438091deg; A_cross=2.5deg; G=0.10mT; L=45mm; duration=0.0375s",
                        "** S4_HEADFORWARD_F60_G0P15_TRUECEL; f=60Hz; A_main=14.343111711438091deg; A_cross=2.5deg; G=0.15mT; L=45mm; duration=0.033333333333s", 1)
    deck = deck.replace("** NON-CEL; ROBOT_LOCAL_ELLIPTIC_ROCKING; FIXED ROUTEA GAUGE; DT=2e-7s; LOW_PRECISION_SCREENING",
                        "** STRICT_CEL_FLUID_PRESENT; ROBOT_LOCAL_ELLIPTIC_ROCKING; FIXED ROUTEA GAUGE; AUTOMATIC EXPLICIT DT", 1)
    deck = deck.replace("** PRESCRIBED BACKGROUND FLOW ONLY; U_FLOW=+10mm/s ALONG CANONICAL +s; DIAGNOSTIC_SCREENING_FLOW_ONLY",
                        "** INITIALIZED_UNIFORM_FLOW_TRUE_CEL; actual Eulerian fluid U0=+10mm/s along canonical +s", 1)
    deck = deck.replace("** Reduced-Hydro: CEL fluid part, material, initialization and CEL contact were removed.",
                        "** TRUE CEL: local fluid part, material initialization, and CEL contact are active.", 1)
    deck = deck.replace("** REDUCED-HYDRO BASELINE; CEL REMOVED; production external Magpylib retained.",
                        "** TRUE-CEL BASELINE DERIVED FROM F60; production external Magpylib retained.", 1)
    deck = deck.replace("** REDUCED-HYDRO LOAD INTERFACE: magnetic Socket loads plus dissipative body-following hydro loads.",
                        "** MAGPYLIB_SOCKET_LOADS = ON; REDUCEDHYDRO_EXTERNAL_LOADS = OFF; CEL_SOLID_COUPLING = ON.", 1)
    deck = deck.replace("** SOCKET_* remain production Magpylib loads; HYDRO_* use RP V_REL with prescribed +s background flow.",
                        "** SOCKET_* are the only external loads; actual CEL fluid supplies all hydrodynamic forces.", 1)
    deck = deck.replace("*Output, field, time interval=5.0e-4, time marks=NO\n*Node Output\nU",
                        "*Output, field, time interval=2.5e-4, time marks=NO\n*Node Output\nU, V\n*Element Output, elset=FLUID_CEL_ALL\nEVF", 1)

    forbidden = ("name=HYDRO_FX", "name=HYDRO_FY", "name=HYDRO_FZ", "name=HYDRO_MX", "name=HYDRO_MY", "name=HYDRO_MZ")
    if any(token in deck for token in forbidden):
        raise RuntimeError("ReducedHydro amplitude survived deck preparation")
    required = ("*Element, type=EC3D8R", "*Eulerian Section", "*Initial Conditions, type=VOLUME FRACTION",
                "*Initial Conditions, type=VELOCITY", "ROBOT_SOLID-1.ROBOT_SOLID_SURF, FLUID_CEL_SURF",
                "*Amplitude, name=SOCKET_FX", "*Element Output, elset=FLUID_CEL_ALL\nEVF")
    if not all(token in deck for token in required):
        raise RuntimeError("required true-CEL marker missing")

    CASE.mkdir(parents=True, exist_ok=True)
    deck_path = CASE / (JOB + ".inp")
    deck_path.write_text(deck, encoding="latin1")
    shutil.copy2(PARENT_CASE / "straight_control_centerline.dxf", CASE / "straight_control_centerline.dxf")
    shutil.copy2(HERE / "vuforc_magnetic_only.f", CASE / "vuforc_magnetic_only.f")

    delta_s = parent["sign_sanity_check"]["gradient_gate"]["delta_s_mm"]
    gl = 45.0
    d_b_ds = -(delta_s / gl ** 2) * math.exp(-0.5 * (delta_s / gl) ** 2) * 0.15e-3 * 1000.0
    f_gradient = -parent["robot_moment_Am2"] * d_b_ds
    times = [0.0, 0.0041667, 0.0083333, 0.0125]
    samples = []
    for time_s in times:
        phase = 2.0 * math.pi * 60.0 * time_s
        a_main = math.radians(parent["rocking_main_amplitude_deg"]) * math.sin(phase)
        a_cross = math.radians(parent["rocking_cross_amplitude_deg"]) * math.cos(phase)
        direction = c - math.tan(a_main) * n - math.tan(a_cross) * b
        direction = unit(direction)
        samples.append({"t_s": time_s, "parent_unit_B": direction.tolist(), "candidate_unit_B": direction.tolist(),
                        "difference_norm": 0.0, "alpha_main_deg": math.degrees(a_main), "alpha_cross_deg": math.degrees(a_cross)})
    magnetic_gate = {
        "parent_case": "S4_HEADFORWARD_F60_G0P10_FLOW", "candidate_case": JOB,
        "same_routeA_basis": True, "same_field_frame_mode": parent["field_frame_mode"] == "ROBOT_LOCAL_ELLIPTIC_ROCKING",
        "same_frequency_Hz": 60.0, "same_A_main_deg": parent["rocking_main_amplitude_deg"],
        "same_A_cross_deg": parent["rocking_cross_amplitude_deg"], "same_B0_mT": parent["B0_mT"],
        "same_phase_deg": 0.0, "only_magnetic_change": "gradient_mT: 0.10 -> 0.15",
        "samples": samples, "F_gradient_dot_canonical_s_N": f_gradient,
        "F_gradient_dot_canonical_s_positive": f_gradient > 0.0,
        "passed": f_gradient > 0.0 and all(item["difference_norm"] == 0.0 for item in samples),
    }
    if not magnetic_gate["passed"]:
        raise RuntimeError("magnetic baseline gate failed")

    source_audit = {
        "classification": "AUTHORITATIVE_LEGACY_ROUTEA_CEL_REUSE",
        "legacy_input": str(LEGACY_CEL), "legacy_run_status_file": str(LEGACY_RUN_STA),
        "legacy_run_evidence": "real FLUID_EULERIAN critical element; stable dt about 1.06e-7 s through at least 55.5 ms",
        "parameters": [
            {"parameter": "fluid_density", "value": RHO_TONNE_MM3, "units": "tonne/mm^3", "source": str(LEGACY_CEL), "action": "reused"},
            {"parameter": "dynamic_viscosity", "value": MU_N_S_MM2, "units": "N*s/mm^2", "source": str(LEGACY_CEL), "action": "reused"},
            {"parameter": "EOS", "value": {"type": "USUP", "c0_mm_s": C0_MM_S, "s": 0.0, "Gamma0": 0.0}, "source": str(LEGACY_CEL), "action": "reused"},
            {"parameter": "Eulerian_element", "value": "EC3D8R", "source": str(LEGACY_CEL), "action": "reused"},
            {"parameter": "Eulerian_section", "value": "MAT_FLUID_CEL, WATER", "source": str(LEGACY_CEL), "action": "reused"},
            {"parameter": "initial_volume_fraction", "value": 1.0, "source": str(LEGACY_CEL), "action": "reused on new straight-lumen mask"},
            {"parameter": "CEL_contact", "value": "Eulerian material surface with general contact", "source": str(LEGACY_CEL), "action": "reused; parent zeta=0.50 retained"},
            {"parameter": "initial_fluid_velocity", "value": flow.tolist(), "units": "mm/s", "source": "task requirement", "action": "new actual Eulerian nodal initial velocity"},
        ],
    }
    presence_gate = {
        "classification": "STRICT_CEL_FLUID_PRESENT", "passed": True,
        "Eulerian_geometric_volume_mm3": domain_volume, "fluid_filled_volume_mm3": filled_volume,
        "initial_material_volume_fraction": 1.0, "initial_total_fluid_mass_tonne": fluid_mass_tonne,
        "initial_total_fluid_mass_mg": fluid_mass_mg, "Eulerian_element_count": n_elem,
        "fluid_containing_element_count": len(filled), "Eulerian_node_count": n_node,
        "cell_size_mm": CELL_MM, "initial_velocity_mm_s": flow.tolist(),
        "checks": {"fluid_filled_volume_gt_0": filled_volume > 0.0, "fluid_mass_gt_0": fluid_mass_tonne > 0.0,
                   "fluid_containing_elements_gt_0": len(filled) > 0},
    }
    setup_audit = {
        "MAGPYLIB_SOCKET_LOADS": "ON", "REDUCEDHYDRO_EXTERNAL_LOADS": "OFF",
        "CEL_FLUID_PRESENT": True, "CEL_SOLID_COUPLING": "ON",
        "flow_classification": "INITIALIZED_UNIFORM_FLOW_TRUE_CEL",
        "structural_contact": {"general_contact": True, "mu": parent["mu"], "zeta": parent["zeta"]},
        "domain": {"axial_length_mm": S_MAX - S_MIN, "cross_section_mm": [Q_MAX - Q_MIN, Q_MAX - Q_MIN],
                   "cell_size_mm": CELL_MM, "open_axial_ends": True},
        "parent_commit": "a6dbd1422e898e708833bd427e4faa8860ee700e",
    }
    identity = {
        "case_id": JOB, "status": "PREPARED", "classification": "STRICT_CEL_FLUID_PRESENT",
        "flow_classification": "INITIALIZED_UNIFORM_FLOW_TRUE_CEL", "source_case": parent["case_id"],
        "source_commit": "a6dbd1422e898e708833bd427e4faa8860ee700e", "frequency_Hz": 60.0,
        "rocking_main_amplitude_deg": parent["rocking_main_amplitude_deg"],
        "rocking_cross_amplitude_deg": parent["rocking_cross_amplitude_deg"], "B0_mT": 10.0,
        "gradient_mT": 0.15, "gradient_length_mm": 45.0, "duration_s": DURATION_S,
        "time_step_control": "Abaqus/Explicit automatic CEL stability", "U_fluid_initial_mm_s": U_FLOW_MM_S,
        "flow_vector_aba_mm_s": flow.tolist(), "fluid_density_tonne_mm3": RHO_TONNE_MM3,
        "fluid_viscosity_N_s_mm2": MU_N_S_MM2, "fluid_EOS_c0_mm_s": C0_MM_S,
        "fluid_filled_volume_mm3": filled_volume, "fluid_mass_mg": fluid_mass_mg,
        "fluid_element_count": len(filled), "Eulerian_element_count": n_elem,
        "field_frame_mode": parent["field_frame_mode"], "initial_center_aba_mm": parent["initial_center_aba_mm"],
        "canonical_plus_s_axis_aba": parent["canonical_plus_s_axis_aba"], "n_routeA_aba": parent["n_routeA_aba"],
        "b_routeA_aba": parent["b_routeA_aba"], "head_tail_axis_aba": parent["head_tail_axis_aba"],
        "initial_magnetic_moment_axis_aba": parent["initial_magnetic_moment_axis_aba"],
        "routeA_gauge_for_flipped_body_deg": parent["routeA_gauge_for_flipped_body_deg"],
        "robot_moment_Am2": parent["robot_moment_Am2"], "robot_mass_mg": parent["robot_mass_mg"],
        "robot_length_mm": parent["robot_length_mm"], "robot_diameter_mm": parent["robot_diameter_mm"],
        "mu": parent["mu"], "zeta": parent["zeta"], "dynamics_run_count": 0,
        "failed_technical_attempt_count": 1,
        "failed_technical_attempt_reason": "attempt 1 terminated at 8.25 ms: centroid-only fluid mask overlapped free rigid robot and released nonphysical CEL offset energy",
        "source_input_sha256": sha256(PARENT_INP), "input_sha256": sha256(deck_path),
        "fortran_sha256": sha256(CASE / "vuforc_magnetic_only.f"),
    }
    for path, payload in ((OUT / "CEL_Fluid_Source_Audit.json", source_audit),
                          (OUT / "CEL_Setup_Audit.json", setup_audit),
                          (OUT / "magnetic_baseline_gate.json", magnetic_gate),
                          (OUT / "fluid_presence_gate.json", presence_gate),
                          (CASE / "case_identity.json", identity)):
        path.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
    print(json.dumps(presence_gate, indent=2))
    print("PREPARED {}".format(JOB))


if __name__ == "__main__":
    main()
