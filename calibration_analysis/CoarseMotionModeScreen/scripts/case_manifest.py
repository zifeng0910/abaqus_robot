"""Authoritative manifest for the human-in-the-loop coarse motion screen."""

from copy import deepcopy


BASE = dict(
    length_mm=2.30, B0_mT=10.0, gradient_mT=6.0, cone_deg=30.0,
    frequency_Hz=30.0, mu=0.03, zeta=0.50,
    cparallel_scale=1.0, kwobble_scale=1.0,
)


def case(case_id, family, **updates):
    value = deepcopy(BASE)
    value.update(case_id=case_id, family=family, **updates)
    value["duration_s"] = 0.2 / value["frequency_Hz"]
    value["job_name"] = "SCR_%s_Wobble" % case_id.replace("_", "")
    return value


CASES = [
    case("GEO_200", "geometry", length_mm=2.00),
    case("GEO_210", "geometry", length_mm=2.10),
    case("GEO_220", "geometry", length_mm=2.20),
    case("GEO_230", "geometry", length_mm=2.30),
    case("GEO_240", "geometry", length_mm=2.40),
    case("GRAD_0", "gradient", gradient_mT=0.0),
    case("GRAD_2", "gradient", gradient_mT=2.0),
    case("GRAD_4", "gradient", gradient_mT=4.0),
    case("GRAD_8", "gradient", gradient_mT=8.0),
    case("CONE_20", "cone", cone_deg=20.0),
    case("CONE_40", "cone", cone_deg=40.0),
    case("CONE_50", "cone", cone_deg=50.0),
    case("B_6", "B", B0_mT=6.0),
    case("B_14", "B", B0_mT=14.0),
    case("FREQ_20", "frequency", frequency_Hz=20.0),
    case("FREQ_40", "frequency", frequency_Hz=40.0),
    case("MU_000", "friction", mu=0.00),
    case("MU_080", "friction", mu=0.08),
    case("MU_150", "friction", mu=0.15),
    case("ZETA_015", "collision_damping", zeta=0.15),
    case("ZETA_030", "collision_damping", zeta=0.30),
    case("ZETA_070", "collision_damping", zeta=0.70),
    case("CPAR_X025", "hydro_translation", cparallel_scale=0.25),
    case("CPAR_X4", "hydro_translation", cparallel_scale=4.0),
    case("KWOB_X025", "hydro_rotation", kwobble_scale=0.25),
    case("KWOB_X4", "hydro_rotation", kwobble_scale=4.0),
    case("COMBO_A", "combo", length_mm=2.10, gradient_mT=0.0, cone_deg=40.0, mu=0.03, zeta=0.30),
    case("COMBO_B", "combo", length_mm=2.10, gradient_mT=2.0, cone_deg=40.0, mu=0.08, zeta=0.30),
    case("COMBO_C", "combo", length_mm=2.20, gradient_mT=0.0, cone_deg=40.0, mu=0.03, zeta=0.30),
    case("COMBO_D", "combo", length_mm=2.20, gradient_mT=2.0, cone_deg=40.0, mu=0.08, zeta=0.30),
    case("COMBO_E", "combo", length_mm=2.20, gradient_mT=2.0, cone_deg=50.0, frequency_Hz=20.0, mu=0.03, zeta=0.30),
    case("COMBO_F", "combo", length_mm=2.20, gradient_mT=2.0, B0_mT=14.0, cone_deg=40.0, mu=0.03, zeta=0.30),
]


ALIASES = {
    "GRAD_6": "GEO_230", "CONE_30": "GEO_230", "B_10": "GEO_230",
    "FREQ_30": "GEO_230", "MU_030": "GEO_230", "ZETA_050": "GEO_230",
    "CPAR_X1": "GEO_230", "KWOB_X1": "GEO_230",
}


def by_id(case_id):
    resolved = ALIASES.get(case_id, case_id)
    return next(item for item in CASES if item["case_id"] == resolved)

