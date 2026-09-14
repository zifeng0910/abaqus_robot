"""Build the unique HeadClearance CAD without altering the authoritative model."""
from pathlib import Path
import importlib.util
import json

import FreeCAD as App
import Part


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SOURCE = REPO / "cad" / "freecad_parametric_robot" / "build_parametric_robot.py"
STEM = "Robot_parametric_L1p800_D0p815_HeadClearance"


def load_builder():
    spec = importlib.util.spec_from_file_location("robot_builder", SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    fit = json.loads((HERE / "head_clearance_fit.json").read_text())
    if not fit["fit_gate_pass"]:
        raise RuntimeError("Head fit gate failed")
    builder = load_builder()
    ratio = fit["R_head_over_R_body"]
    params = builder.profile_definition(1.8, 0.815, ratio)
    shape = builder.make_robot(1.8, 0.815, ratio)
    gate = builder.validate_shape(shape, 1.8, 0.815, STEM)

    fcstd = HERE / f"{STEM}.FCStd"
    step = HERE / f"{STEM}.step"
    brep = HERE / f"{STEM}.brep"
    doc = App.newDocument(STEM)
    obj = doc.addObject("PartDesign::Feature", "Robot")
    obj.Label = "L1800 D0815 HeadClearance Robot"
    obj.Shape = shape
    builder._add_model_properties(obj, params)
    doc.recompute()
    doc.saveAs(str(fcstd))
    Part.export([obj], str(step))
    shape.exportBrep(str(brep))
    App.closeDocument(doc.Name)
    exports = builder.validate_exports({"fcstd": fcstd, "step": step, "brep": brep}, 1.8, 0.815)

    geometry = {
        "model": STEM,
        "geometry_family": "one convex circular-arc head, G1 tangent cylinder, flat tail",
        "axis": "+X from HEAD nose to flat TAIL",
        "L_total_mm": 1.8,
        "D_body_mm": 0.815,
        "R_body_mm": params["r_body_mm"],
        "R_head_design_mm": params["r_head_mm"],
        "R_head_over_R_body": params["head_radius_ratio"],
        "head_axial_extent_mm": params["x_head_mm"],
        "straight_body_length_mm": params["straight_length_mm"],
        "volume_mm3": shape.Volume,
        "volume_ratio_vs_current_exact_CAD": shape.Volume / 0.866162685064448,
        "shape_gate": gate,
        "export_reimport_gate": exports,
        "fit_source": fit,
        "freecad_version": ".".join(App.Version()),
    }
    mass = builder.mass_properties(shape)
    mass.update({"model": STEM, "constant_magnetization_rule": True})
    (HERE / f"{STEM}_geometry.json").write_text(json.dumps(geometry, indent=2) + "\n")
    (HERE / f"{STEM}_mass_properties.json").write_text(json.dumps(mass, indent=2) + "\n")
    print("BUILD_RESULT=" + json.dumps(geometry, sort_keys=True))


if __name__ == "__main__":
    main()
