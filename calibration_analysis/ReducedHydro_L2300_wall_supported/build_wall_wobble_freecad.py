"""Build the independent L2.300 wall-supported-wobble FreeCAD family."""
from pathlib import Path
import importlib.util
import json

import FreeCAD as App
import Part


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BUILDER = REPO / "cad" / "freecad_parametric_robot" / "build_parametric_robot.py"
STEM = "Robot_parametric_L2p300_D0p815_WallWobble"
L_TOTAL = 2.300
D_BODY = 0.815


def load_builder():
    spec = importlib.util.spec_from_file_location("robot_builder", BUILDER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    fit = json.loads((HERE / "L2300_head_fit.json").read_text())
    if not fit["fit_gate_pass"]:
        raise RuntimeError("Head fit gate failed")
    poles = [App.Vector(float(x), float(r), 0.0) for x, r in fit["poles_mm"]]
    curve = Part.BezierCurve()
    curve.setPoles(poles)
    head = curve.toShape()
    xh = float(fit["head_axial_extent_mm"])
    rb = D_BODY / 2.0
    cylinder = Part.makeLine(App.Vector(xh, rb, 0), App.Vector(L_TOTAL, rb, 0))
    tail = Part.makeLine(App.Vector(L_TOTAL, rb, 0), App.Vector(L_TOTAL, 0, 0))
    axis = Part.makeLine(App.Vector(L_TOTAL, 0, 0), App.Vector(0, 0, 0))
    wire = Part.Wire([head, cylinder, tail, axis])
    if not wire.isClosed():
        raise RuntimeError("Profile wire is not closed")
    shape = Part.Face(wire).revolve(App.Vector(0, 0, 0), App.Vector(1, 0, 0), 360).removeSplitter()
    builder = load_builder()
    gate = builder.validate_shape(shape, L_TOTAL, D_BODY, STEM)
    paths = {key: HERE / f"{STEM}.{ext}" for key, ext in
             (("fcstd", "FCStd"), ("step", "step"), ("brep", "brep"))}
    doc = App.newDocument(STEM)
    obj = doc.addObject("PartDesign::Feature", "Robot")
    obj.Label = "L2300 D0815 Wall-supported Wobble Robot"
    obj.Shape = shape
    for name, kind, value in (
        ("LTotal", "Length", L_TOTAL), ("DBody", "Length", D_BODY),
        ("HeadAxialExtent", "Length", xh),
        ("StraightBodyLength", "Length", L_TOTAL - xh),
        ("BezierDegree", "Integer", int(fit["degree"])),
    ):
        obj.addProperty(f"App::Property{kind}", name, "Robot Geometry")
        setattr(obj, name, value)
    doc.recompute(); doc.saveAs(str(paths["fcstd"])); Part.export([obj], str(paths["step"]))
    shape.exportBrep(str(paths["brep"])); App.closeDocument(doc.Name)
    exports = builder.validate_exports(paths, L_TOTAL, D_BODY)
    props = builder.mass_properties(shape)
    props.update({"model": STEM, "axis": "+X from nose to flat tail",
                  "constant_magnetization_rule": True})
    geometry = {
        "model": STEM,
        "geometry_family": "axisymmetric degree-5 Bezier nose, G2 tangent cylinder, flat tail",
        "authoritative_geometry": "FreeCAD Part/OpenCASCADE revolved BRep",
        "L_total_mm": L_TOTAL, "D_body_mm": D_BODY, "R_body_mm": rb,
        "head_axial_extent_mm": xh, "straight_body_length_mm": L_TOTAL - xh,
        "length_change_vs_old_L1800_mm": 0.500,
        "length_change_rule": "head and tail unchanged; add exactly 0.500 mm straight cylinder",
        "volume_mm3": shape.Volume, "surface_area_mm2": shape.Area,
        "fit_source": fit, "shape_gate": gate, "export_reimport_gate": exports,
        "freecad_version": ".".join(App.Version()),
    }
    (HERE / f"{STEM}_geometry.json").write_text(json.dumps(geometry, indent=2) + "\n")
    (HERE / f"{STEM}_mass_properties.json").write_text(json.dumps(props, indent=2) + "\n")
    print("BUILD_RESULT=" + json.dumps(geometry, sort_keys=True))


if __name__ == "__main__":
    main()
