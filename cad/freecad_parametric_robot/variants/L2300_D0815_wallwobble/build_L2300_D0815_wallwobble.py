"""Build and validate the authoritative L2300 x D0815 FreeCAD robot CAD."""

from pathlib import Path
import csv
import hashlib
import json
import math

import FreeCAD as App
import Part


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
STEM = "Robot_L2300_D0815_WallWobble"
L_TOTAL_MM = 2.300
D_BODY_MM = 0.815
R_BODY_MM = D_BODY_MM / 2.0
HEAD_AXIAL_EXTENT_MM = 0.26149477
DENSITY_TONNE_PER_MM3 = 7.80906654321e-9
LINEAR_TOLERANCE_MM = 1.0e-6

# Frozen degree-5 Bezier definition. L changes only the straight cylinder.
HEAD_POLES_MM = (
    (0.0, 0.0),
    (0.0, 0.09664525),
    (0.005, 0.18497527),
    (0.17271875, R_BODY_MM),
    (0.25000000, R_BODY_MM),
    (HEAD_AXIAL_EXTENT_MM, R_BODY_MM),
)

# Constant-magnetization identity inherited from the validated L1800 target.
REFERENCE_VOLUME_MM3 = 0.7889904055827334
REFERENCE_MAGNETIC_MOMENT_AM2 = 0.0007196393110133208


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def bezier_point(t):
    degree = len(HEAD_POLES_MM) - 1
    x = 0.0
    radius = 0.0
    for i, (px, pr) in enumerate(HEAD_POLES_MM):
        weight = math.comb(degree, i) * (1.0-t)**(degree-i) * t**i
        x += weight * px
        radius += weight * pr
    return x, radius


def make_profile_wire(length_mm=L_TOTAL_MM):
    if length_mm <= HEAD_AXIAL_EXTENT_MM:
        raise ValueError("Total length must exceed the frozen head extent")
    curve = Part.BezierCurve()
    curve.setPoles([App.Vector(x, radius, 0.0) for x, radius in HEAD_POLES_MM])
    head = curve.toShape()
    cylinder = Part.makeLine(
        App.Vector(HEAD_AXIAL_EXTENT_MM, R_BODY_MM, 0.0),
        App.Vector(length_mm, R_BODY_MM, 0.0),
    )
    tail = Part.makeLine(
        App.Vector(length_mm, R_BODY_MM, 0.0), App.Vector(length_mm, 0.0, 0.0)
    )
    axis = Part.makeLine(App.Vector(length_mm, 0.0, 0.0), App.Vector(0.0, 0.0, 0.0))
    wire = Part.Wire([head, cylinder, tail, axis])
    if not wire.isClosed():
        raise RuntimeError("Profile wire is open")
    return wire


def make_robot(length_mm=L_TOTAL_MM):
    face = Part.Face(make_profile_wire(length_mm))
    return face.revolve(
        App.Vector(0.0, 0.0, 0.0), App.Vector(1.0, 0.0, 0.0), 360.0
    ).removeSplitter()


def make_head():
    curve = Part.BezierCurve()
    curve.setPoles([App.Vector(x, radius, 0.0) for x, radius in HEAD_POLES_MM])
    edge = curve.toShape()
    radial = Part.makeLine(
        App.Vector(HEAD_AXIAL_EXTENT_MM, R_BODY_MM, 0.0),
        App.Vector(HEAD_AXIAL_EXTENT_MM, 0.0, 0.0),
    )
    axis = Part.makeLine(
        App.Vector(HEAD_AXIAL_EXTENT_MM, 0.0, 0.0), App.Vector(0.0, 0.0, 0.0)
    )
    return Part.Face(Part.Wire([edge, radial, axis])).revolve(
        App.Vector(), App.Vector(1.0, 0.0, 0.0), 360.0
    ).removeSplitter()


def bbox(shape):
    box = shape.BoundBox
    return {"x": box.XLength, "y": box.YLength, "z": box.ZLength}


def validate_shape(shape, expected_length=L_TOTAL_MM, label="shape"):
    problems = []
    if shape.isNull():
        problems.append("null shape")
    if not shape.isValid():
        problems.append("BRepCheck invalid")
    try:
        shape.check(True)
        bop_check = "PASS"
    except Exception as exc:
        bop_check = "FAIL: %s" % exc
        problems.append("self-intersection/BOP check failed")
    if len(shape.Solids) != 1:
        problems.append("solid count %d" % len(shape.Solids))
    if len(shape.Shells) != 1:
        problems.append("shell count %d" % len(shape.Shells))
    if not shape.isClosed() or any(not shell.isClosed() for shell in shape.Shells):
        problems.append("open solid or seam")
    if shape.Volume <= 0.0:
        problems.append("non-positive volume")
    dimensions = bbox(shape)
    targets = {"x": expected_length, "y": D_BODY_MM, "z": D_BODY_MM}
    for axis, target in targets.items():
        if abs(dimensions[axis]-target) >= LINEAR_TOLERANCE_MM:
            problems.append("%s extent %.12g != %.12g" % (axis, dimensions[axis], target))
    if problems:
        raise RuntimeError("%s invalid: %s" % (label, "; ".join(problems)))
    return {
        "valid": True,
        "closed": True,
        "solid_count": len(shape.Solids),
        "shell_count": len(shape.Shells),
        "self_intersection_check": bop_check,
        "bbox_mm": dimensions,
        "volume_mm3": shape.Volume,
        "surface_area_mm2": shape.Area,
    }


def vector_values(vector):
    return [float(vector.x), float(vector.y), float(vector.z)]


def matrix_values(matrix):
    return [
        [float(matrix.A11), float(matrix.A12), float(matrix.A13)],
        [float(matrix.A21), float(matrix.A22), float(matrix.A23)],
        [float(matrix.A31), float(matrix.A32), float(matrix.A33)],
    ]


def mass_properties(shape):
    tensor = matrix_values(shape.MatrixOfInertia)
    principal = shape.PrincipalProperties
    moments = [float(value) for value in principal["Moments"]]
    com = shape.CenterOfMass
    if not shape.isInside(com, 1.0e-9, True):
        raise RuntimeError("Center of mass is not inside the solid")
    if min(moments) <= 0.0:
        raise RuntimeError("Inertia is not positive definite")
    mass_tonne = shape.Volume * DENSITY_TONNE_PER_MM3
    return {
        "density": {"value": DENSITY_TONNE_PER_MM3, "unit": "tonne/mm^3"},
        "volume": {"value": shape.Volume, "unit": "mm^3"},
        "surface_area": {"value": shape.Area, "unit": "mm^2"},
        "mass": {"tonne": mass_tonne, "mg": mass_tonne*1.0e9},
        "center_of_mass": {"value": vector_values(com), "unit": "mm", "inside_solid": True},
        "mass_inertia_tensor": {
            "value_tonne_mm2": [[x*DENSITY_TONNE_PER_MM3 for x in row] for row in tensor],
            "value_mg_mm2": [[x*DENSITY_TONNE_PER_MM3*1.0e9 for x in row] for row in tensor],
        },
        "principal_mass_moments": {
            "value_tonne_mm2": [x*DENSITY_TONNE_PER_MM3 for x in moments],
            "value_mg_mm2": [x*DENSITY_TONNE_PER_MM3*1.0e9 for x in moments],
            "positive_definite": True,
        },
        "principal_axes": {
            "first": vector_values(principal["FirstAxisOfInertia"]),
            "second": vector_values(principal["SecondAxisOfInertia"]),
            "third": vector_values(principal["ThirdAxisOfInertia"]),
        },
        "coordinate_system": "global +X from nose tip to flat tail",
    }


def export_and_roundtrip(shape):
    paths = {
        "fcstd": HERE / (STEM + ".FCStd"),
        "step": HERE / (STEM + ".step"),
        "brep": HERE / (STEM + ".brep"),
    }
    doc = App.newDocument(STEM)
    obj = doc.addObject("PartDesign::Feature", "Robot")
    obj.Label = "Authoritative L2300 D0815 WallWobble Robot"
    obj.Shape = shape
    properties = (
        ("LTotal", "Length", L_TOTAL_MM),
        ("DBody", "Length", D_BODY_MM),
        ("HeadAxialExtent", "Length", HEAD_AXIAL_EXTENT_MM),
        ("StraightBodyLength", "Length", L_TOTAL_MM-HEAD_AXIAL_EXTENT_MM),
        ("BezierDegree", "Integer", 5),
    )
    for name, kind, value in properties:
        obj.addProperty("App::Property" + kind, name, "Frozen Geometry")
        setattr(obj, name, value)
    doc.recompute()
    doc.saveAs(str(paths["fcstd"]))
    Part.export([obj], str(paths["step"]))
    shape.exportBrep(str(paths["brep"]))
    App.closeDocument(doc.Name)

    step = Part.read(str(paths["step"]))
    brep = Part.Shape()
    brep.read(str(paths["brep"]))
    doc = App.openDocument(str(paths["fcstd"]))
    try:
        objects = [obj for obj in doc.Objects if hasattr(obj, "Shape") and not obj.Shape.isNull()]
        if len(objects) != 1:
            raise RuntimeError("FCStd contains %d shape objects" % len(objects))
        fcstd_gate = validate_shape(objects[0].Shape, label="FCStd reopen")
    finally:
        App.closeDocument(doc.Name)
    gates = {
        "step": validate_shape(step, label="STEP round-trip"),
        "brep": validate_shape(brep, label="BREP round-trip"),
        "fcstd": fcstd_gate,
    }
    return paths, gates


def parameterization_checks(shape):
    short = make_robot(1.800)
    head_a = make_head()
    head_b = make_head()
    head_difference = head_a.cut(head_b).Volume + head_b.cut(head_a).Volume
    added_volume = shape.Volume-short.Volume
    expected_added = math.pi*R_BODY_MM**2*0.500
    result = {
        "head_independent_of_total_length": head_difference <= 1.0e-12,
        "head_symmetric_difference_mm3": head_difference,
        "L1800_to_L2300_total_length_delta_mm": L_TOTAL_MM-1.800,
        "straight_cylinder_delta_mm": (L_TOTAL_MM-HEAD_AXIAL_EXTENT_MM)-(1.800-HEAD_AXIAL_EXTENT_MM),
        "added_volume_mm3": added_volume,
        "expected_added_cylinder_volume_mm3": expected_added,
        "added_volume_error_mm3": abs(added_volume-expected_added),
    }
    if not result["head_independent_of_total_length"] or result["added_volume_error_mm3"] > 1.0e-10:
        raise RuntimeError("Length parameterization gate failed")
    return result


def write_profile_csv():
    path = HERE / (STEM + "_profile.csv")
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["region", "parameter", "x_mm", "radius_mm"])
        for i in range(2001):
            t = i/2000.0
            x, radius = bezier_point(t)
            writer.writerow(["HEAD", t, x, radius])
        writer.writerow(["BODY", 0.0, HEAD_AXIAL_EXTENT_MM, R_BODY_MM])
        writer.writerow(["BODY", 1.0, L_TOTAL_MM, R_BODY_MM])
        writer.writerow(["TAIL", 0.0, L_TOTAL_MM, R_BODY_MM])
        writer.writerow(["TAIL", 1.0, L_TOTAL_MM, 0.0])
    return path


def compare_previous(shape, props):
    old_dir = REPO / "calibration_analysis" / "ReducedHydro_L2300_wall_supported"
    old_step_path = old_dir / "Robot_parametric_L2p300_D0p815_WallWobble.step"
    old_props_path = old_dir / "Robot_parametric_L2p300_D0p815_WallWobble_mass_properties.json"
    old_shape = Part.read(str(old_step_path))
    old_props = json.loads(old_props_path.read_text())
    new_moments = props["principal_mass_moments"]["value_tonne_mm2"]
    old_moments = old_props["principal_mass_moments"]["value_tonne_mm2"]
    result = {
        "comparison_role": "previous generated reference only; never a construction input",
        "old_step": str(old_step_path.relative_to(REPO)).replace("\\", "/"),
        "bbox_difference_mm": {k: bbox(shape)[k]-bbox(old_shape)[k] for k in ("x", "y", "z")},
        "volume_difference_mm3": shape.Volume-old_shape.Volume,
        "surface_area_difference_mm2": shape.Area-old_shape.Area,
        "COM_difference_mm": [a-b for a, b in zip(
            props["center_of_mass"]["value"], old_props["center_of_mass"]["value"])],
        "principal_inertia_difference_tonne_mm2": [a-b for a, b in zip(new_moments, old_moments)],
        "head_profile_definition": "identical frozen degree-5 Bezier poles",
    }
    result["geometrically_equivalent_within_roundoff"] = bool(
        max(abs(x) for x in result["bbox_difference_mm"].values()) < 1.0e-12
        and abs(result["volume_difference_mm3"]) < 1.0e-12
        and abs(result["surface_area_difference_mm2"]) < 1.0e-8
        and max(abs(x) for x in result["COM_difference_mm"]) < 1.0e-12
        and max(abs(x) for x in result["principal_inertia_difference_tonne_mm2"]) < 1.0e-20
    )
    if not result["geometrically_equivalent_within_roundoff"]:
        raise RuntimeError("New canonical CAD differs from previous generated reference")
    return result


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    shape = make_robot()
    shape_gate = validate_shape(shape, label="master BRep")
    paths, roundtrip = export_and_roundtrip(shape)
    props = mass_properties(shape)
    parameterization = parameterization_checks(shape)
    profile_path = write_profile_csv()
    comparison = compare_previous(shape, props)
    magnetic_density = REFERENCE_MAGNETIC_MOMENT_AM2/REFERENCE_VOLUME_MM3
    magnetic = {
        "rule": "constant magnetization; magnetic moment is proportional to authoritative CAD volume",
        "reference_volume_mm3": REFERENCE_VOLUME_MM3,
        "reference_magnetic_moment_Am2": REFERENCE_MAGNETIC_MOMENT_AM2,
        "moment_density_Am2_per_mm3": magnetic_density,
        "authoritative_CAD_volume_mm3": shape.Volume,
        "physical_magnetic_moment_target_Am2": magnetic_density*shape.Volume,
        "server_started": False,
    }
    geometry = {
        "classification": "L2300_D0815_FREECAD_CAD_FROZEN",
        "model": STEM,
        "ownership": "FreeCAD Part/OpenCASCADE BRep",
        "axis": "global +X; nose tip X=0; flat tail X=L",
        "L_total_mm": L_TOTAL_MM,
        "D_body_mm": D_BODY_MM,
        "R_body_mm": R_BODY_MM,
        "head": {
            "family": "axisymmetric degree-5 Bezier",
            "degree": 5,
            "poles_mm": HEAD_POLES_MM,
            "axial_extent_mm": HEAD_AXIAL_EXTENT_MM,
            "monotone": True,
            "convex": True,
            "transition_to_cylinder": "G2",
            "previous_fit_max_outward_excess_um": 7.5606686254239275,
        },
        "straight_body_length_mm": L_TOTAL_MM-HEAD_AXIAL_EXTENT_MM,
        "length_parameterization": parameterization,
        "master_shape_gate": shape_gate,
        "roundtrip_gates": roundtrip,
        "freecad_version": ".".join(App.Version()),
    }
    geometry_path = HERE / (STEM + "_geometry.json")
    props_path = HERE / (STEM + "_mass_properties.json")
    magnetic_path = HERE / (STEM + "_magnetic_identity.json")
    comparison_path = HERE / "previous_generated_CAD_regression.json"
    write_json(geometry_path, geometry)
    write_json(props_path, props)
    write_json(magnetic_path, magnetic)
    write_json(comparison_path, comparison)
    identity_files = [paths["step"], paths["brep"], paths["fcstd"], geometry_path]
    manifest = {
        "classification": "L2300_D0815_FREECAD_CAD_FROZEN",
        "hash_algorithm": "SHA256",
        "files": {path.name: sha256(path) for path in identity_files},
        "future_Abaqus_STEP": "cad/freecad_parametric_robot/variants/L2300_D0815_wallwobble/" + paths["step"].name,
    }
    write_json(HERE / "cad_identity_manifest.json", manifest)
    print("BUILD_RESULT=" + json.dumps({
        "classification": manifest["classification"],
        "L_mm": shape.BoundBox.XLength,
        "D_mm": shape.BoundBox.YLength,
        "volume_mm3": shape.Volume,
        "mass_mg": props["mass"]["mg"],
        "COM_mm": props["center_of_mass"]["value"],
        "magnetic_moment_Am2": magnetic["physical_magnetic_moment_target_Am2"],
        "profile_csv": profile_path.name,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
