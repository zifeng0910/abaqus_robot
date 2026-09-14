"""Build and validate the authoritative parametric FreeCAD robot geometry."""

import json
import math
from pathlib import Path

import FreeCAD as App
import Part


# User parameters: these are the only values normally changed for a new robot.
L_TOTAL_MM = 1.800
D_BODY_MM = 0.815
R_BODY_REF_MM = 0.600
R_HEAD_REF_MM = 0.610
HEAD_RADIUS_RATIO = R_HEAD_REF_MM / R_BODY_REF_MM


DENSITY_TONNE_PER_MM3 = 7.80906654321e-9
REFERENCE_LENGTH_MM = 3.050
REFERENCE_DIAMETER_MM = 1.200
LINEAR_TOLERANCE_MM = 1.0e-6
STL_LINEAR_DEFLECTION_MM = 0.005
STL_ANGULAR_DEFLECTION_DEG = 5.0
OUTPUT_DIR = Path(__file__).resolve().parent
TARGET_STEM = "Robot_parametric_L1p800_D0p815"
REFERENCE_STEM = "Robot_reference_L3p050_D1p200"


def geometry_parameters(length_mm, diameter_mm, head_radius_ratio=HEAD_RADIUS_RATIO):
    """Calculate the complete design definition without using length to size the head."""
    length_mm = float(length_mm)
    diameter_mm = float(diameter_mm)
    head_radius_ratio = float(head_radius_ratio)
    if length_mm <= 0.0 or diameter_mm <= 0.0 or head_radius_ratio <= 0.0:
        raise ValueError("Length, diameter, and head radius ratio must be positive")

    r_body = diameter_mm / 2.0
    r_head = r_body * head_radius_ratio
    radicand = 2.0 * r_head * r_body - r_body * r_body
    if radicand <= 0.0:
        raise ValueError("Head radius ratio cannot form the requested tangent circular arc")
    x_head = math.sqrt(radicand)
    if length_mm <= x_head:
        raise ValueError(
            "L_total must exceed the diameter-controlled head axial extent "
            f"({x_head:.12g} mm)"
        )
    return {
        "length_mm": length_mm,
        "diameter_mm": diameter_mm,
        "r_body_mm": r_body,
        "r_head_mm": r_head,
        "head_radius_ratio": head_radius_ratio,
        "x_head_mm": x_head,
        "straight_length_mm": length_mm - x_head,
    }


def profile_definition(length_mm, diameter_mm, head_radius_ratio=HEAD_RADIUS_RATIO):
    """Return the exact upper-profile definition used to create the BRep."""
    p = geometry_parameters(length_mm, diameter_mm, head_radius_ratio)
    r_body = p["r_body_mm"]
    r_head = p["r_head_mm"]
    x_head = p["x_head_mm"]
    center = (x_head, r_body - r_head)
    theta_tip = math.atan2(-center[1], -center[0])
    theta_tangent = math.pi / 2.0
    theta_mid = (theta_tip + theta_tangent) / 2.0
    arc_mid = (
        center[0] + r_head * math.cos(theta_mid),
        center[1] + r_head * math.sin(theta_mid),
    )
    p.update(
        {
            "center": center,
            "nose_tip": (0.0, 0.0),
            "arc_mid": arc_mid,
            "tangent_point": (x_head, r_body),
            "tail_outer": (p["length_mm"], r_body),
            "tail_axis": (p["length_mm"], 0.0),
            "theta_tip_rad": theta_tip,
            "theta_tangent_rad": theta_tangent,
        }
    )
    return p


def _v(point):
    return App.Vector(point[0], point[1], 0.0)


def make_profile_wire(length_mm, diameter_mm, head_radius_ratio=HEAD_RADIUS_RATIO):
    """Construct the exact closed axial profile from one circular arc and three lines."""
    p = profile_definition(length_mm, diameter_mm, head_radius_ratio)
    arc = Part.Arc(_v(p["nose_tip"]), _v(p["arc_mid"]), _v(p["tangent_point"])).toShape()
    cylinder = Part.makeLine(_v(p["tangent_point"]), _v(p["tail_outer"]))
    tail = Part.makeLine(_v(p["tail_outer"]), _v(p["tail_axis"]))
    axis = Part.makeLine(_v(p["tail_axis"]), _v(p["nose_tip"]))
    wire = Part.Wire([arc, cylinder, tail, axis])
    if not wire.isClosed():
        raise RuntimeError("ABORT: axial profile wire is not closed")
    return wire


def make_robot(length_mm, diameter_mm, head_radius_ratio=HEAD_RADIUS_RATIO):
    """Return one exact Part solid revolved 360 degrees around global X."""
    wire = make_profile_wire(length_mm, diameter_mm, head_radius_ratio)
    face = Part.Face(wire)
    shape = face.revolve(App.Vector(0, 0, 0), App.Vector(1, 0, 0), 360.0)
    return shape.removeSplitter()


def make_head(diameter_mm, head_radius_ratio=HEAD_RADIUS_RATIO):
    """Build the actual nose solid, independently of total robot length."""
    p = profile_definition(10.0, diameter_mm, head_radius_ratio)
    tip = _v(p["nose_tip"])
    tangent = _v(p["tangent_point"])
    tangent_axis = App.Vector(p["x_head_mm"], 0.0, 0.0)
    arc = Part.Arc(tip, _v(p["arc_mid"]), tangent).toShape()
    radial = Part.makeLine(tangent, tangent_axis)
    axis = Part.makeLine(tangent_axis, tip)
    wire = Part.Wire([arc, radial, axis])
    head = Part.Face(wire).revolve(
        App.Vector(0, 0, 0), App.Vector(1, 0, 0), 360.0
    )
    return head.removeSplitter()


def _bbox_dimensions(shape):
    box = shape.BoundBox
    return {"x": box.XLength, "y": box.YLength, "z": box.ZLength}


def validate_shape(shape, expected_length_mm, expected_diameter_mm, label="shape"):
    """Apply hard BRep validity and measured-dimension gates."""
    problems = []
    if shape.isNull():
        problems.append("null shape")
    if not shape.isValid():
        problems.append("OpenCASCADE BRepCheck reports invalid")
    try:
        shape.check(True)
    except Exception as exc:
        problems.append(f"BOP/self-intersection check failed: {exc}")
    if len(shape.Solids) != 1:
        problems.append(f"expected one solid, found {len(shape.Solids)}")
    if len(shape.Shells) != 1:
        problems.append(f"expected one shell, found {len(shape.Shells)}")
    if not shape.isClosed() or any(not shell.isClosed() for shell in shape.Shells):
        problems.append("solid or shell is not closed")
    if shape.Volume <= 0.0:
        problems.append(f"non-positive volume {shape.Volume!r}")
    dims = _bbox_dimensions(shape)
    expected = {
        "x": float(expected_length_mm),
        "y": float(expected_diameter_mm),
        "z": float(expected_diameter_mm),
    }
    for axis, target in expected.items():
        if abs(dims[axis] - target) >= LINEAR_TOLERANCE_MM:
            problems.append(
                f"bbox {axis}={dims[axis]:.12g} mm, expected {target:.12g} mm"
            )
    if problems:
        raise RuntimeError(f"ABORT: {label}: " + "; ".join(problems))
    return {
        "valid": True,
        "closed": True,
        "solid_count": len(shape.Solids),
        "shell_count": len(shape.Shells),
        "bbox_mm": dims,
        "volume_mm3": shape.Volume,
        "surface_area_mm2": shape.Area,
    }


def _vector_values(vector):
    return [float(vector.x), float(vector.y), float(vector.z)]


def _matrix3_values(matrix):
    return [
        [float(matrix.A11), float(matrix.A12), float(matrix.A13)],
        [float(matrix.A21), float(matrix.A22), float(matrix.A23)],
        [float(matrix.A31), float(matrix.A32), float(matrix.A33)],
    ]


def mass_properties(shape):
    """Return volume and uniform-density inertias about the shape center of mass."""
    geometric_tensor = _matrix3_values(shape.MatrixOfInertia)
    principal = shape.PrincipalProperties
    geometric_moments = [float(value) for value in principal["Moments"]]
    mass_tonne = shape.Volume * DENSITY_TONNE_PER_MM3
    return {
        "density": {
            "value": DENSITY_TONNE_PER_MM3,
            "unit": "tonne/mm^3",
            "assumption": "uniform density",
        },
        "volume": {"value": shape.Volume, "unit": "mm^3"},
        "mass": {
            "tonne": mass_tonne,
            "mg": mass_tonne * 1.0e9,
        },
        "center_of_mass": {
            "value": _vector_values(shape.CenterOfMass),
            "unit": "mm",
            "coordinate_system": "global XYZ; robot axis is +X",
        },
        "reference_point": "center of mass",
        "geometric_inertia_tensor": {
            "value": geometric_tensor,
            "unit": "mm^5",
            "definition": "integral of squared distance times dV",
            "coordinate_system": "global XYZ",
        },
        "mass_inertia_tensor": {
            "value_tonne_mm2": [
                [value * DENSITY_TONNE_PER_MM3 for value in row]
                for row in geometric_tensor
            ],
            "value_mg_mm2": [
                [value * DENSITY_TONNE_PER_MM3 * 1.0e9 for value in row]
                for row in geometric_tensor
            ],
            "coordinate_system": "global XYZ",
        },
        "principal_geometric_moments": {
            "value": geometric_moments,
            "unit": "mm^5",
        },
        "principal_mass_moments": {
            "value_tonne_mm2": [
                value * DENSITY_TONNE_PER_MM3 for value in geometric_moments
            ],
            "value_mg_mm2": [
                value * DENSITY_TONNE_PER_MM3 * 1.0e9
                for value in geometric_moments
            ],
        },
        "principal_axes": {
            "first": _vector_values(principal["FirstAxisOfInertia"]),
            "second": _vector_values(principal["SecondAxisOfInertia"]),
            "third": _vector_values(principal["ThirdAxisOfInertia"]),
            "coordinate_system": "global XYZ",
        },
    }


def verify_parameterization():
    """Verify length invariance and diameter-controlled scaling on actual BReps."""
    shape_a = make_robot(1.8, 0.815)
    shape_b = make_robot(2.2, 0.815)
    shape_c = make_robot(1.8, 0.9)
    head_a = make_head(0.815)
    head_b = make_head(0.815)
    head_c = make_head(0.9)
    p_a = geometry_parameters(1.8, 0.815)
    p_b = geometry_parameters(2.2, 0.815)
    p_c = geometry_parameters(1.8, 0.9)

    validate_shape(shape_a, 1.8, 0.815, "case A")
    validate_shape(shape_b, 2.2, 0.815, "case B")
    validate_shape(shape_c, 1.8, 0.9, "case C")
    validate_shape(head_a, p_a["x_head_mm"], 0.815, "case A head")
    validate_shape(head_b, p_b["x_head_mm"], 0.815, "case B head")
    validate_shape(head_c, p_c["x_head_mm"], 0.9, "case C head")

    head_difference = head_a.cut(head_b).Volume + head_b.cut(head_a).Volume
    head_embedding_error = (
        abs(shape_a.common(head_a).Volume - head_a.Volume)
        + abs(shape_b.common(head_b).Volume - head_b.Volume)
    )
    volume_delta = shape_b.Volume - shape_a.Volume
    expected_delta = math.pi * p_a["r_body_mm"] ** 2 * (2.2 - 1.8)
    diameter_scale = 0.9 / 0.815
    checks = {
        "same_head_parameters": all(
            abs(p_a[key] - p_b[key]) < 1.0e-12
            for key in ("r_body_mm", "r_head_mm", "x_head_mm")
        ),
        "actual_head_symmetric_difference_mm3": head_difference,
        "actual_head_embedding_volume_error_mm3": head_embedding_error,
        "volume_delta_actual_mm3": volume_delta,
        "volume_delta_expected_mm3": expected_delta,
        "volume_delta_error_mm3": abs(volume_delta - expected_delta),
        "diameter_scale": diameter_scale,
        "head_radius_scale_actual": p_c["r_head_mm"] / p_a["r_head_mm"],
        "head_axial_scale_actual": p_c["x_head_mm"] / p_a["x_head_mm"],
        "head_volume_scale_actual": head_c.Volume / head_a.Volume,
        "head_volume_scale_expected": diameter_scale**3,
    }
    if not checks["same_head_parameters"]:
        raise RuntimeError("FAIL: changing L altered head parameters")
    if head_difference > 1.0e-12 or head_embedding_error > 1.0e-10:
        raise RuntimeError("FAIL: actual A/B head BReps are not geometrically identical")
    if checks["volume_delta_error_mm3"] > 1.0e-10:
        raise RuntimeError("FAIL: L-only volume change is not the added cylinder volume")
    if abs(checks["head_radius_scale_actual"] - diameter_scale) > 1.0e-12:
        raise RuntimeError("FAIL: D does not scale head radius proportionally")
    if abs(checks["head_axial_scale_actual"] - diameter_scale) > 1.0e-12:
        raise RuntimeError("FAIL: D does not scale head axial extent proportionally")
    if abs(checks["head_volume_scale_actual"] - diameter_scale**3) > 1.0e-10:
        raise RuntimeError("FAIL: actual head solid does not scale geometrically with D")
    return checks


def _add_model_properties(obj, p):
    properties = (
        ("LTotal", "Length", p["length_mm"]),
        ("DBody", "Length", p["diameter_mm"]),
        ("RBody", "Length", p["r_body_mm"]),
        ("RHead", "Length", p["r_head_mm"]),
        ("HeadAxialExtent", "Length", p["x_head_mm"]),
        ("StraightBodyLength", "Length", p["straight_length_mm"]),
        ("HeadRadiusRatio", "Float", p["head_radius_ratio"]),
    )
    for name, kind, value in properties:
        obj.addProperty(f"App::Property{kind}", name, "Robot Geometry")
        setattr(obj, name, value)


def export_model(shape, p, stem):
    """Write FCStd, STEP, BREP, and visualization-only STL."""
    validate_shape(shape, p["length_mm"], p["diameter_mm"], stem)
    paths = {
        "fcstd": OUTPUT_DIR / f"{stem}.FCStd",
        "step": OUTPUT_DIR / f"{stem}.step",
        "brep": OUTPUT_DIR / f"{stem}.brep",
        "stl": OUTPUT_DIR / f"{stem}.stl",
    }
    doc = App.newDocument(stem)
    obj = doc.addObject("PartDesign::Feature", "Robot")
    obj.Label = "Parametric Robot"
    obj.Shape = shape
    _add_model_properties(obj, p)
    doc.recompute()
    doc.saveAs(str(paths["fcstd"]))
    Part.export([obj], str(paths["step"]))
    shape.exportBrep(str(paths["brep"]))
    try:
        import MeshPart

        mesh = MeshPart.meshFromShape(
            Shape=shape,
            LinearDeflection=STL_LINEAR_DEFLECTION_MM,
            AngularDeflection=math.radians(STL_ANGULAR_DEFLECTION_DEG),
            Relative=False,
        )
        mesh.write(str(paths["stl"]))
    finally:
        App.closeDocument(doc.Name)
    return paths


def validate_exports(paths, expected_length_mm, expected_diameter_mm):
    """Reimport authoritative outputs and repeat validity/dimension gates."""
    step_shape = Part.read(str(paths["step"]))
    step_result = validate_shape(
        step_shape, expected_length_mm, expected_diameter_mm, "STEP reimport"
    )
    brep_shape = Part.Shape()
    brep_shape.read(str(paths["brep"]))
    brep_result = validate_shape(
        brep_shape, expected_length_mm, expected_diameter_mm, "BREP reimport"
    )

    doc = App.openDocument(str(paths["fcstd"]))
    try:
        objects = [
            obj
            for obj in doc.Objects
            if hasattr(obj, "Shape") and not obj.Shape.isNull()
        ]
        if len(objects) != 1:
            raise RuntimeError(
                f"ABORT: FCStd reopen expected one shape object, found {len(objects)}"
            )
        fcstd_result = validate_shape(
            objects[0].Shape,
            expected_length_mm,
            expected_diameter_mm,
            "FCStd reopen",
        )
    finally:
        App.closeDocument(doc.Name)
    return {"step": step_result, "brep": brep_result, "fcstd": fcstd_result}


def _sample_upper_profile(p, count=500):
    import numpy as np

    theta = np.linspace(p["theta_tip_rad"], p["theta_tangent_rad"], count)
    x_arc = p["center"][0] + p["r_head_mm"] * np.cos(theta)
    y_arc = p["center"][1] + p["r_head_mm"] * np.sin(theta)
    x = np.concatenate((x_arc, [p["length_mm"]]))
    y = np.concatenate((y_arc, [p["r_body_mm"]]))
    return x, y


def create_profile_plots(target_p, reference_p):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams["font.family"] = "DejaVu Sans"
    fig, ax = plt.subplots(figsize=(12, 5.5), dpi=180)
    for p, color, label, text_position in (
        (
            reference_p,
            "#31688e",
            "Reference 3.050 x 1.200 mm",
            (1.88, 0.39),
        ),
        (target_p, "#d1495b", "Target 1.800 x 0.815 mm", (0.72, 0.12)),
    ):
        x, y = _sample_upper_profile(p)
        ax.plot(x, y, color=color, linewidth=2.0, label=label)
        ax.plot(x, -y, color=color, linewidth=2.0)
        ax.plot([p["length_mm"], p["length_mm"]], [-p["r_body_mm"], p["r_body_mm"]], color=color, linewidth=2.0)
        ax.scatter([0, p["x_head_mm"]], [0, p["r_body_mm"]], color=color, s=22)
        ax.text(
            text_position[0],
            text_position[1],
            f"L={p['length_mm']:.3f}, D={p['diameter_mm']:.3f}\n"
            f"Rbody={p['r_body_mm']:.6f}, Rhead={p['r_head_mm']:.6f}\n"
            f"xhead={p['x_head_mm']:.6f} mm",
            color=color,
            fontsize=8,
            ha="center",
            va="center",
            bbox={"facecolor": "white", "edgecolor": color, "alpha": 0.82, "pad": 3},
        )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Global X / axial position (mm)")
    ax.set_ylabel("Radial coordinate (mm)")
    ax.set_title("Exact circular-arc axial profiles", pad=12)
    ax.grid(True, color="#d7d7d7", linewidth=0.6)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "robot_profile_reference_vs_target.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6), dpi=180)
    for p, color, label in (
        (reference_p, "#31688e", "Reference nose"),
        (target_p, "#d1495b", "Target nose"),
    ):
        theta = np.linspace(p["theta_tip_rad"], p["theta_tangent_rad"], 600)
        x = p["center"][0] + p["r_head_mm"] * np.cos(theta)
        y = p["center"][1] + p["r_head_mm"] * np.sin(theta)
        ax.plot(x / p["r_body_mm"], y / p["r_body_mm"], color=color, linewidth=2.4, label=label)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x / R_body")
    ax.set_ylabel("y / R_body")
    ax.set_title("Diameter-normalized nose overlay")
    ax.grid(True, color="#d7d7d7", linewidth=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "robot_head_normalized_overlay.png", dpi=300)
    plt.close(fig)


def create_brep_preview(shape):
    """Render an orthographic isometric preview from FreeCAD's BRep tessellation."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    vertices, facets = shape.tessellate(0.0025)
    triangles = [
        [
            (vertices[index].x, vertices[index].y, vertices[index].z)
            for index in facet
        ]
        for facet in facets
    ]
    fig = plt.figure(figsize=(10, 6), dpi=180, facecolor="#f3f4f5")
    ax = fig.add_subplot(111, projection="3d", facecolor="#f3f4f5")
    collection = Poly3DCollection(
        triangles,
        facecolors="#b9c5c9",
        linewidth=0.0,
        alpha=1.0,
        shade=True,
        lightsource=LightSource(azdeg=300, altdeg=45),
    )
    ax.add_collection3d(collection)
    box = shape.BoundBox
    ax.set_xlim(box.XMin, box.XMax)
    ax.set_ylim(box.YMin, box.YMax)
    ax.set_zlim(box.ZMin, box.ZMax)
    ax.set_box_aspect((box.XLength, box.YLength, box.ZLength))
    ax.set_proj_type("ortho")
    ax.view_init(elev=24, azim=-58)
    ax.set_axis_off()
    fig.tight_layout(pad=0)
    fig.savefig(
        OUTPUT_DIR / f"{TARGET_STEM}_preview.png",
        dpi=300,
        facecolor=fig.get_facecolor(),
        bbox_inches="tight",
        pad_inches=0.08,
    )
    plt.close(fig)


def build_geometry_json(shape, p, reference_shape, export_validation, parameter_checks):
    dims = _bbox_dimensions(shape)
    return {
        "geometry_family": "axisymmetric circular-arc nose with straight cylindrical body and flat tail",
        "authoritative_geometry": "FreeCAD Part/OpenCASCADE BRep revolved from exact 2D edges",
        "axis": "+X",
        "length_definition": "nose tip to flat tail",
        "diameter_definition": "straight cylindrical outer diameter",
        "L_total_mm": p["length_mm"],
        "D_body_mm": p["diameter_mm"],
        "R_body_mm": p["r_body_mm"],
        "head_radius_ratio": p["head_radius_ratio"],
        "R_head_mm": p["r_head_mm"],
        "head_axial_extent_mm": p["x_head_mm"],
        "straight_body_length_mm": p["straight_length_mm"],
        "volume_mm3": shape.Volume,
        "surface_area_mm2": shape.Area,
        "bbox_x_mm": dims["x"],
        "bbox_y_mm": dims["y"],
        "bbox_z_mm": dims["z"],
        "solid_valid": shape.isValid(),
        "closed": shape.isClosed(),
        "solid_count": len(shape.Solids),
        "shell_count": len(shape.Shells),
        "volume_ratio_to_SW_reference": shape.Volume / reference_shape.Volume,
        "density_tonne_per_mm3": DENSITY_TONNE_PER_MM3,
        "estimated_mass_mg": shape.Volume * DENSITY_TONNE_PER_MM3 * 1.0e9,
        "stl_role": "visualization only; not a dimensional source",
        "stl_linear_deflection_mm": STL_LINEAR_DEFLECTION_MM,
        "stl_angular_deflection_deg": STL_ANGULAR_DEFLECTION_DEG,
        "export_reimport_validation": export_validation,
        "parameterization_validation": parameter_checks,
        "freecad_version": ".".join(App.Version()),
    }


def _write_json(path, value):
    with open(path, "w", encoding="ascii", newline="\n") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    target_p = profile_definition(L_TOTAL_MM, D_BODY_MM)
    reference_p = profile_definition(REFERENCE_LENGTH_MM, REFERENCE_DIAMETER_MM)
    target = make_robot(L_TOTAL_MM, D_BODY_MM)
    reference = make_robot(REFERENCE_LENGTH_MM, REFERENCE_DIAMETER_MM)
    target_gate = validate_shape(target, L_TOTAL_MM, D_BODY_MM, "target")
    reference_gate = validate_shape(
        reference, REFERENCE_LENGTH_MM, REFERENCE_DIAMETER_MM, "reference"
    )
    parameter_checks = verify_parameterization()
    target_paths = export_model(target, target_p, TARGET_STEM)
    reference_paths = export_model(reference, reference_p, REFERENCE_STEM)
    target_exports = validate_exports(target_paths, L_TOTAL_MM, D_BODY_MM)
    reference_exports = validate_exports(
        reference_paths, REFERENCE_LENGTH_MM, REFERENCE_DIAMETER_MM
    )

    # Meshing populates a discretized bounding-box cache in FreeCAD. Rebuild fresh
    # BReps so authoritative JSON and mass properties remain exact BRep measurements.
    report_target = make_robot(L_TOTAL_MM, D_BODY_MM)
    report_reference = make_robot(REFERENCE_LENGTH_MM, REFERENCE_DIAMETER_MM)
    validate_shape(report_target, L_TOTAL_MM, D_BODY_MM, "target report BRep")
    validate_shape(
        report_reference,
        REFERENCE_LENGTH_MM,
        REFERENCE_DIAMETER_MM,
        "reference report BRep",
    )

    geometry = build_geometry_json(
        report_target,
        target_p,
        report_reference,
        target_exports,
        parameter_checks,
    )
    geometry["reference_geometry"] = {
        "L_total_mm": reference_p["length_mm"],
        "D_body_mm": reference_p["diameter_mm"],
        "R_body_mm": reference_p["r_body_mm"],
        "R_head_mm": reference_p["r_head_mm"],
        "volume_mm3": report_reference.Volume,
        "validation": reference_gate,
        "export_reimport_validation": reference_exports,
    }
    _write_json(OUTPUT_DIR / f"{TARGET_STEM}_geometry.json", geometry)
    props = mass_properties(report_target)
    props["model"] = TARGET_STEM
    props["axis"] = "+X"
    props["note"] = "CAD theoretical properties; audit again after later Abaqus import"
    _write_json(OUTPUT_DIR / f"{TARGET_STEM}_mass_properties.json", props)
    create_profile_plots(target_p, reference_p)
    create_brep_preview(report_target)

    summary = {
        "freecad_version": ".".join(App.Version()),
        "target_gate": target_gate,
        "target_parameters": target_p,
        "mass_properties": props,
        "volume_ratio_to_SW_reference": report_target.Volume / report_reference.Volume,
        "parameterization_validation": parameter_checks,
        "outputs": {key: str(path) for key, path in target_paths.items()},
    }
    print("BUILD_RESULT=" + json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
