"""FreeCAD/OpenCASCADE tests for the parametric robot."""

import math
import sys
import unittest
from pathlib import Path

import FreeCAD as App
import Part

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_parametric_robot as robot


class ParametricRobotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p = robot.geometry_parameters(1.8, 0.815)
        cls.shape = robot.make_robot(1.8, 0.815)

    def test_target_length(self):
        self.assertAlmostEqual(self.shape.BoundBox.XLength, 1.8, places=12)

    def test_target_diameter(self):
        self.assertAlmostEqual(self.shape.BoundBox.YLength, 0.815, places=12)
        self.assertAlmostEqual(self.shape.BoundBox.ZLength, 0.815, places=12)

    def test_body_radius(self):
        self.assertAlmostEqual(self.p["r_body_mm"], 0.4075, places=12)

    def test_head_radius_and_ratio(self):
        self.assertAlmostEqual(
            self.p["head_radius_ratio"], 0.61 / 0.60, places=14
        )
        self.assertAlmostEqual(
            self.p["r_head_mm"], 0.4075 * 0.61 / 0.60, places=14
        )

    def test_shape_valid_closed_single_solid(self):
        result = robot.validate_shape(self.shape, 1.8, 0.815, "unit target")
        self.assertTrue(result["valid"])
        self.assertTrue(result["closed"])
        self.assertEqual(result["solid_count"], 1)
        self.assertEqual(result["shell_count"], 1)
        self.assertGreater(self.shape.Volume, 0.0)

    def test_length_does_not_change_actual_head_geometry(self):
        shape_a = robot.make_robot(1.8, 0.815)
        shape_b = robot.make_robot(2.2, 0.815)
        head_a = robot.make_head(0.815)
        head_b = robot.make_head(0.815)
        symmetric_difference = head_a.cut(head_b).Volume + head_b.cut(head_a).Volume
        self.assertLessEqual(symmetric_difference, 1.0e-12)
        self.assertAlmostEqual(shape_a.common(head_a).Volume, head_a.Volume, places=12)
        self.assertAlmostEqual(shape_b.common(head_b).Volume, head_b.Volume, places=12)

    def test_length_only_volume_change_is_cylindrical(self):
        shape_a = robot.make_robot(1.8, 0.815)
        shape_b = robot.make_robot(2.2, 0.815)
        expected = math.pi * (0.815 / 2.0) ** 2 * (2.2 - 1.8)
        self.assertAlmostEqual(shape_b.Volume - shape_a.Volume, expected, places=12)

    def test_diameter_scales_actual_head(self):
        p_a = robot.geometry_parameters(1.8, 0.815)
        p_c = robot.geometry_parameters(1.8, 0.9)
        head_a = robot.make_head(0.815)
        head_c = robot.make_head(0.9)
        scale = 0.9 / 0.815
        self.assertAlmostEqual(p_c["r_head_mm"] / p_a["r_head_mm"], scale, places=13)
        self.assertAlmostEqual(p_c["x_head_mm"] / p_a["x_head_mm"], scale, places=13)
        self.assertAlmostEqual(head_c.Volume / head_a.Volume, scale**3, places=12)
        self.assertAlmostEqual(head_c.Area / head_a.Area, scale**2, places=12)

    def test_invalid_short_body_is_rejected(self):
        x_head = self.p["x_head_mm"]
        with self.assertRaises(ValueError):
            robot.geometry_parameters(x_head - 1.0e-5, 0.815)

    def test_exported_step_reimport(self):
        path = Path(__file__).resolve().parent / f"{robot.TARGET_STEM}.step"
        self.assertTrue(path.is_file())
        imported = Part.read(str(path))
        result = robot.validate_shape(imported, 1.8, 0.815, "test STEP")
        self.assertEqual(result["solid_count"], 1)

    def test_saved_fcstd_reopen(self):
        path = Path(__file__).resolve().parent / f"{robot.TARGET_STEM}.FCStd"
        self.assertTrue(path.is_file())
        doc = App.openDocument(str(path))
        try:
            objects = [
                obj
                for obj in doc.Objects
                if hasattr(obj, "Shape") and not obj.Shape.isNull()
            ]
            self.assertEqual(len(objects), 1)
            robot.validate_shape(objects[0].Shape, 1.8, 0.815, "test FCStd")
        finally:
            App.closeDocument(doc.Name)


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]], verbosity=2)
