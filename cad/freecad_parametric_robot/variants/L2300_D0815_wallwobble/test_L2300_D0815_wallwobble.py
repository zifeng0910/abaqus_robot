"""FreeCAD/OpenCASCADE tests for the frozen L2300 D0815 CAD."""

from pathlib import Path
import json
import math
import sys
import unittest

import FreeCAD as App
import Part

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_L2300_D0815_wallwobble as cad


class FrozenCADTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.shape = cad.make_robot()
        cls.short = cad.make_robot(1.800)

    def test_01_dimensions(self):
        self.assertAlmostEqual(self.shape.BoundBox.XLength, 2.300, places=12)
        self.assertAlmostEqual(self.shape.BoundBox.YLength, 0.815, places=12)
        self.assertAlmostEqual(self.shape.BoundBox.ZLength, 0.815, places=12)

    def test_02_valid_closed_single_solid(self):
        gate = cad.validate_shape(self.shape)
        self.assertTrue(gate["valid"] and gate["closed"])
        self.assertEqual(gate["solid_count"], 1)
        self.assertEqual(gate["shell_count"], 1)
        self.assertGreater(self.shape.Volume, 0.0)

    def test_03_frozen_head_profile(self):
        self.assertEqual(len(cad.HEAD_POLES_MM)-1, 5)
        self.assertEqual(cad.HEAD_POLES_MM[-1], (cad.HEAD_AXIAL_EXTENT_MM, cad.R_BODY_MM))
        samples = [cad.bezier_point(i/1000.0) for i in range(1001)]
        self.assertTrue(all(b[0] >= a[0] for a, b in zip(samples, samples[1:])))
        self.assertTrue(all(b[1] >= a[1] for a, b in zip(samples, samples[1:])))

    def test_04_head_independent_of_length(self):
        head = cad.make_head()
        self.assertAlmostEqual(self.shape.common(head).Volume, head.Volume, places=12)
        self.assertAlmostEqual(self.short.common(head).Volume, head.Volume, places=12)

    def test_05_length_change_is_exact_cylinder(self):
        self.assertAlmostEqual(self.shape.BoundBox.XLength-self.short.BoundBox.XLength, 0.500, places=12)
        expected = math.pi*cad.R_BODY_MM**2*0.500
        self.assertAlmostEqual(self.shape.Volume-self.short.Volume, expected, places=12)

    def test_06_step_roundtrip(self):
        shape = Part.read(str(cad.HERE/(cad.STEM+".step")))
        cad.validate_shape(shape, label="test STEP")

    def test_07_brep_roundtrip(self):
        shape = Part.Shape()
        shape.read(str(cad.HERE/(cad.STEM+".brep")))
        cad.validate_shape(shape, label="test BREP")

    def test_08_fcstd_reopen(self):
        doc = App.openDocument(str(cad.HERE/(cad.STEM+".FCStd")))
        try:
            objects = [x for x in doc.Objects if hasattr(x, "Shape") and not x.Shape.isNull()]
            self.assertEqual(len(objects), 1)
            cad.validate_shape(objects[0].Shape, label="test FCStd")
        finally:
            App.closeDocument(doc.Name)

    def test_09_mass_COM_inertia(self):
        props = cad.mass_properties(self.shape)
        self.assertTrue(props["center_of_mass"]["inside_solid"])
        self.assertTrue(props["principal_mass_moments"]["positive_definite"])

    def test_10_magnetic_volume_scaling(self):
        data = json.loads((cad.HERE/(cad.STEM+"_magnetic_identity.json")).read_text())
        expected = (data["reference_magnetic_moment_Am2"]
                    * data["authoritative_CAD_volume_mm3"]/data["reference_volume_mm3"])
        self.assertAlmostEqual(data["physical_magnetic_moment_target_Am2"], expected, places=15)
        self.assertFalse(data["server_started"])

    def test_11_previous_CAD_regression(self):
        data = json.loads((cad.HERE/"previous_generated_CAD_regression.json").read_text())
        self.assertTrue(data["geometrically_equivalent_within_roundoff"])


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]], verbosity=2)
