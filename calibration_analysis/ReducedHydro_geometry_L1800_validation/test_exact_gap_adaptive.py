"""Regression: adaptive triangle candidates must equal exhaustive search."""
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent / "ReducedHydro_hidden_impact_audit"
sys.path.insert(0, str(AUDIT))
from exact_gap_audit import Wall, tests


def main():
    assert tests().startswith("PASS")
    wall = Wall()
    rng = np.random.default_rng(7)
    points = wall.centers[rng.integers(0, len(wall.centers), 100)] + rng.normal(0, .1, (100, 3))
    gap, triangle, closest = wall.query(points)
    reference_gap, reference_triangle, reference_closest = wall.query(points, all_faces=True)
    assert np.allclose(gap, reference_gap, rtol=0, atol=1e-12)
    assert np.array_equal(triangle, reference_triangle)
    assert np.allclose(closest, reference_closest, rtol=0, atol=1e-12)
    print("PASS: adaptive exact-gap search equals exhaustive triangle search")


if __name__ == "__main__":
    main()
