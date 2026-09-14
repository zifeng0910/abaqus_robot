"""Focused numerical checks for L2300 post-processing conventions."""
import numpy as np

from analyze_l2300_dynamic import intervals, local_frames, unit


def main():
    assert np.allclose(unit([[3., 0., 0.]])[0], [1., 0., 0.])
    tangent = np.tile([1., 0., 0.], (4, 1))
    e1, e2 = local_frames(tangent)
    assert np.allclose(np.einsum("ij,ij->i", tangent, e1), 0)
    assert np.allclose(np.cross(tangent, e1), e2)
    t = np.arange(6)*1e-7
    found = intervals(t, np.array([0, 1, 1, 0, 1, 0], bool))
    assert found[0][:2] == (1, 2) and abs(found[0][4]-2e-7) < 1e-15
    assert found[1][:2] == (4, 4)
    print("PASS: dynamic-analysis numerical conventions")


if __name__ == "__main__":
    main()
