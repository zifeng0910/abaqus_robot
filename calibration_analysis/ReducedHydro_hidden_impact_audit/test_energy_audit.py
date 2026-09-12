"""Deterministic tests for tetra inertia, energy units and exported identities."""
import numpy as np
import pandas as pd

from analyze_energy_audit import HERE, exact_tetra_inertia


# Unit-tetra continuum inertia is checked against its closed-form moments
# (rho=1, mass=1/6).
verts = np.array([[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]])
mesh = {'verts': verts, 'vol': np.array([1. / 6.]), 'mass': 1. / 6., 'com': np.array([.25, .25, .25])}
inertia = exact_tetra_inertia(mesh)
assert np.allclose(np.diag(inertia), [1. / 80.] * 3)
assert np.allclose(inertia[np.triu_indices(3, 1)], [1. / 480.] * 3)

# tonne*(mm/s)^2 is N mm, and 1 N mm = 1e-3 J.
assert np.isclose(0.5 * 1e-8 * 1000.0**2 * 1e-3, 5e-6)

w = pd.read_csv(HERE / 'first_impact_energy_windows.csv')
assert np.max(np.abs(w.energy_identity_error_J)) < 1e-18
a = pd.read_csv(HERE / 'first_impact_abaqus_energy_crosscheck.csv').iloc[0]
assert abs(a.identity_error_J) < 1e-11
assert abs(a.dense_contact_dissipation_J - a.contact_dissipation_effective_J) < 1e-11
print('PASS: exact tetra inertia, N-mm/J conversion, dense work identity, Abaqus energy identity')
