"""Small independent tests: projection, tetra mass and impulse units/factors."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from exact_gap_audit import closest, tests, Wall

assert tests().startswith('PASS')
t=np.linspace(0,.001,101);F=np.broadcast_to([1.,-2.,3.],(len(t),3))
J=cumulative_trapezoid(F,t,axis=0,initial=0)
assert np.allclose(J[-1],[.001,-.002,.003],atol=1e-15)
m_tonne=1e-8; v=J[-1]/m_tonne
assert np.allclose(1e-5*(v*.001),J[-1]) # SI momentum equals N s
tetra=np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1.]])
assert abs(np.linalg.det(tetra[1:]-tetra[0])/6-1/6)<1e-15
assert np.allclose(tetra.mean(axis=0),[.25,.25,.25])
wall=Wall()
# Deterministic close-to-wall samples, candidate exclusion vs all triangles.
p=wall.centers[::23]+wall.normals[::23]*.001
g,i,q=wall.query(p);gg,ii,qq=wall.query(p,all_faces=True)
assert np.max(np.abs(g-gg))<1e-12
assert np.max(np.linalg.norm(q-qq,axis=1))<1e-12
print('PASS: triangle interior/edge/vertex, exact candidate bound, trapezoid scale, SI/Abaqus momentum, tetra mass/COM')
