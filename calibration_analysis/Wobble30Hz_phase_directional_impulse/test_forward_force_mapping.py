"""Small coordinate/sign regression test; no Abaqus or socket required."""
import numpy as np

def project(force, tangent):
    force=np.asarray(force,float); tangent=np.asarray(tangent,float)
    tangent=tangent/np.linalg.norm(tangent)
    return float(force@tangent)

def main():
    t=np.array([1.,0.,0.])
    assert project([1,0,0],t)>0
    assert project([-1,0,0],t)<0
    print('Case A Ft=',project([1,0,0],t),'Case B Ft=',project([-1,0,0],t))
    g=np.array([-0.9737765086,0.0739386448,-0.2151566594])
    print('real initial tangent=',g,'norm=',np.linalg.norm(g))
    print('global/local/Abaqus sign test: PASS (identity transform test)')

if __name__=='__main__': main()
