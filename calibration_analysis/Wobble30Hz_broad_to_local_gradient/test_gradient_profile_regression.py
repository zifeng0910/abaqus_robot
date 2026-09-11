"""Small offline regression for broad-to-local profile timing."""
import math

def weight(t, start, duration):
    x=max(0.0,min(1.0,(t-start)/duration))
    return x*x*(3.0-2.0*x)

def force(gb, gl, lb, ll, delta, t, start=0.0, duration=0.0005, mode='legacy'):
    eb=math.exp(-0.5*(delta/lb)**2); fb=-(delta/(lb*lb))*eb*gb
    if mode=='legacy': return fb
    el=math.exp(-0.5*(delta/ll)**2); fl=-(delta/(ll*ll))*el*gl
    w=weight(t,start,duration)
    return (1-w)*fb+w*fl

# Legacy mode must remain exactly the broad profile, independent of gated args.
for t in (0.0, 0.0002, 0.001, 0.01):
    a=force(0.006,0.006,45.0,3.75,5.1,t,mode='legacy')
    b=-(5.1/(45.0*45.0))*math.exp(-0.5*(5.1/45.0)**2)*0.006
    assert a == b, (a,b)
# Cross-fade endpoints and monotonic bounded weight.
assert weight(0.0017,0.0017,0.0005) == 0.0
assert weight(0.0022,0.0017,0.0005) == 1.0
assert 0.0 < weight(0.00195,0.0017,0.0005) < 1.0
for t in [0.0017,0.0018,0.0019,0.0020,0.0022]:
    assert force(0.006,0.006,45.0,3.75,5.1,t,start=0.0017) < 0.0
print('PASS: legacy profile exact and broad-to-local endpoints/force sign valid')
