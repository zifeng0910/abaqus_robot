from odbAccess import openOdb
import sys
o=openOdb(sys.argv[1]+'.odb',readOnly=True)
s=list(o.steps.values())[-1]
print(len(s.frames), sorted(s.frames[-1].fieldOutputs.keys()))
o.close()
