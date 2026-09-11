import sys
from odbAccess import openOdb
for job in sys.argv[1:]:
    odb=openOdb(job+'.odb',readOnly=True)
    print('JOB',job)
    step=odb.steps[list(odb.steps.keys())[0]]
    print('history regions',len(step.historyRegions))
    for name,reg in step.historyRegions.items():
        keys=list(reg.historyOutputs.keys())
        if any(k in keys for k in ('V1','V2','V3','U1','U2','U3','VR1','VR2','VR3')):
            print(name,keys)
    odb.close()
