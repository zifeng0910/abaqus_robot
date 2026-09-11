import csv, sys
from odbAccess import openOdb
for job in sys.argv[1:]:
    odb=openOdb(job+'.odb',readOnly=True); step=odb.steps[list(odb.steps.keys())[0]]
    reg=None
    for name,r in step.historyRegions.items():
        if 'Node ASSEMBLY.2' in name: reg=r; break
    if reg is None: raise RuntimeError('RP history region not found')
    out=job+'_rp_velocity_history.csv'
    keys=['V1','V2','V3','U1','U2','U3','VR1','VR2','VR3']
    data={k:dict(reg.historyOutputs[k].data) for k in keys}
    times=sorted(set().union(*[set(v) for v in data.values()]))
    with open(out,'w',newline='') as fh:
        w=csv.writer(fh); w.writerow(['time_s']+keys)
        for t in times: w.writerow([t]+[data[k].get(t,'') for k in keys])
    print('WROTE',out,len(times)); odb.close()
