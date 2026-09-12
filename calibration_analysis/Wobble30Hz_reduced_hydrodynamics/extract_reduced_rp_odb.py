from odbAccess import openOdb
import csv, sys
job=sys.argv[1]; out=sys.argv[2]
o=openOdb(job+'.odb',readOnly=True); st=list(o.steps.values())[-1]
rp=o.rootAssembly.nodeSets['RP_ROBOT']; rows=[]
for fr in st.frames:
    row={'time_s':float(fr.frameValue)}
    for name,var,n in [('U','U',3),('UR','UR',3),('V','V',3),('VR','VR',3),('A','A',3),('AR','AR',3)]:
        if var not in fr.fieldOutputs: continue
        vals=fr.fieldOutputs[var].getSubset(region=rp).values
        if not vals: continue
        dat=vals[0].data
        for i in range(min(n,len(dat))): row[name+str(i+1)]=float(dat[i])
    rows.append(row)
keys=['time_s']+sum(([nm+str(i) for i in range(1,4)] for nm in ('U','UR','V','VR','A','AR')),[])
with open(out,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=keys); w.writeheader()
    for r in rows: w.writerow({k:r.get(k,'') for k in keys})
print('WROTE',out,'frames',len(rows),'fields',sorted(set(k.rstrip('123') for r in rows for k in r if k!='time_s')))
o.close()
