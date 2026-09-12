from odbAccess import openOdb
import csv,sys,math
job=sys.argv[1]
o=openOdb(job+'.odb',readOnly=True); st=list(o.steps.values())[-1]
rows=[]
for fr in st.frames:
    vals=[]
    for key in fr.fieldOutputs.keys():
        if key.startswith('CPRESS'):
            try:
                for v in fr.fieldOutputs[key].values:
                    if getattr(getattr(v,'instance',None),'name','')=='ROBOT_SOLID-1':
                        try: vals.append(float(v.data))
                        except: pass
            except: pass
    rows.append((float(fr.frameValue), max(vals) if vals else 0.0, sum(x>1e-12 for x in vals)))
out=job+'_cpress_timeseries.csv'
with open(out,'w',newline='') as f:
    w=csv.writer(f); w.writerow(['time_s','CPRESS_max_MPa','active_robot_values']); w.writerows(rows)
print('WROTE',out,'frames',len(rows),'max',max(r[1] for r in rows))
o.close()
