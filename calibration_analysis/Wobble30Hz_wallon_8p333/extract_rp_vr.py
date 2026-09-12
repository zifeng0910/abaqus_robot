from odbAccess import openOdb
import csv,sys
job=sys.argv[1]; out=sys.argv[2]
o=openOdb(job+'.odb',readOnly=True); st=list(o.steps.values())[-1]; inst=o.rootAssembly.instances['ROBOT_SOLID-1']; rp=o.rootAssembly.nodeSets['RP_ROBOT']
with open(out,'w',newline='') as f:
 w=csv.writer(f); w.writerow(['time_s','VR1','VR2','VR3'])
 for fr in st.frames:
  vals=fr.fieldOutputs['VR'].getSubset(region=rp).values
  v=vals[0].data if vals else None
  if v is not None: w.writerow([float(fr.frameValue),float(v[0]),float(v[1]),float(v[2])])
o.close(); print('WROTE',out)
