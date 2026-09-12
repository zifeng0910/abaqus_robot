from odbAccess import openOdb
import csv,sys
job=sys.argv[1]; out=sys.argv[2]
o=openOdb(job+'.odb',readOnly=True); st=list(o.steps.values())[-1]
with open(out,'w',newline='') as f:
 w=csv.writer(f); w.writerow(['time_s','CPRESS_max_MPa','CPRESS_active_values','CNORMF_max_N','CSHEARF_max_N'])
 for fr in st.frames:
  def vals(prefix):
   a=[]
   for k,fo in fr.fieldOutputs.items():
    if k.startswith(prefix):
     for v in fo.values:
      try:
       d=v.data; a.append(float(d if not hasattr(d,'__len__') else max(abs(float(x)) for x in d)))
      except: pass
   return a
  cp=vals('CPRESS'); nf=vals('CNORMF'); sf=vals('CSHEARF')
  w.writerow([float(fr.frameValue),max(cp) if cp else 0.,sum(x>1e-12 for x in cp),max(nf) if nf else 0.,max(sf) if sf else 0.])
o.close(); print('WROTE',out)
