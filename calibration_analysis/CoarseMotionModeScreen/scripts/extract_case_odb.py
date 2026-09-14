"""Abaqus Python read-only extraction of dense RP/contact history to private NPZ."""

from odbAccess import openOdb
from pathlib import Path
import json
import numpy as np
import sys


job, destination = sys.argv[1:3]
out=Path(destination); out.mkdir(parents=True,exist_ok=True)
odb=openOdb(job+".odb",readOnly=True); step=list(odb.steps.values())[-1]
rp={}; contact={}; inventory=[]
for region_name,region in step.historyRegions.items():
    for variable,history in region.historyOutputs.items():
        array=np.asarray(history.data,dtype=float)
        inventory.append((region_name,variable,len(array)))
        if variable in [p+str(i) for p in ("U","UR","V","VR") for i in (1,2,3)]:
            current=rp.get(variable)
            if current is None or len(array)>len(current): rp[variable]=array
        if variable.startswith(("CFN","CFS","CFT")): contact[region_name+"|"+variable]=array
np.savez_compressed(out/"rp_history_private.npz",**rp)
np.savez_compressed(out/"contact_history_private.npz",**contact)
(out/"odb_inventory.json").write_text(json.dumps({"history":inventory,"frames":len(step.frames)},indent=2))
odb.close()
print("EXTRACTED",job,"RP",sorted(rp),"CONTACT",len(contact))

