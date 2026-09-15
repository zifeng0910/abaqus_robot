"""Run with Abaqus Python to extract dense replay reaction histories."""
from odbAccess import openOdb
from pathlib import Path
import numpy as np
import sys

job, destination = sys.argv[1:3]
out = Path(destination)
odb = openOdb(job + ".odb", readOnly=True)
step = list(odb.steps.values())[-1]
wanted = {prefix + str(i) for prefix in ("RM", "RF", "CM", "CF", "A", "AR") for i in (1, 2, 3)}
result = {}
for region in step.historyRegions.values():
    for name, history in region.historyOutputs.items():
        if name in wanted:
            array = np.asarray(history.data, dtype=float)
            if name not in result or len(array) > len(result[name]):
                result[name] = array
np.savez_compressed(out / "reaction_history_private.npz", **result)
odb.close()
print("EXTRACTED", sorted(result))
