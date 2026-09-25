"""Extract prescribed RP reactions and pair-resolved contact history from Abaqus."""
from __future__ import print_function

import json
import sys
from pathlib import Path

import numpy as np
from odbAccess import openOdb

job, destination = sys.argv[1:3]
out = Path(destination)
out.mkdir(parents=True, exist_ok=True)
odb = openOdb(job + ".odb", readOnly=True)
try:
    step = list(odb.steps.values())[-1]
    rp, contact, energy, inventory = {}, {}, {}, []
    rp_names = {stem + str(axis) for stem in ("U", "UR", "V", "VR", "A", "AR", "RF", "RM", "CF")
                for axis in (1, 2, 3)}
    for region_name, region in step.historyRegions.items():
        for variable, history in region.historyOutputs.items():
            data = np.asarray(history.data, dtype=float)
            inventory.append((region_name, variable, len(data)))
            if variable in rp_names:
                current = rp.get(variable)
                if current is None or len(data) > len(current):
                    rp[variable] = data
            if variable.startswith(("CFN", "CFS", "CFT")):
                contact[region_name + "|" + variable] = data
            if variable.startswith(("ALL", "ETOTAL")):
                current = energy.get(variable)
                if current is None or len(data) > len(current):
                    energy[variable] = data
    np.savez_compressed(out / "rp_replay_private.npz", **rp)
    np.savez_compressed(out / "contact_replay_private.npz", **contact)
    np.savez_compressed(out / "energy_replay_private.npz", **energy)
    (out / "replay_odb_inventory.json").write_text(json.dumps({
        "job": job, "history": inventory, "frames": len(step.frames)
    }, indent=2) + "\n")
    print("EXTRACTED", job, "RP", sorted(rp), "CONTACT", len(contact), "ENERGY", sorted(energy))
finally:
    odb.close()
