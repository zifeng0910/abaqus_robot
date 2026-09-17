from __future__ import print_function

import sys
from odbAccess import openOdb


odb = openOdb(sys.argv[1], readOnly=True)
assembly = odb.rootAssembly
for name in sorted(assembly.elementSets.keys()):
    if "INITPENET" in name.upper() or "WARN" in name.upper():
        region = assembly.elementSets[name]
        print(name, [(block[0].instanceName, len(block)) for block in region.elements if block])
for instance_name, instance in assembly.instances.items():
    for name in sorted(instance.elementSets.keys()):
        if "INITPENET" in name.upper() or "WARN" in name.upper():
            print(instance_name, name, len(instance.elementSets[name].elements))
odb.close()
