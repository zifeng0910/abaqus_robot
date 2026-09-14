"""Abaqus/CAE noGUI script: create global-0.10 mm C3D4 screening meshes."""

from abaqus import mdb
from abaqusConstants import C3D4, DEFORMABLE_BODY, FREE, OFF, STANDARD, TET, THREE_D
from caeModules import *
from mesh import ElemType
import csv
import json
import os
import time


def main():
    repo = os.path.abspath(os.getcwd())
    cad = os.path.join(repo, "cad", "freecad_parametric_robot", "screening", "coarse_motion_mode")
    out = os.path.join(repo, "calibration_analysis", "CoarseMotionModeScreen", "cases", "_meshes")
    if not os.path.isdir(out): os.makedirs(out)
    results = []
    opener = getattr(mdb, "openStep", None) or getattr(mdb, "openSTEP", None)
    for length in (2.00, 2.10, 2.20, 2.30, 2.40):
        code = "%04d" % round(length*1000)
        stem = "Robot_SCREENING_L%s_D0815" % code
        started = time.time()
        model_name = "SCREENING_L%s" % code
        model = mdb.Model(name=model_name)
        geometry = opener(os.path.join(cad, stem+".step"), scaleFromFile=OFF)
        part = model.PartFromGeometryFile(name="Robot_SOLID", geometryFile=geometry,
            combine=False, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        part.seedPart(size=.100, deviationFactor=.05, minSizeFactor=.50)
        part.setMeshControls(regions=part.cells, elemShape=TET, technique=FREE)
        part.setElementType(regions=(part.cells,), elemTypes=(ElemType(elemCode=C3D4, elemLibrary=STANDARD),))
        part.generateMesh()
        with open(os.path.join(out, stem+"_nodes.csv"), "w", newline="") as handle:
            writer=csv.writer(handle); writer.writerow(("node","x_mm","y_mm","z_mm"))
            for node in part.nodes: writer.writerow([node.label]+list(node.coordinates))
        with open(os.path.join(out, stem+"_elements.csv"), "w", newline="") as handle:
            writer=csv.writer(handle); writer.writerow(("element","n1","n2","n3","n4"))
            for element in part.elements:
                writer.writerow([element.label]+[part.nodes[i].label for i in element.connectivity])
        result={"length_mm":length,"mesh_size_mm":.100,"nodes":len(part.nodes),
                "elements":len(part.elements),"element_type":"C3D4",
                "wallclock_s":time.time()-started}
        if not (5000 <= len(part.elements) < 20000):
            raise RuntimeError("Element-count gate failed: %r" % result)
        results.append(result)
        del mdb.models[model_name]
        print(json.dumps(result), flush=True)
    with open(os.path.join(out,"coarse_mesh_manifest.json"),"w") as handle:
        json.dump(results,handle,indent=2)


if __name__ == "__main__": main()
