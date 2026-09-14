"""Import a specified STEP in Abaqus/CAE and export one C3D4 mesh."""
from abaqus import mdb
from abaqusConstants import C3D4, DEFORMABLE_BODY, FREE, OFF, STANDARD, TET, THREE_D
from caeModules import *
from mesh import ElemType
import csv
import json
import os
import sys
import time


def main():
    step = os.path.abspath(sys.argv[-3])
    size = float(sys.argv[-2])
    prefix = sys.argv[-1]
    here = os.path.join(os.getcwd(), "calibration_analysis",
                        "ReducedHydro_L1800_head_clearance")
    started = time.time()
    model = mdb.Model(name="HeadClearanceMesh")
    opener = getattr(mdb, "openStep", None) or getattr(mdb, "openSTEP", None)
    geometry = opener(step, scaleFromFile=OFF)
    part = model.PartFromGeometryFile(name="Robot_SOLID", geometryFile=geometry,
                                      combine=False, dimensionality=THREE_D,
                                      type=DEFORMABLE_BODY)
    bbox = part.cells.getBoundingBox(); lo, hi = bbox["low"], bbox["high"]
    extents = [hi[i] - lo[i] for i in range(3)]
    if abs(extents[0] - 1.8) > 1e-6 or max(abs(extents[i] - .815) for i in (1, 2)) > 1e-6:
        raise RuntimeError("STEP extent gate failed: %r" % (extents,))
    part.seedPart(size=size, deviationFactor=0.05, minSizeFactor=0.05)
    part.setMeshControls(regions=part.cells, elemShape=TET, technique=FREE)
    part.setElementType(regions=(part.cells,), elemTypes=(ElemType(elemCode=C3D4, elemLibrary=STANDARD),))
    part.generateMesh()
    with open(os.path.join(here, prefix + "_nodes_local.csv"), "w", newline="") as stream:
        writer = csv.writer(stream); writer.writerow(["node", "x_mm", "y_mm", "z_mm"])
        for node in part.nodes: writer.writerow([node.label] + list(node.coordinates))
    with open(os.path.join(here, prefix + "_elements.csv"), "w", newline="") as stream:
        writer = csv.writer(stream); writer.writerow(["element", "n1", "n2", "n3", "n4"])
        for element in part.elements:
            writer.writerow([element.label] + [part.nodes[i].label for i in element.connectivity])
    result = {"source_step": step, "mesh_size_mm": size, "bbox_low_mm": list(lo),
              "bbox_high_mm": list(hi), "bbox_extents_mm": extents,
              "node_count": len(part.nodes), "element_count": len(part.elements),
              "preprocessing_wall_clock_s": time.time() - started,
              "abaqus_import": "openStep + PartFromGeometryFile + C3D4 FREE TET"}
    with open(os.path.join(here, prefix + "_raw.json"), "w") as stream:
        json.dump(result, stream, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
