"""Import the L2300 STEP in Abaqus/CAE and export the sole 0.060 mm C3D4 mesh."""
from abaqus import mdb
from abaqusConstants import C3D4, DEFORMABLE_BODY, FREE, OFF, STANDARD, TET, THREE_D
from caeModules import *
from mesh import ElemType
import csv
import json
import os
import time


def main():
    here = os.path.join(os.getcwd(), "abaqus_robot", "calibration_analysis",
                        "ReducedHydro_L2300_wall_supported")
    step = os.path.join(here, "Robot_parametric_L2p300_D0p815_WallWobble.step")
    started = time.time(); model = mdb.Model(name="L2300WallWobbleMesh")
    opener = getattr(mdb, "openStep", None) or getattr(mdb, "openSTEP", None)
    geometry = opener(step, scaleFromFile=OFF)
    part = model.PartFromGeometryFile(name="Robot_SOLID", geometryFile=geometry,
        combine=False, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    bbox = part.cells.getBoundingBox(); lo, hi = bbox["low"], bbox["high"]
    extents = [hi[i] - lo[i] for i in range(3)]
    if abs(extents[0] - 2.3) > 1e-6 or max(abs(extents[i] - .815) for i in (1, 2)) > 1e-6:
        raise RuntimeError("STEP extent gate failed: %r" % (extents,))
    # Keep the authorized nominal size while preventing the Bezier nose from
    # creating a few-micron local seed that dominates a rigid-body volume mesh.
    part.seedPart(size=.060, deviationFactor=.005, minSizeFactor=.05)
    part.setMeshControls(regions=part.cells, elemShape=TET, technique=FREE)
    part.setElementType(regions=(part.cells,), elemTypes=(ElemType(elemCode=C3D4, elemLibrary=STANDARD),))
    part.generateMesh()
    with open(os.path.join(here, "L2300_mesh060_nodes_local.csv"), "w", newline="") as stream:
        writer = csv.writer(stream); writer.writerow(["node", "x_mm", "y_mm", "z_mm"])
        for node in part.nodes: writer.writerow([node.label] + list(node.coordinates))
    with open(os.path.join(here, "L2300_mesh060_elements.csv"), "w", newline="") as stream:
        writer = csv.writer(stream); writer.writerow(["element", "n1", "n2", "n3", "n4"])
        for element in part.elements:
            writer.writerow([element.label] + [part.nodes[i].label for i in element.connectivity])
    result = {"source_step": step, "mesh_size_mm": .060, "bbox_low_mm": list(lo),
              "bbox_high_mm": list(hi), "bbox_extents_mm": extents,
              "node_count": len(part.nodes), "element_count": len(part.elements),
              "preprocessing_wall_clock_s": time.time() - started,
              "deviation_factor": .005, "min_size_factor": .05,
              "abaqus_import": "openStep + PartFromGeometryFile + C3D4 FREE TET"}
    with open(os.path.join(here, "L2300_mesh060_raw.json"), "w") as stream:
        json.dump(result, stream, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
