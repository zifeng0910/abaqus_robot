"""Import only the frozen canonical STEP and create the global-0.060 C3D4 mesh."""

from abaqus import mdb
from abaqusConstants import C3D4, DEFORMABLE_BODY, FREE, OFF, STANDARD, TET, THREE_D
from caeModules import *
from mesh import ElemType
import csv
import hashlib
import json
import os
import time


EXPECTED_SHA256 = "461a36dabd0cc1f94390b3e3dc24b2e059b8a74a98d3080ca67694c3a5c0d0e2"


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        while True:
            block = stream.read(1024*1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def main():
    repo = os.path.abspath(os.getcwd())
    here = os.path.join(repo, "calibration_analysis", "ReducedHydro_L2300_frozenCAD_validation")
    step = os.path.join(repo, "cad", "freecad_parametric_robot", "variants",
                        "L2300_D0815_wallwobble", "Robot_L2300_D0815_WallWobble.step")
    actual_sha = sha256(step)
    if actual_sha != EXPECTED_SHA256:
        raise RuntimeError("ABORT: frozen STEP SHA256 mismatch: %s" % actual_sha)
    started = time.time()
    model = mdb.Model(name="FrozenCAD_L2300_Global060")
    opener = getattr(mdb, "openStep", None) or getattr(mdb, "openSTEP", None)
    geometry = opener(step, scaleFromFile=OFF)
    part = model.PartFromGeometryFile(name="Robot_SOLID", geometryFile=geometry,
        combine=False, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    bbox = part.cells.getBoundingBox()
    low, high = bbox["low"], bbox["high"]
    extents = [high[i]-low[i] for i in range(3)]
    if abs(extents[0]-2.300) >= 1e-6 or max(abs(extents[i]-.815) for i in (1, 2)) >= 1e-6:
        raise RuntimeError("ABORT: imported STEP dimension mismatch: %r" % (extents,))
    part.seedPart(size=.060, deviationFactor=.05, minSizeFactor=.50)
    part.setMeshControls(regions=part.cells, elemShape=TET, technique=FREE)
    part.setElementType(regions=(part.cells,), elemTypes=(
        ElemType(elemCode=C3D4, elemLibrary=STANDARD),))
    part.generateMesh()
    nodes_path = os.path.join(here, "L2300_solver_mesh_nodes_local.csv")
    elements_path = os.path.join(here, "L2300_solver_mesh_elements.csv")
    with open(nodes_path, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["node", "x_mm", "y_mm", "z_mm"])
        for node in part.nodes:
            writer.writerow([node.label]+list(node.coordinates))
    with open(elements_path, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["element", "n1", "n2", "n3", "n4"])
        for element in part.elements:
            writer.writerow([element.label]+[part.nodes[i].label for i in element.connectivity])
    result = {
        "source_STEP": step,
        "source_STEP_SHA256": actual_sha,
        "scaleFromFile": False,
        "coordinate_scaling": False,
        "bbox_low_mm": list(low), "bbox_high_mm": list(high),
        "bbox_extents_mm": extents,
        "global_nominal_size_mm": .060,
        "HEAD_local_refinement": False,
        "deviation_factor": .05, "min_size_factor": .50,
        "node_count": len(part.nodes), "element_count": len(part.elements),
        "preprocessing_wall_clock_s": time.time()-started,
        "element_type": "C3D4", "mesh_technique": "FREE TET",
    }
    with open(os.path.join(here, "L2300_solver_mesh_raw.json"), "w") as stream:
        json.dump(result, stream, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
