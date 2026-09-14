"""Import the authoritative STEP in Abaqus/CAE and export its C3D4 mesh."""
from abaqus import mdb
from abaqusConstants import C3D4, DEFORMABLE_BODY, FREE, OFF, STANDARD, TET, THREE_D
from caeModules import *
import csv
import json
import os
import sys


HERE = os.path.join(os.getcwd(), 'abaqus_robot', 'calibration_analysis',
                    'ReducedHydro_FreeCAD_L1800_validation')
STEP = os.path.abspath(os.path.join(HERE, '..', '..', 'cad', 'freecad_parametric_robot',
                                    'Robot_parametric_L1p800_D0p815.step'))
SIZE = float(sys.argv[-1]) if sys.argv[-1].replace('.', '', 1).isdigit() else 0.035


def main():
    model = mdb.Model(name='FreeCADImport')
    opener = getattr(mdb, 'openStep', None) or getattr(mdb, 'openSTEP', None)
    if opener is None:
        print('MDB_IMPORT_METHODS', [x for x in dir(mdb) if 'open' in x.lower() or 'step' in x.lower()])
        raise RuntimeError('No STEP importer exposed by this Abaqus kernel')
    geometry = opener(STEP, scaleFromFile=OFF)
    part = model.PartFromGeometryFile(name='Robot_SOLID', geometryFile=geometry,
                                      combine=False, dimensionality=THREE_D,
                                      type=DEFORMABLE_BODY)
    bbox = part.cells.getBoundingBox()
    lo, hi = bbox['low'], bbox['high']
    extents = [hi[i] - lo[i] for i in range(3)]
    if abs(extents[0] - 1.8) > 1e-6 or max(abs(extents[i] - .815) for i in (1, 2)) > 1e-6:
        raise RuntimeError('STEP unit/extent gate failed: %r' % (extents,))
    part.seedPart(size=SIZE, deviationFactor=0.05, minSizeFactor=0.05)
    part.setMeshControls(regions=part.cells, elemShape=TET, technique=FREE)
    part.setElementType(regions=(part.cells,), elemTypes=(ElemType(elemCode=C3D4, elemLibrary=STANDARD),))
    part.generateMesh()
    with open(os.path.join(HERE, 'abaqus_step_nodes_local.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['node', 'x_mm', 'y_mm', 'z_mm'])
        for n in part.nodes:
            w.writerow([n.label] + list(n.coordinates))
    with open(os.path.join(HERE, 'abaqus_step_elements.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['element', 'n1', 'n2', 'n3', 'n4'])
        for e in part.elements:
            w.writerow([e.label] + [part.nodes[i].label for i in e.connectivity])
    identity = {'source_step': STEP, 'mesh_size_mm': SIZE, 'bbox_low_mm': list(lo),
                'bbox_high_mm': list(hi), 'bbox_extents_mm': extents,
                'node_count': len(part.nodes), 'element_count': len(part.elements),
                'abaqus_import': 'mdb.openStep(scaleFromFile=OFF) + PartFromGeometryFile'}
    with open(os.path.join(HERE, 'abaqus_step_import_raw.json'), 'w') as f:
        json.dump(identity, f, indent=2)
    print(json.dumps(identity, indent=2))


if __name__ == '__main__':
    from mesh import ElemType
    main()
