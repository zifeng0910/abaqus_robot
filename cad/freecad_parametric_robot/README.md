# Parametric FreeCAD Robot

This directory contains the authoritative, parameter-driven CAD for the
axisymmetric robot. The master geometry is an exact FreeCAD Part/OpenCASCADE
BRep, not a mesh and not a scaled Abaqus model.

## Change the dimensions

Edit only these values near the top of `build_parametric_robot.py`:

```python
L_TOTAL_MM = 1.800
D_BODY_MM = 0.815
```

Then run the generator with the FreeCAD command-line executable:

```powershell
& 'I:\Program Files\FreeCAD 1.1\bin\freecadcmd.exe' `
  -c "import sys; sys.path.insert(0,r'J:\abaqusfangzhen\abaqus_robot\cad\freecad_parametric_robot'); import build_parametric_robot as b; b.main()"
```

Run the tests after generation:

```powershell
& 'I:\Program Files\FreeCAD 1.1\bin\freecadcmd.exe' `
  -c "import runpy; runpy.run_path(r'J:\abaqusfangzhen\abaqus_robot\cad\freecad_parametric_robot\test_parametric_robot.py',run_name='__main__')"
```

The generator aborts if the solid is invalid, open, non-positive in volume,
not exactly one solid, or outside the 1e-6 mm bounding-box tolerance. It also
reimports STEP and BREP and reopens FCStd before reporting success.

## Geometry rules

### Rule 1

D controls robot radial scale and head scale.

### Rule 2

L only changes the straight cylindrical section.

### Rule 3

Changing L must NEVER stretch the head.

### Rule 4

R_head/R_body = 0.61/0.60 unless the user explicitly redesigns the nose family.

The robot axis is global +X. The nose tip is at X=0 and the flat tail is at
X=L_total. The upper axial profile uses one exact circular arc, tangent to the
straight line at Y=R_body, followed by the flat tail and axis closure. That
profile is revolved 360 degrees around global X.

## Outputs

`Robot_parametric_L1p800_D0p815.step` is the authoritative interchange CAD.
The matching `.FCStd` and `.brep` files retain native and direct BRep forms.
The `.stl` file is visualization-only and must never be used as the dimensional
source. Its linear and angular deflections are recorded in the geometry JSON.

The geometry and mass-property JSON files are intended as the shared identity
record for later Abaqus, mass/inertia, and magnetic-volume work. The present
generator does not run or modify Abaqus, SolidWorks, or Magpylib.
