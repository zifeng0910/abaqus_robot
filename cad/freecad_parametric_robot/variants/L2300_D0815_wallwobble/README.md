# L2300 D0815 WallWobble FreeCAD CAD

This directory is the sole geometry source of truth for the `2.300 x 0.815 mm`
wall-wobble robot variant.

## Geometry ownership

FreeCAD/OpenCASCADE owns the total length, body diameter, degree-5 Bezier HEAD,
flat TAIL, volume, center of mass, and inertia. The robot axis is global `+X`,
with the nose tip at `X=0` and the flat tail at `X=L`.

Abaqus owns only mesh discretization, rigid-body definition, contact, and
dynamics. Scaling Abaqus nodes is forbidden as a way to modify CAD dimensions.

The frozen HEAD occupies `0.26149477 mm`. Changing `L` must leave that HEAD
unchanged and alter only the straight cylindrical body. Therefore the L1.800 to
L2.300 change adds exactly `0.500 mm` of straight cylinder.

## Build

Run the builder headlessly:

```powershell
& 'I:\Program Files\FreeCAD 1.1\bin\freecadcmd.exe' -c "import sys; sys.path.insert(0,r'J:\abaqusfangzhen\abaqus_robot\cad\freecad_parametric_robot\variants\L2300_D0815_wallwobble'); import build_L2300_D0815_wallwobble as b; b.main()"
```

Then run the FreeCAD tests:

```powershell
& 'I:\Program Files\FreeCAD 1.1\bin\freecadcmd.exe' -c "import runpy; runpy.run_path(r'J:\abaqusfangzhen\abaqus_robot\cad\freecad_parametric_robot\variants\L2300_D0815_wallwobble\test_L2300_D0815_wallwobble.py',run_name='__main__')"
```

The builder creates exact FCStd, STEP, and BREP files, reopens each format, and
applies validity, closure, single-solid, single-shell, dimension, volume, COM,
and inertia gates. The PNGs are previews only and never geometry inputs.

## Future Abaqus input

Future Abaqus scripts must import only:

`cad/freecad_parametric_robot/variants/L2300_D0815_wallwobble/Robot_L2300_D0815_WallWobble.step`

They must not import a STEP from `calibration_analysis/`. The earlier files in
that directory remain historical references and are not authoritative.

`cad_identity_manifest.json` records the exact SHA256 identity of the STEP,
BREP, FCStd, and geometry JSON used by later work.
