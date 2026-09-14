# Authoritative FreeCAD L2300 D0815 WallWobble Build

Final classification: **`L2300_D0815_FREECAD_CAD_FROZEN`**

## 1. Why CAD and simulation directories are separated

This canonical directory owns the robot geometry. The previous CAD files under
`calibration_analysis/` are retained only as historical generated references.
Abaqus may later import the frozen STEP, but it may not redefine length or
diameter by scaling nodes.

## 2. Authoritative dimensions

The master BRep measures `L = 2.300000 mm` from nose tip to flat tail and
`D = 0.815000 mm` across the straight cylindrical body. The axis is global `+X`.

## 3. HEAD definition

The HEAD is the frozen monotone, convex, axisymmetric degree-5 Bezier profile
already selected for this design. Its axial extent is `0.26149477 mm`; its final
three radial control coordinates equal `0.4075 mm`, imposing the intended G2
transition to the cylinder. The earlier envelope comparison gave a maximum
outward excess of `7.5607 um`. No new fit or alternate nose family was used.

## 4. Length parameterization

The straight cylindrical extent is `2.03850523 mm`. Relative to L1.800, L2.300
adds exactly `0.500000 mm` of cylinder. The computed volume increase is
`0.260840547541 mm3`, agreeing with `pi R^2 (0.500 mm)` within
`2.22e-16 mm3`. The HEAD solid symmetric-difference volume is zero.

## 5. FreeCAD build

FreeCAD `1.1.3` generated one valid, closed OpenCASCADE solid with one shell,
positive volume, no open seam, and a passing BOP/self-intersection check.

## 6. STEP/BREP/FCStd round-trip

STEP, BREP, and FCStd were independently reopened. Every format remained valid,
closed, one-solid/one-shell, and measured L/D within the `1e-6 mm` hard tolerance.
The STEP volume differs from the master by only `1.55e-15 mm3`.

## 7. Volume and surface area

- Volume: `1.140268022399 mm3`
- Surface area: `6.476763231800 mm2`

## 8. Mass / COM / inertia

At `7.80906654321e-9 tonne/mm3`, the mass is `8.904428864 mg`. The COM is
`(1.206018693716, 0, 0) mm` to numerical precision and lies inside the solid.
Principal inertias are positive: `3.929384741`, `3.929384741`, and
`0.725781069 mg mm2`.

## 9. Magnetic-volume identity

Constant magnetization gives a physical target moment of
`0.001040040143 A m2` from the authoritative CAD volume. This is an identity
calculation only; no Magpylib server or simulation was started.

## 10. Comparison to previous calibration-directory CAD

The rebuilt CAD is geometrically equivalent to the previous generated CAD
within roundoff: L difference `0`, transverse bounding difference
`1.11e-16 mm`, volume difference `1.55e-15 mm3`, COM difference `0`, and
principal-inertia difference `0`. The HEAD control points are identical.

## 11. SHA256 manifest

`cad_identity_manifest.json` records the SHA256 hashes of the authoritative
STEP, BREP, FCStd, and geometry JSON. Later Abaqus reports must record the exact
STEP hash from this manifest.

## 12. Final CAD decision

All master-shape, parameterization, round-trip, mass-property, magnetic-volume,
and regression gates pass. The canonical CAD is frozen as
`L2300_D0815_FREECAD_CAD_FROZEN`.

## 13. Exactly one next step

In a separate future simulation commit, import the frozen canonical STEP
directly into Abaqus and perform the simplified rigid solver-mesh audit.
