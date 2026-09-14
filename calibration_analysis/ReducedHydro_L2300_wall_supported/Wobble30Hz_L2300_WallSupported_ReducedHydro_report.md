# Wobble30Hz L2300 wall-supported Reduced-Hydrodynamics report

## 1. Why L1.800 is rejected as a final design

The prior L1.800 result crosses 90 degrees and reverses HEAD/TAIL axial polarity. It is a tumble case, not the requested wall-supported wobble, and is not treated as a successful final design.

## 2. Desired wall-supported wobble regime

The target is repeated wall approach/contact with continued phase progression, no 90-degree tumble, no polarity reversal, and no persistent opposing-wall jam. Forward translation is a separate late-time gate.

## 3. Why L2.300 was selected

The prescribed 2.300 mm length places the straight-section estimate near the 30-32 degree wall-support range. This estimate was not used as a substitute for the required exact SmoothWall114 audit.

## 4. Old successful head-envelope reconstruction

The head reference was reconstructed from the actual exterior nodes of the successful old scaled L1800 C3D4 mesh using the continuous upper axial-radial convex hull. The final smooth degree-5 Bezier profile is monotone and convex, has a smooth nose and a G2 cylinder transition. Maximum outward excess is **7.561 um**, below the 10 um hard limit. RMS difference is **5.142 um**, narrowly above the preferred 5 um target.

## 5. FreeCAD source-of-truth and length parameterization

The independent CAD family is valid, closed, and reimports consistently from STEP, BREP, and FCStd. Its dimensions are 2.300 x 0.815 mm. The head and flat tail definition are fixed; the L1800-to-L2300 change adds exactly 0.500 mm of straight cylinder.

## 6. CAD mass properties and magnetic moment

CAD volume is **1.140268022 mm3**, uniform-density mass is **8.904428864 mg**, and the constant-magnetization rule gives **0.0010400401427 A m2**. No magnetic, contact, damping, or Reduced-Hydrodynamics coefficient was changed.

## 7. Nominal 0.060 mm mesh audit

The final curvature-refined nominal 0.060 mm mesh has 29,302 nodes and 145,677 C3D4 elements. Length and diameter are 2.299999952 and 0.815000034 mm. Mass error is 0.0326%, COM error 0.0612 um, and maximum principal-inertia error 0.0616%. Overall surface-normal P95 is 0.6450 degrees.

The HEAD surface-normal P95 is **2.3948 degrees**, above the mandatory **2.0 degree** limit. Coarser variants reached the expected 25k-35k element range but had still larger HEAD P95 errors. Curvature refinement increased the final mesh to 145,677 elements without clearing the hard gate.

## 8. Exact SmoothWall114 static pose audit

Not run. The workflow requires all CAD and mesh gates to pass before the tilt-by-azimuth scan. Running it after the failed HEAD-normal gate would violate the prescribed order.

## 9. Dynamic eligibility and solver budget

The exact job `Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L2300_D0815_WallSupported_0083` was **not run**. No ODB was created and the one-job dynamic budget remains unconsumed.

## 10. Classification

**L2300_MESH_GATE_FAILED_STOP**

This is a preprocessing failure, not evidence that L2.300 tumbles or jams. The requested dynamic classifications cannot be assigned without an eligible solve.

## 11. Evidence boundary and next decision

The CAD shape passes its hard envelope and solid-identity gates, while its RMS fit misses the preferred target by 0.142 um. The nominal 0.060 mm discretization cannot simultaneously meet the expected mesh scale and the HEAD P95 gate for this short, blunt smooth head. The next step requires an explicit choice: permit targeted head-surface refinement beyond the nominal strategy, or revise the smooth-head representation while retaining L2.300. No automatic length sweep or physics tuning is justified.
