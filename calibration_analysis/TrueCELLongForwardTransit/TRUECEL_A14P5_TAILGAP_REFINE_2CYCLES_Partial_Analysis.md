# TRUECEL_A14P5_TAILGAP_REFINE_2CYCLES: terminated partial diagnostic

This is NOT a completed two-cycle convergence result. Abaqus was terminated at 12.720 / 16.667 ms on user request; last field frame 12.700 ms.

## Existing evidence

- Common 8.0-12.5 ms inferred robot-fluid axial impulse: coarse `-6.7774488e-07` N s; refined `-3.2697940e-04` N s; ratio `482.4520`.
- Refined momentum closure error in the common window: `0.02%`.
- Available displacement through 12.700 ms: `-0.542796` mm (not a two-cycle delta).
- Refined ETOTAL max absolute drift to stop: `1322.57` N mm; minimum reported stable increment `1.07e-08` s.
- Abaqus warned that General Contact nodes penetrated tracked faces by over 50% of the typical 0.14599-mm element dimension (InfoNodeDeepPenetFirst).
- Robot-fluid force is inferred as whole-robot General Contact minus direct robot-wall contact; it is not a native pair-isolated robot-fluid measurement.

## Decision boundary

Mesh classification and second-cycle recoil: INDETERMINATE because the run stopped before 2T. The 8.0-12.5 ms impulse is measured over a common window, but deep penetration and large energy drift make a physical mesh-convergence claim unsafe. NUMERICS: WORSE. This was a fixed-17,600-element local radial r-refinement, not added-cell h-refinement.

All new GIFs are labeled PARTIAL/TERMINATED, use existing ODB frames only, and are not a completed scientific case.
