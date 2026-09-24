# Magnetic-only rotational damping screen

**MAGNETIC_ROCKING_LIMIT_CYCLE_NOT_FOUND**

The passive wall has no contact interaction. Its crossing count is a geometric overlap diagnostic, not a simulated collision count.

| c/(I 2pi f) | c (N mm s) | bulk-water multiple | amplitude spread C2-C5 | phase spread (deg) | virtual wall events | min gap (mm) | longest transverse hover (ms) | gates |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0.1 | 3.501e-07 | 13.4 | 5.5% | 5.9 | 11 | -0.2037 | 0.000 | FAIL |
| 0.3 | 1.050e-06 | 40.2 | 2.1% | 1.3 | 10 | -0.1811 | 0.000 | FAIL |
| 1 | 3.501e-06 | 134.0 | 0.2% | 0.2 | 10 | -0.1187 | 0.000 | FAIL |
| 3 | 1.050e-05 | 401.9 | 0.0% | 0.0 | 10 | -0.0555 | 0.000 | FAIL |

Undamped baseline: 13 virtual wall crossings; minimum signed gap -0.2605 mm.

C2-C5 after 1-ms start-up ramp: amplitude spread <10%, phase spread <10 deg, harmonic residual <20%; no geometric wall crossing; no >=0.5-ms transverse hover (|theta|>=70 deg and |omega|<=20 rad/s)

Candidate c values are dynamic screening scales, 13-402 times the bulk-water slender-rod estimate. A passing dynamics gate alone does not calibrate the drag or prove hydrodynamic plausibility.
The input-mesh mass integral is 10.0188 mg, whereas the inherited case metadata records 9.2078 mg. The inertia and screening scale use the actual input mesh and density; the metadata discrepancy needs resolution before any calibrated hydrodynamic claim.

The GIF shows TRUECEL_MAGNETIC_ONLY_DAMPING_3P0X; it is the first passing candidate, or the smallest-wall-crossing diagnostic candidate if none passes.

No translational damping, propulsion, position constraint, magnetic parameter change, fluid domain, or active wall contact was added. All cases start at t=0 and use one uninterrupted Explicit step. The magnetic-only torque log records magnetic torque before rotational damping is subtracted in VUAMP.
