# Straight-pipe control of the production local-frame wobble case

## Decision

**THE CORE LIMITATION IS NOT SPECIFIC TO THE BEND.**

The straight tube supports a clear rotating robot-axis response, but it does not produce a clean sequence of isolated wall impacts and separations. After the initial transient, the robot rotates at approximately 15 deg directed tilt while remaining in an opposing two-end geometric bridge for almost the entire cycle. Curvature strongly worsens axial translation and orbit irregularity, but the persistent contact/geometry limitation is already present in the straight tube.

## Why a dedicated control was required

The existing `PROD_LOCAL30_G0_3CYCLE` trajectory starts at centerline arc length 13.821574 mm, beyond the configured bend start at 13.49 mm. Its complete visited interval, 13.821245--15.960064 mm, has non-negligible curvature (mean 0.205689 mm^-1; minimum 0.118971 mm^-1). There is no pre-bend straight interval to extract. Consequently, the proposed `StraightSection_FromCurrentCase_*` animations are not applicable and were deliberately not generated.

## Controlled model change

The solved control is `PROD_LOCAL30_G0_STRAIGHT_CTRL`, duration 33.333333 ms.

- Unchanged: Robot part and its 1588-node/7302-C3D4 mesh, initial robot pose, contact law (`mu=0.03`, damping fraction `0.5`), Reduced-Hydro implementation, direct step `dt=1e-7 s`, `B=10 mT`, 30 Hz, 30 deg cone, Bias=0, G=0, 1 ms ramp, phase 248 deg, production Fortran, and production socket server.
- Changed: the rigid wall helper and its magnetic centerline were replaced by a straight circular tube with effective lumen radius 0.66734524 mm and 24 mm physical length.
- Initialization: the straight tube axis passes through the unchanged Robot COM, as required for placement in a straight-tube centerline reference frame. This gives 187.2 um initial analytic clearance and Abaqus reported no initial overclosures or nodal adjustments.
- Phase control: a remote, nonphysical DXF prefix ends 40 mm behind the Robot and transports the validated transverse gauge into the straight reachable segment. The Robot's initial projection, tangent, and transverse basis match the curved production case to numerical precision; the prefix cannot be reached or selected by the Robot projection.
- Facet control: the 24-sided wall-node radius is `R/cos(pi/24)`, so the contact-face apothem, rather than the vertex radius, equals the target lumen radius.

The first two trial initializations were stopped at 3.3 ms and 1.1 ms because geometry gates detected, respectively, a faceting-induced radius error and an incompatible curved-projection offset. They are excluded from all results. The final control started without overclosure and completed all 333,334 increments successfully in 2858 s of Explicit wall time (3134 s including setup and extraction).

## One-cycle comparison

| Metric | Curved case, cycle 1 | Straight control | Interpretation |
|---|---:|---:|---|
| Local B winding | 0.999897 | 0.999898 | Identical commanded rotation |
| Robot-axis winding | 1.407365 | 1.411182 | Clear rotation in both; startup phase is included |
| Directed tilt, mean / max | 26.938 / 38.525 deg | 14.312 / 15.079 deg | Straight tube gives a compact, bounded orbit |
| Axial displacement | +2.014010 mm | +0.066080 mm | Curvature strongly amplifies translation |
| Contact events | 25 | 353 | Straight response contains extensive contact chatter |
| Longest contact | 2.012 us | 1164.263 us | Straight contact is not an isolated-impact regime |
| True separations after events | 25 | 126 | Separations occur, but not between most events |
| Contact sectors visited | irregular sequence | all 12 sectors | The bridge/contact region rotates around the tube |
| Raw contact-time fraction | not reported | 65.874% | Contact occupies most of the cycle |
| Both-end 20 um support fraction | 3.892% | 93.163% | Opposing support dominates the straight control |
| Longest 20 um bridge | 0.900 ms | 31.043 ms | Continuous from 2.300 ms to cycle end |
| Peak contact force | not used for baseline table | 4.4267 N | Short high-force events remain present |
| Minimum analytic side gap | n/a for curved wall | -8.065 um | Small contact penetration; no end interaction |
| Minimum tube-end margin | n/a | 10.720 mm | End effects are excluded |

The post-transient straight tilt is 14.829 deg on average over 5--33.333 ms (range 14.369--15.079 deg). The field magnitude remains 9.999999--10.000000 mT. The final COM radial offset is only 0.0434 mm, so the long bridge is caused by the tilted finite-length body contacting opposing ends, not by bulk COM drift to one wall.

## Interpretation boundary

One cycle is sufficient to reject the narrow hypothesis that bridge/stick behavior appears only in the bend: the straight control is already bridged for 31.04 of 33.33 ms. A second cycle is therefore not justified for this decision. One cycle is not sufficient to claim a multi-cycle straight-tube limit cycle or repeatability.

The straight response is rotationally acceptable but contact-mode unacceptable. The bend remains the dominant cause of large axial drift and irregular global motion, while the robot/tube clearance, finite body length, contact regularization, and Reduced-Hydro/contact interaction remain the core causes of prolonged two-end support.

## Best next step

Run one geometry-focused straight-tube diagnostic that reduces the finite-body bridge ratio while keeping the validated magnetic frame fixed. The most informative single variable is clearance-to-length ratio (tube inner diameter or robot effective length), evaluated first with a short one-cycle case and the same dense contact history. Do not return to magnetic-frame debugging.

## Reproducibility

- Source baseline commit: `ba8b5071ce1707912b8682097f687e1ef33b83f3`
- Production server SHA-256: `776736f93268deeef7194a9ee7ffd81ad59e8b27ac51c46f344bba60ce9b0e6d`
- Control input SHA-256: `4133af73a4b911f7c931f230e50c01f3fb74782fbef15d60e75198ad7143c16d`
- Control Fortran SHA-256: `940cd6817b091239db97e6ef40adf54f77fb107b3a7645898ff324d3db0e8f1d`
- Quantitative sources: `straight_pose.csv`, `straight_exact_gap_timeseries.csv`, `straight_contact_events.csv`, and `straight_metrics.json`

