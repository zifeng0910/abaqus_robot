# Wobble30Hz geometric jam audit

Classification: **PERSISTENT_GEOMETRIC_JAM**

## Actual robot geometry
```text
 n_exterior_nodes  axial_length_mm  max_diameter_mm  p95_diameter_mm  median_diameter_mm  head_diameter_mm  tail_diameter_mm   rp_x_mm   rp_y_mm   rp_z_mm  head_x_mm  head_y_mm  head_z_mm  tail_x_mm  tail_y_mm  tail_z_mm
              334         2.931415         0.903075         0.861107            0.721622          0.875326          0.903075 -7.481362 -3.681165 -9.554105  -6.135879    -3.8667  -9.293127  -8.827986   -3.54069  -9.838923
```

## Local SmoothWall114 section
```text
 station  station_arc_mm  center_x_mm  center_y_mm  center_z_mm  tangent_x  tangent_y  tangent_z       n_x      n_y      n_z      b_x      b_y          b_z  slab_halfwidth_mm  n_span_mm  b_span_mm  radial_median_mm  radial_min_mm  wall_nodes_in_slab
       2       13.873035    -7.438677    -3.227893     -9.80332   -0.96753   0.100573  -0.231886 -0.230643 0.023975 0.972743 0.103391 0.994641 3.469447e-18               0.45   2.044327   2.058291           1.17512       0.677412                  19
```

The wall section uses validated Abaqus-frame SmoothWall114 R3D4 triangles (764 quads split into 1,528 triangles) and a ±0.45 mm tangent slab. The orientation map uses the full 334-node exterior mesh, exact point-to-triangle nearest distances, and a ±0.24 mm transverse COM-offset grid; it is a coarse θ=2°, χ=10° map with explicit geometric/5 µm/20 µm thresholds in `orientation_feasibility_map.csv`.

Cylinder sanity: using the conservative smaller local span, asin((D-d)/L) gives theta_max = **22.91 deg**, while the dynamic cone30 path reaches 35.00 deg. This formula is only a sanity check; the mesh/triangle map is the acceptance test.

## Dynamic overlay
- first opposing-wall bridge: `0.0022500304039567` s
- maximum near-20-µm exterior nodes: 11
- maximum penetrating exterior nodes: 3

The dynamic path is overlaid without altering the physical model. A cone25 job is **not** submitted by this audit; it is permitted only because the geometry evidence now identifies a tilt-and-jam window.
