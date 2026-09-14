# Figure QA

Backend: Python/matplotlib only. Source data are the complete 360-pose static and CAD-vs-solver CSV tables; no rows were excluded. These are deterministic geometry audits, so replicate statistics and uncertainty intervals do not apply.

| Figure | Unique claim | Summary | Labels / legend | Collision check | Pass |
|---|---|---|---|---|---|
| `L2300_static_clearance_map` | Near-wall geometry emerges at intermediate tilt and depends on azimuth. | All exact pose gaps; 0 and 20 um contours. | Axes, color bar, and contour key are visible. | Contours and key do not obscure threshold structure. | Yes |
| `L2300_static_wall_support_fraction_vs_tilt` | Wall support grows through 28-35 deg while opposing geometry is delayed to 38 deg. | Fraction of all 36 azimuths at each prescribed tilt. | Direct legend; distinct markers and colors. | Legend and curves do not overlap. | Yes |
| `L2300_CAD_vs_solver_surface` | Positional error is small, but threshold-adjacent classification mismatches block dynamics. | All 360 paired CAD/mesh poses; identity and 0/20 um thresholds. | Mismatches use color plus larger outlined marks. | Threshold lines and mismatch marks remain distinguishable. | Yes |

Automated source preflight: 20 pass, 0 warn, 0 fail. PDF text audit: all three PDFs pass the 5 pt floor (minimum detected 6-7 pt). SVG text remains editable; PDF/SVG and 600 dpi TIFF/PNG exports were generated. Final-size visual inspection found no clipped labels, incoherent overlaps, or unreadable legends.
