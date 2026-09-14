# Figure QA

## Dynamic validation figures

All dynamic panels use the saved Python backend. The core visual claim is that L2.3 remains below the tumble boundary, maintains threshold-robust wall support, changes supported wall sector with phase, and advances forward without the phase/posture lock required for jam. Quantitative panels display deterministic time histories or complete pose sets, so replicate/spread statistics do not apply. Every panel is sourced from the public CSV files in this directory; no samples are excluded. The exact-gap line plots use the complete 2 us three-region records, while the 0.1 us single-region rows remain available for event minima and are not joined across missing regions.

| Figure/panel | Unique claim | Summary/spread | Labels and collisions | Pass |
|---|---|---|---|---|
| Surface regression | Four mismatches are uncertainty-bound | All 360 poses; 5 um band | Ambiguous poses use color plus edge | Yes |
| Directed tilt | No 90 deg crossing | Full 83,331-point history | Boundary remains distinct | Yes |
| Support fractions | 15/20/25 um robustness | Complete 2 us timeline | Grouped bars and legend clear | Yes |
| Sector versus phase | Supported sectors move | All supported 2 us samples | HEAD/TAIL encoded by color and labels | Yes |
| Gap/contact/bridge | Geometry and solver contact are separate | Complete gap plus every event | No panel overlap; initial gaps clipped at 100 um | Yes |
| Phase progression | No sustained plateau/reversal | Full history; 0.5 ms window | Threshold and traces clear | Yes |
| Torque | Magnetic drive persists | Full history; logarithmic norm | Both torque scales visible | Yes |
| Translation | Both forward gates pass | Full history | 2 ms gate is unobstructed | Yes |
| Energy | Bounded energy and dissipative hydro | Full history | KE and ALLKE remain distinguishable | Yes |
| Length comparison | Bridge alone is insufficient for jam | Three complete 8.333 ms cases | Four panels use consistent colors | Yes |

Both GIFs were inspected at five distributed frames including the first and last. The fixed camera, zoom, wall opacity, physical time, and vector scale remain constant. The single-case GIF shows moderate bounded tilt and translation; the three-case GIF visibly separates long-case constraint, L1.8 tumble, and the L2.3 intermediate wall-supported trajectory.

The source validator's TIFF/600 dpi warnings are accepted because repository policy explicitly excludes TIFF; 300 dpi PNG previews plus editable PDF/SVG are delivered. The rotated x tick warning applies only to the three-case labels, which were visually checked and do not collide. All PDF text runs pass the 5 pt floor after replacing logarithmic mathtext ticks with plain scientific notation.
Backend: Python/matplotlib only. Source data are the complete 360-pose static and CAD-vs-solver CSV tables; no rows were excluded. These are deterministic geometry audits, so replicate statistics and uncertainty intervals do not apply.

| Figure | Unique claim | Summary | Labels / legend | Collision check | Pass |
|---|---|---|---|---|---|
| `L2300_static_clearance_map` | Near-wall geometry emerges at intermediate tilt and depends on azimuth. | All exact pose gaps; 0 and 20 um contours. | Axes, color bar, and contour key are visible. | Contours and key do not obscure threshold structure. | Yes |
| `L2300_static_wall_support_fraction_vs_tilt` | Wall support grows through 28-35 deg while opposing geometry is delayed to 38 deg. | Fraction of all 36 azimuths at each prescribed tilt. | Direct legend; distinct markers and colors. | Legend and curves do not overlap. | Yes |
| `L2300_CAD_vs_solver_surface` | Positional error is small, but threshold-adjacent classification mismatches block dynamics. | All 360 paired CAD/mesh poses; identity and 0/20 um thresholds. | Mismatches use color plus larger outlined marks. | Threshold lines and mismatch marks remain distinguishable. | Yes |

Automated source preflight: 20 pass, 0 warn, 0 fail. PDF text audit: all three PDFs pass the 5 pt floor (minimum detected 6-7 pt). SVG text remains editable; PDF/SVG and 600 dpi TIFF/PNG exports were generated. Final-size visual inspection found no clipped labels, incoherent overlaps, or unreadable legends.
