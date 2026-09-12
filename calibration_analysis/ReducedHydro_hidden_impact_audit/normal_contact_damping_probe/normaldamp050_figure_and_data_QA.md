# NormalDamp050 figure and data QA

## Figure contract

- Core conclusion: increasing only normal critical damping from 0.20 to 0.50 lowers first-impact contact-point restitution into the provisional 0.30-0.65 window without sticking, unacceptable penetration, energy creation, or momentum imbalance.
- Evidence chain: restitution establishes the primary gate; gap establishes penetration and reopening; wall force establishes event timing and loading; energy establishes dissipation; COM recoil is a secondary diagnostic; the fixed-camera animation shows the same contact geometry and vectors.
- Archetype: quantitative comparison series plus one fixed-camera mixed-modality animation.
- Backend: Python/matplotlib only.
- Export: 89 mm static plots and 183 mm poster, editable PDF/SVG, 600 dpi PNG/TIFF, minimum rendered glyph size 5 pt.

## Data integrity

- Reference and candidate RP histories each contain 13,001 samples at 0.1 µs spacing.
- Contact force histories use every explicit increment; no impact samples were dropped.
- Candidate exact gap was recomputed for all robot surface nodes at all 13,001 increments, not only at field frames or around detected events.
- Static impact plots display 1.000-1.010 ms for legibility. This is an explicit axis window, not removal from numerical analysis.
- The animation uses the same prescribed frame times for both cases and a fixed shared camera. It does not auto-fit per frame.
- Curves are deterministic single-run outputs; no replicate distribution or uncertainty interval exists, so error bars are not applicable.

## Rendered panel audit

| Artifact | Unique claim | Source | Labels/legend | Collision check | Pass |
|---|---|---|---|---|---|
| Restitution | ζ=0.50 enters the provisional window | `normaldamp020_vs_050_summary.csv` | direct values and shaded window | clear | yes |
| Gap | rebound opening slope decreases while penetration remains limited | first-impact history CSVs | shared legend | clear | yes |
| Wall force | force pulse remains short and finite | contact-force CSVs | shared legend | clear | yes |
| Energy | ζ=0.50 removes more total kinetic energy | energy-history CSVs | shared legend | clear | yes |
| COM recoil | recoil reduction is modest and secondary | summary CSV | values and diagnostic note | clear after annotation adjustment | yes |
| Poster/GIF | same geometry, contact node, normal and velocity definitions | private pose arrays plus public force CSVs | per-panel status text | clear | yes |

## Automated checks

- Static source preflight: 20 pass, 0 warn, 0 fail.
- Animation/poster source preflight: 20 pass, 0 warn, 0 fail.
- PDF text audits: all six PDFs pass; minimum glyph size is 6 pt for static plots and 5 pt for the poster.
- Numerical regression: decision gates, momentum, energy, gap, event count, output inventory, and 15 report sections pass.
- Python compilation: all new analysis, plotting, animation, report, and test sources compile.

TIFF and SVG files are retained locally for QA but excluded from Git as requested. ODB, private NPZ arrays, solver binaries, and raw telemetry are also excluded.
