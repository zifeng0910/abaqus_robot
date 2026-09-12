# Figure and data QA

## Figure contract

- Core conclusion: damping fraction 0.20 with tangent fraction 0 lowers restitution but does not adequately suppress first-impact recoil.
- Evidence chain: exact gap establishes penetration; contact-point normal velocity establishes restitution; wall-force history establishes pulse shape; kinetic-energy history establishes dissipation; the recoil comparison applies the selection gate.
- Archetype: quantitative comparison grid plus a fixed-camera geometry animation.
- Backend: Python/matplotlib only.
- Export: 89 mm static panels as editable SVG/PDF and 600 dpi PNG; 183 mm impact poster and GIF. Local TIFF QA exports are intentionally excluded from Git because the requested PNG files are the delivery raster format.

## Data integrity

- Baseline and candidate use every 0.1 µs history sample in the contact and momentum calculations.
- Gap is an exact nearest-point query over all 217 exterior robot nodes and actual wall triangles in every solver-active event plus 5 µs shoulders; no interpolation crosses an impact.
- The plotted first-impact window is 1.000-1.010 ms for both cases. Cropping changes only the displayed time range, not any calculated result.
- No stochastic statistics, uncertainty interval, or replicate aggregation applies to this deterministic Abaqus comparison.
- CFS combines friction and tangential contact damping in the available output; those mechanisms are not falsely separated.

## Rendered panel audit

| Panel | Unique claim | Source | Labels / legend | Collision check | Pass |
|---|---|---|---|---|---|
| Gap | Penetration remains below 1 µm | current/candidate contact histories | Shared case colors and units | Curves and zero band clear | yes |
| Normal velocity | Restitution decreases only to 0.826 | current/candidate contact histories | Direct restitution annotation | Annotation clear of curves | yes |
| Wall force | Pulse remains microsecond-scale | baseline/candidate force histories | Shared case colors and units | Peaks fully visible | yes |
| Energy | Candidate dissipates more total KE | baseline/candidate energy histories | Delta K referenced to 1.0035 ms | Legend clear of curves | yes |
| COM recoil | Reduction is only 6.1% | comparison summary | Values shown above bars | Labels inside canvas | yes |
| GIF/poster | Actual geometry and contact vectors remain inspectable | private ODB extraction transformed by public scripts | Fixed camera, time, gap, force, normal velocity | Two views do not overlap | yes |

## Automated QA

- `test_probe_numerics.py`: pass.
- Python compilation: pass.
- Nature figure static preflight for `plot_normaldamp020.py`: 20 pass, 0 warning, 0 fail.
- Nature figure static preflight for `make_first_impact_gif.py`: 20 pass, 0 warning, 0 fail.
- PDF glyph audit: all six PDFs pass the 5 pt floor; minimum detected glyph is 5 pt.
- GIF audit: 85 frames, fixed camera, deterministic frame schedule with dense 0.2 µs sampling around first impact.
