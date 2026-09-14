# Figure QA

- Core conclusion: shortening to L=1.800 mm removes the long-model persistent bridge, while the exact circular head retains a short late bridge.
- Archetype: synchronized 3D image plate plus quantitative validation.
- Backend: Python/matplotlib only.
- Three-way GIF: synchronized physical time, fixed center/camera/zoom, common wall opacity and vector scales; event windows are oversampled.
- Slow-motion GIF: 5.8-7.0 ms only; near-wall nodes and corresponding wall triangles are highlighted without modifying geometry.
- Static source preflight: 20 pass, 0 warnings, 0 failures.
- Animation source preflight: 17 pass, 3 justified warnings, 0 failures. The warnings reflect GIF-oriented PNG output at 300 dpi and a 274 mm three-panel screen canvas rather than a journal raster panel.
- PDF text audits: all four PDFs pass the 5 pt floor; observed minima are 5-8 pt.
- Visual inspection: all static panels and four frames from each GIF were inspected; no clipping, blank panels, incoherent overlap, or camera-scale drift was observed.
- Data integrity: all requested physical time intervals are retained; nonuniform frame density adds frames and does not exclude any interval.
- Statistics: not applicable; panels show deterministic mesh and trajectory metrics without replicate aggregation or inferential statistics.
