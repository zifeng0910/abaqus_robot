# Figure and data QA

- Backend: Python/matplotlib only.
- Data: one deterministic 100 ms Abaqus run; no observations were excluded and
  no smoothing, stochastic aggregation or invented uncertainty was used.
- Static figure contract: quantitative comparison grid; the claim is that the
  Cycle 2/3 orientation orbits nearly repeat while the contact/translation
  metrics do not.
- Source preflight: 19 PASS, 0 FAIL, 1 accepted warning. The warning is the
  intentional absence of TIFF, which is excluded from this repository.
- PDF text audit: minimum rendered text 6 pt, required minimum 5 pt, PASS.
- Visual QA: all three orbit panels and the assembled figure were inspected at
  final size. The global/local dual-view preview and initial, middle and final
  GIF frames were inspected. Axes are fixed, labels are readable, and no data,
  annotation or panel collisions were observed.
- Exports: editable SVG and PDF plus 600 dpi PNG preview; GIFs use fixed limits.
- Traceability: plotted coordinates come from `three_cycle_pose.csv`; summary
  and contact statements come from the cycle, comparison, event and exact-gap
  CSV files generated from the completed ODB and telemetry.
