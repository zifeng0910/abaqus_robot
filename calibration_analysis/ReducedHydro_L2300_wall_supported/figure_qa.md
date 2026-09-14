# L2300 preflight figure QA

Core conclusion: the L2300 CAD passes the old-envelope hard gate, but no dynamic solve is eligible because the nominal 0.060 mm mesh misses the HEAD-normal P95 gate.

Figure archetype: quantitative grid.

Target/output: double-column scientific report figure, 182.9 mm wide, Python/matplotlib, PNG/PDF/SVG.

Evidence map:

- `L1800_vs_L2300_geometry`: old-envelope agreement and the isolated 0.500 mm straight-body extension.
- `L2300_mesh_gate_failure`: overall versus HEAD normal error and the element-count cost of refinement.

Source data: `L2300_head_profile_fit.csv` and `L2300_mesh_trial_summary.csv`; no rows were excluded and no stochastic summary or uncertainty interval applies.

Automated QA: source syntax, editable vector text, font configuration, 600 dpi PNG, journal width, and backend exclusivity pass. Both PDF files pass the 5 pt glyph-floor audit with a 6 pt minimum. TIFF was intentionally omitted because the requested repository deliverable is PNG and prior TIFF exports are excluded from version control.

Visual QA: both panels were inspected at final size. Titles, axes, legends, annotations, and curves do not overlap or clip. The old envelope remains neutral, the L2300 profile is the primary signal, and the failed hard gate is redundantly encoded by position, color, and a dashed threshold rather than color alone.
