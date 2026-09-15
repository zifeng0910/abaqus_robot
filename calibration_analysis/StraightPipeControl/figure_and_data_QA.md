# Figure and data QA

## Figure contract

- Core conclusion: the straight tube permits bounded rotation but still produces prolonged opposing two-end bridge/contact, so the core limitation is not bend-specific.
- Evidence chain: synchronized motion views establish geometry and timing; the local orbit establishes bounded rotation; the phase-matched comparison establishes the curved/straight contrast; the static summary quantifies tilt, winding, event count, and bridge duration.
- Archetype: image plate plus quantitative comparison.
- Backend: Python/matplotlib only.
- Export: GIF for requested motion products; editable PDF/SVG and 600 dpi PNG for the static comparison. TIFF was intentionally excluded by project instruction.

## Data integrity

- The Abaqus solve completed 333,334 increments at direct `dt=1e-7 s`.
- RP displacement, rotation, velocity, angular velocity, surface contact force, and global energy were extracted from the completed ODB.
- Contact events use every increment. Gaps and animation pose use every 100th increment (10 us) plus the final increment; counts are 333,335 dense records and 3335 pose/gap records.
- No rows were excluded. GIF frames are evenly selected from the full pose table only for display rate.
- The straight analytical side gap uses the actual 1588 Robot nodes reconstructed by the ODB rigid-body pose. The target wall radius is 0.66734524 mm.
- Contact truth comes from surface-force history; the 20 um bridge is a geometric proximity diagnostic and is not substituted for contact force.
- Phase interpolation asserts strictly increasing source time.

## Rendered audit

| Product/panel | Unique claim | Data summary | Visual QA | Pass |
|---|---|---|---|---|
| Dual-view GIF, global | Robot remains inside the straight control tube | Raw pose and fixed tube geometry | Initial/middle/final frames inspected; fixed limits; no clipping | yes |
| Dual-view GIF, local | Body, B, axis, sector and state are synchronized | Raw pose plus nearest-sector gap | Initial/middle/final frames inspected; labels clear | yes |
| Local-orbit GIF | The axis follows a compact rotating orbit | Full q1/q2 trail | Initial/middle/final frames inspected; stable dimensions | yes |
| Comparison GIF | Straight orbit is compact while curved cycle 1 is irregular | Phase-matched raw cycle data | Initial/middle/final frames inspected; identical axes | yes |
| Static panel a | Straight tilt rapidly settles near 15 deg | Raw 10 us pose samples | Curves and legend clear | yes |
| Static panel b | Straight orbit is bounded and compact | Raw q1/q2 samples | Identical coordinate meaning and equal aspect | yes |
| Static panel c | Straight greatly increases event/bridge metrics | Ratio to curved cycle 1, log axis | Reference ratio shown; labels fit | yes |

## Automated checks

- `validate_figure.py`: 19 PASS, 0 FAIL; the sole warning is absence of TIFF, which is intentional.
- `audit_pdf_text.py --min-pt 5`: run on the final PDF before commit.
- GIF dimensions/frame counts: comparison 720x335/151; dual view 720x345/181; local orbit 420x380/151.
- GIF initial, midpoint, and final frames were visually inspected for every requested animation.

No stochastic aggregate or inferential statistic is shown; uncertainty intervals and replicate counts are therefore not applicable.
