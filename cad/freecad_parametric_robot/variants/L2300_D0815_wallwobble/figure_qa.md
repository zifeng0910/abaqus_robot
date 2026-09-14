# Figure QA

Core conclusion: the L2.300 CAD retains the frozen Bezier HEAD and differs from
L1.800 only by `0.500 mm` of straight cylindrical body.

Archetype: quantitative geometry identity views. The profile is the primary
dimensional view; the 3D preview and length-parameterization overlay are
supporting views.

All figures were generated with the saved Python backend from the canonical
profile CSV and geometry/mass JSON. No simulation output or synthetic geometry
was used. PNG, editable-text SVG/PDF, and 600 dpi TIFF exports are included.

| Figure | Unique claim | Labels | Collision check | Result |
|---|---|---|---|---|
| Authoritative profile | L, D, HEAD, body, TAIL, and COM identity | Direct | Clear at final size | PASS |
| BRep preview | Overall 3D form and identity values | Two fixed header rows | Clear of model and axes | PASS |
| L1800 vs L2300 | HEAD is identical; only 0.500 mm body is added | Direct arrow and legend | Clear at final size | PASS |

The strict source validator reports 20 passes, zero warnings, and zero failures.
All three PDF glyph audits pass with a minimum detected text size of 7 pt against
the 5 pt floor. Final rendered PNGs were inspected individually; geometry,
dimensions, labels, axes, and titles are visible without overlap.
