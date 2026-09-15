# Figure QA notes

## Contract

- Core claim: the production magnetic frame does not drive tube-axis winding, whereas the controlled local-tangent field does; a 5 ms ramp reduces the initial force peak without changing that conclusion.
- Evidence chain: field-only panels establish the frame mismatch; robot-axis orbit panels distinguish slanted lock from winding; synchronized GIFs show the motion sequence; the winding plot provides the numerical time-history comparison.
- Archetype: quantitative grid with synchronized mechanism animations.
- Backend: Python/Matplotlib only.
- Export contract: PNG at 600 dpi plus editable-text SVG/PDF for static plots; GIF for animations; minimum PDF glyph size 5 pt.
- TIFF exception: the audit directory intentionally excludes TIFF from version control; the source preflight's TIFF warning is accepted because equivalent 600 dpi PNG previews and vector SVG/PDF exports are delivered.

## Rendered checks

| Artifact/panel | Unique claim | Visual check | Result |
|---|---|---|---|
| Field diagnostics a | Production field remains far from the tube axis | Axes, legend, and curves clear | Pass |
| Field diagnostics b | Production local azimuth does not complete one turn | Shared time axis and line hierarchy clear | Pass |
| Field-tip orbits a-c | Production and local cone centers differ | Start/end markers and all three panels visible | Pass |
| Dynamic winding a | Only local-tangent cases produce full robot-axis winding | Four traces distinguishable; ramp onset delay visible | Pass |
| Dynamic winding b | Local-tangent B winds once while production B does not | Four traces and legend clear | Pass |
| Four robot-axis orbit plots | Current is a short arc; local cases form complete orbits | Start/end/contact markers and 20/30/40 deg guides clear | Pass |
| Three-way comparison GIF | Frame definition changes the motion mode | First, middle, and final frames checked; fixed axes and labels clear | Pass |
| Four dual-view GIFs | Global and local motion remain interpretable through the cycle | Persistent title/status objects prevent frame accumulation; ramp first/middle/final frames checked | Pass |

No observations or cases were excluded. Contact markers use a window around each plotted sample so microsecond impacts are not silently lost during temporal downsampling. The animations are diagnostic visualization, while numerical metrics remain sourced from dense ODB histories.
