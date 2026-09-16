# S4 Head-Forward Low-G Screen

This directory contains exactly four new short non-CEL Explicit dynamics runs:

- `S4_HEADFORWARD_G0`
- `S4_HEADFORWARD_G0P25`
- `S4_HEADFORWARD_G0P5`
- `S4_HEADFORWARD_G1P0`

The robot mesh is physically rotated 180 degrees about fixed `n_routeA` so the
semantic HEAD is forward in canonical increasing centerline arclength. The
UR=0 magnetic-moment and ReducedHydro body-axis references are transformed by
the same rigid rotation. Geometry, tube, contact, damping, hydrodynamics,
field amplitude, rocking mode, RouteA gauge, and time step are frozen.

The required deliverables are the four individual GIFs,
`S4_HEADFORWARD_LowG_4Way.gif`, `S4_HEADFORWARD_LowG_ContactSheet.png`, the
screening report, and the setup/forward-gate audit. The old S4 G=6 mT result is
comparison-only. No automatic ranking is assigned.

**USER VISUAL SELECTION REQUIRED**
