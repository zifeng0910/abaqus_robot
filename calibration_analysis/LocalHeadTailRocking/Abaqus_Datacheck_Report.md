# Abaqus datacheck: PROD_LOCAL_ROCK10_5HZ_G0_STRAIGHT

- Abaqus version: 2025, double precision, one CPU
- User subroutine: `vuforc_production_local.f`, compiled and linked successfully with Intel ifx
- Step: direct explicit, `dt=1e-7 s`, duration `0.200 s`
- Field interval: `5e-4 s` (401 intended field frames including the initial frame)
- Initial node-face or edge-edge overclosures: **none**
- Initial nodal position adjustments: **none**
- Unresolved initial overclosures: **none**
- Datacheck completion: **successful**

The 13 preprocessing warnings comprise the existing USER-amplitude dependency warnings and the four pre-existing distorted elements inherited unchanged from the validated source deck. The General Contact output warning also predates this case and states that only the portion of the named robot surface in the explicitly included robot-wall pair participates. No new geometry/contact error was reported.
