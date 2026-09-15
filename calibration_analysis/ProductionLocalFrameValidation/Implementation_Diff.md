# Production magnetic-frame candidate

## Scope

The candidate is a vendored copy of the live production socket server. The live
file at `J:\magpy\magpylib_socket_server.py` was deliberately left unchanged
until all offline gates pass.

## Candidate changes

- Names the historical behavior `LEGACY_DRIVER_FRAME` and keeps it as default.
- Adds the explicit `ROBOT_LOCAL_TANGENT` analytic field mode.
- Projects the current robot position onto DXF segments independently of
  `adaptive_lead_mm`, with previous-arc tie continuity and endpoint extrapolation.
- Builds a polarity-fixed tangent and parallel-transported transverse frame.
- Rejects physical-source drives and the legacy global-Z cone bias in local mode.
- Uses the requested local-cone equation and startup ramp without endpoint taper.
- Preserves the historical `sense=-1` equation in legacy mode.
- Reports the field-frame mode and the actual segment-projected robot arc.

## Deployment status

**NOT DEPLOYED.** The moving-trajectory reproduction gate failed because the
completed screening wrapper used nearest-resampled-vertex projection, whereas
the requested production mode uses continuous segment projection. See
`Production_Local_Tangent_Frame_3Cycle_Validation.md`.
