# F60/G0.15 Strict CEL Screen

`S4_HEADFORWARD_F60_G0P15_TRUECEL` is the only new case and is classified `STRICT_CEL_FLUID_PRESENT` / `INITIALIZED_UNIFORM_FLOW_TRUE_CEL`. It uses the authoritative F60 parent magnetic path unchanged except `G=0.15 mT`; ReducedHydro is absent. Actual ODB `EVF` and `V` fields show 1536 initially filled cells and an early near-robot velocity disturbance of 357.282569 mm/s, so the CEL coupling gate passed.

The fixed view is canonical `+s` left-to-right, TAIL left/blue and HEAD right/red. The fluid snapshot is derived from actual ODB volume fraction and velocity fields, not a schematic.

| metric | value |
| --- | ---: |
| actual rocking frequency (Hz) | 60.621450 |
| fitted rocking amplitude (deg) | 6.373434 |
| rocking min/max (deg) | -9.656845 / +13.422439 |
| phase lag (deg) | +0.342824 |
| delta_s (mm) | -0.893101 |
| mean/final/max axial velocity (mm/s) | -34.547481 / +782.284401 / +782.284401 |
| fraction v_s > 0 | 0.668271 |
| axial reversals | 22 |
| max radial COM displacement (mm) | 0.260611 |
| HEAD/TAIL/BODY contact episodes | 8 / 10 / 11 |
| longest contact / BOTH bridge (ms) | 22.712417 / 0.000000 |
| fluid-filled volume / mass | 5.184000 mm3 / 5.187629 mg |
| CEL dt min/typical (s) | 4.265310e-08 / 4.675145e-08 |

The CEL force series is a robot momentum-balance estimate after subtracting Magpylib and therefore includes simultaneous robot-wall contact; it is not falsely presented as a pair-isolated Abaqus contact output. Comparison is only against `S4_HEADFORWARD_F60_G0P10_FLOW`. This run was stopped at 25.851 ms before the requested 33.333 ms endpoint, so all metrics and GIF frames are partial-run screening evidence.

## Screen conclusion

The true-CEL coupling gate passed, but this partial run is not an acceptable natural head-forward candidate: `delta_s` is negative, axial reversals are numerous, the final axial speed is nonphysical, and the ODB fluid velocity field contains extreme late-time values. Treat the GIF and metrics as a technical CEL-coupling diagnostic, not as a successful motion result.
