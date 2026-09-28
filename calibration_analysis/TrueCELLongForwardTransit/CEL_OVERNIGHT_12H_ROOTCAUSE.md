# CEL Overnight 12 h Root-Cause Investigation

## Scope and controls

This report closes the INITFIX20 fluid-solid tangential-friction isolation campaign. The S6/S4 outlet correction removed most of the original axial-boundary artifact; INITFIX restored near-complete Eulerian occupancy and amplified reversal. All four runs use the same 0–20 ms step, 801 frames, EVF initialization, magnetic forcing, damping, mesh, outlets, fluid normal `HARD` contact, robot-wall `mu=0.03`, and unchanged Fortran physics. Only the named tangential fluid pair(s) are changed to frictionless.

## Results

| Case | robot-fluid μ | wall-fluid μ | C1 Δs (mm) | C2 Δs (mm) | min v (mm/s) | C2 reverse (mm) | ALLFD20 (N·mm) | ALLPW20 (N·mm) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `INITFIX20` reference | 0.03 | 0.03 | -0.05618 | -0.85748 | -205.13 | 0.85748 | 0.076247 | 0.276250 |
| `FLUIDFRIC0` both fluid pairs | 0 | 0 | -3.06459 | -5.21363 | -734.12 | 5.21363 | 0.000417 | 36.308872 |
| `RF_FRIC0` robot-fluid only | 0 | 0.03 | -0.79488 | -0.00327 | -418.71 | 0.00327* | 5.436177 | 25.305626 |
| `WF_FRIC0` wall-fluid only | 0.03 | 0 | -0.04439 | -1.46925 | -226.79 | 1.46925 | 0.063645 | 0.238181 |

`*` RF-only has a large first-cycle excursion (0.79846 mm) followed by near-zero C2 increment; the C2 value alone therefore understates its early reversal.

The simultaneous frictionless case reverses at 1.269673 ms and reaches -734.12 mm/s at 7.129842 ms. Relative to INITFIX20, its first 0.01 mm displacement divergence occurs at 1.892970 ms; maximum/RMS displacement differences are 7.36456/3.98755 mm and velocity differences are 727.660/415.096 mm/s. RF-only first diverges at 1.875716 ms (0.79626 mm maximum displacement difference), while WF-only first diverges at 6.303658 ms (0.59998 mm maximum difference; only 0.01771 mm through 8 ms).

## Causal interpretation

The three perturbations are strongly non-additive. RF-only nearly removes the C2 continuing reversal, WF-only does not, and changing both pairs frictionless amplifies the reversal and fluid work dramatically. The correct root-cause classification is therefore **combined/nonlinear tangential sensitivity**. These data show that the tangential fluid-solid formulation materially controls the trajectory, but they do not validate `mu=0.03` as a water-viscosity model or frictionless contact as a no-slip boundary. The ALLFD/ALLPW changes are consequences of the altered contact response, not an energy-based proof of a unique physical mechanism.

## Contact and solver audit

`PROP_CEL_HARD`, normal `pressure-overclosure=HARD`, EVF, S6/S4, field/output requests, mesh, magnetic field, damping, and step duration are unchanged in all child inputs. Contact assignments remain robot-wall to `PROP_CEL_HARD` and both Eulerian fluid pairs to `PROP_CEL_FLUID_HARD`; the RF/WF manifests and contact audits record the single pair change. All three child datachecks pass with the expected INITFIX warning family plus the empty-friction notice, and all three dynamics analyses complete successfully with 801 frames. The unchanged Fortran source is byte-identical after child output-directory normalization.

## Decision and next experiment

The four-run overnight cap is reached. No fourth dynamics run and no penalty-stiffness test are launched. The next allowed experiment is one normally parameterized contact/interface reconstruction validation: preserve the successful INITFIX EVF and outlet setup, use physically justified fluid-solid tangential behavior, and audit the RF/WF pair assignments and interface geometry before interpreting trajectory differences. This is the smallest test that can distinguish contact semantics from the combined numerical sensitivity found here.

## Artifacts and limits

The comparison GIF is `CEL_FLUID_CONTACT_FRICTION_comparison.gif`; the three single-case GIFs, summary JSON/CSV files, input identity checks, contact audits, manifests, and solver `.sta/.dat/.log` files are included in the result package. Large ODB/SIM/RES/restart files and private NPZ histories remain local and are intentionally excluded from Git and the ZIP. Unrelated working-tree changes remain untouched.
