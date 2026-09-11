# Wobble30Hz cone-frame azimuth report

## 1. Why the directional pulse branch is not continued

The previous single Maxwell-consistent directional pulse reduced reverse speed only about 10%. It is therefore recorded as `SINGLE_PHASE_DIRECTIONAL_PULSE_DOES_NOT_RECOVER_FORWARD_MOTION`; its `No forward phase window found` wording applies only to that one candidate, not to all phases. This report isolates a different spatial degree of freedom.

## 2. Current production cone frame

The production analytic field is 30 Hz, 10 mT, cone half-angle 30°, Bias40, phase 248°, sense +1, G6/L45 (`G=3 mT = 0.003 T`, `L=45 mm`), robot arc following, and SmoothWall114 single-wall General Contact. The field is evaluated in the authoritative Abaqus global frame.

## 3. Definition of cone-frame azimuth

At each local pipe tangent `t`, the projected +Z normal is `n`, `c0=cos(40°)t+sin(40°)n`, `e10=normalize(n-(n·c0)c0)`, and `e20=c0×e10`. The new `χ` applies one Rodrigues rotation about `t` to all three vectors (`c0,e10,e20`), then the unchanged production field formula is used. χ is therefore a cross-sectional spatial orientation, not a phase change.

## 4. Legacy regression

`test_cone_frame_azimuth_regression.py` compares the archived production server with the modified server at χ=0 over multiple poses. B, F and T maximum absolute difference is 0.0; regression PASS. The frame definition file reports unit norms, orthogonality, handedness +1, 40° bias and 30° cone invariants.

## 5. Baseline wall-contact geometry

The G6/L45 baseline exact-wall audit reproduces the established −0.532947 μm minimum signed gap in the 0.7–2.0 ms danger window. The worst baseline contact is HEAD-side (node 74 at the minimum-gap frame; first strong aggregate CPRESS is node 57 at 1.300 ms). Aggregate CPRESS remains explicitly a General Contact quantity, not wall-only.

## 6. Zero-cost azimuth scan

The scan used all 31 G6/L45 rigid-body frames, the 217 true exterior robot nodes, and the 430-face/764-facet authoritative SmoothWall triangle mesh. For each χ=0…355° in 5° steps, and then 125…145° in 1° steps, the baseline robot COM trajectory was held fixed and the complete rigid robot orbit was rotated about the local tangent. No Abaqus solve was used.

## 7. Scan objective and selection

The danger-window objective used minimum signed gap, 5th percentile, mean of the lowest 10%, negative-gap exposure, gap<20 μm exposure, and maximum wall-normal closing proxy, with separate HEAD/TAIL/SIDE values. χ=134° is the local optimum for minimum clearance: `min gap=+0.286997 mm`, an improvement of `287.530 μm` over χ=0, with zero negative-gap frames and no endpoint becoming worse. It passes Gate A. The closing-speed proxy does not improve (it increases by about 9.7%), so the selection is a clearance/geometry result, not a claim of reduced impact speed.

## 8. HEAD/TAIL clearance

At χ=134° in the proxy, danger-window minima are HEAD +0.287134 mm, TAIL +0.286997 mm and SIDE +0.338406 mm; all are positive. The worst contact exposure is therefore removed in the proxy without transferring a negative gap to the opposite endpoint. χ=134° is the only candidate selected for dynamic testing.

## 9. Experimental interpretation

χ rotates the entire biased rotating-field cone in the pipe cross-section while preserving 30 Hz, 10 mT, 30° cone angle, Bias40, phase law, gradient, mass, inertia and contact. It is an open-loop spatial reorientation intended to change which HEAD/TAIL surface approaches the wall.

## 10. 3 ms dynamic probe

`Wobble_F30_ConeFrameAzimuth_Chi134_WallOn_Free_003` completed normally after datacheck. It used the unchanged G6/L45 deck and no directional tensor pulse. Stable time increment remained about 1.056×10⁻⁷ s; 31 frames were written.

## 11. First-impact momentum

The candidate first aggregate CPRESS threshold event occurs near 1.1 ms, and its first-impact signed tangential velocity change is −0.783 mm/s versus −3.821 mm/s for G6 in the same coarse window, a reduction of about 79.5%. This is a velocity proxy, not a replacement for a surface-resolved wall impulse. The resolved magnetic/nonmagnetic momentum accounting remains the required source of impulse interpretation; the retired cumulative `Jcontact/Jmag` ratio is not used.

## 12. Exact wall clearance

The dynamic χ=134° minimum exact Pipe_WALL_HELPER gap is **+8.661 μm** at 2.700 ms (node 141, wall element 441), versus baseline −0.533 μm at 1.500 ms. This passes the −1 μm dynamic target and improves the baseline clearance. Aggregate CPRESS nevertheless reaches 7.427 MPa at 2.000 ms (node 74), because it includes CEL-fluid contact as well as robot-wall contact.

## 13. Contact-location change

The baseline strong aggregate peak is node 57 at 1.300 ms; χ=134° moves the aggregate peak to node 74 at 2.000 ms and the exact-wall minimum to node 141 at 2.700 ms. Thus the spatial orbit was changed, but the dynamic contact load was not eliminated and the aggregate peak became larger. This is a real geometry change, not a visualization-only effect.

## 14. Wobble preservation

The true long-axis proxy (`axis_tangent_angle_deg`) spans 26.566° in G6 and 18.224° for χ=134°, a 31.4% reduction, just beyond the 30% preservation gate. UR1 peak-to-peak falls from 0.847 rad to 0.239 rad (−71.7%). Therefore the 3 ms candidate does not preserve the user-approved high-amplitude wobble.

## 15. Canonical translation

The candidate remains net reverse: Δs −0.147898 mm and final Vt −132.809 mm/s, compared with −0.316767 mm and −251.077 mm/s for G6. Forward fraction is only 3.2%; this is not forward recovery. The improvement in reverse speed is meaningful but does not satisfy a forward-motion gate.

## 16. Dynamic comparison

| metric | G6/L45 χ=0 | χ=134° | change |
|---|---:|---:|---:|
| Δs (mm) | −0.316767 | −0.147898 | 53.3% less reverse magnitude |
| final Vt (mm/s) | −251.077 | −132.809 | 47.1% less reverse magnitude |
| first-impact ΔVt (mm/s) | −3.821 | −0.783 | 79.5% lower magnitude |
| aggregate CPRESS max (MPa) | 0.450 | 7.427 | increased; not wall-only |
| exact minimum gap (μm) | −0.533 | +8.661 | improved |
| UR1 pp (rad) | 0.847 | 0.239 | −71.7% |
| true-axis range (deg) | 26.566 | 18.224 | −31.4% |

## 17. Decision

The offline χ geometry gate passed, and the dynamic probe confirmed exact-wall clearance and reduced reverse displacement. However, the candidate failed the wobble-preservation gate and did not produce forward motion; aggregate CPRESS also increased. The correct decision is **`RIGID_CONE_FRAME_AZIMUTH_INSUFFICIENT`**. Do not continue to 8.333 ms and do not run χ refinements under this branch.

## 18. Exactly one next step

The next permitted architecture is an offline Maxwell/vector audit of **elliptically polarized wobble**: keep 30 Hz, cone/field scale and production wall/contact fixed, but independently vary the two transverse field amplitudes to reduce the dangerous wall-normal component while retaining one high-amplitude wobble direction. No new Abaqus job is submitted under the present report.

## Reproducibility and files

Scripts: `scan_cone_frame_azimuth_geometry.py`, `test_cone_frame_azimuth_regression.py`, `run_cone_frame_azimuth_chi134_3ms.ps1`, `analyze_cone_frame_azimuth_probe.py`. Outputs include `cone_frame_current_definition.csv`, `cone_frame_azimuth_geometry_scan.csv`, `cone_frame_azimuth_rank.csv`, `chi_candidate_exact_gap_history.csv`, `chi_candidate_head_tail_clearance.csv`, `chi_candidate_closing_speed_proxy.csv`, dynamic CSVs, plots and the fixed-camera GIF. No ODB or large raw telemetry is part of the GitHub bundle.
