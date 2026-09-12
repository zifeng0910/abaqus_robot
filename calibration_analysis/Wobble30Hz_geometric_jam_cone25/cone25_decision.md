

## Cone25 decision addendum

The authorized cone25 run completed at 8.333 ms and was post-processed with the same 50 µs telemetry/ODB cadence and exact SmoothWall114 triangle audit. The comparison is in `cone25_vs_cone30_summary.csv`.

- Cone30 baseline robot arc change: **-2.0000 mm**; cone25: **-2.1000 mm** (both negative, so the motion remains reverse in the canonical arc convention).
- Mean tangential speed proxy: cone30 0.0039 mm/s; cone25 0.0037 mm/s.
- True-axis UR ranges (cone25): UR1 2.578 rad, UR2 0.726 rad, UR3 1.628 rad.
- Exact minimum gap: cone30 -2.059 µm; cone25 -2.318 µm. Both remain below zero; cone25 is not a no-penetration pass.
- Aggregate CPRESS: cone30 6.868 MPa; cone25 6.524 MPa.
- Cone25 late post-impact window: phase-equivalent rate 5.042 Hz, R=0.168, wobble RMS 128.284 rad/s; classification remains **B_SUSTAINED_POST_IMPACT_NONMAGNETIC_ROTATIONAL_STALL**.

Conclusion: reducing cone30 to cone25 did not remove the geometric wedge or restore forward motion. It slightly increased the reverse arc excursion and made the exact minimum gap more negative, although CPRESS decreased modestly. The zero-cost map's cylinder sanity limit (~23°) explains why 25° can still enter the constrained region once translation and transverse offset are included. No further cone-angle or contact tuning is authorized by this stage; the next action should be selected from the mechanism classification, not from a claim that cone25 solved the jam.
