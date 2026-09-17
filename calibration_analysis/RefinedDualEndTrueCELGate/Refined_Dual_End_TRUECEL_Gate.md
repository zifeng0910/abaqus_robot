# Refined Dual-End TRUE-CEL Gate

## Superseded Geometry Decision

The original formal geometry decision below was superseded by explicit user
acceptance of the bridge-free near-optimal candidate. Its correct historical
classification was `GEOMETRY_GATE_FORMALLY_MISSED / CEL_NOT_RUN`, not a CEL
numerical failure. The accepted geometry was subsequently run once at
`A_main=14.8 deg`; the authoritative result is
`REFINED_DUALEND_F60_A14P8_G0P15_TRUECEL_10MS_Gate_Report.md`.

The asymmetry is shape-driven rather than pose-driven. In the exact source
mesh, the pipe center and robot COM are essentially coincident, while the HEAD
and TAIL maximum radial envelopes are approximately 0.389125 mm and 0.407500
mm, respectively.

## Bounded refinement

The audit retained the requested size bounds and tested only small changes near
GEO-C. In addition to axial length and a small `n_routeA` centering correction,
it tested a smooth one-sided HEAD profile correction over the final 25% of the
robot length. This correction addresses the measured local shape asymmetry; it
does not exchange HEAD/TAIL semantics.

The bridge-free near-candidate later accepted for dynamics was:

| Quantity | Value |
|---|---:|
| L | 2.600 mm |
| nominal D | 0.815 mm |
| L/D | 3.190184 |
| radial COM offset | +0.015 mm along `n_routeA` |
| one-sided HEAD profile correction | 0.0425 mm along `-n_routeA` |
| HEAD touch | 14.387838 deg |
| TAIL touch | 13.423313 deg |
| threshold difference | 0.964525 deg |
| initial minimum gap | 0.245067 mm |
| same-phase static BOTH bridge | no |

This candidate misses the requested HEAD upper limit by 0.087838 deg. Increasing
the same correction to 0.04375 mm reduces HEAD touch to 14.340641 deg, but a
same-phase static BOTH bridge then appears. At 0.045 mm, HEAD touch reaches
14.293452 deg and the threshold difference reaches 0.870139 deg, but the bridge
remains.

The transition is geometric: for centered rigid rocking in the straight lumen,
bringing both endpoint thresholds below the same 14.343111711438091 deg command
also brings the opposite endpoints against opposite walls at an extremum. The
requested reachability and no-static-bridge constraints therefore do not
overlap within this tightly bounded refinement family.

## Current Decision

See the authoritative 10 ms gate report. No FWD/REV pair was run.
