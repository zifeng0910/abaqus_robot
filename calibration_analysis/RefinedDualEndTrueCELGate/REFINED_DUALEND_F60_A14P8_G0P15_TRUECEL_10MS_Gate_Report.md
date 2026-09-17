# Refined Dual-End 10 ms TRUE-CEL Gate

This was one strict-CEL dynamics run. It completed 10.0 ms with real water EVF,
initial `+10 mm/s` flow, Magpylib socket loading, General Contact, and no
ReducedHydro load.

The accepted physical mesh change is a smooth one-sided 0.0425 mm HEAD profile
increase over the final 25% of the body, plus a +0.015 mm radial pose offset.
It is an actual node-coordinate geometry change, not a label or measurement hack.

## Gate result

Final status: `TRUECEL_GATE_FAILED`

Failed hard gates: `energy_drift_below_0p01_Nmm, near_robot_fluid_velocity_below_10000_mm_s`.

| quantity | value |
|---|---:|
| ETOTAL max absolute drift | 0.0109777253 N mm |
| near-robot fluid velocity max | 75573.6256 mm/s |
| robot total General Contact resultant max | 0.00119379731 N |
| minimum stable dt | 4.63521e-08 s |
| rocking min / max | -6.638598 / 12.738038 deg |
| HEAD / TAIL / BODY minimum gap | 0.013468 / 0.001122 / 0.001443 mm |
| HEAD / TAIL / BODY total-General-Contact episodes | 2 / 1 / 2 |
| longest simultaneous HEAD+TAIL wall proximity | 0.000000 ms |
| delta_s | +0.00920900187 mm |
| Jmag_s | +1.95028265e-08 N s |
| J robot-total-contact_s | +7.14629682e-08 N s |

Abaqus directly provides total robot General Contact fields (`CNORMF`,
`CSHEARF`, `COPEN`). It does not pair-isolate robot-wall from robot-fluid in
this ODB, so no momentum residual or geometry-filtered proxy is labeled as a
direct wall-only force. The energy hard gate failed at the endpoint; therefore
no propulsion interpretation and no FWD/REV pair are permitted.
