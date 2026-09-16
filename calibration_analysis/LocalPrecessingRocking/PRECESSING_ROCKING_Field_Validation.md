# ROBOT_LOCAL_PRECESSING_ROCKING field validation

Status: **PASS**

## Command

- `A_rock = 14.343111711438091 deg`, reused exactly from the solved LIGHTCONTACT case.
- `f_rock = 10 Hz`, `f_prec = 2.5 Hz`, `B0 = 10 mT`, `G = 0`.
- RouteA gauge: `-61.37284757596327 deg` from the production PT `e1` axis.
- `psi0 = -83.87284757596327 deg`, so `psi(25 ms)` equals the RouteA gauge at the first positive rocking peak.
- `B = B0 [cos(alpha)c_hat + sin(alpha)n_p]`, where `alpha = A_rock sin(2 pi f_rock t)` and `n_p = cos(psi)e1 + sin(psi)e2`.

## Peak directions

| Time (ms) | alpha (deg) | Plane azimuth from e1 (deg) | Effective signed direction from e1 (deg) | From RouteA (deg) |
|---:|---:|---:|---:|---:|
| 25 | +14.343111711438091 | -61.372847576 | -61.372847576 | 0 |
| 75 | -14.343111711438091 | -16.372847576 | +163.627152424 | -135 |
| 125 | +14.343111711438091 | +28.627152424 | +28.627152424 | +90 |
| 175 | -14.343111711438091 | +73.627152424 | -106.372847576 | -45 |

The rocking sign is included in the effective direction. The four extrema therefore do not collapse onto a fixed line.

## Gates

- Maximum `|B|-10 mT` error: `3.469446951953614e-18 T`.
- Command extrema: exactly `-14.343111711438091 .. +14.343111711438091 deg` on the sampled grid.
- Precession advance over 200 ms: `180 deg`, equivalent to `2.5 Hz`.
- Minimum adjacent PT-frame dot product over the reachable straight-case neighborhood: `0.9999999999999999`.
- Maximum adjacent 0.1 ms field step: `1.572887693878714e-05 T`; no discontinuity was found.
- Transverse path singular-value ratio: `0.940525312247557`, rejecting a straight-line path.
- Transverse radius range: approximately `0 .. 2.477280702 mT`, rejecting a constant-radius cone.
- Existing `LEGACY_DRIVER_FRAME`, `ROBOT_LOCAL_TANGENT`, and `ROBOT_LOCAL_ROCKING` regressions pass. The live old server and changed vendored server agree exactly for legacy `B/F/T`.

The implementation gate is passed. Abaqus must still be run only for the single authorized case after the vendored server is committed, pushed, backed up, and deployed with exact SHA equality.
