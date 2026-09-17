# FAST screen performance

- Magnetic table precompute: 1.324737 s (reused validated table)
- Abaqus wallclock: 1455.775 s
- Physical simulated time: 16.666668 ms
- Wallclock per physical millisecond: 87.3465 s/ms
- CPUs: 1
- Completed Explicit increments / magnetic lookups: 83,334 / 83,334
- Socket calls: 0

The 10 us contact-field output generated 1,668 ODB frames and dominated this run's I/O. Therefore its wallclock is not directly comparable with the sparse-output 0.5 ms backend benchmark (19 s live socket versus 5 s precomputed-table solver time). Even with dense output, this 16.667 ms FAST run was substantially cheaper than the historical strict-CEL context (5959.387 s for 10 ms; estimated 9932.3 s for 16.667 ms).

The increment CSV retained only its final row because the solved subroutine treated `VEXTERNALDB lOp=3` as final close. The ODB contains 83,335 strictly increasing increment-history samples and is authoritative for the 83,334 completed lookup count. The maintained subroutine now closes only at `lOp=6`; this diagnostic-only fix did not alter or rerun dynamics.
