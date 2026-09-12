# Abaqus 2025 contact-damping keyword audit

Official source: [Abaqus 2025 `*CONTACT DAMPING`](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-contactdamping.htm).

- `DEFINITION=CRITICAL DAMPING FRACTION` makes the data-line value the dimensionless normal damping fraction.
- `TANGENT FRACTION` is tangential damping divided by normal damping.
- Abaqus/Explicit default `TANGENT FRACTION` is **1.0**; Abaqus/Standard default is 0.0.
- `TANGENT FRACTION=0.0` requests no tangential contact damping.

Candidate syntax, subsequently required to pass Abaqus datacheck:

```text
*Contact Damping, definition=CRITICAL DAMPING FRACTION, tangent fraction=0.0
0.20
```
