# TAIL No-Rebound zeta=0.80, 8 ms TRUE-CEL Gate

Final classification: `TRUECEL_CONTACT_GATE_FAILED_REBOUND`

The sole dynamics run completed successfully. TAIL contact was diagnosed from
the actual transformed TAIL mesh geometry, wall-normal gap velocity, and direct
robot-total General Contact fields. The nominal 0.005 mm plane-gap threshold
misses the first directly loaded event because the wall is faceted; the event
gate therefore requires direct contact force and a documented 0.010 mm
faceted-wall neighborhood. The ODB does not pair-isolate robot-wall from
robot-fluid contact.

- First touch: 2.500029979273677 ms
- First release: 2.750028856098652 ms
- Second touch in the same half-cycle: 5.250026471912861 ms
- Same-half-cycle TAIL contacts: 2
- Impact / rebound gap velocity: -215.24936499399405 / 94.93977238027678 mm/s
- R_v: 0.44106876869497763
- Failed CEL numerical gates: near_robot_fluid_velocity_below_10000_mm_s

No FWD/REV interpretation is made.
