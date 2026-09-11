# COM-fixed calibration contact-scope audit

The 40 Hz calibration deck is explicitly marked `Wall=OFF; CEL=ON`.
Its only contact inclusions are:

```text
ROBOT_SOLID-1.ROBOT_SOLID_SURF, FLUID_CEL_SURF
Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF, FLUID_CEL_SURF
```

There is no `ROBOT_SOLID ↔ PIPE_WALL_HELPER` contact inclusion.  Therefore
the robot is free to overlap the rendered pipe-wall geometry in the GIF; this
is not evidence of a missed wall-contact event.  In the extracted 40 Hz ODB,
all nonzero CPRESS values are reported on `ROBOT_SOLID-1` and the pipe-side
CPRESS maximum is zero.  The CPRESS values are robot–Eulerian-fluid interface
pressure, not robot–wall pressure.

This calibration scene is suitable only for rotational field-lock diagnosis.
It cannot pass a no-penetration or wall-guidance acceptance test.  Such a test
must be performed later with the validated Wall-On production architecture,
after the COM-fixed frequency calibration decision.
