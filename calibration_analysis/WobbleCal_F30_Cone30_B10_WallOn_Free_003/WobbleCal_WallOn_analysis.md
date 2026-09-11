# Wall-on diagnostic status

## COM-fixed Wall-on probe

`WobbleCal_COMFixed_F30_Cone30_B10_WallOn_003Rerun` completed 3 ms with
robot--wall contact enabled. The robot translations were still fixed, so this
was a wall-reaction diagnostic, not a transport test.

- No initial node-face or edge-edge overclosure.
- Robot--wall contact nodes: 0 for every frame.
- Maximum robot-fluid CPRESS: 0.182 MPa at about 2.000 ms.
- UR ranges: UR1 0 to 0.551 rad, UR2 -0.405 to 0 rad, UR3 -0.113 to 0.324 rad.
- Therefore the initial COM-fixed pose did not reach the helper wall within
  3 ms; adding Wall=ON alone did not yet alter the rotation response.

The earlier `...WallOn_003` attempt stopped at about 0.9 ms due to a terminal
window-close event and is treated as an invalid partial run, not as a physical
failure.

## Free-translation Wall-on probe

`WobbleCal_F30_Cone30_B10_WallOn_Free_003` is the next probe. It removes the
`RP_ROBOT, 1, 3` translational constraint while retaining the same 30-Hz,
10-mT analytic field, 3-mT gradient (`0.003 T`), CSF, and Wall-on contact.
It is running for 3 ms and is the first test in this branch that can reveal
whether wall reaction reverses the robot's net transport direction.

## Free-translation result

The free-translation probe completed 3 ms and generated the fixed-camera GIF.
The robot moved backward along the authoritative centerline:

- `s_robot`: 13.8142 to 13.4981 mm (`Delta s=-0.3161 mm`)
- magnetic driver arc: +0.0174 mm over the same window
- peak COM speed: about 511.5 mm/s; final speed about 294.2 mm/s
- maximum CPRESS: 8.99 MPa at 2.900 ms
- exact minimum wall gap: -0.002727 mm (-2.727 um), node 55 / wall element 442
- wall-contact nodes: maximum 4, active near the 2.9 ms impact

The reconstructed tangential impulse ratio is approximately
`|J_wall,t|/|J_mag,t| = 1.00e3`; this is a diagnostic ratio until a full
momentum-balance audit is applied, but its timing and sign are consistent with
the observed post-impact reversal. The centerline in the GIF was plotted after
the validated Abaqus-frame transform, so the backward motion is not a
centerline visualization artifact.
