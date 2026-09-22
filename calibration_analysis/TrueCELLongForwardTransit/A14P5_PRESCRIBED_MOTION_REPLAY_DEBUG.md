# A14P5 prescribed-motion replay debug

**UR_COMPONENT_FINITE_ROTATION_REPLAY_ERROR**

**KINEMATIC_REPLAY_VALIDATED**

The invalid refined replay remains labelled `FIXED_TRAJECTORY_REPLAY_INVALID` and is excluded from CEL-resolution inference. Its restart was reproducible; the failure was the rotational boundary implementation.

## Why the old replay failed

The submitted deck used six unit-reference `*BOUNDARY` cards with `TYPE=DISPLACEMENT`, default `OP=MOD`, and `STEP TIME` tabular amplitudes. Thus every imposed value was exactly `1.0 x amplitude`; scaling and time basis were correct. Translation reproduced correctly, while the three independently prescribed rotational displacement components did not reproduce the target finite-rotation path. Abaqus/Explicit composes finite rotation increments; the ODB `UR` is the equivalent total axis-angle vector, not three independently enforceable absolute orientation coordinates for a time-varying rotation axis.

Componentwise 5-us rotation-vector interpolation versus quaternion SLERP differed by only 6.196e-07 deg, so sampling/interpolation cannot explain the old 3.211751-deg error. No duplicate rotational BC, non-unit scale, time reset, initial-pose error, or RP-selection error was found.

## Target UR/VR consistency and frame

Log-map differentiation of target orientation gives {"spatial_global": {"rms_rad_s": 0.23356812254880918, "p99_rad_s": 0.9579078315085312, "max_rad_s": 1.747370905472966}, "body_material": {"rms_rad_s": 23.877012827908384, "p99_rad_s": 63.65792515852893, "max_rad_s": 64.50825438473866}, "selected": "GLOBAL_SPATIAL"}. The spatial/global convention is authoritative; the body/material candidate is strongly inconsistent.

## Replacement replay

One robot-only Abaqus/Explicit harness prescribed target global `V(t)` and spatial/global `VR(t)` at the same RP. It contained the exact robot part/RP rigid-body definition but no CEL, wall, contact, fluid, magnetic load, or external load. Initial gate: {"target_U": [0.0, 0.0, 0.0], "replay_U": [0.0, 0.0, 0.0], "position_error_mm": 0.0, "target_UR": [0.0, 0.0, 0.0], "replay_UR": [0.0, 0.0, 0.0], "orientation_error_deg": 0.0}.

Validation metrics: {"translation_vector_mm": {"full": {"rms": 1.4567116191125543e-08, "p95": 3.470364848712539e-08, "p99": 5.604556967450068e-08, "max": 7.631342676065967e-08, "max_time_ms": 9.371000342071056}, "critical_8p333333_to_10p682021_ms": {"rms": 2.5557004190589047e-08, "p95": 5.365044459155729e-08, "p99": 6.759481675265382e-08, "max": 7.631342676065967e-08, "max_time_ms": 9.371000342071056}}, "orientation_geodesic_deg": {"full": {"rms": 2.3515912312848086e-06, "p95": 4.971119543369419e-06, "p99": 9.63270391494992e-06, "max": 1.5209011955989004e-05, "max_time_ms": 8.92999954521656}, "critical_8p333333_to_10p682021_ms": {"rms": 4.26897047156563e-06, "p95": 9.581093713024603e-06, "p99": 1.3016094413089953e-05, "max": 1.5209011955989004e-05, "max_time_ms": 8.92999954521656}}, "axial_velocity_mm_s": {"full": {"rms": 0.0016877603551075458, "p95": 0.002398821662189789, "p99": 0.006449217276179526, "max": 0.0294321687912659, "max_time_ms": 9.166000410914421}, "critical_8p333333_to_10p682021_ms": {"rms": 0.0032283485044688532, "p95": 0.0060215499558408546, "p99": 0.012621390395576927, "max": 0.0294321687912659, "max_time_ms": 9.166000410914421}}, "rocking_angular_velocity_rad_s": {"full": {"rms": 0.0461562491404919, "p95": 0.007934234913851842, "p99": 0.03287306771075891, "max": 1.990709215159248, "max_time_ms": 1.8360000103712082}, "critical_8p333333_to_10p682021_ms": {"rms": 0.014796667648980903, "p95": 0.012264608027860304, "p99": 0.026967707916496844, "max": 0.24521264817016453, "max_time_ms": 9.305999614298344}}}

Critical former force-event timestamps: [{"time_ms": 9.305999614298344, "event": "peak_negative_robot_CEL", "position_error_mm": 3.3763259349847916e-08, "orientation_error_deg": 2.194175049737399e-06, "v_s_error_mm_s": 0.0017681605229629643, "omega_rock_error_rad_s": -0.24521264817016453}, {"time_ms": 9.331000037491322, "event": "major_force_change", "position_error_mm": 2.4723174240625943e-09, "orientation_error_deg": 2.1790272152928247e-06, "v_s_error_mm_s": -0.0064077107615584565, "omega_rock_error_rad_s": 0.025447168497958424}, {"time_ms": 9.320000186562538, "event": "major_force_change", "position_error_mm": 1.2289211906054413e-08, "orientation_error_deg": 5.4778903772835034e-06, "v_s_error_mm_s": -0.0064433760004915935, "omega_rock_error_rad_s": 0.02632958643212601}, {"time_ms": 9.320000186562538, "event": "major_force_change", "position_error_mm": 1.2289211906054413e-08, "orientation_error_deg": 5.4778903772835034e-06, "v_s_error_mm_s": -0.0064433760004915935, "omega_rock_error_rad_s": 0.02632958643212601}]

The required gates pass: position max <1 micron; geodesic orientation max <0.02 deg; critical-window axial-velocity p99 <0.05 mm/s; critical-window rocking-angular-velocity p99 <0.5 rad/s. Maximum instantaneous angular-velocity error is reported and occurs at initialization/output staggering; it does not control the critical-window p99 gate.

No CEL model was launched. The validated driver is global `V + VR`; rebuilding coarse and refined CEL replays is a separate next task.
