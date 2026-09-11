# 30 Hz local gradient length-scale report

Production equation uses G as localized field amplitude in T and L as a fixed Gaussian length scale in mm; dB/ds is in T/m. No force multiplication shortcut was used.

Real G6 trajectory `delta_s = s_driver - s_robot`: min 5.0857, median 5.1013, mean 5.1755, max 5.4199 mm. Instantaneous Lopt=|delta_s|/sqrt(2): min 3.5962, median 3.6072, mean 3.6596, max 3.8325 mm.

Trajectory-specific fixed optimum by post-impact impulse is L=3.750 mm. At G=6 mT this gives mean/post-impact Ft 394.81/458.01 uN versus L45 7.10/8.42 uN (mean enhancement 55.9x, post-impact impulse enhancement 54.4x). Torque is unchanged in the production force-only gradient implementation.

## G/L budget

 G_mT  L_mm  mean_Ft_uN  postimpact_mean_Ft_uN  Ft_impulse_Ns  positive_fraction
  6.0  3.75  394.805692             458.006611       0.000001                1.0
 10.0  3.75  658.009487             763.344352       0.000002                1.0
 15.0  3.75  987.014231            1145.016528       0.000003                1.0
 20.0  3.75 1316.018975            1526.688704       0.000004                1.0

## Recovery estimate

 recovery_time_ms  target_Vt_mm_s  Vt_reference_mm_s  J_required_Ns  mean_force_required_uN  estimated_G_mT_at_Lbest  within_20mT                          label
            8.333             0.0        -250.938694       0.000003              367.245271                 4.811004         True MINIMUM_IMPULSE_SCALE_ESTIMATE
            8.333             5.0        -250.938694       0.000003              374.562701                 4.906864         True MINIMUM_IMPULSE_SCALE_ESTIMATE
           16.667             0.0        -250.938694       0.000003              165.450448                 2.167442         True MINIMUM_IMPULSE_SCALE_ESTIMATE
           16.667             5.0        -250.938694       0.000003              168.747078                 2.210629         True MINIMUM_IMPULSE_SCALE_ESTIMATE
           33.333             0.0        -250.938694       0.000003               78.829735                 1.032689         True MINIMUM_IMPULSE_SCALE_ESTIMATE
           33.333             5.0        -250.938694       0.000003               80.400432                 1.053266         True MINIMUM_IMPULSE_SCALE_ESTIMATE

At the optimized fixed L, the minimum-impulse estimates are within the 20 mT guard (about 4.8 mT for neutral recovery by 8.333 ms, decreasing for longer windows). Classification: `SHORT_LENGTH_SCALE_GRADIENT_CAN_REACH_FORWARD_IMPULSE_SCALE`. This remains a minimum scale estimate and requires one dynamic 8.333 ms validation. Changing L represents a different coil/magnet spacing or localized gradient hardware, not a software knob in the experiment.



## 6. Dynamic 8.333 ms validation (G=6 mT, L=3.75 mm)

The production Gaussian gradient was tested in the Wall-ON CEL model with no other physical parameter change.

|case|job|duration_ms|s_start_mm|s_end_mm|delta_s_mm|Vt_mean_mm_s|Vt_final_mm_s|Vt_min_mm_s|Vt_max_mm_s|positive_Vt_fraction|precontact_mean_Vt_mm_s|Ft_mean_uN|Ft_min_uN|Ft_max_uN|Jmag_t_Ns|contact_onset_ms|CPRESS_max_MPa|CPRESS_peak_time_ms|min_gap_um|gap_time_ms|UR1_pp_rad|UR2_pp_rad|UR3_pp_rad|max_abs_axis_angle_deg|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|G6_L45|WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003|3.000000026077|13.814214119926389|13.497447083509444|-0.3167670364169446|-110.29240931029734|-250.93869388542558|-250.93869388542558|0.0|0.0|-2.8722134361426983|5.738406866386683|3.616350086833513|7.871011356556878|1.6581801681953173e-08|1.3000426115468|0.4501559138298034|1.3000426115468|-0.5329466905074001|1.5000704443082|0.8471474647521973|0.4959879219532013|0.7033452689647675|29.42167574843427|
|Local_G6_L3p75|Wobble_F30_LocalForwardGradient_WallOn_Free_0083|8.3330003544688|13.814214119926389|11.938640317218582|-1.875573802707807|-224.1481768277101|-283.1547073272204|-324.4070521390255|7.980212371044799|0.12941176470588237|-0.7484514734909072|183.23985770596957|32.10319796431655|452.3954267879|1.5114930546666127e-06|1.4000082155689|13.97852897644043|1.5000303974375|-5.4048796443144|1.5000303974375|2.1803243160247803|0.688329815864563|1.4151552468538284|32.302971453604385|

The local profile ends at **Δs = -1.8756 mm** and **Vt = -283.2 mm/s**, while the G6/L45 baseline ends at **Δs = -0.3168 mm** and **Vt = -250.9 mm/s**. The candidate CPRESS peak is 13.98 MPa and exact minimum signed gap -5.405 µm.

### Stage-wise canonical velocity

|case|window|start_ms|end_ms|mean_Vt_mm_s|median_Vt_mm_s|positive_fraction|
|---|---|---|---|---|---|---|
|G6_L45|precontact|0.0|1.3000426115468|-2.8722134361426983|-1.0466513660465988|0.0|
|G6_L45|impact|1.1000426115467998|1.5000426115468|-25.21902374102859|-23.251920646516794|0.0|
|G6_L45|postimpact|1.5000426115468|2.8000426115468002|-202.2510189622076|-197.44337650144138|0.0|
|G6_L45|late|2.8000426115468002|3.000000027077|-232.39603326634972|-223.4151374125038|0.0|
|Local_G6_L3p75|precontact|0.0|1.4000082155689|-0.7484514734909072|0.9933123337637098|0.7857142857142857|
|Local_G6_L3p75|impact|1.2000082155689|1.6000082155689002|-46.65373283335789|-20.951234636303163|0.0|
|Local_G6_L3p75|postimpact|1.6000082155689002|2.9000082155689|-184.59752758090113|-189.14906668879527|0.0|
|Local_G6_L3p75|late|2.9000082155689|8.333000355468801|-294.43167730811064|-300.51296948351074|0.0|

The offline force increase is real, but concentrating it near the measured positive Δs (~5.1 mm) injects a large force during the contact-sensitive initial state. In the coupled solve the first-contact CPRESS and penetration worsen, and no post-impact forward recovery appears. UR1 peak-to-peak remains 2.180 rad, so this is not loss of magnetic torque; it is a translation/contact failure.

**Decision: `LOCAL_GRADIENT_PROFILE_DOES_NOT_RECOVER_FORWARD_MOTION`.** Do not submit a 16.667 ms continuation or increase G / further shrink L in this scalar Gaussian architecture. Preserve G6 as the softer-contact wobble reference and move to impulse-cancellation or another spatial/temporal field architecture.
