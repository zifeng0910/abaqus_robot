# 30 Hz broad-to-local gradient report



## 10. 8.333 ms dynamic result

The production broad-to-local profile was tested in the Wall-ON CEL model with no other physical parameter change.

|case|job|duration_ms|s_start_mm|s_end_mm|delta_s_mm|Vt_mean_mm_s|Vt_final_mm_s|Vt_min_mm_s|Vt_max_mm_s|positive_Vt_fraction|precontact_mean_Vt_mm_s|Ft_mean_uN|Ft_min_uN|Ft_max_uN|Jmag_t_Ns|contact_onset_ms|CPRESS_max_MPa|CPRESS_peak_time_ms|min_gap_um|gap_time_ms|UR1_pp_rad|UR2_pp_rad|UR3_pp_rad|max_abs_axis_angle_deg|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|G6_L45|WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003|3.000000026077|13.814214119926389|13.497447083509444|-0.3167670364169446|-110.29240931029734|-250.93869388542558|-250.93869388542558|0.0|0.0|-2.8722134361426983|5.738406866386683|3.616350086833513|7.871011356556878|1.6581801681953173e-08|1.3000426115468|0.4501559138298034|1.3000426115468|-0.5329466905074001|1.5000704443082|0.8471474647521973|0.4959879219532013|0.7033452689647675|29.42167574843427|
|Gated_G6_L45_to_L3p75|Wobble_F30_BroadToLocalGradient_WallOn_Free_0083_R2|8.3330003544688|13.814214119926389|11.884982149075752|-1.929231970850637|-231.39096063422505|-293.810316114757|-336.2013289450801|0.0|0.0|-2.8722134361426983|97.58846997709044|3.616350086833513|300.0149679416|8.157842440635248e-07|1.3000426115468|6.216233730316162|2.9000353533773997|-1.9200641654272|6.3000731170177|2.141020536422729|0.6895099878311157|1.3923374116420746|32.899590290931386|

The gated profile ends at **Δs = -1.9292 mm** and **Vt = -293.8 mm/s**, while G6/L45 ends at **Δs = -0.3168 mm** and **Vt = -250.9 mm/s**. Candidate CPRESS peak is 6.22 MPa and exact minimum signed gap -1.920 µm.

### Stage-wise canonical velocity

|case|window|start_ms|end_ms|mean_Vt_mm_s|median_Vt_mm_s|positive_fraction|
|---|---|---|---|---|---|---|
|G6_L45|precontact|0.0|1.3000426115468|-2.8722134361426983|-1.0466513660465988|0.0|
|G6_L45|impact|1.1000426115467998|1.5000426115468|-25.21902374102859|-23.251920646516794|0.0|
|G6_L45|postimpact|1.5000426115468|2.8000426115468002|-202.2510189622076|-197.44337650144138|0.0|
|G6_L45|late|2.8000426115468002|3.000000027077|-232.39603326634972|-223.4151374125038|0.0|
|Gated_G6_L45_to_L3p75|precontact|0.0|1.3000426115468|-2.8722134361426983|-1.0466513660465988|0.0|
|Gated_G6_L45_to_L3p75|impact|1.1000426115467998|1.5000426115468|-25.21902374102859|-23.251920646516794|0.0|
|Gated_G6_L45_to_L3p75|postimpact|1.5000426115468|2.8000426115468002|-201.1537489303685|-197.75656707195543|0.0|
|Gated_G6_L45_to_L3p75|late|2.8000426115468002|8.333000355468801|-299.27617547450893|-309.98033239844574|0.0|

The gated profile delays the strong local force until after the broad-contact phase. UR1 peak-to-peak is 2.141 rad.

**Decision:** see the quantitative gates above; no continuation is automated before review.

## 1–9. Timing, architecture, continuity and impulse audit

The G6/L45 reference reaches its first meaningful contact at **1.300 ms**
(CPRESS peak 0.450 MPa at the same output frame; the exact minimum gap is
−0.533 µm at 1.500 ms).  Using the conservative criterion CPRESS ≤0.10 MPa,
positive exact gap, and persistence for at least 0.2 ms, the first contact-clear
interval begins at **1.600 ms**.  A 0.200 ms safety delay gives
`switch_start = 1.800 ms`; the C1 smoothstep cross-fade lasts **0.500 ms** and
ends at **2.300 ms**.  The broad and local profiles are cross-faded, not added;
both use G=6 mT and the rotating 30 Hz field phase is continuous.

The offline replay has `Ft_gated > 0` for 100% of its valid samples.  The largest
sample-to-sample force change is **132.0 µN per 0.1 ms telemetry interval**;
there is no step discontinuity.  Its estimated forward impulse is
`0.421 µN·s` from switch start and `0.304 µN·s` after the transition, or only
about **0.17 / 0.12** of the G6 negative tangential momentum magnitude
(`2.509 µN·s`).  Thus the offline estimate is below the preferred 0.5 ratio.

## 11–15. Event timeline and post-switch recovery

The corrected R2 Wall-ON job completed all 8.333 ms.  It reaches its most
negative canonical velocity, **−336.2 mm/s at 6.100 ms**, and has no
`Vt=0` crossing.  After the local channel opens, `Vt` does not trend toward
zero; the late-window mean is approximately **−299.3 mm/s**.  The robot therefore
does not show a reverse-to-forward inflection in this window.

## 12–14. Early contact and wobble preservation

Early contact is substantially softer than the always-on local probe
(13.98 MPa, −5.405 µm): the gated case peaks at **6.216 MPa** and reaches
**−1.920 µm** exact gap.  This is an improvement but does not meet the
preferred `<2 MPa` / `>−1 µm` gates.  Wobble is preserved and even larger in
the short window (UR1 peak-to-peak **2.141 rad**, UR2 0.690 rad, UR3 1.392 rad),
so the field gate did not suppress the magnetic rotation; it failed to produce
forward translation without renewed wall interaction.

## 15–17. Final decision

The timed architecture is clearly better than always-on local for early
contact severity (CPRESS and penetration reduced by roughly 55–65%), but it
does not recover forward motion and still develops a post-switch contact peak
above 6 MPa.  The measured magnetic impulse is only `0.816 µN·s`, about 0.33
of the reference negative momentum, and the resolved momentum residual remains
large because contact/CEL/geometry terms dominate.

**Decision: `TIME_GATED_GAUSSIAN_TRANSLATION_ARCHITECTURE_CLOSED`.** Do not
continue to 16.667 ms and do not tune switch time, G, or L within this scalar
Gaussian translation architecture.  The single next step is to design a
phase-synchronized directional impulse-cancellation field while retaining the
current 30 Hz wobble operating point.

Generated artifacts include `broad_to_local_event_timeline.csv`,
`event_timeline_force_contact_velocity.png`, `broad_to_local_force_continuity.csv`,
and the fixed-camera validation GIF.  ODB and raw telemetry are intentionally
not part of the reproducibility bundle.
