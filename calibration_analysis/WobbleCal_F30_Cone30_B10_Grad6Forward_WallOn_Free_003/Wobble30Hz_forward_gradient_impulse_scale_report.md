# 30 Hz forward-gradient impulse budget

## 1. Why the G6 3 ms result is not a gradient failure

The 3 ms window is only 0.09 of a 30 Hz cycle and is dominated by the initial rotation/contact transient. G6 does not cancel that initial reverse transient; it does, however, preserve the wobble and substantially soften contact.

## 2. Canonical tangential velocity and stage split

Forward is increasing canonical centerline arclength. Vt is computed by projecting the ODB RP global velocity (V1,V2,V3) onto the continuously interpolated local tangent; continuous `ds/dt` is retained as a cross-check. No legacy arc steps or absolute-value force are used.

|job|delta_s_mm|Vt_final_mm_s|Vt_min_mm_s|Vt_max_mm_s|positive_Vt_fraction|Ft_mean_uN|Jmag_t_Ns|contact_onset_ms|CPRESS_max_MPa|min_gap_um|
|---|---|---|---|---|---|---|---|---|---|---|
|WobbleCal_F30_Cone30_B10_WallOn_Free_003|-0.3160894526006839|-221.23341332764943|-226.87018764166615|0.0|0.0|2.872549570291945|8.299452226482937e-09|1.5000409912317|8.98662281036377|-2.7268386971292|
|WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003|-0.3167670364169446|-250.93869388542558|-250.93869388542558|0.0|0.0|5.738406866386683|1.6581801681953173e-08|1.3000426115468|0.4501559138298034|-0.5329466905074001|

|case|job|window|start_ms|end_ms|mean_Vt_mm_s|median_Vt_mm_s|
|---|---|---|---|---|---|---|
|G3|WobbleCal_F30_Cone30_B10_WallOn_Free_003|W0_precontact|0.0|1.5000409912317|-6.2520493136808275|-1.9501835669349055|
|G3|WobbleCal_F30_Cone30_B10_WallOn_Free_003|W1_first_contact|1.3000409912317|1.7000409912317|-155.61911314632914|-196.65966472136316|
|G3|WobbleCal_F30_Cone30_B10_WallOn_Free_003|W2_postimpact|1.7000409912317|2.5|-195.10407026803347|-186.43632116248241|
|G3|WobbleCal_F30_Cone30_B10_WallOn_Free_003|W3_late|2.5|3.000000027077|-223.31525788223655|-224.83654332834416|
|G6|WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003|W0_precontact|0.0|1.3000426115468|-2.8722134361426983|-1.0466513660465988|
|G6|WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003|W1_first_contact|1.1000426115467998|1.5000426115468|-25.21902374102859|-23.251920646516794|
|G6|WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003|W2_postimpact|1.5000426115468|2.5|-195.77237378476565|-194.45106497867337|
|G6|WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003|W3_late|2.5|3.000000027077|-228.121268076682|-223.9104779008374|

At 3 ms both cases are still moving rapidly toward −s (G3 -221.2 mm/s, G6 -250.9 mm/s); neither is in case C (instantaneous Vt already positive). The reverse motion is therefore an initial inertia/contact-dominated transient.

## 3. Magnetic impulse and sensitivity

Measured magnetic tangential impulse is Jmag,t=8.299e-09 N·s (G3) and 1.658e-08 N·s (G6). The extra G6 impulse is 8.282e-09 N·s, while the G6 3 ms tangential momentum magnitude is about 2.509e-06 N·s; the added impulse is only 0.33% of that magnitude.

Two-point measured sensitivity is K_F=0.955 µN per mT, with Ft≈0.955·G(mT)+0.007 µN. This is a short-range estimate only. CPRESS falling 95% does not prove tangential wall impulse stayed constant; the available data only support the resolved nonmagnetic remainder interpretation.

## 4. Minimum impulse-scale gradient estimate

|target_total_time_ms|available_postimpact_time_ms|Vt_reference_mm_s|target_Vt_mm_s|J_required_Ns|mean_force_required_uN|estimated_gradient_amplitude_mT|within_existing_20mT_guard|interpretation|
|---|---|---|---|---|---|---|---|---|
|8.333|6.8329590087683|-250.93869388542558|0.0|2.509386938854256e-06|367.24747443005583|384.43028826622907|False|minimum impulse scale estimate; excludes CEL drag and repeated wall impulses|
|8.333|6.8329590087683|-250.93869388542558|5.0|2.559386938854256e-06|374.56494844619414|392.0902726200659|False|minimum impulse scale estimate; excludes CEL drag and repeated wall impulses|
|16.667|15.1669590087683|-250.93869388542558|0.0|2.509386938854256e-06|165.45089476430527|173.18817937887886|False|minimum impulse scale estimate; excludes CEL drag and repeated wall impulses|
|16.667|15.1669590087683|-250.93869388542558|5.0|2.559386938854256e-06|168.7475345172771|176.6391255485965|False|minimum impulse scale estimate; excludes CEL drag and repeated wall impulses|
|33.333|31.8329590087683|-250.93869388542558|0.0|2.509386938854256e-06|78.82983602507868|82.51263298241608|False|minimum impulse scale estimate; excludes CEL drag and repeated wall impulses|
|33.333|31.8329590087683|-250.93869388542558|5.0|2.559386938854256e-06|80.4005351230239|84.15685207882983|False|minimum impulse scale estimate; excludes CEL drag and repeated wall impulses|

These are `MINIMUM_IMPULSE_SCALE_ESTIMATE` values: they ignore CEL drag, changing orientation, and repeated contacts, so they are not a full trajectory prediction. Every neutral/+5 mm/s estimate is above the existing 20 mT guard (approximately 83–392 mT over 8.333–33.333 ms).

## 5. Decision

`FORWARD_GRADIENT_WITHIN_CURRENT_LIMIT_TOO_WEAK`. No ≤20 mT candidate has a credible impulse scale to recover the observed −s momentum within a quarter, half, or full 30 Hz cycle. Therefore no new 8.333 ms Abaqus job is submitted in this stage. The next physically meaningful work is a hardware/field-architecture change (or an explicitly approved relaxation of the 20 mT guard), while preserving the 30 Hz wobble operating point.
