# Coarse Motion Mode Screen Report

## Scope and decision boundary

This is a cheap, broad, human-in-the-loop screen of 32 unique rigid-robot Abaqus/Explicit cases. Every case uses direct `dt = 1e-7 s`, a global 0.10 mm C3D4 mesh, dense RP/contact history, and 0.20 magnetic cycle. The outputs are motion-mode evidence, not calibrated predictions. No composite score was computed and no winner is selected.

The 0.10 mm L2.3 baseline passed the macro-fidelity gate against the existing 0.060 mm result over 6.667 ms: maximum directed tilt 41.24 vs 41.10 deg, robot phase advance -54.82 vs -60.05 deg, and axial displacement +1.920 vs +1.953 mm. The coarse baseline remains translation-dominated wall sliding, so no 0.08 mm rerun was authorized or needed.

## Screen integrity

- Unique solved cases: 32; requested baseline labels are deduplicated through `GEO_230`.
- Mesh range: 6,118-7,302 C3D4, below the 20,000 hard limit.
- Durations: 20 Hz = 10.000 ms, 30 Hz = 6.667 ms, 40 Hz = 5.000 ms.
- Case-level Socket lifecycles were sequential and ports were confirmed closed before advancing.
- Proximity is a lightweight 20 us sampled, coarse-mesh diagnostic; it is not finalist-grade exact-CAD reconstruction.
- STEP round-trip and coarse tessellation introduce a small geometry/mass approximation; these variants are `SCREENING ONLY`.

## Motion boundaries observed

- 90 deg crossing / tumble-like candidates: none in this short screen.
- Sustained low-motion contact / stuck-like candidates: MU_150.
- Translation-heavy relative to L2.3: GEO_230, GEO_240, GRAD_0, GRAD_4, CONE_20, B_14, FREQ_20, MU_000, ZETA_015, ZETA_030, CPAR_X4, KWOB_X025, KWOB_X4, COMBO_C, COMBO_E, COMBO_F.
- Rotation/translation index above L2.3: GEO_200, GEO_210, GEO_240, GRAD_8, CONE_40, CONE_50, B_6, FREQ_40, MU_080, MU_150, ZETA_030, ZETA_070, CPAR_X025, KWOB_X4, COMBO_A, COMBO_B, COMBO_D.
- Both-end 20 um support below 75% of L2.3: ZETA_015, ZETA_030, COMBO_C, COMBO_D.

## Required questions

1. **Does reducing G reduce translation?** No clear monotonic reduction appears in this short screen: GRAD_0=+1.95, GRAD_2=+1.88, GRAD_4=+1.94, GEO_230=+1.92, GRAD_8=+1.88. Gradient alone did not remove the translation-dominated response.
2. **Does G=0 retain rotation?** `GRAD_0` phase advance is -51.79 deg with maximum tilt 41.20 deg.
3. **Does cone angle enhance wall-impact wobble?** The response is non-monotonic. The 40 deg case has the largest impact fraction (0.07308 versus 0.03118, 0.01737, and 0.01683 at 20, 30, and 50 deg), making `CONE_40` important for visual impact/separation review; occupancy alone does not establish convincing wobble.
4. **Does B increase rotation or tumble?** Increasing B raises the absolute phase advance from 49.0 deg (`B_6`) through 54.8 deg (`B_10`) to 59.7 deg (`B_14`), without a 90 deg crossing in any of the three short cases. Thus the screened range increases rotation modestly but does not produce short-window tumble.
5. **Which frequency follows more readily?** The 20 Hz case has the largest absolute phase advance (59.6 deg versus 54.8 deg at 30 Hz and 51.5 deg at 40 Hz) over the same 0.20 command cycle. This is the strongest quantitative following response in this frequency family, subject to GIF review.
6. **Does low zeta increase impact separation?** No. Impact fraction rises from 0.00783 at zeta=0.15 to 0.02995 at zeta=0.70; low zeta did not create more detected impact/separation under the coarse event rule. `IMPACT` requires short contact followed by at least 20 us separation.
7. **Does high zeta favor wall sliding?** Yes in this screen: sliding fraction rises from 0 at zeta=0.15 to 0.142, 0.412, and 0.552 at zeta=0.30, 0.50, and 0.70.
8. **Does friction promote rotation or sticking?** Moderate friction increases absolute phase advance (44.0 deg at mu=0, 54.8 deg at 0.03, and 61.6 deg at 0.08), whereas mu=0.15 adds a 0.35 stuck fraction without further phase gain (60.9 deg). Friction therefore promotes rotation up to the screened intermediate level, then introduces sticking at the highest level.
9. **How hydro-sensitive is the mode?** Sensitivity is modest over these multipliers: hydro cases span 52.55-62.65 deg absolute phase advance and +1.895 to +1.952 mm displacement, without a categorical mode change in the automated diagnostics. This remains model sensitivity, not a design recommendation.
10. **Which lengths still tumble?** none in this short screen.
11. **Which lengths tend toward sliding/jam?** Geometry sliding fractions are GEO_200=+0.23, GEO_210=+0.253, GEO_220=+0.253, GEO_230=+0.412, GEO_240=+0.399; stuck-like geometry cases (>5% simulated time) are none in this short screen.
12. **Which combos are rotation-rich?** Combos above the L2.3 rotation/translation index are COMBO_A, COMBO_B, COMBO_D. This is a review set, not a ranking.
13. **Which cases reduce both-end support?** ZETA_015, ZETA_030, COMBO_C, COMBO_D.
14. **Which cases clearly exceed the L2.3 rotation/translation index?** At a descriptive 25% increase (not a success gate): GEO_200, CONE_50, B_6, FREQ_40, MU_080, MU_150, COMBO_A, COMBO_B, COMBO_D.
15. **Which GIFs are worth focused human review?** GEO_200, GEO_210, GEO_220, GEO_240, GRAD_0, GRAD_2, GRAD_8, CONE_20, CONE_40, CONE_50, B_6, B_14, FREQ_40, MU_000, MU_080, MU_150, ZETA_015, ZETA_030, ZETA_070, CPAR_X025, CPAR_X4, KWOB_X025, KWOB_X4, COMBO_A, COMBO_B, COMBO_C, COMBO_D, COMBO_F. They are listed in manifest/family order and are not ranked.

## Human review

Use `gallery/motion_screening_gallery.html`, the nine family mosaics, and `human_review_template.csv`. Judge the target sequence directly: rotate, approach, impact, separate, continue rotating, move contact sector, impact again. Lean-and-slide, tumble, and jam remain boundary cases rather than automatic successes.

AWAITING HUMAN GIF REVIEW
