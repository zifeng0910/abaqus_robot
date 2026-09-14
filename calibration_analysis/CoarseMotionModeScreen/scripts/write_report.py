"""Write the no-winner screening report from lightweight metrics."""
from pathlib import Path
import pandas as pd

HERE=Path(__file__).resolve().parent; SCREEN=HERE.parent; REPO=SCREEN.parents[1]

def values(df,ids,column):
    x=df.set_index("case_id"); return ", ".join("%s=%+.3g"%(i,x.loc[i,column]) for i in ids)

def ids_where(df,mask): return ", ".join(df.loc[mask,"case_id"].tolist()) or "none in this short screen"

def main():
    df=pd.read_csv(SCREEN/"metrics/all_cases_metrics.csv"); x=df.set_index("case_id"); b=x.loc["GEO_230"]
    geo=df[df.family=="geometry"]; grad=df[df.family=="gradient"]; cone=df[df.family=="cone"]; mag=df[df.family=="B"]
    freq=df[df.family=="frequency"]; mu=df[df.family=="friction"]; zeta=df[df.family=="collision_damping"]
    hydro=df[df.family.str.startswith("hydro_")]; combo=df[df.family=="combo"]
    aliases=pd.DataFrame([
        ("GRAD_6","GEO_230"),("CONE_30","GEO_230"),("B_10","GEO_230"),("FREQ_30","GEO_230"),
        ("MU_030","GEO_230"),("ZETA_050","GEO_230"),("CPAR_X1","GEO_230"),("KWOB_X1","GEO_230")],columns=["label","case"])
    review_mask=(df.rotation_translation_index>b.rotation_translation_index)|(df.impact_fraction>df.impact_fraction.median())|(df.both_end_support_fraction<b.both_end_support_fraction*.75)
    worth=df.loc[review_mask,"case_id"].tolist()
    text=f"""# Coarse Motion Mode Screen Report

## Scope and decision boundary

This is a cheap, broad, human-in-the-loop screen of 32 unique rigid-robot Abaqus/Explicit cases. Every case uses direct `dt = 1e-7 s`, a global 0.10 mm C3D4 mesh, dense RP/contact history, and 0.20 magnetic cycle. The outputs are motion-mode evidence, not calibrated predictions. No composite score was computed and no winner is selected.

The 0.10 mm L2.3 baseline passed the macro-fidelity gate against the existing 0.060 mm result over 6.667 ms: maximum directed tilt 41.24 vs 41.10 deg, robot phase advance -54.82 vs -60.05 deg, and axial displacement +1.920 vs +1.953 mm. The coarse baseline remains translation-dominated wall sliding, so no 0.08 mm rerun was authorized or needed.

## Screen integrity

- Unique solved cases: {len(df)}; requested baseline labels are deduplicated through `GEO_230`.
- Mesh range: {int(df.element_count.min()):,}-{int(df.element_count.max()):,} C3D4, below the 20,000 hard limit.
- Durations: 20 Hz = 10.000 ms, 30 Hz = 6.667 ms, 40 Hz = 5.000 ms.
- Case-level Socket lifecycles were sequential and ports were confirmed closed before advancing.
- Proximity is a lightweight 20 us sampled, coarse-mesh diagnostic; it is not finalist-grade exact-CAD reconstruction.
- STEP round-trip and coarse tessellation introduce a small geometry/mass approximation; these variants are `SCREENING ONLY`.

## Motion boundaries observed

- 90 deg crossing / tumble-like candidates: {ids_where(df,df.crossing_90deg.astype(bool))}.
- Sustained low-motion contact / stuck-like candidates: {ids_where(df,df.stuck_fraction>0.05)}.
- Translation-heavy relative to L2.3: {ids_where(df,abs(df.delta_s_mm)>=abs(b.delta_s_mm))}.
- Rotation/translation index above L2.3: {ids_where(df,df.rotation_translation_index>b.rotation_translation_index)}.
- Both-end 20 um support below 75% of L2.3: {ids_where(df,df.both_end_support_fraction<b.both_end_support_fraction*.75)}.

## Required questions

1. **Does reducing G reduce translation?** No clear monotonic reduction appears in this short screen: {values(df,["GRAD_0","GRAD_2","GRAD_4","GEO_230","GRAD_8"],"delta_s_mm")}. Gradient alone did not remove the translation-dominated response.
2. **Does G=0 retain rotation?** `GRAD_0` phase advance is {x.loc['GRAD_0','robot_phase_advance_deg']:+.2f} deg with maximum tilt {x.loc['GRAD_0','max_directed_tilt_deg']:.2f} deg.
3. **Does cone angle enhance wall-impact wobble?** The response is non-monotonic. The 40 deg case has the largest impact fraction ({x.loc['CONE_40','impact_fraction']:.4g} versus {x.loc['CONE_20','impact_fraction']:.4g}, {b.impact_fraction:.4g}, and {x.loc['CONE_50','impact_fraction']:.4g} at 20, 30, and 50 deg), making `CONE_40` important for visual impact/separation review; occupancy alone does not establish convincing wobble.
4. **Does B increase rotation or tumble?** Increasing B raises the absolute phase advance from {abs(x.loc['B_6','robot_phase_advance_deg']):.1f} deg (`B_6`) through {abs(b.robot_phase_advance_deg):.1f} deg (`B_10`) to {abs(x.loc['B_14','robot_phase_advance_deg']):.1f} deg (`B_14`), without a 90 deg crossing in any of the three short cases. Thus the screened range increases rotation modestly but does not produce short-window tumble.
5. **Which frequency follows more readily?** The 20 Hz case has the largest absolute phase advance ({abs(x.loc['FREQ_20','robot_phase_advance_deg']):.1f} deg versus {abs(b.robot_phase_advance_deg):.1f} deg at 30 Hz and {abs(x.loc['FREQ_40','robot_phase_advance_deg']):.1f} deg at 40 Hz) over the same 0.20 command cycle. This is the strongest quantitative following response in this frequency family, subject to GIF review.
6. **Does low zeta increase impact separation?** No. Impact fraction rises from {x.loc['ZETA_015','impact_fraction']:.4g} at zeta=0.15 to {x.loc['ZETA_070','impact_fraction']:.4g} at zeta=0.70; low zeta did not create more detected impact/separation under the coarse event rule. `IMPACT` requires short contact followed by at least 20 us separation.
7. **Does high zeta favor wall sliding?** Yes in this screen: sliding fraction rises from {x.loc['ZETA_015','sliding_fraction']:.3g} at zeta=0.15 to {x.loc['ZETA_030','sliding_fraction']:.3g}, {b.sliding_fraction:.3g}, and {x.loc['ZETA_070','sliding_fraction']:.3g} at zeta=0.30, 0.50, and 0.70.
8. **Does friction promote rotation or sticking?** Moderate friction increases absolute phase advance ({abs(x.loc['MU_000','robot_phase_advance_deg']):.1f} deg at mu=0, {abs(b.robot_phase_advance_deg):.1f} deg at 0.03, and {abs(x.loc['MU_080','robot_phase_advance_deg']):.1f} deg at 0.08), whereas mu=0.15 adds a {x.loc['MU_150','stuck_fraction']:.2f} stuck fraction without further phase gain ({abs(x.loc['MU_150','robot_phase_advance_deg']):.1f} deg). Friction therefore promotes rotation up to the screened intermediate level, then introduces sticking at the highest level.
9. **How hydro-sensitive is the mode?** Sensitivity is modest over these multipliers: hydro cases span {abs(hydro.robot_phase_advance_deg).min():.2f}-{abs(hydro.robot_phase_advance_deg).max():.2f} deg absolute phase advance and {hydro.delta_s_mm.min():+.3f} to {hydro.delta_s_mm.max():+.3f} mm displacement, without a categorical mode change in the automated diagnostics. This remains model sensitivity, not a design recommendation.
10. **Which lengths still tumble?** {ids_where(geo,geo.crossing_90deg.astype(bool))}.
11. **Which lengths tend toward sliding/jam?** Geometry sliding fractions are {values(df,geo.case_id.tolist(),"sliding_fraction")}; stuck-like geometry cases (>5% simulated time) are {ids_where(geo,geo.stuck_fraction>0.05)}.
12. **Which combos are rotation-rich?** Combos above the L2.3 rotation/translation index are {ids_where(combo,combo.rotation_translation_index>b.rotation_translation_index)}. This is a review set, not a ranking.
13. **Which cases reduce both-end support?** {ids_where(df,df.both_end_support_fraction<b.both_end_support_fraction*.75)}.
14. **Which cases clearly exceed the L2.3 rotation/translation index?** At a descriptive 25% increase (not a success gate): {ids_where(df,df.rotation_translation_index>1.25*b.rotation_translation_index)}.
15. **Which GIFs are worth focused human review?** {', '.join(worth)}. They are listed in manifest/family order and are not ranked.

## Human review

Use `gallery/motion_screening_gallery.html`, the nine family mosaics, and `human_review_template.csv`. Judge the target sequence directly: rotate, approach, impact, separate, continue rotating, move contact sector, impact again. Lean-and-slide, tumble, and jam remain boundary cases rather than automatic successes.

AWAITING HUMAN GIF REVIEW
"""
    (SCREEN/"Coarse_Motion_Mode_Screen_Report.md").write_text(text,encoding="utf-8")
    reassess="""# L2300 motion-mode reassessment

The previous automatic quantitative gate labelled the L2.3 case `L2300_WALL_SUPPORTED_FORWARD_WOBBLE`. The result and historical report are preserved unchanged. Subsequent user GIF review identifies that motion as a **translation-dominated wall sliding reference**, not an automatically selected success or winner.

Future project motion-mode classification therefore requires both quantitative metrics and human GIF assessment. Wall support, phase advance, and positive axial displacement alone cannot distinguish convincing repeated impact/separation wobble from sustained lean-and-slide motion.
"""
    (SCREEN/"L2300_motion_mode_reassessment.md").write_text(reassess,encoding="utf-8")

if __name__=="__main__":main()
