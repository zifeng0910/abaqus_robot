# Precomputed Magnetic + FAST_SURROGATE Forward Short Screen

Selected result: `FAST_PRECOMP_F120_G1P20_FLOW_1CYCLE`

This is a low-cost non-CEL screen. It uses the validated offline Magpylib-derived field/gradient table during dynamics and the existing relative-velocity `FAST_SURROGATE_FLUID` model. No live Python, TCP, socket, CEL, CFD, or FSI coupling was used.

Frozen settings: rigid refined robot, straight tube, current General Contact law, `mu=0.03`, `zeta=1.0`, `B0=10 mT`, `f=120 Hz`, `A_main=14.8 deg`, `A_cross=2.5 deg`, initial pose, `U_flow=+10 mm/s`, `Cparallel=4e-9 N s/mm`, `Cperp=1.2e-8 N s/mm`, `Kspin=1e-9 N mm s`, and `Kwobble=3e-9 N mm s`.

| Case | G (mT) | delta_s (mm) | fraction v_s > 0 | mean v_s (mm/s) | final v_s (mm/s) | rocking min/max (deg) | tumble |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Parent first cycle | 0.15 | -0.03240515 | 0.102309 | - | - | - | no |
| Initial short trial | 0.30 | -0.02641529 | 0.226577 | -3.16982 | -5.01358 | -12.3524 / +13.3278 | no |
| Selected correction | 1.20 | +0.00957765 | 0.792311 | +1.14934 | +4.07839 | -12.3614 / +13.3278 | no |

The `0.30 mT` trial improved directionality but did not reverse the net displacement. Linear extrapolation of the measured `0.15` and `0.30 mT` first-cycle results placed the zero crossing near `0.96 mT`; the evidence-based `1.20 mT` correction then produced net left-to-right progression while preserving the prior rocking range.

Classification: `PASSED_FORWARD_SHORT_SCREEN`

This only establishes a one-cycle qualitative forward screen. The current contact still has three TAIL episodes and should not be interpreted as final no-rebound contact validation.
