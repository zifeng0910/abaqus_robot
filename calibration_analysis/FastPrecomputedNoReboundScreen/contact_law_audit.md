# General Contact audit

Old: General Contact penalty, HARD pressure-overclosure, mu=0.03, normal zeta=1.0.

New (single modification): `SCALE FACTOR` with `r=5%`, geometric stiffness multiplier `10`, and initial stiffness scale `0.01`; damping and friction are unchanged. For the measured 0.0449331796 mm minimum robot tetra edge, the first transition is about 0.00224666 mm. This targets a weak initial impulse while progressively recovering nonpenetration. It does not add adhesion, no-separation, a velocity clamp, Contact Pair, or SplitWall.
