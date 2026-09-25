# F100_G2P20_NOFLUID_REALWALL50

NOFLUID_REALWALL_CHATTER_FAIL

WALL_CONTACT_PERTURBS_BUT_DOES_NOT_CAUSE_SUSTAINED_RECOIL

Completed 50 ms: True. Net displacement +2.573047 mm; max backtrack 0.000074 mm; forward-time fraction 0.995.

| Cycle | delta_s mm | mean v_s mm/s | end v_s | minimum v_s | max backtrack mm | rocking amp deg | fundamental deg | lag deg | peak omega rad/s |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.07920137341790508 | 7.920137341790507 | 19.2053078050083 | -0.6418596751978327 | 7.396943792856853e-05 | 12.770736335127763 | 12.823222996887266 | 17.266974504169696 | 301.09575763068943 |
| 2 | 0.31771205015297005 | 31.771205015297003 | 43.590712243292835 | 19.2053078050083 | 0.0 | 11.52780603242714 | 13.255702432401296 | 13.378917473768695 | 324.835355389428 |
| 3 | 0.5386434030139107 | 53.86434030139107 | 64.71186010681046 | 43.59071224329283 | 0.0 | 12.774583267749616 | 13.110095438666866 | 8.278442539559393 | 294.39293028783743 |
| 4 | 0.7319108310339485 | 73.19108310339485 | 82.74034565064902 | 64.71186010681046 | 0.0 | 13.268271123967722 | 13.560561585979555 | 5.981734308790749 | 247.2494152309099 |
| 5 | 0.9055796209743481 | 90.55796209743481 | 99.95310921929673 | 82.67743093228975 | 0.0 | 13.052855636788408 | 13.419059617432481 | 6.572512168617038 | 228.85956276094007 |

## Three-way comparison

| Cycle | no fluid, no wall delta_s mm | no fluid, real wall delta_s mm | clean FULL CEL delta_s mm |
|---:|---:|---:|---:|
| 1 | +0.127606 | +0.079201 | +0.092337 |
| 2 | +0.409383 | +0.317712 | +0.139908 |
| 3 | +0.694075 | +0.538643 | +0.013300 |
| 4 | +0.981784 | +0.731911 | -0.564697 |
| 5 | +1.271742 | +0.905580 | None |

No-wall C1-C5 rocking half-range is 21.18-25.67 deg; real-wall C1-C5 is 11.53-13.27 deg. Real-wall lag changes from 17.27 to 6.57 deg without a sustained positional reversal. FULL CEL C4 is -0.565 mm with 0.565 mm backtrack; C5 was not completed. No-wall has no active wall contact. The FULL CEL parent lacks 25-us end-resolved contact fields, so its same-end episode count cannot be compared at this resolution.

## Contact

Same-end recontacts under broad COPEN/gap gate: {'HEAD': 0, 'TAIL': 27}. TAIL episodes with negative gap and peak wall force >0.001 N: 21.

The repeated-contact stop condition was established only in offline end-resolved fields; the solver therefore reached the full 50 ms.

Contact requires robot-surface COPEN, signed cylindrical gap and direct wall-only CFNM; 25-us fields, 0.15-ms episode merge. End-specific force is attributed only where that end contacts alone.

| End | start ms | end ms | duration ms | minimum gap mm | peak normal N | tangent impulse N s | attribution | same-end recontact |
|---|---:|---:|---:|---:|---:|---:|---|---|
| HEAD | 8.124999701976776 | 8.124999701976776 | 0.0 | 0.005408611961508325 | 0.008898177184164524 | None | exclusive end | False |
| TAIL | 10.924999602138996 | 10.924999602138996 | 0.0 | 0.0003160676759965675 | 0.0687878429889679 | None | exclusive end | False |
| TAIL | 13.074999675154686 | 13.650000095367432 | 0.5750004202127457 | -0.000938816613796245 | 0.00037027665530331433 | 3.784010728034439e-09 | exclusive end | True |
| TAIL | 15.87500050663948 | 15.87500050663948 | 0.0 | -0.0003035880543333924 | 0.4321447014808655 | None | exclusive end | True |
| TAIL | 17.224999144673347 | 17.224999144673347 | 0.0 | -0.0007990915531211895 | 0.026527272537350655 | None | exclusive end | True |
| TAIL | 17.625000327825546 | 18.60000006854534 | 0.9749997407197952 | -0.0013318843026728366 | 0.00543658621609211 | 1.0261591658274343e-08 | exclusive end | True |
| HEAD | 21.22499980032444 | 21.22499980032444 | 0.0 | 0.009760478288057484 | None | None | overlap; exclusive samples only | False |
| TAIL | 21.22499980032444 | 21.22499980032444 | 0.0 | -0.0013277460755274362 | None | None | overlap; exclusive samples only | False |
| TAIL | 22.199999541044235 | 22.199999541044235 | 0.0 | -0.000692613586876778 | 0.03490365669131279 | None | exclusive end | True |
| TAIL | 22.749999538064003 | 23.05000089108944 | 0.3000013530254364 | -0.003667957845550185 | 0.3049054443836212 | 7.225260883568627e-08 | exclusive end | True |
| TAIL | 23.900000378489494 | 23.97499978542328 | 0.07499940693378448 | -0.0015323245541585262 | 0.0005422941758297384 | 6.106933706880637e-10 | exclusive end | True |
| TAIL | 24.399999529123306 | 24.399999529123306 | 0.0 | -0.004076257035966369 | 0.003321112133562565 | None | exclusive end | True |
| TAIL | 24.824999272823334 | 24.824999272823334 | 0.0 | -0.0021976297363393016 | 0.01581496372818947 | None | exclusive end | True |
| TAIL | 25.499999523162842 | 25.499999523162842 | 0.0 | -0.0009166131069728856 | 0.107785664498806 | None | exclusive end | True |
| TAIL | 25.824999436736107 | 25.824999436736107 | 0.0 | -0.0025169126016036936 | 0.03979974985122681 | None | exclusive end | True |
| TAIL | 26.17499977350235 | 26.225000619888306 | 0.05000084638595581 | -0.004657615606078247 | 0.000915126409381628 | 7.288905554056751e-10 | exclusive end | True |
| TAIL | 26.42500028014183 | 27.574999257922173 | 1.149998977780342 | -0.002438956217589494 | 0.034463945776224136 | 2.4351120894512813e-08 | exclusive end | True |
| TAIL | 30.424999073147774 | 30.424999073147774 | 0.0 | -0.0020343498020277684 | 0.012234887108206749 | None | exclusive end | True |
| TAIL | 31.575001776218414 | 31.599998474121094 | 0.024996697902679443 | -0.002565128609070655 | 0.0003151801647618413 | 1.8462864272875774e-10 | exclusive end | True |
| TAIL | 31.800001859664917 | 31.925000250339508 | 0.12499839067459106 | -0.00197306922444751 | 0.0012478722492232919 | 2.7043567456855775e-09 | exclusive end | True |
| TAIL | 32.12499991059303 | 32.15000033378601 | 0.025000423192977905 | -0.003332644239618432 | 0.00015725780394859612 | 1.1479720487397936e-10 | exclusive end | True |
| HEAD | 32.40000084042549 | 32.40000084042549 | 0.0 | 0.011140721307632884 | None | None | overlap; exclusive samples only | False |
| TAIL | 32.40000084042549 | 32.40000084042549 | 0.0 | 0.0005174544061379294 | None | None | overlap; exclusive samples only | False |
| TAIL | 34.974999725818634 | 34.974999725818634 | 0.0 | -0.000229058626595946 | 0.008521017618477345 | None | exclusive end | True |
| TAIL | 35.349998623132706 | 35.54999828338623 | 0.19999966025352478 | -0.0035720081056193065 | 0.023371465504169464 | 1.4890203599884607e-08 | exclusive end | True |
| TAIL | 35.82499921321869 | 35.82499921321869 | 0.0 | -0.004551226543465292 | 0.0243939571082592 | None | exclusive end | True |
| TAIL | 36.04999929666519 | 36.32500022649765 | 0.2750009298324585 | -0.0010778376772451193 | 0.0007140971720218658 | 3.656173220993234e-09 | exclusive end | True |
| TAIL | 36.4999994635582 | 37.324998527765274 | 0.824999064207077 | -0.0011814684247999363 | 0.0006499161827377975 | 1.3430741035172378e-08 | overlap; exclusive samples only | True |
| HEAD | 37.174999713897705 | 37.324998527765274 | 0.14999881386756897 | 0.009838955644336123 | 0.0 | 0.0 | overlap; exclusive samples only | False |
| TAIL | 39.07499834895134 | 39.07499834895134 | 0.0 | -0.00024082617992160493 | 0.03044973313808441 | None | exclusive end | False |
| TAIL | 40.699999779462814 | 40.800001472234726 | 0.10000169277191162 | -0.005614879284443264 | 0.04790526255965233 | 3.153501577511826e-08 | exclusive end | True |
| TAIL | 41.12499952316284 | 42.27500036358833 | 1.1500008404254913 | -0.002836868376150159 | 0.02287614718079567 | 1.7998131794996086e-08 | exclusive end | True |
| HEAD | 42.47500002384186 | 42.55000129342079 | 0.07500126957893372 | 0.002789377358058176 | None | None | overlap; exclusive samples only | False |
| TAIL | 42.47500002384186 | 42.55000129342079 | 0.07500126957893372 | -0.002000244845761334 | None | None | overlap; exclusive samples only | False |
| TAIL | 45.57500034570694 | 45.57500034570694 | 0.0 | -0.0025774712833291913 | 0.0011132939253002405 | None | exclusive end | True |
| TAIL | 45.875001698732376 | 45.875001698732376 | 0.0 | -0.0009389492541177491 | 0.02385125495493412 | None | exclusive end | True |
| TAIL | 46.12499848008156 | 47.35000059008598 | 1.225002110004425 | -0.0018841095188348866 | 0.00247932318598032 | 2.151756685983766e-08 | overlap; exclusive samples only | True |
| HEAD | 47.325000166893005 | 47.35000059008598 | 0.025000423192977905 | 0.010733374911158 | None | None | overlap; exclusive samples only | False |

## Magnetic and comparison

Magnetic 10-us trace and 0.2-ms collision windows are in the CSV/JSON. |B| remains 10.99997-11.00000 mT; +s magnetic force remains positive. VUAMP uses 100 Hz; table JSON's 120 Hz is stale metadata.

Maximum absolute change in theta-command error across a 0.4-ms contact-centered window: 2.350 deg. This and the cycle phase lags show no catastrophic loss of magnetic phase lock, though the lag drifts over five cycles.

The 21 stronger TAIL episodes establish contact chatter, but wall contact alone does not reproduce FULL CEL's C4 sustained recoil. Fluid/CEL coupling remains required for that observed recoil in this matched control; the exact force pathway is not resolved here.

The clean FULL CEL baseline ended partially at 40.755 ms; no C5 comparison is inferred.

Three-way cycle metrics are in the JSON. This one-case isolation changes only removal of fluid physics and retains parent robot-wall law.
