# Legacy RouteA Motion Target Report

## Authoritative reconstruction

The paired generator and report define `beta(t) = 10 deg * r(t) * sin(2*pi*666.667*t)` for the 0.03 s display probe, where `r(t)` is a C1 smoothstep from 0 to 1 over the first 2 ms. The 666.667 Hz value produces 20 compressed display cycles. The physical reference is **5 Hz**, so the magnetic validation uses `alpha_B(t) = 10 deg * sin(2*pi*5*t)` over 200 ms.

The generator constructs `Q=[T,N,B]`, with `B=T x N`, and applies `Q Rz(beta) Q0^T`. Therefore `c_hat=T`, `n_rock=N`, and `b_rock=B`. Legacy `N` is initialized from the robot transverse PCA axis, while production `e1` is initialized from projected global Z, so the names cannot be equated. Projection at their common inlet establishes the fixed transported gauge **`chi=-61.372847575963 deg`**: `n_rock=cos(chi)e1+sin(chi)e2 = 0.479107880 e1 - 0.877756025 e2`, and `b_rock=c_hat x n_rock = 0.877756025 e1 + 0.479107880 e2`. Positive rocking gives `a=cos(theta)c_hat+sin(theta)n_rock`. The robot axis `a` is HEAD-to-TAIL; HEAD is the low axial end and TAIL the high axial end.

The legacy initial rocking angle is 0 deg (the ramped sine starts at zero). Its RP follows all three prescribed centerline path translations; it has no separately imposed global-Y/Z wobble, but radial translation is not dynamically free. This and the reported CEL deep-penetration warning limit the GIF to motion topology, plane, phase convention, and approximate amplitude.

## Geometry gate

The actual straight Abaqus wall nodes/faces measure an inscribed solver radius of **0.667345240 mm** (ID **1.334690480 mm**, face-to-face spread 7.392e-11 mm). For `L=2.40 mm`, `D=0.815 mm`, `W(10 deg)=1.219373945 mm`, leaving **0.115316535 mm** diametral analytic margin. Exact faceted-wall checks at -10, -7.5, -5, 0, +5, +7.5, +10 deg give a minimum surface gap of **69.198 um** at centered COM. All requested poses are penetration-free; geometry is unchanged and the Abaqus gate is open.

The centered pose does not create both-end wall support at +/-10 deg. That is acceptable for a clean Level-1 motion-mode test and is not altered by artificial penetration.

## Secondary HighEndStop reference

The centerline-distance audit stays below 0.25 mm through frame 27. Frame 28 is the first suspicious frame (0.3197 mm), and frame 32 is the first obviously invalid frame (0.5153 mm for this robot/lumen scale). `Legacy_HighEndStop_TrustedEarlyWindow.gif` therefore contains frames 0-27 only. This secondary reference does not override RouteA.
