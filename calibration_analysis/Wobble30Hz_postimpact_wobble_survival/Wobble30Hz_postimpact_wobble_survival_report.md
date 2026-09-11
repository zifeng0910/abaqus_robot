# 30 Hz post-impact wobble survival audit

## Decision
`D_POST_IMPACT_NONMAGNETIC_ROTATIONAL_LOADING_DOMINATES`. This turn used completed telemetry/ODB exports only; no new Abaqus job was submitted.

## Runtime and replay identity
The paired Wall-ON reference records G6/L45 exactly: 6 mT = 0.006 T, L=45 mm, B0=10 mT, f=30 Hz, cone30°, Bias40°, phase248°, sense+1, χ=0. Server SHA256: `d6e8cd46521a6921e58d2cd7628cef42ae999dedf414ebb6a71165b307de2b43`. The 55.5 ms Wall-OFF COM-fixed record logs **G3** (3 mT), so it is not mislabeled G6. The archived production evaluator has now been replayed exactly; B/F/T regression passes at machine precision. The earlier mismatch came from an approximate flattened-CSV/frame replay, not the production server.

## Wall-OFF true-axis result
`walloff_f30_true_axis_phase.csv` uses the directed HEAD→TAIL axis and separates body spin (ω·a) from wobble (ω−(ω·a)a). `walloff_f30_lock_windows.csv` reports field rate, robot phase rate, R_lock, wobble rate and true-axis range over startup, cycle 1 and cycle 2 windows.

## Wall-ON post-impact result
The 3 ms paired run is split into pre-contact, contact, immediate post-impact and latest windows. Exact geometric gap and aggregate CPRESS are kept separate in `wallon_contact_timing.csv`; rates and recorded applied torque are in `wallon_pre_postimpact_wobble.csv` and `postimpact_torque_budget.csv`. No GIF was regenerated because no dynamic case was authorized in this diagnostic turn.

## B estimate
Offline-only candidates 12.5/15/20 mT are listed in `required_B_candidate_estimate.csv`; no B increase is selected because the post-impact torque remains strong (1.301× pre-contact). These are not dynamic predictions.

## Final interpretation
Classification `D_POST_IMPACT_NONMAGNETIC_ROTATIONAL_LOADING_DOMINATES` is based on Wall-OFF lock and Wall-ON post-impact phase/torque ratios, not on a single UR1 component. Wall-OFF cycle 1/2 are synchronized, while the Wall-ON latest window loses phase rate but retains torque. Because that latest window contains only five telemetry rows, the next evidence should be a longer paired Wall-ON record; no B-amplitude increase is justified by this audit alone.

## Runtime regression correction
The exact archived evaluator was re-run in the configured SIMULIA Python environment (`runtime_exact_evaluator_replay.csv`). B, F and T maximum absolute errors are below 2e-12, so the earlier 0.31 mT / 0.0116 N·mm discrepancy was a replay-frame/curve approximation error, not a production-server mismatch. The numerical field follows `driver_arc` in the archived telemetry even though the banner says `ROBOT_ARC`; this metadata discrepancy is retained for future audits.

## Quantitative classification caveat
Wall-OFF cycle-1/2 directed phase rates are approximately 32.16 and 29.58 Hz (absolute R_lock 1.072, 0.986), supporting B10 synchronization without the wall. Wall-ON latest-window R_phase is 0.134, while recorded |Tmag| remains 1.301× pre-contact; this supports `D_POST_IMPACT_NONMAGNETIC_ROTATIONAL_LOADING_DOMINATES` rather than magnetic torque collapse. The latest window contains only five output rows, so a longer wall-on wobble-survival probe is required before claiming a sustained stall.

Angular impulse bookkeeping is in `angular_momentum_audit.csv`; ΔH is intentionally marked unavailable because the archived exports contain no independent rigid-body inertia/history channel.
