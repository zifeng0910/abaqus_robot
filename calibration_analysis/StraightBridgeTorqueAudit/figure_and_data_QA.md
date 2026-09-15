# Figure and data QA

- Core conclusion: the nominal 30 deg prescribed trajectory remains oppositely supported and activates penalty-scale reaction; it does not demonstrate clean wobble.
- Archetype: quantitative grid with synchronized torque, angular velocity, gaps, and bridge state; fixed-camera dual-view animations provide geometric validation.
- Backend: Python/matplotlib exclusively.
- Source-data sampling: dense ODB histories were used for all metrics; public plotting CSVs retain every 100th increment (10 us) plus the final increment. Event impulses use the full 0.1 us histories.
- Figure exports: PNG and PDF; all rendered text is at least 5 pt by source definition.
- GIF integrity: fixed limits and camera, transparent wall, complete robot, HEAD/TAIL, local cross-section, B direction, torque direction, gaps, gap-rate normal-speed proxy, and bridge state. No frame-wise auto-fit.
- Contact limitation: baseline contact torque is balance-inferred. ODB has no node-level contact force/contact-point coordinate pair, so normal and tangential moment contributions cannot be separated.
- Replay limitation: deep penetration makes reaction magnitude penalty-dominated. Ratio values are diagnostic and not actuator-sizing data.
