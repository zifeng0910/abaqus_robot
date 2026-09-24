# F80 / G2.20 coarse TRUE-CEL result

Classification: **F80_RECOIL_FAIL**.

The fresh-start solve halted after the second complete 12.5-ms cycle. C1 advanced +0.166813 mm; C2 reversed -0.102779 mm. MAX_BACKTRACK reached 0.263828 mm. The solver completed its controlled stop, so this is a physical-screen failure rather than numerical invalidity.

F100 / G2.20 had C1/C2 positive and first negative cycle C3. Lowering to F80 moved reversal to C2; it did not eliminate late-cycle reversal. No further frequency tuning is authorized by this task.

The private contact history provides whole robot and pipe General Contact resultants, not a pair-isolated direct robot-wall resultant; these are not interchangeable.

Pre-run F100 rocking check: `ROCKING_REMAINS_PHASE_COHERENT_THROUGH_CYCLE3`. Fundamental rocking amplitudes C1/C2/C3 were 12.937°, 13.856°, 14.477°; C2−C1 and C3−C1 phase lags were +3.03° and +3.87°. RMS phase-matched alpha errors were 2.563° and 2.363°; p95 absolute rocking-rate differences were 259.06 and 225.73 rad/s. At F100 C3 negative-velocity onset, signed TAIL/HEAD wall gaps were +0.00218/+0.04441 mm; at maximum backtrack, −0.000082/+0.04476 mm. The F100 saved history lacks a pair-isolated robot-wall resultant, so force attribution remains unavailable.

Screen conditions: reduced-sound-speed coarse CEL (`c0=100000 mm/s`, 44×20×20, Explicit scale factor 0.4). The F80 run stopped at 25 ms by the complete-cycle failure gate. The planned third and five-cycle continuations were not run.
