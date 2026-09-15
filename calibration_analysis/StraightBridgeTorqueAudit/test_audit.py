import json
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
zero = json.loads((HERE / "zero_solve_summary.json").read_text())
replay = json.loads((HERE / "prescribed_replay_summary.json").read_text())
balance = pd.read_csv(HERE / "straight_bridge_torque_balance.csv")
reaction = pd.read_csv(HERE / "prescribed_reaction_vs_magnetic.csv")

assert zero["direct_contact_torque_available"] is False
assert 0.99 <= zero["bridge_fraction"] <= 1.0
assert zero["zeta020_authorized"] is False
assert replay["primary_classification"] == "PRESCRIBED_MOTION_MASKED_GEOMETRIC_BRIDGE"
assert replay["both_end_support_fraction"] == 1.0
assert replay["contact_active_fraction"] == 1.0
assert replay["minimum_gap_um"] < -200
assert replay["case_A_supported"] is False
assert replay["zeta020_magnetic_run_authorized"] is False
assert len(balance) == 3335
assert len(reaction) == 1668
assert balance.time_s.is_monotonic_increasing
assert reaction.time_s.is_monotonic_increasing
print("PASS: straight bridge torque audit numerical gates")
