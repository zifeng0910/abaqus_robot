"""Stage A: current contact-point restitution and normal/tangential work."""
import json

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))
from audit_stage_a import mesh_properties
from exact_gap_audit import Wall
from contact_probe_common import AUDIT, HERE, ROOT, dense_rp, summarize_case


JOB = 'Wobble_F30_G6L45_ReducedHydroFixed_ContactAudit_0012'


def main():
    folder = AUDIT / 'diagnostic'
    t, data = dense_rp(folder)
    meshes, rp, _ = mesh_properties((AUDIT / (JOB + '.inp')).read_text())
    wall = Wall()
    gaps = pd.read_csv(AUDIT / 'first_impact_gap_history.csv')
    # The published gap file covers 0.98--1.08 ms; retain that complete window.
    i0 = int(round(gaps.time_s.iloc[0] / 1e-7))
    i1 = i0 + len(gaps)
    summary, history, events, windows = summarize_case(
        'baseline_zeta0p055_default_tangent1', t[i0:i1],
        {k: v[i0:i1] for k, v in data.items()}, meshes['Robot_SOLID'], rp,
        wall, gaps, folder, .055, 1.0)
    pd.DataFrame([summary]).to_csv(HERE / 'current_contact_restitution_audit.csv', index=False)
    history.to_csv(HERE / 'current_contact_normal_tangential_work.csv', index=False)
    events.to_csv(HERE / 'baseline_contact_event_catalog.csv', index=False)
    (HERE / 'current_contact_restitution_summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
