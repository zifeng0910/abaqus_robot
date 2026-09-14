from pathlib import Path
import json

import pandas as pd

HERE=Path(__file__).resolve().parent


def main():
    s=pd.read_csv(HERE/"cross_replay_bridge_summary.csv"); identity=json.loads((HERE/"cross_replay_identity.json").read_text())
    case=lambda g,t:s[(s.geometry==g)&(s.trajectory==t)].iloc[0]
    a,b,c,d=case("OLD","OLD"),case("NEW","OLD"),case("OLD","NEW"),case("NEW","NEW")
    assert abs(a.longest20_ms-.2588)<.002
    assert abs(d.longest20_ms-.547)<.002
    assert b.longest20_ms>.5 and c.longest20_ms<.5
    assert all(getattr(b,f"longest{x}_ms")>.5 for x in (20,10,5,0))
    assert d.longest10_ms<.5 and d.longest5_ms==0 and d.longest0_ms==0
    assert abs(a.penetration_um+1.3575)<.05
    assert abs(d.penetration_um+.7916)<.2
    assert c.penetration_um>0
    assert identity["wall"]=="SmoothWall114 exact triangles"
    assert identity["pose_map"]=="COM(t)+R(t)*(x_ref-COM_ref)"
    assert identity["causal_classification"]=="EXACT_HEAD_GEOMETRY_DOMINATES_BRIDGE_REINTRODUCTION"
    head=pd.read_csv(HERE/"head_profile_difference_summary.csv").iloc[0]
    assert head.head_max_outward_difference_um>70 and head.head_max_difference_x_mm<-.8
    poses=pd.read_csv(HERE/"bridge_pose_geometry.csv")
    assert len(poses)==3 and set(poses.relative_to_first_flip)=={"after"}
    print("PASS: four-way exact-wall bridge causality gates")


if __name__=="__main__":main()
