"""Independent U/UR/V/VR kinematic gate for a prescribed V+VR CEL replay."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation

ROOT=Path(__file__).resolve().parents[1]
TARGET="TRUECEL_A14P5_DIAG_2CYCLES"
LO,HI=8.333333e-3,10.682021e-3

def interp(z,stem,t):
    return np.column_stack([np.interp(t,z[f"{stem}1"][:,0],z[f"{stem}{i}"][:,1]) for i in (1,2,3)])
def stat(x,t):
    a=np.abs(np.asarray(x));i=int(np.argmax(a))
    return {"rms":float(np.sqrt(np.mean(np.asarray(x)**2))),"p99":float(np.percentile(a,99)),"max":float(a[i]),"max_time_ms":float(t[i]*1e3)}

def main():
    if len(sys.argv)!=2: raise SystemExit("usage: analyze_a14p5_vvr_case.py JOB")
    job=sys.argv[1]; case=ROOT/"case"/job; private=case/"private"
    actual=np.load(private/"rp_history_private.npz"); target=np.load(ROOT/"case"/TARGET/"private"/"rp_history_private.npz")
    ident=json.loads((ROOT/"case"/TARGET/"case_identity.json").read_text())
    t=actual["U1"][:,0]; A={s:interp(actual,s,t) for s in ("U","UR","V","VR")}; D={s:interp(target,s,t) for s in ("U","UR","V","VR")}
    pos=np.linalg.norm(A["U"]-D["U"],axis=1)
    ori=np.rad2deg((Rotation.from_rotvec(A["UR"]).inv()*Rotation.from_rotvec(D["UR"])).magnitude())
    s=np.asarray(ident["canonical_plus_s_axis_aba"]); b=np.asarray(ident["b_routeA_aba"])
    vs=(A["V"]-D["V"])@s; wr=(A["VR"]-D["VR"])@b; crit=(t>=LO)&(t<=HI)
    metrics={"translation_error_mm":{"full":stat(pos,t),"critical":stat(pos[crit],t[crit])},"orientation_geodesic_deg":{"full":stat(ori,t),"critical":stat(ori[crit],t[crit])},"v_s_error_mm_s":{"full":stat(vs,t),"critical":stat(vs[crit],t[crit])},"omega_rock_error_rad_s":{"full":stat(wr,t),"critical":stat(wr[crit],t[crit])}}
    passed=pos.max()<1e-3 and ori.max()<.02 and metrics["v_s_error_mm_s"]["critical"]["p99"]<.05 and metrics["omega_rock_error_rad_s"]["critical"]["p99"]<.5
    out={"job":job,"method":"rotation-matrix/geodesic comparison of target vs actual global V + spatial/global VR replay","passed":bool(passed),"thresholds":{"max_translation_error_mm":.001,"max_orientation_error_deg":.02,"critical_v_s_p99_mm_s":.05,"critical_omega_rock_p99_rad_s":.5},"metrics":metrics,"initial_pose":{"target_U":D["U"][0].tolist(),"actual_U":A["U"][0].tolist(),"target_UR":D["UR"][0].tolist(),"actual_UR":A["UR"][0].tolist()},"window_ms":[8.333333,10.682021],"classification":"COARSE_VVR_REPLAY_VALID" if passed else "COARSE_VVR_REPLAY_INVALID"}
    (ROOT/(job+"_Kinematic_Gate.json")).write_text(json.dumps(out,indent=2)+"\n",encoding="ascii")
    print(json.dumps(out,indent=2)); raise SystemExit(0 if passed else 2)
if __name__=="__main__":main()
