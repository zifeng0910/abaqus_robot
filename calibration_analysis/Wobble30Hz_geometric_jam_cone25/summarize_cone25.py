from pathlib import Path
import pandas as pd, numpy as np
out=Path(__file__).resolve().parent; root=out.parents[2]
def metrics(job,adir,cone):
 t=pd.read_csv(root/(job+'_telemetry.csv')); rp=pd.DataFrame(__import__('json').loads((root/(job+'_rp_fields.json')).read_text())['rows']); ur=np.vstack(rp.UR.to_numpy()); gap=pd.read_csv(root/(job+'_exact_wall_penetration.csv')); cp=pd.read_csv(root/(job+'_cpress_timeseries.csv')); av=pd.read_csv(adir/'wallon_8p333_angular_velocity.csv');
 return dict(job=job,cone_deg=cone,duration_ms=t.t_s.max()*1000,driver_delta_mm=t.driver_arc_mm.iloc[-1]-t.driver_arc_mm.iloc[0],robot_delta_mm=t.robot_arc_mm.iloc[-1]-t.robot_arc_mm.iloc[0],mean_Vt_mm_s=t.force_tangent_N.mean()*1000,max_Vt_mm_s=t.force_tangent_N.max()*1000,UR1_range_rad=np.ptp(ur[:,0]),UR2_range_rad=np.ptp(ur[:,1]),UR3_range_rad=np.ptp(ur[:,2]),wobble_rms_rad_s=av.omega_wobble_rad_s.mean(),min_exact_gap_um=gap.min_signed_gap_mm.min()*1000,max_CPRESS_MPa=cp.CPRESS_max_MPa.max())
cone25=metrics('Wobble_F30_G6L45_WallOn_Free_0083_Cone25_WobbleSurvival',out,25)
base=metrics('Wobble_F30_G6L45_WallOn_Free_0083_WobbleSurvival',out.parent/'Wobble30Hz_wallon_8p333',30)
pd.DataFrame([base,cone25]).to_csv(out/'cone25_vs_cone30_summary.csv',index=False)
print(pd.DataFrame([base,cone25]).to_string(index=False))
