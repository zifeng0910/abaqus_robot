"""Build ONE instrumentation-only case after the explicit Stage A gate."""
from pathlib import Path
import json,hashlib
import pandas as pd
from audit_stage_a import HERE,ROOT,JOB

NEW='Wobble_F30_G6L45_ReducedHydroFixed_ContactAudit_0012'
s=json.loads((HERE/'stage_a_summary.json').read_text())
g=pd.read_csv(HERE/'verified_sparse_wall_gap.csv')
windows=pd.read_csv(HERE/'preimpact_momentum_windows.csv')
hist=pd.read_csv(HERE/'baseline'/'existing_odb_history_inventory.csv')
gate=dict(clean_momentum_closes=s['clean_residual']<.1,nearwall_at_1ms=abs(g.iloc[(g.time_s-.001).abs().argmin()].gap_um)<5,large_missing_impulse=windows.iloc[1].missing_norm_Ns>100*windows.iloc[0].missing_norm_Ns,no_contact_history=not hist.variable.str.match(r'CF[NST]|CNORMF|CSHEARF').any())
gate={k:bool(v) for k,v in gate.items()}
assert all(gate.values()),gate
original=(ROOT/(JOB+'.inp')).read_text()
split=original.index('*Output, field')
before=original[:split].replace('1.0e-7, 0.008333','1.0e-7, 0.001300')
assert before.replace('1.0e-7, 0.001300','1.0e-7, 0.008333')==original[:split]
newoutput='''*Output, field, time interval=5.0e-7, time marks=NO
*Node Output
U
*Node Output, nset=RP_ROBOT
U, UR, V, VR, A, AR, CF, RF, RM
*Contact Output, general contact
CSTRESS, CFORCE
*Contact Output, surface=ROBOT_SOLID-1.ROBOT_SOLID_SURF
CDISP
*Output, history, frequency=1
*Node Output, nset=RP_ROBOT
COORD, U, UR, V, VR, A, AR, CF, RF, RM
*Contact Output, surface=ROBOT_SOLID-1.ROBOT_SOLID_SURF
CFN, CFS, CFT
*Contact Output, surface=Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF
CFN, CFS, CFT
*Output, history, time interval=5.0e-5
*Energy Output, variable=ALL
'''
sensors=original[original.index('*Output, history, frequency=1, sensor'):]
dest=ROOT/(NEW+'.inp')
if (ROOT/(NEW+'.odb')).exists():
    assert dest.read_text()==before+newoutput+sensors,'Existing ODB/input must not be overwritten'
else:
    dest.write_text(before+newoutput+sensors)
# Exactly the existing production runner and bridge; only job/output/log cadence.
runner=(HERE.parent/'Wobble30Hz_reduced_hydrodynamics'/'run_reduced_hydro_0083.ps1').read_text()
runner=runner.replace(JOB,NEW).replace('[int]$Port=65504','[int]$Port=65505')
runner=runner.replace("'--telemetry-interval-s','0.00005'","'--telemetry-interval-s','0.0'")
runner=runner.replace('duration=8.333ms','duration=1.3ms; CONTACT_HISTORY=every_increment; FIELD=0.5us')
runner=runner.replace('Write-Host "Completed $Job"',"Copy-Item -LiteralPath 'J:\\abaqusfangzhen\\reduced_hydro_runtime_load.csv' -Destination (Join-Path $WorkDir ($Job+'_hydro_increment.csv'))\nWrite-Host \"Completed $Job\"")
(HERE/'run_contact_audit_0012.ps1').write_text(runner)
# All physics up to output block byte-identical except end time.
identity=dict(gates=gate,source_job=JOB,diagnostic_job=NEW,physics_unchanged=True,duration_s=.0013,increment_s=1e-7,field_interval_s=5e-7,history_interval_s=1e-7,source_sha256=hashlib.sha256(original.encode()).hexdigest(),new_sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),bridge_sha256=hashlib.sha256((HERE.parent/'Wobble30Hz_reduced_hydrodynamics'/'vuforc_socket_bridge_reduced_hydro.f').read_bytes()).hexdigest())
(HERE/'instrumentation_identity.json').write_text(json.dumps(identity,indent=2));print(identity)
