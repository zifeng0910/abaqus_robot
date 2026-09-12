"""Build the single approved 1.3-ms damping probe; never launches Abaqus."""
from pathlib import Path
import hashlib
import json


HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent
ROOT = AUDIT.parents[2]
SOURCE = 'Wobble_F30_G6L45_ReducedHydroFixed_ContactAudit_0012'
JOB = 'Wobble_F30_G6L45_ReducedHydroFixed_NormalDamp020_ContactAudit_0013'
OLD = '*Contact Damping, definition=CRITICAL DAMPING FRACTION\n0.055'
NEW = '*Contact Damping, definition=CRITICAL DAMPING FRACTION, tangent fraction=0.0\n0.20'


def main():
    source = (AUDIT / (SOURCE + '.inp')).read_text()
    assert source.count(OLD) == 1
    candidate = source.replace(OLD, NEW)
    restored = candidate.replace(NEW, OLD)
    assert restored == source
    # Frozen parameters must remain explicitly present.
    for token in ('0.03,', '1.0e-7, 0.001300', 'phase=248.0 deg',
                  'gradient=6.000 mT=0.006000 T'):
        assert token in candidate, token
    bridge = (AUDIT.parent / 'Wobble30Hz_reduced_hydrodynamics' / 'vuforc_socket_bridge_reduced_hydro.f').read_text()
    for token in ('cpar /4.0D-9/', 'cperp /1.2D-8/', 'kspin /1.0D-9/', 'kwob /3.0D-9/'):
        assert token in bridge, token
    local = HERE / (JOB + '.inp')
    root = ROOT / (JOB + '.inp')
    local.write_text(candidate)
    root.write_text(candidate)

    runner = (AUDIT / 'run_contact_audit_0012.ps1').read_text()
    runner = runner.replace(SOURCE, JOB).replace('[int]$Port=65505', '[int]$Port=65506')
    runner = runner.replace('PORT=$Port;', 'PORT=$Port; NORMAL_DAMPING=0.20; TANGENT_FRACTION=0.0;')
    (HERE / 'run_normaldamp020_contact_audit_0013.ps1').write_text(runner)
    datacheck = """$ErrorActionPreference='Stop'
Set-Location -LiteralPath 'J:\\abaqusfangzhen'
. '.\\launch_abaqus_with_env.ps1' -SkipVerify
$job='Wobble_F30_G6L45_ReducedHydroFixed_NormalDamp020_ContactAudit_0013'
if(Test-Path "$job.odb"){throw 'Existing ODB; do not repeat the candidate'}
& 'I:\\SIMULIA\\Commands\\abq2025.bat' job=$job input="$job.inp" user='J:\\abaqusfangzhen\\abaqus_robot\\calibration_analysis\\Wobble30Hz_reduced_hydrodynamics\\vuforc_socket_bridge_reduced_hydro.f' datacheck double=both cpus=1 ask_delete=OFF interactive
if($LASTEXITCODE -ne 0){throw "Datacheck failed: $LASTEXITCODE"}
if(Select-String -LiteralPath "$job.dat" -Pattern '\\*\\*\\*ERROR') {throw 'Datacheck ERROR found'}
Write-Host 'DATACHECK PASS: zeta=0.20, tangent fraction=0.0 accepted by Abaqus 2025'
"""
    (HERE / 'datacheck_normaldamp020.ps1').write_text(datacheck)
    identity = {
        'label': 'PROVISIONAL_CONTACT_DISSIPATION_PROBE',
        'source_job': SOURCE, 'candidate_job': JOB,
        'source_sha256': hashlib.sha256((AUDIT / (SOURCE + '.inp')).read_bytes()).hexdigest(),
        'candidate_sha256': hashlib.sha256(local.read_bytes()).hexdigest(),
        'only_physics_change': {
            'critical_damping_fraction': [0.055, 0.20],
            'tangent_fraction': ['Abaqus Explicit default 1.0', 0.0],
        },
        'frozen': ['HARD pressure-overclosure', 'General Contact', 'SmoothWall114',
                   'friction 0.03', 'surfaces and penalty formulation', 'dt 1e-7 s',
                   'magnetic architecture', 'robot geometry/mass/inertia/pose', 'Reduced-Hydro'],
    }
    (HERE / 'normaldamp020_input_identity.json').write_text(json.dumps(identity, indent=2))
    print(json.dumps(identity, indent=2))


if __name__ == '__main__':
    main()
