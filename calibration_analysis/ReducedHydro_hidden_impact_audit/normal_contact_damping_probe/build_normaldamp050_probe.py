"""Build the one approved 0.50 normal-only probe; never launches Abaqus."""
from pathlib import Path
import hashlib
import json


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
WORKDIR = REPO.parent
SOURCE = 'Wobble_F30_G6L45_ReducedHydroFixed_NormalDamp020_ContactAudit_0013'
JOB = 'Wobble_F30_G6L45_ReducedHydroFixed_NormalDamp050_ContactAudit_0013'
OLD = '*Contact Damping, definition=CRITICAL DAMPING FRACTION, tangent fraction=0.0\n0.20'
NEW = '*Contact Damping, definition=CRITICAL DAMPING FRACTION, tangent fraction=0.0\n0.50'


def main():
    source_path = HERE / (SOURCE + '.inp')
    source = source_path.read_text()
    assert source.count(OLD) == 1
    candidate = source.replace(OLD, NEW)
    assert candidate.replace(NEW, OLD) == source
    for token in ('0.03,', '1.0e-7, 0.001300', 'phase=248.0 deg',
                  'gradient=6.000 mT=0.006000 T'):
        assert token in candidate, token
    bridge = REPO / 'calibration_analysis' / 'Wobble30Hz_reduced_hydrodynamics' / 'vuforc_socket_bridge_reduced_hydro.f'
    bridge_text = bridge.read_text()
    for token in ('cpar /4.0D-9/', 'cperp /1.2D-8/', 'kspin /1.0D-9/', 'kwob /3.0D-9/'):
        assert token in bridge_text, token

    local = HERE / (JOB + '.inp')
    work = WORKDIR / (JOB + '.inp')
    local.write_text(candidate)
    work.write_text(candidate)

    source_runner = (HERE / 'run_normaldamp020_contact_audit_0013.ps1').read_text()
    runner = source_runner.replace(SOURCE, JOB).replace('[int]$Port=65506', '[int]$Port=65507')
    runner = runner.replace('NORMAL_DAMPING=0.20', 'NORMAL_DAMPING=0.50')
    runner = runner.replace("if(-not(Test-Path \"$Job.inp\")){throw \"Missing input deck $Job.inp\"}",
                            "if(Test-Path \"$Job.odb\"){throw \"Existing ODB; refusing repeat run\"}; if(-not(Test-Path \"$Job.inp\")){throw \"Missing input deck $Job.inp\"}")
    (HERE / 'run_normaldamp050_contact_audit_0013.ps1').write_text(runner)

    datacheck = """$ErrorActionPreference='Stop'
Set-Location -LiteralPath 'J:\\abaqusfangzhen'
. '.\\launch_abaqus_with_env.ps1' -SkipVerify
$job='Wobble_F30_G6L45_ReducedHydroFixed_NormalDamp050_ContactAudit_0013_Datacheck'
$input='Wobble_F30_G6L45_ReducedHydroFixed_NormalDamp050_ContactAudit_0013.inp'
if(Test-Path "$job.odb"){throw 'Existing datacheck ODB; refusing repeat'}
& 'I:\\SIMULIA\\Commands\\abq2025.bat' job=$job input=$input user='J:\\abaqusfangzhen\\abaqus_robot\\calibration_analysis\\Wobble30Hz_reduced_hydrodynamics\\vuforc_socket_bridge_reduced_hydro.f' datacheck double=both cpus=1 ask_delete=OFF interactive
if($LASTEXITCODE -ne 0){throw "Datacheck failed: $LASTEXITCODE"}
if(Select-String -LiteralPath "$job.dat" -Pattern '\\*\\*\\*ERROR'){throw 'Datacheck ERROR found'}
Write-Host 'DATACHECK PASS: zeta=0.50, tangent fraction=0.0 accepted by Abaqus 2025'
"""
    (HERE / 'datacheck_normaldamp050.ps1').write_text(datacheck)

    identity = {
        'label': 'PROVISIONAL_NORMAL_DAMPING_PROBE',
        'reference_job': SOURCE,
        'candidate_job': JOB,
        'reference_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
        'candidate_sha256': hashlib.sha256(local.read_bytes()).hexdigest(),
        'only_physics_change': {'critical_damping_fraction': [0.20, 0.50]},
        'unchanged_contact_parameter': {'tangent_fraction': 0.0},
        'selection_basis': {
            'empirical_de_dzeta_055_to_020': -0.81535,
            'rough_extrapolated_e_at_050': 0.5813,
            'qualification': 'candidate selection only; not theory, experiment, or Abaqus mapping',
        },
        'frozen': ['HARD pressure-overclosure', 'General Contact', 'SmoothWall114',
                   'friction 0.03', 'surfaces and penalty formulation', 'dt 1e-7 s',
                   'magnetic architecture', 'robot geometry/mass/inertia/pose', 'Reduced-Hydro'],
    }
    (HERE / 'normaldamp050_input_identity.json').write_text(json.dumps(identity, indent=2))
    print(json.dumps(identity, indent=2))


if __name__ == '__main__':
    main()
