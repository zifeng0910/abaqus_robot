[CmdletBinding()]
param()
$ErrorActionPreference = 'Stop'
$repo = 'J:\abaqusfangzhen\abaqus_robot'
$out = Join-Path $repo 'calibration_analysis\TrueCELLongForwardTransit'
$job = 'TRUECEL_A14P5_TAILGAP_REFINE_2CYCLES'
$case = Join-Path $out ('case\' + $job)
$identityPath = Join-Path $case 'case_identity.json'
$identity = Get-Content -Raw $identityPath | ConvertFrom-Json
$setup = Get-Content -Raw (Join-Path $out ($job + '_Setup_Audit.json')) | ConvertFrom-Json
$audit = Get-Content -Raw (Join-Path $out ($job + '_PreSolve_Initialization_Audit.json')) | ConvertFrom-Json

if ($identity.status -ne 'PREPARED' -or $identity.dynamics_run_count -ne 0) {
    throw 'One-run gate refused: candidate is not PREPARED with dynamics_run_count=0'
}
if (-not $audit.initialization_gate_passed) { throw 'Initialization gate failed' }
if ($audit.max_EVF_geometrically_inside_robot -gt 1e-6) { throw 'Robot interior EVF gate failed' }
if ($audit.estimated_water_overlap_volume_upper_bound_mm3 -gt 1e-12) { throw 'Water overlap gate failed' }
if ($audit.robot_wall_initial_penetration_mm -gt 1e-9) { throw 'Initial wall penetration gate failed' }
if ($setup.minimum_cell_size_mm -lt 0.0370 -or $setup.minimum_cell_size_mm -gt 0.0400) {
    throw 'Local mesh size is outside the accepted 0.0375-0.040 mm band'
}
if ($setup.target_corridor_refined_elements -le 0) { throw 'No effective TAIL-corridor elements were refined' }
if (-not $setup.physics_frozen) { throw 'Physics freeze audit failed' }

$deck = Get-Content -Raw (Join-Path $case ($job + '.inp'))
foreach ($required in @(
    '*Dynamic, Explicit, SCALE FACTOR=0.4',
    '*Element, type=EC3D8R',
    'FLUID_CEL_INIT_ALL',
    'PROP_CEL_FLUID_HARD',
    'SECOND SURFACE=Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF',
    ', 0.016666666667'
)) {
    if (-not $deck.Contains($required)) { throw "Missing frozen-model token: $required" }
}
foreach ($forbidden in @('name=HYDRO_', '*Fixed Mass Scaling', '*Variable Mass Scaling', '*Contact Pair')) {
    if ($deck.Contains($forbidden)) { throw "Forbidden token: $forbidden" }
}

. 'J:\abaqusfangzhen\launch_abaqus_with_env.ps1' -SkipVerify
Set-Location -LiteralPath $case
$sub = Join-Path $case 'vuamp_precomputed_truecel.f90'
$gateJob = $job + '_DATACHECK'
$gateLog = Join-Path $out ($gateJob + '.log')
& 'I:\SIMULIA\Commands\abq2025.bat' job=$gateJob input="$job.inp" user=$sub double=both datacheck ask_delete=ON interactive *> $gateLog
$gateExit = $LASTEXITCODE
$gateSta = Join-Path $case ($gateJob + '.sta')
$gateDat = Join-Path $case ($gateJob + '.dat')
$staText = if (Test-Path $gateSta) { Get-Content -Raw $gateSta } else { '' }
$datText = if (Test-Path $gateDat) { Get-Content -Raw $gateDat } else { '' }
$gatePassed = ($gateExit -eq 0 -and $staText -match 'ANALYSIS HAS COMPLETED SUCCESSFULLY' -and
               $datText -notmatch '(?i)\*\*\*ERROR|EMPTY_EULERIAN|NEGATIVE VOLUME|ZERO OR NEGATIVE')
[ordered]@{
    job = $gateJob
    passed = $gatePassed
    exit_code = $gateExit
    dynamics_launched_by_gate = $false
    local_mesh_size_mm = $setup.minimum_cell_size_mm
    target_corridor_refined_elements = $setup.target_corridor_refined_elements
} | ConvertTo-Json | Set-Content (Join-Path $out ($job + '_datacheck_gate.json')) -Encoding ascii
if (-not $gatePassed) { throw 'Datacheck failed; dynamics was not launched' }

$identity.status = 'RUNNING'
$identity.dynamics_run_count = 1
$identity | ConvertTo-Json -Depth 20 | Set-Content $identityPath -Encoding ascii
$started = Get-Date
try {
    $solverLog = Join-Path $out ($job + '_solver.log')
    & 'I:\SIMULIA\Commands\abq2025.bat' job=$job input="$job.inp" user=$sub double=both cpus=1 ask_delete=ON interactive *> $solverLog
    if ($LASTEXITCODE -ne 0) { throw "Abaqus dynamics exit code $LASTEXITCODE" }
    $sta = Get-Content -Raw (Join-Path $case ($job + '.sta'))
    if ($sta -notmatch 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY') { throw 'Dynamics incomplete' }

    $private = Join-Path $case 'private'
    New-Item -ItemType Directory -Path $private -Force | Out-Null
    $extractLog = Join-Path $out ($job + '_extract.log')
    $extractors = @(
        (Join-Path $repo 'calibration_analysis\ProductionLocalFrameValidation\scripts\extract_dynamic_odb.py'),
        (Join-Path $repo 'calibration_analysis\RefinedDualEndTrueCELGate\scripts\extract_robot_contact_fields.py'),
        (Join-Path $repo 'calibration_analysis\TrueCELShort3msTransfer\scripts\extract_truecel_fields.py'),
        (Join-Path $out 'scripts\extract_force_flux_histories.py')
    )
    foreach ($extractor in $extractors) {
        & 'I:\SIMULIA\Commands\abq2025.bat' python $extractor $job $private *>> $extractLog
        if ($LASTEXITCODE -ne 0) { throw "Extraction failed: $extractor" }
    }
    $identity.status = 'SOLVED'
}
catch {
    $identity.status = 'FAILED'
    $identity | Add-Member failure $_.Exception.Message -Force
    throw
}
finally {
    $identity | Add-Member wallclock_s (((Get-Date) - $started).TotalSeconds) -Force
    $identity | Add-Member cpus 1 -Force
    $identity | Add-Member socket_calls 0 -Force
    $identity | ConvertTo-Json -Depth 20 | Set-Content $identityPath -Encoding ascii
}
Write-Host ("COMPLETE {0}; wallclock={1:N1}s" -f $job, $identity.wallclock_s)
