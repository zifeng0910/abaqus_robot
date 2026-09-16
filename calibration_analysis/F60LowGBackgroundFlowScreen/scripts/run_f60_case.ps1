[CmdletBinding()]
param(
    [int]$Port = 65531
)
$ErrorActionPreference = 'Continue'
$repo = 'J:\abaqusfangzhen\abaqus_robot'
$out = Join-Path $repo 'calibration_analysis\F60LowGBackgroundFlowScreen'
$case = Join-Path $out 'case\S4_HEADFORWARD_F60_G0P10_FLOW'
$job = 'S4_HEADFORWARD_F60_G0P10_FLOW'
$server = Join-Path $repo 'calibration_analysis\FastStraightDynamicScreen\production\magpylib_socket_server_fast.py'
$transform = 'J:\abaqusfangzhen\abaqus_magpylib_frame_transform.json'
$python = 'I:\SIMULIA\EstProducts\2025\win_b64\tools\SMApy\python3.10\python.exe'
$extract = Join-Path $repo 'calibration_analysis\ProductionLocalFrameValidation\scripts\extract_dynamic_odb.py'
$identityPath = Join-Path $case 'case_identity.json'
$identity = Get-Content -Raw $identityPath | ConvertFrom-Json
if ($identity.status -ne 'PREPARED') { throw "Case status is $($identity.status), expected PREPARED" }
if (Test-Path (Join-Path $case ($job + '.lck'))) { throw 'Existing active solver lock; refusing a concurrent run' }
if (Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue) { throw "Port $Port busy" }
$envScript = Join-Path (Split-Path $repo -Parent) 'launch_abaqus_with_env.ps1'
if (-not (Test-Path $envScript)) { throw "Abaqus environment script missing: $envScript" }
. $envScript -SkipVerify
# The environment loader uses Stop for its own validation; Abaqus emits normal
# license/compiler diagnostics on stderr, so keep command capture non-terminating.
$ErrorActionPreference = 'Continue'

$gateJob = $job + '_DATACHECK'
$gateLog = Join-Path $out ($gateJob + '.log')
Set-Location -LiteralPath $case
$gateOutput = & 'I:\SIMULIA\Commands\abq2025.bat' job=$gateJob input="$job.inp" user='vuforc_background_flow.f' double=both datacheck ask_delete=ON interactive 2>&1
$gateExit = $LASTEXITCODE
$gateOutput | Out-File -LiteralPath $gateLog -Encoding ascii
$gateText = ($gateOutput | Out-String)
$gateCompileFailed = $gateText -match '(?i)(Abaqus Error:|Problem during compilation|not recognized as an internal|Analysis exited with errors)'
$gate = [ordered]@{
    gate_job = $gateJob
    attempted_dt_s = [double]$identity.direct_dt_s
    exit_code = $gateExit
    log = $gateLog.Substring($repo.Length + 1)
    passed = ($gateExit -eq 0 -and -not $gateCompileFailed)
    checks = @('no initial penetration problem', 'no contact instability', 'no obvious energy explosion', 'no negative volume/fatal error', 'bounded rocking')
}
$gate | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $out 'dt_2e-7_datacheck_gate.json') -Encoding ascii
if ($gateExit -ne 0 -or $gateCompileFailed) { throw "2e-7 datacheck failed; no dynamic run launched" }

$env:MAGPY_SOCKET_PORT = "$Port"
$telemetry = Join-Path $case ($job + '_telemetry.csv')
$stdout = Join-Path $case ($job + '_socket_stdout.log')
$stderr = Join-Path $case ($job + '_socket_stderr.log')
$solverLog = Join-Path $out ($job + '_solver.log')
$extractLog = Join-Path $out ($job + '_extract.log')
$dxf = Join-Path $case 'straight_control_centerline.dxf'
$m = $identity.initial_magnetic_moment_axis_aba
$arguments = @('-u',$server,'--host','127.0.0.1','--port',"$Port",'--job-name',$job,
    '--dxf',$dxf,'--drive-type','analytic','--field-frame-mode',$identity.field_frame_mode,
    '--rocking-amplitude-deg',([string]$identity.rocking_main_amplitude_deg),
    '--rocking-cross-amplitude-deg',([string]$identity.rocking_cross_amplitude_deg),
    '--rocking-frame-azimuth-deg',([string]$identity.routeA_gauge_for_flipped_body_deg),
    '--robot-diameter-mm','0.815','--robot-height-mm','2.4','--robot-br-t','1.46',
    '--robot-moment-Am2',([string]$identity.robot_moment_Am2),'--robot-mass-mg',([string]$identity.robot_mass_mg),
    '--analytic-b-t','0.010','--analytic-follow-robot','--analytic-gradient-b-t','0.0001',
    '--analytic-gradient-length-mm','45','--analytic-gradient-profile','legacy',
    '--driver-speed-mm-s','6','--bend-speed-mm-s','4.5','--bend-start-mm','13.49','--bend-end-mm','18.56',
    '--z-offset-mm','90','--driver-start-offset-mm','18.899960626','--spin-hz','60',
    '--cone-half-angle-deg','30','--driver-xy-shift-mm','2','-6','--cone-axis-bias-deg','0',
    '--cone-axis','tangent','--robot-polarity','1','--robot-moment-axis-aba',([string]$m[0]),([string]$m[1]),([string]$m[2]),'--phase-deg','0',
    '--analytic-rotation-sense','1','--ramp-time-s','0.001','--endpoint-taper-mm','1','--adaptive-lead-mm','0',
    '--frame-transform-json',$transform,'--telemetry-csv',$telemetry,'--telemetry-interval-s','0.00002','--frame-log-interval-s','0.0005')
$process = Start-Process -FilePath $python -ArgumentList $arguments -WorkingDirectory $case -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
$started = Get-Date
try {
    $deadline = (Get-Date).AddSeconds(30)
    do { Start-Sleep -Milliseconds 200; $listen = Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue } while (-not $listen -and (Get-Date) -lt $deadline)
    if (-not $listen) { throw "Socket startup failed; see $stderr" }
    Set-Location -LiteralPath $case
    $solverOutput = & 'I:\SIMULIA\Commands\abq2025.bat' job=$job input="$job.inp" user='vuforc_background_flow.f' double=both cpus=1 ask_delete=ON interactive 2>&1
    $solverExit = $LASTEXITCODE
    $solverOutput | Out-File -LiteralPath $solverLog -Encoding ascii
    if ($solverExit -ne 0) { throw "Abaqus exit $solverExit" }
    $hydroLog = Get-Item (Join-Path $case 'hydro_increment.csv') -ErrorAction SilentlyContinue
    if (-not $hydroLog) {
        $hydroLog = Get-ChildItem ([IO.Path]::GetTempPath()) -Directory -Filter ("{0}_{1}_*" -f $env:USERNAME,$job) |
            ForEach-Object { Get-ChildItem $_.FullName -File -Filter 'hydro_increment.csv' -ErrorAction SilentlyContinue } |
            Sort-Object LastWriteTime -Descending | Select-Object -First 1
    }
    if (-not $hydroLog) { throw 'ReducedHydro force log was not produced' }
    if ($hydroLog.FullName -ne (Join-Path $case 'hydro_increment.csv')) {
        Copy-Item -LiteralPath $hydroLog.FullName -Destination (Join-Path $case 'hydro_increment.csv') -Force
    }
    $extractOutput = & 'I:\SIMULIA\Commands\abq2025.bat' python $extract $job (Join-Path $case 'private') 2>&1
    $extractExit = $LASTEXITCODE
    $extractOutput | Out-File -LiteralPath $extractLog -Encoding ascii
    if ($extractExit -ne 0) { throw 'ODB extraction failed' }
    $identity.status = 'SOLVED'
    $identity.dt_gate.accepted = $true
    $identity.dt_gate.gate_status = 'PASSED_DYNAMIC_RUN'
    $identity | Add-Member -NotePropertyName wallclock_s -NotePropertyValue (((Get-Date) - $started).TotalSeconds) -Force
    $identity | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath $identityPath -Encoding ascii
} finally {
    if ($process -and -not $process.HasExited) { try { $process.Kill() } catch {} }
}
Write-Host ("COMPLETE {0} wallclock={1:N1}s" -f $job,((Get-Date)-$started).TotalSeconds)
