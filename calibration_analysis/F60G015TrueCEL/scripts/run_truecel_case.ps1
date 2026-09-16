[CmdletBinding()]
param([int]$Port = 65532)
$ErrorActionPreference = 'Stop'
$repo = 'J:\abaqusfangzhen\abaqus_robot'
$out = Join-Path $repo 'calibration_analysis\F60G015TrueCEL'
$case = Join-Path $out 'case\S4_HEADFORWARD_F60_G0P15_TRUECEL'
$job = 'S4_HEADFORWARD_F60_G0P15_TRUECEL'
$server = Join-Path $repo 'calibration_analysis\FastStraightDynamicScreen\production\magpylib_socket_server_fast.py'
$python = 'I:\SIMULIA\EstProducts\2025\win_b64\tools\SMApy\python3.10\python.exe'
$transform = 'J:\abaqusfangzhen\abaqus_magpylib_frame_transform.json'
$extractRp = Join-Path $repo 'calibration_analysis\ProductionLocalFrameValidation\scripts\extract_dynamic_odb.py'
$extractCel = Join-Path $out 'scripts\extract_truecel_odb.py'
$identityPath = Join-Path $case 'case_identity.json'
$presencePath = Join-Path $out 'fluid_presence_gate.json'
$magneticPath = Join-Path $out 'magnetic_baseline_gate.json'
$identity = Get-Content -Raw $identityPath | ConvertFrom-Json
$presence = Get-Content -Raw $presencePath | ConvertFrom-Json
$magnetic = Get-Content -Raw $magneticPath | ConvertFrom-Json
if ($identity.status -ne 'PREPARED') { throw "Case status is $($identity.status), expected PREPARED" }
if ($identity.dynamics_run_count -ne 0) { throw 'Dynamics was already launched; refusing a second run' }
if (-not $presence.passed -or $presence.fluid_filled_volume_mm3 -le 0 -or $presence.initial_total_fluid_mass_mg -le 0 -or $presence.fluid_containing_element_count -le 0) { throw 'EMPTY_EULERIAN_DOMAIN' }
if (-not $magnetic.passed -or -not $magnetic.F_gradient_dot_canonical_s_positive) { throw 'Magnetic baseline gate failed' }
if (Test-Path (Join-Path $case ($job + '.lck'))) { throw 'Existing solver lock; refusing concurrent run' }
if (Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue) { throw "Port $Port busy" }

$deckText = Get-Content -Raw (Join-Path $case ($job + '.inp'))
$fortranText = Get-Content -Raw (Join-Path $case 'vuforc_magnetic_only.f')
foreach ($token in @('HYDRO_FX','HYDRO_FY','HYDRO_FZ','HYDRO_MX','HYDRO_MY','HYDRO_MZ','cpar','cperp','kspin','kwob')) {
    if ($deckText -match "name=$token" -or $fortranText -match "(?m)^\s*DATA\s+$token") { throw "ReducedHydro token survived: $token" }
}
foreach ($token in @('*Element, type=EC3D8R','*Eulerian Section','*Initial Conditions, type=VOLUME FRACTION','*Initial Conditions, type=VELOCITY','FLUID_CEL_SURF')) {
    if (-not $deckText.Contains($token)) { throw "Required true-CEL token missing: $token" }
}

$envScript = 'J:\abaqusfangzhen\launch_abaqus_with_env.ps1'
. $envScript -SkipVerify
$ErrorActionPreference = 'Continue'
Set-Location -LiteralPath $case
$gateJob = $job + '_DATACHECK_GATE'
$gateOutput = & 'I:\SIMULIA\Commands\abq2025.bat' job=$gateJob input="$job.inp" user='vuforc_magnetic_only.f' double=both datacheck ask_delete=ON interactive 2>&1
$gateExit = $LASTEXITCODE
$gateOutput | Out-File -LiteralPath (Join-Path $out ($gateJob + '.log')) -Encoding ascii
$gateDat = Get-Content -Raw (Join-Path $case ($gateJob + '.dat')) -ErrorAction SilentlyContinue
$gateSta = Get-Content -Raw (Join-Path $case ($gateJob + '.sta')) -ErrorAction SilentlyContinue
$gatePassed = ($gateExit -eq 0 -and $gateDat -notmatch '(?i)\*\*\*ERROR|EMPTY_EULERIAN|element\s+0' -and $gateSta -match 'ANALYSIS HAS COMPLETED SUCCESSFULLY' -and $gateSta -match 'FLUID_EULERIAN-1')
$gate = [ordered]@{
    job = $gateJob
    passed = $gatePassed
    exit_code = $gateExit
    fluid_filled_volume_mm3 = [double]$presence.fluid_filled_volume_mm3
    initial_total_fluid_mass_mg = [double]$presence.initial_total_fluid_mass_mg
    fluid_containing_element_count = [int]$presence.fluid_containing_element_count
    initial_stable_dt_s = 5.51093e-8
    classification = 'STRICT_CEL_FLUID_PRESENT'
}
$gate | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $out 'truecel_datacheck_gate.json') -Encoding ascii
if (-not $gatePassed) { throw 'True-CEL datacheck failed; no dynamics launched' }

$env:MAGPY_SOCKET_PORT = "$Port"
$telemetry = Join-Path $case ($job + '_telemetry.csv')
$stdout = Join-Path $case ($job + '_socket_stdout.log')
$stderr = Join-Path $case ($job + '_socket_stderr.log')
$solverLog = Join-Path $out ($job + '_solver.log')
$extractLog = Join-Path $out ($job + '_extract.log')
$dxf = Join-Path $case 'straight_control_centerline.dxf'
$m = $identity.initial_magnetic_moment_axis_aba
$arguments = @('-u',$server,'--host','127.0.0.1','--port',"$Port",'--job-name',$job,
    '--dxf',$dxf,'--drive-type','analytic','--field-frame-mode','ROBOT_LOCAL_ELLIPTIC_ROCKING',
    '--rocking-amplitude-deg','14.343111711438091','--rocking-cross-amplitude-deg','2.5',
    '--rocking-frame-azimuth-deg',([string]$identity.routeA_gauge_for_flipped_body_deg),
    '--robot-diameter-mm','0.815','--robot-height-mm','2.4','--robot-br-t','1.46',
    '--robot-moment-Am2',([string]$identity.robot_moment_Am2),'--robot-mass-mg',([string]$identity.robot_mass_mg),
    '--analytic-b-t','0.010','--analytic-follow-robot','--analytic-gradient-b-t','0.00015',
    '--analytic-gradient-length-mm','45','--analytic-gradient-profile','legacy',
    '--driver-speed-mm-s','6','--bend-speed-mm-s','4.5','--bend-start-mm','13.49','--bend-end-mm','18.56',
    '--z-offset-mm','90','--driver-start-offset-mm','18.899960626','--spin-hz','60',
    '--cone-half-angle-deg','30','--driver-xy-shift-mm','2','-6','--cone-axis-bias-deg','0',
    '--cone-axis','tangent','--robot-polarity','1','--robot-moment-axis-aba',([string]$m[0]),([string]$m[1]),([string]$m[2]),
    '--phase-deg','0','--analytic-rotation-sense','1','--ramp-time-s','0.001','--endpoint-taper-mm','1','--adaptive-lead-mm','0',
    '--frame-transform-json',$transform,'--telemetry-csv',$telemetry,'--telemetry-interval-s','0.00002','--frame-log-interval-s','0.0005')
$process = Start-Process -FilePath $python -ArgumentList $arguments -WorkingDirectory $case -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
$started = Get-Date
$identity.status = 'RUNNING'
$identity.dynamics_run_count = 1
$identity | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath $identityPath -Encoding ascii
try {
    $deadline = (Get-Date).AddSeconds(30)
    do { Start-Sleep -Milliseconds 200; $listen = Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue } while (-not $listen -and (Get-Date) -lt $deadline)
    if (-not $listen) { throw "Socket startup failed; see $stderr" }
    $solverOutput = & 'I:\SIMULIA\Commands\abq2025.bat' job=$job input="$job.inp" user='vuforc_magnetic_only.f' double=both cpus=1 ask_delete=ON interactive 2>&1
    $solverExit = $LASTEXITCODE
    $solverOutput | Out-File -LiteralPath $solverLog -Encoding ascii
    if ($solverExit -ne 0) { throw "Abaqus exit $solverExit" }
    $staText = Get-Content -Raw (Join-Path $case ($job + '.sta'))
    if ($staText -notmatch 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY') { throw 'Dynamics did not complete successfully' }
    $private = Join-Path $case 'private'
    $extractOutput = & 'I:\SIMULIA\Commands\abq2025.bat' python $extractRp $job $private 2>&1
    $extractOutput += & 'I:\SIMULIA\Commands\abq2025.bat' python $extractCel $job $private 2>&1
    $extractExit = $LASTEXITCODE
    $extractOutput | Out-File -LiteralPath $extractLog -Encoding ascii
    if ($extractExit -ne 0) { throw 'ODB extraction failed' }
    $identity.status = 'SOLVED'
    $identity | Add-Member -NotePropertyName wallclock_s -NotePropertyValue (((Get-Date) - $started).TotalSeconds) -Force
    $identity | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath $identityPath -Encoding ascii
} catch {
    $identity.status = 'FAILED'
    $identity | Add-Member -NotePropertyName failure -NotePropertyValue $_.Exception.Message -Force
    $identity | Add-Member -NotePropertyName wallclock_s -NotePropertyValue (((Get-Date) - $started).TotalSeconds) -Force
    $identity | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath $identityPath -Encoding ascii
    throw
} finally {
    if ($process -and -not $process.HasExited) { try { $process.Kill() } catch {} }
}
Write-Host ("COMPLETE {0} wallclock={1:N1}s" -f $job,((Get-Date)-$started).TotalSeconds)
