[CmdletBinding()]
param([int]$Port = 65534)
$ErrorActionPreference = 'Stop'
$repo = 'J:\abaqusfangzhen\abaqus_robot'
$out = Join-Path $repo 'calibration_analysis\RefinedDualEndTrueCELGate'
$job = 'TAIL_STICK_F120_ZETA100_8P333MS_TRUECEL'
$case = Join-Path $out ('case\' + $job)
$identityPath = Join-Path $case 'case_identity.json'
$identity = Get-Content -Raw $identityPath | ConvertFrom-Json
if ($identity.status -ne 'PREPARED' -or $identity.dynamics_run_count -ne 0) { throw 'One-run gate refused by case identity' }
if (Test-Path (Join-Path $case ($job + '.lck'))) { throw 'Existing solver lock' }
if (Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue) { throw "Port $Port busy" }

$server = Join-Path $repo 'calibration_analysis\FastStraightDynamicScreen\production\magpylib_socket_server_fast.py'
$python = 'I:\SIMULIA\EstProducts\2025\win_b64\tools\SMApy\python3.10\python.exe'
$transform = 'J:\abaqusfangzhen\abaqus_magpylib_frame_transform.json'
$telemetry = Join-Path $case ($job + '_telemetry.csv')
$stdout = Join-Path $case ($job + '_socket_stdout.log')
$stderr = Join-Path $case ($job + '_socket_stderr.log')
$m = $identity.initial_magnetic_moment_axis_aba
$args = @('-u',$server,'--host','127.0.0.1','--port',"$Port",'--job-name',$job,
    '--dxf',(Join-Path $case 'straight_control_centerline.dxf'),'--drive-type','analytic',
    '--field-frame-mode','ROBOT_LOCAL_ELLIPTIC_ROCKING','--rocking-amplitude-deg','14.8',
    '--rocking-cross-amplitude-deg','2.5','--rocking-frame-azimuth-deg',([string]$identity.routeA_gauge_for_flipped_body_deg),
    '--robot-diameter-mm','0.815','--robot-height-mm','2.6','--robot-br-t','1.46',
    '--robot-moment-Am2',([string]$identity.robot_moment_Am2),'--robot-mass-mg',([string]$identity.robot_mass_mg),
    '--analytic-b-t','0.010','--analytic-follow-robot','--analytic-gradient-b-t','0.00015',
    '--analytic-gradient-length-mm','45','--analytic-gradient-profile','legacy',
    '--driver-speed-mm-s','6','--bend-speed-mm-s','4.5','--bend-start-mm','13.49','--bend-end-mm','18.56',
    '--z-offset-mm','90','--driver-start-offset-mm','18.899960626','--spin-hz','120',
    '--cone-half-angle-deg','30','--driver-xy-shift-mm','2','-6','--cone-axis-bias-deg','0',
    '--cone-axis','tangent','--robot-polarity','1','--robot-moment-axis-aba',([string]$m[0]),([string]$m[1]),([string]$m[2]),
    '--phase-deg','0','--analytic-rotation-sense','1','--ramp-time-s','0.001','--endpoint-taper-mm','1',
    '--adaptive-lead-mm','0','--frame-transform-json',$transform,'--telemetry-csv',$telemetry,
    '--telemetry-interval-s','0.00002','--frame-log-interval-s','0.0005')

. 'J:\abaqusfangzhen\launch_abaqus_with_env.ps1' -SkipVerify
$env:MAGPY_SOCKET_PORT = "$Port"
Set-Location -LiteralPath $case
$process = Start-Process -FilePath $python -ArgumentList $args -WorkingDirectory $case -WindowStyle Hidden `
    -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
$identity.status = 'RUNNING'; $identity.dynamics_run_count = 1
$identity | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath $identityPath -Encoding ascii
$started = Get-Date
try {
    $deadline = (Get-Date).AddSeconds(30)
    do { Start-Sleep -Milliseconds 200; $listen = Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue } while (-not $listen -and (Get-Date) -lt $deadline)
    if (-not $listen) { throw 'Magpylib socket startup failed' }
    $inputFile = $job + '.inp'
    & 'I:\SIMULIA\Commands\abq2025.bat' job=$job input=$inputFile user='vuforc_magnetic_only.f' double=both cpus=1 ask_delete=ON interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus exit $LASTEXITCODE" }
    $sta = Get-Content -Raw (Join-Path $case ($job + '.sta'))
    if ($sta -notmatch 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY') { throw 'Dynamics incomplete' }
    $identity.status = 'SOLVED'
} catch {
    $identity.status = 'FAILED'
    $identity | Add-Member -NotePropertyName failure -NotePropertyValue $_.Exception.Message -Force
    throw
} finally {
    $identity | Add-Member -NotePropertyName wallclock_s -NotePropertyValue (((Get-Date) - $started).TotalSeconds) -Force
    $identity | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath $identityPath -Encoding ascii
    if ($process -and -not $process.HasExited) { try { $process.Kill() } catch {} }
}
Write-Host ("COMPLETE {0} wallclock={1:N1}s" -f $job,((Get-Date)-$started).TotalSeconds)
