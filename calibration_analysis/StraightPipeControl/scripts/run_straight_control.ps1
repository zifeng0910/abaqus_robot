[CmdletBinding()]
param([int]$Port = 65431)
$ErrorActionPreference = 'Stop'
$repo = 'J:\abaqusfangzhen\abaqus_robot'
$out = Join-Path $repo 'calibration_analysis\StraightPipeControl'
$job = 'PROD_LOCAL30_G0_STRAIGHT_CTRL'
$case = Join-Path $out ('case\' + $job)
$server = 'J:\magpy\magpylib_socket_server.py'
$transform = 'J:\abaqusfangzhen\abaqus_magpylib_frame_transform.json'
$python = 'I:\SIMULIA\EstProducts\2025\win_b64\tools\SMApy\python3.10\python.exe'
$extract = Join-Path $repo 'calibration_analysis\ProductionLocalFrameValidation\scripts\extract_dynamic_odb.py'
$identityPath = Join-Path $case 'case_identity.json'
$identity = Get-Content -Raw $identityPath | ConvertFrom-Json
$vendoredHash = (Get-FileHash (Join-Path $repo 'calibration_analysis\ProductionLocalFrameValidation\production\magpylib_socket_server.py') -Algorithm SHA256).Hash
$liveHash = (Get-FileHash $server -Algorithm SHA256).Hash
if ($vendoredHash -ne $liveHash -or $liveHash -ne $identity.production_server_sha256) { throw 'Production server identity mismatch' }
if ($identity.status -ne 'PREPARED') { throw "Case status is $($identity.status), expected PREPARED" }
if (Test-Path (Join-Path $case ($job + '.odb'))) { throw 'Existing ODB; refusing a second dynamic run' }
if (Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue) { throw "Port $Port busy" }
$env:MAGPY_SOCKET_PORT = "$Port"
$telemetry = Join-Path $case ($job + '_telemetry.csv')
$stdout = Join-Path $case ($job + '_socket_stdout.log')
$stderr = Join-Path $case ($job + '_socket_stderr.log')
$dxf = Join-Path $case 'straight_control_centerline.dxf'
$arguments = @('-u',$server,'--host','127.0.0.1','--port',"$Port",'--job-name',$job,
    '--dxf',$dxf,'--drive-type','analytic','--field-frame-mode','ROBOT_LOCAL_TANGENT',
    '--robot-diameter-mm','0.815','--robot-height-mm','2.4','--robot-br-t','1.46',
    '--robot-moment-Am2','0.0010876227522174417','--robot-mass-mg','9.207793514166587',
    '--analytic-b-t','0.010','--analytic-follow-robot','--analytic-gradient-b-t','0',
    '--analytic-gradient-length-mm','45','--analytic-gradient-profile','legacy',
    '--driver-speed-mm-s','6','--bend-speed-mm-s','4.5','--bend-start-mm','13.49','--bend-end-mm','18.56',
    '--z-offset-mm','90','--driver-start-offset-mm','18.899960626','--spin-hz','30',
    '--cone-half-angle-deg','30','--driver-xy-shift-mm','2','-6','--cone-axis-bias-deg','0',
    '--cone-axis','tangent','--robot-polarity','-1','--robot-axis-tangent','--phase-deg','248',
    '--analytic-rotation-sense','1','--ramp-time-s','0.001','--endpoint-taper-mm','1',
    '--adaptive-lead-mm','0','--frame-transform-json',$transform,'--telemetry-csv',$telemetry,
    '--telemetry-interval-s','0.00001','--frame-log-interval-s','0.0001')
$process = Start-Process -FilePath $python -ArgumentList $arguments -WorkingDirectory 'J:\magpy' -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
$started = Get-Date
try {
    $deadline = (Get-Date).AddSeconds(30)
    do { Start-Sleep -Milliseconds 200; $listen = Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue } while (-not $listen -and (Get-Date) -lt $deadline)
    if (-not $listen) { throw "Socket startup failed; see $stderr" }
    . 'J:\abaqusfangzhen\launch_abaqus_with_env.ps1' -SkipVerify
    Set-Location -LiteralPath $case
    & 'I:\SIMULIA\Commands\abq2025.bat' job=$job input="$job.inp" user='vuforc_production_local.f' double=both cpus=1 ask_delete=OFF interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus exit $LASTEXITCODE" }
    & 'I:\SIMULIA\Commands\abq2025.bat' python $extract $job (Join-Path $case 'private')
    if ($LASTEXITCODE -ne 0) { throw 'ODB extraction failed' }
    if (Test-Path 'hydro_increment.csv') { Move-Item -LiteralPath 'hydro_increment.csv' -Destination (Join-Path $case ($job + '_hydro_increment.csv')) -Force }
    $identity.status = 'SOLVED'
    $identity | Add-Member -NotePropertyName wallclock_s -NotePropertyValue (((Get-Date) - $started).TotalSeconds) -Force
    $identity | ConvertTo-Json -Depth 10 | Set-Content $identityPath
} finally {
    if ($process -and -not $process.HasExited) { try { $process.Kill() } catch {} }
}
Write-Host ("COMPLETE {0} wallclock={1:N1}s" -f $job,((Get-Date)-$started).TotalSeconds)
