[CmdletBinding()]
param(
    [string[]]$CaseId = @('FRAME_CURRENT','FRAME_LOCAL30','FRAME_LOCAL40'),
    [int]$BasePort = 65100
)
$ErrorActionPreference = 'Stop'
$repo = 'J:\abaqusfangzhen\abaqus_robot'
$audit = Join-Path $repo 'calibration_analysis\MagneticDriveFrameAudit'
$productionServer = 'J:\magpy\magpylib_socket_server.py'
$localServer = Join-Path $audit 'scripts\local_tangent_socket_server.py'
$transform = 'J:\abaqusfangzhen\abaqus_magpylib_frame_transform.json'
$python = 'I:\SIMULIA\EstProducts\2025\win_b64\tools\SMApy\python3.10\python.exe'
$extract = Join-Path $repo 'calibration_analysis\CoarseMotionModeScreen\scripts\extract_case_odb.py'
. 'J:\abaqusfangzhen\launch_abaqus_with_env.ps1' -SkipVerify
$index = 0
foreach ($caseName in $CaseId) {
    $folder = Join-Path $audit ('cases\' + $caseName)
    $identityPath = Join-Path $folder 'case_identity.json'
    if (-not (Test-Path $identityPath)) { throw "Unknown case $caseName" }
    $case = Get-Content -Raw $identityPath | ConvertFrom-Json
    $job = $case.job_name
    $rampTime = if ($null -ne $case.ramp_time_s) { [double]$case.ramp_time_s } else { 0.001 }
    if (Test-Path (Join-Path $folder ($job + '.odb'))) { Write-Host "SKIP existing $caseName"; continue }
    $port = $BasePort + $index; $index++
    if (Get-NetTCPConnection -LocalPort $port -State Listen -ErrorAction SilentlyContinue) { throw "Port $port busy" }
    $env:MAGPY_SOCKET_PORT = "$port"
    $server = if ($case.field_mode -eq 'LOCAL_TANGENT_CONE') { $localServer } else { $productionServer }
    $telemetry = Join-Path $folder ($job + '_telemetry.csv')
    $stdout = Join-Path $folder ($job + '_socket_stdout.log')
    $stderr = Join-Path $folder ($job + '_socket_stderr.log')
    $arguments = @('-u',$server,'--host','127.0.0.1','--port',"$port",'--job-name',$job,
        '--dxf','J:\magpy\curvenew_CEL_xyrot56_exact.dxf','--drive-type','analytic',
        '--robot-diameter-mm','0.815','--robot-height-mm','2.4','--robot-br-t','1.46',
        '--robot-moment-Am2',"$($case.magnetic_moment_Am2)",'--robot-mass-mg','9.207793514166587',
        '--analytic-b-t','0.010','--analytic-follow-robot','--analytic-gradient-b-t','0',
        '--analytic-gradient-length-mm','45','--analytic-gradient-profile','legacy',
        '--driver-speed-mm-s','6','--bend-speed-mm-s','4.5','--bend-start-mm','13.49','--bend-end-mm','18.56',
        '--z-offset-mm','90','--driver-start-offset-mm','18.899960626','--spin-hz','30',
        '--cone-half-angle-deg',"$($case.cone_deg)",'--driver-xy-shift-mm','2','-6',
        '--cone-axis-bias-deg',"$($case.cone_axis_bias_deg)",'--cone-axis','tangent',
        '--robot-polarity','-1','--robot-axis-tangent','--phase-deg','248','--analytic-rotation-sense','1',
        '--ramp-time-s',"$rampTime",'--endpoint-taper-mm','1','--adaptive-lead-mm','0',
        '--frame-transform-json',$transform,'--telemetry-csv',$telemetry,
        '--telemetry-interval-s','0.00001','--frame-log-interval-s','0.0001')
    $process = Start-Process -FilePath $python -ArgumentList $arguments -WorkingDirectory 'J:\magpy' -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
    $started = Get-Date
    try {
        $deadline = (Get-Date).AddSeconds(30)
        do { Start-Sleep -Milliseconds 200; $listen = Get-NetTCPConnection -LocalPort $port -State Listen -ErrorAction SilentlyContinue } while (-not $listen -and (Get-Date) -lt $deadline)
        if (-not $listen) { throw "Socket failed for $caseName" }
        Set-Location -LiteralPath $folder
        & 'I:\SIMULIA\Commands\abq2025.bat' job=$job input="$job.inp" user='vuforc_audit.f' double=both cpus=1 ask_delete=OFF interactive
        if ($LASTEXITCODE -ne 0) { throw "Abaqus exit $LASTEXITCODE for $caseName" }
        & 'I:\SIMULIA\Commands\abq2025.bat' python $extract $job (Join-Path $folder 'private')
        if ($LASTEXITCODE -ne 0) { throw "ODB extraction failed for $caseName" }
        $case.status = 'SOLVED'; $case | Add-Member -Force NoteProperty wallclock_s ((Get-Date)-$started).TotalSeconds
        $case | ConvertTo-Json -Depth 8 | Set-Content $identityPath
    } finally {
        if ($process -and -not $process.HasExited) { try { $process.Kill() } catch {} }
        Start-Sleep -Milliseconds 300
        if (Get-NetTCPConnection -LocalPort $port -State Listen -ErrorAction SilentlyContinue) { throw "Socket remained open after $caseName" }
    }
    Write-Host ("COMPLETE {0} wallclock={1:N1}s" -f $caseName,((Get-Date)-$started).TotalSeconds)
}
