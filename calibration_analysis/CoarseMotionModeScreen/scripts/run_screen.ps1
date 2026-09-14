[CmdletBinding()]
param(
    [string[]]$CaseId = @('GEO_230'),
    [int]$BasePort = 64000
)
$ErrorActionPreference = 'Stop'
$repo = 'J:\abaqusfangzhen\abaqus_robot'
$screen = Join-Path $repo 'calibration_analysis\CoarseMotionModeScreen'
$server = 'J:\magpy\magpylib_socket_server.py'
$transform = 'J:\abaqusfangzhen\abaqus_magpylib_frame_transform.json'
$py = 'I:\SIMULIA\EstProducts\2025\win_b64\tools\SMApy\python3.10\python.exe'
$extract = Join-Path $screen 'scripts\extract_case_odb.py'
$manifest = Import-Csv (Join-Path $screen 'metrics\case_manifest.csv')
$aliases = @{}
Import-Csv (Join-Path $screen 'metrics\case_aliases.csv') | ForEach-Object { $aliases[$_.requested_case_id] = $_.executed_as }
. 'J:\abaqusfangzhen\launch_abaqus_with_env.ps1' -SkipVerify
$index = 0
foreach ($requested in $CaseId) {
    $resolved = if ($aliases.ContainsKey($requested)) { $aliases[$requested] } else { $requested }
    $case = $manifest | Where-Object case_id -eq $resolved
    if (-not $case) { throw "Unknown case $requested" }
    $folder = Join-Path $screen ('cases\' + $resolved)
    $job = $case.job_name
    $port = $BasePort + $index
    $index++
    if (Test-Path (Join-Path $folder ($job + '.odb'))) {
        Write-Host "SKIP existing $resolved"
        continue
    }
    if (Get-NetTCPConnection -LocalPort $port -State Listen -ErrorAction SilentlyContinue) { throw "Port $port busy" }
    $env:MAGPY_SOCKET_PORT = "$port"
    $telemetry = Join-Path $folder ($job + '_telemetry.csv')
    $stdout = Join-Path $folder ($job + '_socket_stdout.log')
    $stderr = Join-Path $folder ($job + '_socket_stderr.log')
    $args = @('-u',$server,'--host','127.0.0.1','--port',"$port",'--job-name',$job,
        '--dxf','J:\magpy\curvenew_CEL_xyrot56_exact.dxf','--drive-type','analytic',
        '--robot-diameter-mm','0.815','--robot-height-mm',"$($case.length_mm)",'--robot-br-t','1.46',
        '--robot-moment-Am2',"$($case.magnetic_moment_Am2)",'--robot-mass-mg',"$($case.mesh_mass_mg)",
        '--analytic-b-t',"$([double]$case.B0_mT/1000)",'--analytic-follow-robot',
        '--analytic-gradient-b-t',"$([double]$case.gradient_mT/1000)",'--analytic-gradient-length-mm','45',
        '--analytic-gradient-profile','legacy','--analytic-gradient-switch-start-s','0',
        '--analytic-gradient-transition-s','0.0005','--driver-speed-mm-s','6','--bend-speed-mm-s','4.5',
        '--bend-start-mm','13.49','--bend-end-mm','18.56','--z-offset-mm','90',
        '--driver-start-offset-mm','18.899960626','--spin-hz',"$($case.frequency_Hz)",
        '--cone-half-angle-deg',"$($case.cone_deg)",'--driver-xy-shift-mm','2','-6',
        '--cone-axis-bias-deg','40','--cone-axis','tangent','--robot-polarity','-1',
        '--robot-axis-tangent','--phase-deg','248','--analytic-rotation-sense','1',
        '--ramp-time-s','0.001','--endpoint-taper-mm','1','--adaptive-lead-mm','0',
        '--frame-transform-json',$transform,'--telemetry-csv',$telemetry,
        '--telemetry-interval-s','0.00001','--frame-log-interval-s','0.0001')
    $process = Start-Process -FilePath $py -ArgumentList $args -WorkingDirectory 'J:\magpy' -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
    $started = Get-Date
    try {
        $deadline = (Get-Date).AddSeconds(30)
        do { Start-Sleep -Milliseconds 200; $listen = Get-NetTCPConnection -LocalPort $port -State Listen -ErrorAction SilentlyContinue } while (-not $listen -and (Get-Date) -lt $deadline)
        if (-not $listen) { throw "Socket failed for $resolved" }
        Set-Location -LiteralPath $folder
        & 'I:\SIMULIA\Commands\abq2025.bat' job=$job input="$job.inp" user='vuforc_screen.f' double=both cpus=1 ask_delete=OFF interactive
        if ($LASTEXITCODE -ne 0) { throw "Abaqus exit $LASTEXITCODE for $resolved" }
        & 'I:\SIMULIA\Commands\abq2025.bat' python $extract $job (Join-Path $folder 'private')
        if ($LASTEXITCODE -ne 0) { throw "ODB extraction failed for $resolved" }
        $identityPath = Join-Path $folder 'case_identity.json'
        $identity = Get-Content -Raw $identityPath | ConvertFrom-Json
        $identity.status = 'SOLVED'
        $identity | Add-Member -Force NoteProperty wallclock_s ((Get-Date)-$started).TotalSeconds
        $identity | ConvertTo-Json -Depth 8 | Set-Content $identityPath
    } finally {
        if ($process -and -not $process.HasExited) { try { $process.Kill() } catch {} }
        Start-Sleep -Milliseconds 300
        if (Get-NetTCPConnection -LocalPort $port -State Listen -ErrorAction SilentlyContinue) { throw "Socket remained open after $resolved" }
    }
    Write-Host ("COMPLETE {0} wallclock={1:N1}s" -f $resolved,((Get-Date)-$started).TotalSeconds)
}
