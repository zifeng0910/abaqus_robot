[CmdletBinding()]
param([string]$WorkDir = 'J:\abaqusfangzhen', [int]$Port = 65520)
$ErrorActionPreference = 'Stop'
$Job = 'Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L2300_D0815_WallSupported_0083'
$server = 'J:\magpy\magpylib_socket_server.py'
$bridge = 'J:\abaqusfangzhen\abaqus_robot\calibration_analysis\Wobble30Hz_reduced_hydrodynamics\vuforc_socket_bridge_reduced_hydro.f'
$transform = 'J:\abaqusfangzhen\abaqus_magpylib_frame_transform.json'
$here = 'J:\abaqusfangzhen\abaqus_robot\calibration_analysis\ReducedHydro_L2300_frozenCAD_validation'
$identity = Get-Content -Raw "$here\L2300_deck_identity.json" | ConvertFrom-Json
$datacheck = Import-Csv "$here\L2300_datacheck_identity.csv"
if (-not $identity.ready_for_datacheck) { throw 'Preflight gates not passed' }
if ([int](($datacheck | Where-Object metric -eq 'errors').value) -ne 0) { throw 'Datacheck identity not passed' }
Set-Location -LiteralPath $WorkDir
$moment = [string]::Format([Globalization.CultureInfo]::InvariantCulture, '{0:R}', [double]$identity.magnetic_moment_Am2)
$mass = [string]::Format([Globalization.CultureInfo]::InvariantCulture, '{0:R}', [double]$identity.mesh_mass_mg)
$telemetry = Join-Path $WorkDir ($Job + '_telemetry.csv')
$stdout = Join-Path $WorkDir ($Job + '_socket_stdout.log')
$stderr = Join-Path $WorkDir ($Job + '_socket_stderr.log')
$py = 'I:\SIMULIA\EstProducts\2025\win_b64\tools\SMApy\python3.10\python.exe'
if (Test-Path -LiteralPath "$WorkDir\$Job.odb") { throw 'Existing ODB; refusing repeat run' }
if (-not (Test-Path -LiteralPath "$WorkDir\$Job.inp")) { throw "Missing input deck $Job.inp" }
if (Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue) { throw "Port $Port already listening" }
$env:MAGPY_SOCKET_PORT = "$Port"
$args = @('-u', $server, '--host', '127.0.0.1', '--port', "$Port", '--job-name', $Job,
    '--dxf', 'J:\magpy\curvenew_CEL_xyrot56_exact.dxf', '--drive-type', 'analytic',
    '--robot-diameter-mm', '0.815', '--robot-height-mm', '2.3', '--robot-br-t', '1.46',
    '--robot-moment-Am2', $moment, '--robot-mass-mg', $mass,
    '--analytic-b-t', '0.010', '--analytic-follow-robot', '--analytic-gradient-b-t', '0.006',
    '--analytic-gradient-length-mm', '45', '--analytic-gradient-profile', 'legacy',
    '--analytic-gradient-switch-start-s', '0', '--analytic-gradient-transition-s', '0.0005',
    '--driver-speed-mm-s', '6', '--bend-speed-mm-s', '4.5', '--bend-start-mm', '13.49',
    '--bend-end-mm', '18.56', '--z-offset-mm', '90', '--driver-start-offset-mm', '18.899960626',
    '--spin-hz', '30', '--cone-half-angle-deg', '30', '--driver-xy-shift-mm', '2', '-6',
    '--cone-axis-bias-deg', '40', '--cone-axis', 'tangent', '--robot-polarity', '-1',
    '--robot-axis-tangent', '--phase-deg', '248', '--analytic-rotation-sense', '1',
    '--ramp-time-s', '0.001', '--endpoint-taper-mm', '1', '--adaptive-lead-mm', '0',
    '--frame-transform-json', $transform, '--telemetry-csv', $telemetry,
    '--telemetry-interval-s', '0.0', '--frame-log-interval-s', '0.00005')
$process = Start-Process -FilePath $py -ArgumentList $args -WorkingDirectory 'J:\magpy' -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
try {
    $deadline = (Get-Date).AddSeconds(30)
    do { Start-Sleep -Milliseconds 200; $listen = Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue } while (-not $listen -and (Get-Date) -lt $deadline)
    if (-not $listen) { throw "Socket failed; see $stderr" }
    . "$WorkDir\launch_abaqus_with_env.ps1" -SkipVerify
    & 'I:\SIMULIA\Commands\abq2025.bat' job=$Job input="$Job.inp" user="$bridge" double=both cpus=1 ask_delete=OFF interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus exit $LASTEXITCODE" }
} finally {
    if ($process -and -not $process.HasExited) { try { $process.Kill() } catch {} }
}
Copy-Item -LiteralPath 'J:\abaqusfangzhen\reduced_hydro_runtime_load.csv' -Destination (Join-Path $WorkDir ($Job + '_hydro_increment.csv'))
