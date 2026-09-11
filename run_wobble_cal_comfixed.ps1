[CmdletBinding()]
param(
  [string]$WorkDir = 'J:\abaqusfangzhen',
  [string]$JobName = 'WobbleCal_COMFixed_F40_Cone30_B10',
  [int]$Port = 65491
)
$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $WorkDir
$server = 'J:\magpy\magpylib_socket_server.py'
$bridge = Join-Path $WorkDir 'vuforc_socket_bridge_comfixed.f'
$transform = Join-Path $WorkDir 'abaqus_magpylib_frame_transform.json'
$telemetry = Join-Path $WorkDir ($JobName + '_telemetry.csv')
$stdout = Join-Path $WorkDir ($JobName + '_socket_stdout.log')
$stderr = Join-Path $WorkDir ($JobName + '_socket_stderr.log')
$env:MAGPY_SOCKET_PORT = "$Port"

if (-not (Test-Path -LiteralPath "$JobName.inp")) { throw "Missing input deck: $JobName.inp" }
if (Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue) { throw "Port $Port is already listening" }

# Keep load updates at every Explicit increment. Only file diagnostics are throttled.
$serverArgs = @(
  $server, '--host','127.0.0.1','--port',"$Port", '--job-name',$JobName,
  '--dxf','J:\magpy\curvenew_CEL_xyrot56_exact.dxf', '--drive-type','analytic',
  '--robot-diameter-mm','1.22','--robot-height-mm','2.81',
  '--robot-br-t','1.46','--robot-moment-Am2','0.001168', '--robot-mass-mg','10',
  '--analytic-b-t','0.010','--analytic-follow-robot',
  '--analytic-gradient-b-t','0.003','--analytic-gradient-length-mm','45',
  '--driver-speed-mm-s','6','--bend-speed-mm-s','4.5', '--bend-start-mm','13.49', '--bend-end-mm','18.56',
  '--z-offset-mm','90','--driver-start-offset-mm','18.899960626', '--spin-hz','40',
  '--cone-half-angle-deg','30','--driver-xy-shift-mm','2','-6',
  '--cone-axis-bias-deg','40','--cone-axis','tangent','--robot-polarity','-1',
  '--robot-axis-tangent','--phase-deg','248','--analytic-rotation-sense','1',
  '--ramp-time-s','0.001','--endpoint-taper-mm','1','--adaptive-lead-mm','0',
  '--frame-transform-json',$transform, '--telemetry-csv',$telemetry,
  '--telemetry-interval-s','0.0001','--frame-log-interval-s','0.0001'
)
$py = 'I:\SIMULIA\EstProducts\2025\win_b64\tools\SMApy\python3.10\python.exe'
if (-not (Test-Path -LiteralPath $py)) { $py = 'python' }
Write-Host "Starting Magpylib server on port $Port (F=40 Hz; B=10 mT; gradient=3 mT=0.003 T; COM-fixed calibration)"
$p = Start-Process -FilePath $py -ArgumentList $serverArgs -WorkingDirectory 'J:\magpy' -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
try {
  $deadline = (Get-Date).AddSeconds(30)
  do {
    Start-Sleep -Milliseconds 250
    $listen = Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue
  } while (-not $listen -and (Get-Date) -lt $deadline)
  if (-not $listen) { throw "Socket did not enter LISTENING state. See $stderr" }
  Write-Host "Socket $Port LISTENING; launching Abaqus job $JobName"
  . "$WorkDir\launch_abaqus_with_env.ps1" -SkipVerify
  & 'I:\SIMULIA\Commands\abq2025.bat' job=$JobName input="$JobName.inp" user="$bridge" double=both cpus=1 ask_delete=OFF interactive
  if ($LASTEXITCODE -ne 0) { throw "Abaqus exited with code $LASTEXITCODE" }
}
finally {
  if ($p -and -not $p.HasExited) {
    try { $p.CloseMainWindow() | Out-Null } catch {}
    try { $p.Kill() } catch {}
  }
}
Write-Host "Calibration job finished. Telemetry: $telemetry"
