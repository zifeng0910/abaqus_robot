[CmdletBinding()]
param(
    [Parameter(Mandatory=$true)][string]$Job,
    [Parameter(Mandatory=$true)][int]$Port
)
$ErrorActionPreference = 'Stop'
$repo = 'J:\abaqusfangzhen\abaqus_robot'
$out = Join-Path $repo 'calibration_analysis\FastStraightDynamicScreen'
$case = Join-Path $out ('cases\' + $Job)
$server = Join-Path $out 'production\magpylib_socket_server_fast.py'
$transform = 'J:\abaqusfangzhen\abaqus_magpylib_frame_transform.json'
$python = 'I:\SIMULIA\EstProducts\2025\win_b64\tools\SMApy\python3.10\python.exe'
$extract = Join-Path $repo 'calibration_analysis\ProductionLocalFrameValidation\scripts\extract_dynamic_odb.py'
$identityPath = Join-Path $case 'case_identity.json'
$identity = Get-Content -Raw $identityPath | ConvertFrom-Json
function Get-Sha256([string]$Path) {
    $stream = [System.IO.File]::OpenRead($Path)
    try {
        $sha = [System.Security.Cryptography.SHA256]::Create()
        try { return ([BitConverter]::ToString($sha.ComputeHash($stream))).Replace('-','').ToLower() }
        finally { $sha.Dispose() }
    } finally { $stream.Dispose() }
}
$serverHash = Get-Sha256 $server
if ($serverHash -ne $identity.fast_server_sha256) { throw 'Fast server identity mismatch' }
if ($identity.status -ne 'PREPARED') { throw "Case status is $($identity.status), expected PREPARED" }
if (Test-Path (Join-Path $case ($Job + '.odb'))) { throw 'Existing ODB; refusing a second dynamic run' }
if (Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue) { throw "Port $Port busy" }

$env:MAGPY_SOCKET_PORT = "$Port"
$telemetry = Join-Path $case ($Job + '_telemetry.csv')
$stdout = Join-Path $case ($Job + '_socket_stdout.log')
$stderr = Join-Path $case ($Job + '_socket_stderr.log')
$dxf = Join-Path $case 'straight_control_centerline.dxf'
$arguments = @('-u',$server,'--host','127.0.0.1','--port',"$Port",'--job-name',$Job,
    '--dxf',$dxf,'--drive-type','analytic','--field-frame-mode',$identity.field_frame_mode,
    '--rocking-amplitude-deg',([string]$identity.rocking_main_amplitude_deg),
    '--rocking-cross-amplitude-deg',([string]$identity.rocking_cross_amplitude_deg),
    '--rocking-frame-azimuth-deg','-61.37284757596327','--robot-diameter-mm','0.815','--robot-height-mm','2.4',
    '--robot-br-t','1.46','--robot-moment-Am2','0.0010876227522174417','--robot-mass-mg','9.207793514166587',
    '--analytic-b-t','0.010','--analytic-follow-robot','--analytic-gradient-b-t','0.006',
    '--analytic-gradient-length-mm','45','--analytic-gradient-profile','legacy',
    '--driver-speed-mm-s','6','--bend-speed-mm-s','4.5','--bend-start-mm','13.49','--bend-end-mm','18.56',
    '--z-offset-mm','90','--driver-start-offset-mm','18.899960626','--spin-hz',([string]$identity.frequency_Hz),
    '--cone-half-angle-deg','30','--driver-xy-shift-mm','2','-6','--cone-axis-bias-deg','0',
    '--cone-axis','tangent','--robot-polarity','-1','--robot-moment-axis-aba','0.9762799602464296','-0.0618320247446977','0.2074951564186524','--phase-deg','0',
    '--analytic-rotation-sense','1','--ramp-time-s','0.001','--endpoint-taper-mm','1',
    '--adaptive-lead-mm','0','--frame-transform-json',$transform,'--telemetry-csv',$telemetry,
    '--telemetry-interval-s','0.00002','--frame-log-interval-s','0.0005')
$process = Start-Process -FilePath $python -ArgumentList $arguments -WorkingDirectory $case -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
$started = Get-Date
try {
    $deadline = (Get-Date).AddSeconds(30)
    do {
        Start-Sleep -Milliseconds 200
        $listen = Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue
    } while (-not $listen -and (Get-Date) -lt $deadline)
    if (-not $listen) { throw "Socket startup failed; see $stderr" }
    . 'J:\abaqusfangzhen\launch_abaqus_with_env.ps1' -SkipVerify
    Set-Location -LiteralPath $case
    & 'I:\SIMULIA\Commands\abq2025.bat' job=$Job input="$Job.inp" user='vuforc_production_local.f' double=both cpus=1 ask_delete=OFF interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus exit $LASTEXITCODE" }
    & 'I:\SIMULIA\Commands\abq2025.bat' python $extract $Job (Join-Path $case 'private')
    if ($LASTEXITCODE -ne 0) { throw 'ODB extraction failed' }
    $identity.status = 'SOLVED'
    $identity | Add-Member -NotePropertyName wallclock_s -NotePropertyValue (((Get-Date) - $started).TotalSeconds) -Force
    $identity | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath $identityPath
} finally {
    if ($process -and -not $process.HasExited) { try { $process.Kill() } catch {} }
}
Write-Host ("COMPLETE {0} wallclock={1:N1}s" -f $Job,((Get-Date)-$started).TotalSeconds)
