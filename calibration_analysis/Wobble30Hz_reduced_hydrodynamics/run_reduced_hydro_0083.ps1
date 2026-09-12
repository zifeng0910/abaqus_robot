[CmdletBinding()]
param([string]$WorkDir='J:\abaqusfangzhen',[int]$Port=65504)
$ErrorActionPreference='Stop'; Set-Location -LiteralPath $WorkDir
$Job='Wobble_F30_G6L45_ReducedHydroFixed_WallOn_Free_0083'
$server='J:\magpy\magpylib_socket_server.py'
$bridge='J:\abaqusfangzhen\abaqus_robot\calibration_analysis\Wobble30Hz_reduced_hydrodynamics\vuforc_socket_bridge_reduced_hydro.f'
$transform='J:\abaqusfangzhen\abaqus_magpylib_frame_transform.json'
$telemetry=Join-Path $WorkDir ($Job+'_telemetry.csv')
$stdout=Join-Path $WorkDir ($Job+'_socket_stdout.log'); $stderr=Join-Path $WorkDir ($Job+'_socket_stderr.log')
$py='I:\SIMULIA\EstProducts\2025\win_b64\tools\SMApy\python3.10\python.exe'; if(-not(Test-Path $py)){$py='python'}
if(-not(Test-Path "$Job.inp")){throw "Missing input deck $Job.inp"}
if(Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue){throw "Port $Port already listening"}
$env:MAGPY_SOCKET_PORT="$Port"
$args=@('-u',$server,'--host','127.0.0.1','--port',"$Port",'--job-name',$Job,'--dxf','J:\magpy\curvenew_CEL_xyrot56_exact.dxf','--drive-type','analytic','--robot-diameter-mm','1.22','--robot-height-mm','2.81','--robot-br-t','1.46','--robot-moment-Am2','0.001168','--robot-mass-mg','10','--analytic-b-t','0.010','--analytic-follow-robot','--analytic-gradient-b-t','0.006','--analytic-gradient-length-mm','45','--analytic-gradient-profile','legacy','--analytic-gradient-switch-start-s','0','--analytic-gradient-transition-s','0.0005','--driver-speed-mm-s','6','--bend-speed-mm-s','4.5','--bend-start-mm','13.49','--bend-end-mm','18.56','--z-offset-mm','90','--driver-start-offset-mm','18.899960626','--spin-hz','30','--cone-half-angle-deg','30','--driver-xy-shift-mm','2','-6','--cone-axis-bias-deg','40','--cone-axis','tangent','--robot-polarity','-1','--robot-axis-tangent','--phase-deg','248','--analytic-rotation-sense','1','--ramp-time-s','0.001','--endpoint-taper-mm','1','--adaptive-lead-mm','0','--frame-transform-json',$transform,'--telemetry-csv',$telemetry,'--telemetry-interval-s','0.00005','--frame-log-interval-s','0.00005')
Write-Host 'REDUCED-HYDRO: CEL fluid removed; SmoothWall114 single-wall contact; dissipative RP hydro.'
Write-Host 'FIELD: B0=10mT, f=30Hz, cone=30deg, Bias=40deg, G=6mT=0.006T, L=45mm, phase=248deg, sense=+1'
Write-Host "TRAJECTORY_DIRECTION=FORWARD; PORT=$Port; DIRECT_DT=0.1us; duration=8.333ms"
$p=Start-Process -FilePath $py -ArgumentList $args -WorkingDirectory 'J:\magpy' -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
try{$deadline=(Get-Date).AddSeconds(30); do{Start-Sleep -Milliseconds 200;$listen=Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue}while(-not$listen -and (Get-Date)-lt$deadline);if(-not$listen){throw "Socket failed; see $stderr"};. "$WorkDir\launch_abaqus_with_env.ps1" -SkipVerify;& 'I:\SIMULIA\Commands\abq2025.bat' job=$Job input="$Job.inp" user="$bridge" double=both cpus=1 ask_delete=OFF interactive;if($LASTEXITCODE -ne 0){throw "Abaqus exit $LASTEXITCODE"}}finally{if($p -and -not$p.HasExited){try{$p.Kill()}catch{}}}
Write-Host "Completed $Job"
