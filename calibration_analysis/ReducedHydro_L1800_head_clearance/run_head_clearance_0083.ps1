[CmdletBinding()]
param([string]$WorkDir='J:\abaqusfangzhen',[int]$Port=65523)
$ErrorActionPreference='Stop'
Set-Location -LiteralPath $WorkDir
$Job='Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_HeadClearance_0083'
$Study='J:\abaqusfangzhen\abaqus_robot\calibration_analysis\ReducedHydro_L1800_head_clearance'
$mesh=Get-Content -Raw (Join-Path $Study 'headclearance_mesh060_audit.json')|ConvertFrom-Json
$replay=Get-Content -Raw (Join-Path $Study 'head_clearance_replay_identity.json')|ConvertFrom-Json
$props=Get-Content -Raw (Join-Path $Study 'Robot_parametric_L1p800_D0p815_HeadClearance_mass_properties.json')|ConvertFrom-Json
$base=Get-Content -Raw 'J:\abaqusfangzhen\abaqus_robot\calibration_analysis\ReducedHydro_FreeCAD_L1800_validation\freecad_preflight_identity.json'|ConvertFrom-Json
if(-not $mesh.all_required_geometry_gates_pass){throw 'Simplified candidate mesh gates failed'}
if(-not $replay.all_zero_abaqus_gates_pass){throw 'Zero-Abaqus replay gates failed'}
if(Test-Path "$Job.odb"){throw 'Existing ODB; refusing repeat dynamic run'}
$moment=[double]$base.moment_density_Am2_per_mm3*[double]$props.volume.value
$momentText=[string]::Format([Globalization.CultureInfo]::InvariantCulture,'{0:R}',$moment)
$massText=[string]::Format([Globalization.CultureInfo]::InvariantCulture,'{0:R}',[double]$props.mass.mg)
$server='J:\magpy\magpylib_socket_server.py'
$bridge='J:\abaqusfangzhen\abaqus_robot\calibration_analysis\Wobble30Hz_reduced_hydrodynamics\vuforc_socket_bridge_reduced_hydro.f'
$transform='J:\abaqusfangzhen\abaqus_magpylib_frame_transform.json'
$telemetry=Join-Path $WorkDir ($Job+'_telemetry.csv')
$stdout=Join-Path $WorkDir ($Job+'_socket_stdout.log')
$stderr=Join-Path $WorkDir ($Job+'_socket_stderr.log')
$py='I:\SIMULIA\EstProducts\2025\win_b64\tools\SMApy\python3.10\python.exe'
Copy-Item -LiteralPath (Join-Path $Study ($Job+'.inp')) -Destination (Join-Path $WorkDir ($Job+'.inp'))
if(Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue){throw "Port $Port already listening"}
$env:MAGPY_SOCKET_PORT="$Port"
$arguments=@('-u',$server,'--host','127.0.0.1','--port',"$Port",'--job-name',$Job,
 '--dxf','J:\magpy\curvenew_CEL_xyrot56_exact.dxf','--drive-type','analytic','--robot-diameter-mm','0.815',
 '--robot-height-mm','1.8','--robot-br-t','1.46','--robot-moment-Am2',$momentText,'--robot-mass-mg',$massText,
 '--analytic-b-t','0.010','--analytic-follow-robot','--analytic-gradient-b-t','0.006','--analytic-gradient-length-mm','45',
 '--analytic-gradient-profile','legacy','--analytic-gradient-switch-start-s','0','--analytic-gradient-transition-s','0.0005',
 '--driver-speed-mm-s','6','--bend-speed-mm-s','4.5','--bend-start-mm','13.49','--bend-end-mm','18.56',
 '--z-offset-mm','90','--driver-start-offset-mm','18.899960626','--spin-hz','30','--cone-half-angle-deg','30',
 '--driver-xy-shift-mm','2','-6','--cone-axis-bias-deg','40','--cone-axis','tangent','--robot-polarity','-1',
 '--robot-axis-tangent','--phase-deg','248','--analytic-rotation-sense','1','--ramp-time-s','0.001',
 '--endpoint-taper-mm','1','--adaptive-lead-mm','0','--frame-transform-json',$transform,'--telemetry-csv',$telemetry,
 '--telemetry-interval-s','0.0','--frame-log-interval-s','0.00005')
$serverProcess=Start-Process -FilePath $py -ArgumentList $arguments -WorkingDirectory 'J:\magpy' -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
$started=Get-Date
try{
 $deadline=(Get-Date).AddSeconds(30)
 do{Start-Sleep -Milliseconds 200;$listen=Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue}while(-not $listen -and (Get-Date)-lt $deadline)
 if(-not $listen){throw "Socket failed; see $stderr"}
 . "$WorkDir\launch_abaqus_with_env.ps1" -SkipVerify
 & 'I:\SIMULIA\Commands\abq2025.bat' job=$Job input="$Job.inp" user="$bridge" double=both cpus=1 ask_delete=OFF interactive
 if($LASTEXITCODE -ne 0){throw "Abaqus exit $LASTEXITCODE"}
}finally{if($serverProcess-and-not $serverProcess.HasExited){try{$serverProcess.Kill()}catch{}}}
$elapsed=((Get-Date)-$started).TotalSeconds
Copy-Item -LiteralPath 'J:\abaqusfangzhen\reduced_hydro_runtime_load.csv' -Destination (Join-Path $WorkDir ($Job+'_hydro_increment.csv'))
@{job=$Job;dynamic_wall_clock_s=$elapsed;mesh_elements=$mesh.element_count;mesh_nodes=$mesh.node_count;magnetic_moment_Am2=$moment;mass_mg=$props.mass.mg}|ConvertTo-Json|Set-Content (Join-Path $Study 'headclearance_dynamic_run_identity.json')
