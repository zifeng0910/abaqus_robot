$ErrorActionPreference = 'Stop'
$work = 'J:\abaqusfangzhen'
$job = 'Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L2300_D0815_WallSupported_0083_Datacheck'
$input = 'Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L2300_D0815_WallSupported_0083.inp'
$bridge = 'J:\abaqusfangzhen\abaqus_robot\calibration_analysis\Wobble30Hz_reduced_hydrodynamics\vuforc_socket_bridge_reduced_hydro.f'
Set-Location -LiteralPath $work
. '.\launch_abaqus_with_env.ps1' -SkipVerify
if (Test-Path -LiteralPath "$work\$job.odb") { throw 'Existing datacheck ODB; refusing repeat' }
& 'I:\SIMULIA\Commands\abq2025.bat' job=$job input=$input user=$bridge datacheck double=both cpus=1 ask_delete=OFF interactive
if ($LASTEXITCODE -ne 0) { throw "Datacheck failed: $LASTEXITCODE" }
if (Select-String -LiteralPath "$work\$job.dat" -Pattern '\*\*\*ERROR') { throw 'Datacheck ERROR found' }
Write-Host 'DATACHECK PASS: frozen L2300 accepted rigid mesh'
