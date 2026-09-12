$ErrorActionPreference='Stop'
Set-Location -LiteralPath 'J:\abaqusfangzhen'
. '.\launch_abaqus_with_env.ps1' -SkipVerify
$job='Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083_Datacheck'
$input='Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083.inp'
if(Test-Path "$job.odb"){throw 'Existing datacheck ODB; refusing repeat'}
& 'I:\SIMULIA\Commands\abq2025.bat' job=$job input=$input user='J:\abaqusfangzhen\abaqus_robot\calibration_analysis\Wobble30Hz_reduced_hydrodynamics\vuforc_socket_bridge_reduced_hydro.f' datacheck double=both cpus=1 ask_delete=OFF interactive
if($LASTEXITCODE -ne 0){throw "Datacheck failed: $LASTEXITCODE"}
if(Select-String -LiteralPath "$job.dat" -Pattern '\*\*\*ERROR'){throw 'Datacheck ERROR found'}
Write-Host 'DATACHECK PASS: L=1.800 mm geometry candidate accepted by Abaqus 2025'
