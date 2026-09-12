$ErrorActionPreference='Stop'
Set-Location -LiteralPath 'J:\abaqusfangzhen'
. '.\launch_abaqus_with_env.ps1' -SkipVerify
$job='Wobble_F30_G6L45_ReducedHydroFixed_ContactAudit_0012'
if(Test-Path "$job.sta"){throw 'Existing solver STA; inspect instead of repeating'}
& 'I:\SIMULIA\Commands\abq2025.bat' job=$job input="$job.inp" user='J:\abaqusfangzhen\abaqus_robot\calibration_analysis\Wobble30Hz_reduced_hydrodynamics\vuforc_socket_bridge_reduced_hydro.f' datacheck double=both cpus=1 ask_delete=OFF interactive
if($LASTEXITCODE -ne 0){throw "Datacheck failed: $LASTEXITCODE"}
if(Select-String -LiteralPath "$job.dat" -Pattern '\*\*\*ERROR') {throw 'Datacheck ERROR found'}
Write-Host 'DATACHECK PASS: no ***ERROR entries'
