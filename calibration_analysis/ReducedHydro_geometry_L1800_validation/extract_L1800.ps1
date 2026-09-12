$ErrorActionPreference='Stop'
Set-Location -LiteralPath 'J:\abaqusfangzhen'
$job='Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083'
$dest='J:\abaqusfangzhen\abaqus_robot\calibration_analysis\ReducedHydro_geometry_L1800_validation\candidate_private'
$extract='J:\abaqusfangzhen\abaqus_robot\calibration_analysis\ReducedHydro_hidden_impact_audit\extract_odb_audit.py'
if(-not(Test-Path "$job.odb")){throw "Missing ODB $job.odb"}
& 'I:\SIMULIA\Commands\abq2025.bat' python $extract $job $dest
if($LASTEXITCODE -ne 0){throw "ODB extraction failed: $LASTEXITCODE"}
