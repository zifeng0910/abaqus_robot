$ErrorActionPreference='Stop';Set-Location -LiteralPath 'J:\abaqusfangzhen'
$job='Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_WallOn_Free_0083';$dest='J:\abaqusfangzhen\abaqus_robot\calibration_analysis\ReducedHydro_FreeCAD_L1800_validation\candidate_private'
& 'I:\SIMULIA\Commands\abq2025.bat' python 'J:\abaqusfangzhen\abaqus_robot\calibration_analysis\ReducedHydro_hidden_impact_audit\extract_odb_audit.py' $job $dest
if($LASTEXITCODE -ne 0){throw "ODB extraction failed: $LASTEXITCODE"}
