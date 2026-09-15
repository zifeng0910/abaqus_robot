[CmdletBinding()]
param()
$ErrorActionPreference = 'Stop'
$repo = 'J:\abaqusfangzhen\abaqus_robot'
$job = 'STRAIGHT_PRESCRIBED_WOBBLE_CONTACT_AUDIT'
$case = Join-Path $repo ('calibration_analysis\StraightBridgeTorqueAudit\case\' + $job)
$identityPath = Join-Path $case 'case_identity.json'
$identity = Get-Content -Raw $identityPath | ConvertFrom-Json
if ($identity.status -ne 'PREPARED') { throw "Case status is $($identity.status), expected PREPARED" }
if (Test-Path (Join-Path $case ($job + '.odb'))) { throw 'Existing ODB; refusing a second run' }
. 'J:\abaqusfangzhen\launch_abaqus_with_env.ps1' -SkipVerify
$started = Get-Date
Set-Location -LiteralPath $case
& 'I:\SIMULIA\Commands\abq2025.bat' job=$job input="$job.inp" user='vuforc_prescribed_hydro.f' double=both cpus=1 ask_delete=OFF interactive
if ($LASTEXITCODE -ne 0) { throw "Abaqus exit $LASTEXITCODE" }
$extract = Join-Path $repo 'calibration_analysis\ProductionLocalFrameValidation\scripts\extract_dynamic_odb.py'
& 'I:\SIMULIA\Commands\abq2025.bat' python $extract $job (Join-Path $case 'private')
if ($LASTEXITCODE -ne 0) { throw 'ODB extraction failed' }
$identity.status = 'SOLVED'
$identity | Add-Member -NotePropertyName wallclock_s -NotePropertyValue (((Get-Date) - $started).TotalSeconds) -Force
$identity | ConvertTo-Json -Depth 10 | Set-Content $identityPath
Write-Host ("COMPLETE {0} wallclock={1:N1}s" -f $job,((Get-Date)-$started).TotalSeconds)
