[CmdletBinding()]
param()
$ErrorActionPreference = 'Stop'
$repo = 'J:\abaqusfangzhen\abaqus_robot'
$out = Join-Path $repo 'calibration_analysis\FastPrecomputedNoReboundScreen'
$job = 'FAST_PRECOMP_F120_DUALEND_NOREBOUND'
$case = Join-Path $out ('case\' + $job)
$identityPath = Join-Path $case 'case_identity.json'
$identity = Get-Content -Raw $identityPath | ConvertFrom-Json
if ($identity.status -ne 'PREPARED') { throw "Case status is $($identity.status), expected PREPARED" }
if (Test-Path (Join-Path $case ($job + '.lck'))) { throw 'Existing solver lock; refusing concurrent run' }

$envScript = Join-Path (Split-Path $repo -Parent) 'launch_abaqus_with_env.ps1'
if (-not (Test-Path $envScript)) { throw "Missing Abaqus environment script: $envScript" }
. $envScript -SkipVerify
$ErrorActionPreference = 'Continue'
$userSub = Join-Path $case 'vuamp_precomputed_fast_surrogate.f90'
Set-Location -LiteralPath $case

$gateJob = $job + '_DATACHECK'
$gateStarted = Get-Date
$gateOutput = & 'I:\SIMULIA\Commands\abq2025.bat' job=$gateJob input="$job.inp" user=$userSub double=both datacheck ask_delete=ON interactive 2>&1
$gateExit = $LASTEXITCODE
$gateOutput | Out-File -LiteralPath (Join-Path $out ($gateJob + '.log')) -Encoding ascii
$gateText = $gateOutput | Out-String
$gatePassed = $gateExit -eq 0 -and $gateText -notmatch '(?i)(Abaqus Error:|Problem during compilation|Analysis exited with errors)'
[ordered]@{job=$gateJob; exit_code=$gateExit; passed=$gatePassed; wallclock_s=((Get-Date)-$gateStarted).TotalSeconds} |
    ConvertTo-Json | Set-Content -LiteralPath (Join-Path $out 'datacheck_gate.json') -Encoding ascii
if (-not $gatePassed) { throw 'Datacheck failed; dynamics was not launched' }

$started = Get-Date
$solverOutput = & 'I:\SIMULIA\Commands\abq2025.bat' job=$job input="$job.inp" user=$userSub double=both cpus=1 ask_delete=ON interactive 2>&1
$solverExit = $LASTEXITCODE
$solverOutput | Out-File -LiteralPath (Join-Path $out ($job + '_solver.log')) -Encoding ascii
if ($solverExit -ne 0) { throw "Abaqus dynamics failed with exit code $solverExit" }
$wallclock = ((Get-Date)-$started).TotalSeconds

$private = Join-Path $case 'private'
$extract = Join-Path $out 'scripts\extract_fast_odb.py'
$extractOutput = & 'I:\SIMULIA\Commands\abq2025.bat' python $extract $job $private 2>&1
$extractExit = $LASTEXITCODE
$extractOutput | Out-File -LiteralPath (Join-Path $out ($job + '_extract.log')) -Encoding ascii
if ($extractExit -ne 0) { throw 'ODB extraction failed' }

$lookupLine = Select-String -Path (Join-Path $case ($job + '.msg')),(Join-Path $case ($job + '.sta')) -Pattern 'PRECOMPUTED_TABLE_LOOKUPS' -ErrorAction SilentlyContinue | Select-Object -Last 1
$identity.status = 'SOLVED'
$identity.dynamics_run_count = 1
$identity | Add-Member -NotePropertyName wallclock_s -NotePropertyValue $wallclock -Force
$identity | Add-Member -NotePropertyName cpus -NotePropertyValue 1 -Force
$identity | Add-Member -NotePropertyName socket_calls -NotePropertyValue 0 -Force
$identity | Add-Member -NotePropertyName lookup_log_line -NotePropertyValue ([string]$lookupLine.Line) -Force
$identity | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $identityPath -Encoding ascii
Write-Host ("COMPLETE {0}; wallclock={1:N1}s" -f $job,$wallclock)
