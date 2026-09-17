[CmdletBinding()]
param()
$ErrorActionPreference = 'Stop'
$repo = 'J:\abaqusfangzhen\abaqus_robot'
$out = Join-Path $repo 'calibration_analysis\FastPrecomputedForwardFlowScreen'
$job = 'FAST_PRECOMP_F120_G0P30_FLOW_1CYCLE'
$case = Join-Path $out ('case\' + $job)
$identityPath = Join-Path $case 'case_identity.json'
$identity = Get-Content -Raw $identityPath | ConvertFrom-Json
if ($identity.status -ne 'PREPARED') { throw "Expected PREPARED, got $($identity.status)" }
if (-not $identity.forward_sign_sanity.passed) { throw 'Forward sign sanity did not pass' }
if (Test-Path (Join-Path $case ($job + '.lck'))) { throw 'Existing solver lock' }
. (Join-Path (Split-Path $repo -Parent) 'launch_abaqus_with_env.ps1') -SkipVerify
$ErrorActionPreference = 'Continue'
Set-Location -LiteralPath $case
$sub = Join-Path $case 'vuamp_precomputed_g0p30_flow.f90'
$gateJob = $job + '_DATACHECK'
$gateOutput = & 'I:\SIMULIA\Commands\abq2025.bat' job=$gateJob input="$job.inp" user=$sub double=both datacheck ask_delete=ON interactive 2>&1
$gateExit = $LASTEXITCODE
$gateOutput | Out-File (Join-Path $out ($gateJob + '.log')) -Encoding ascii
$gatePassed = $gateExit -eq 0 -and ($gateOutput | Out-String) -notmatch '(?i)(Abaqus Error:|Problem during compilation|Analysis exited with errors)'
[ordered]@{job=$gateJob; exit_code=$gateExit; passed=$gatePassed} | ConvertTo-Json | Set-Content (Join-Path $out 'datacheck_gate.json') -Encoding ascii
if (-not $gatePassed) { throw 'Datacheck failed; dynamics not launched' }
$started = Get-Date
$solverOutput = & 'I:\SIMULIA\Commands\abq2025.bat' job=$job input="$job.inp" user=$sub double=both cpus=1 ask_delete=ON interactive 2>&1
$solverExit = $LASTEXITCODE
$solverOutput | Out-File (Join-Path $out ($job + '_solver.log')) -Encoding ascii
if ($solverExit -ne 0) { throw "Dynamics failed: $solverExit" }
$wallclock = ((Get-Date)-$started).TotalSeconds
$extract = Join-Path $repo 'calibration_analysis\FastPrecomputedNoReboundScreen\scripts\extract_fast_odb.py'
$extractOutput = & 'I:\SIMULIA\Commands\abq2025.bat' python $extract $job (Join-Path $case 'private') 2>&1
$extractExit = $LASTEXITCODE
$extractOutput | Out-File (Join-Path $out ($job + '_extract.log')) -Encoding ascii
if ($extractExit -ne 0) { throw 'ODB extraction failed' }
$identity.status = 'SOLVED'
$identity.dynamics_run_count = 1
$identity | Add-Member wallclock_s $wallclock -Force
$identity | Add-Member cpus 1 -Force
$identity | Add-Member socket_calls 0 -Force
$identity | ConvertTo-Json -Depth 12 | Set-Content $identityPath -Encoding ascii
Write-Host ("COMPLETE {0}; wallclock={1:N1}s" -f $job,$wallclock)
