[CmdletBinding()]
param()
$ErrorActionPreference = 'Stop'
$repo = 'J:\abaqusfangzhen\abaqus_robot'
$out = Join-Path $repo 'calibration_analysis\FastVMean5Calibration'
$job = 'FAST_VMEAN5_G1P184_UFLOW0P5_2CYCLES'
$case = Join-Path $out ('case\' + $job)
$identityPath = Join-Path $case 'case_identity.json'
$identity = Get-Content -Raw $identityPath | ConvertFrom-Json
if ($identity.status -ne 'PREPARED') { throw "$job status is $($identity.status), expected PREPARED" }
if (Test-Path (Join-Path $case ($job + '.lck'))) { throw "$job has an active lock" }
. (Join-Path (Split-Path $repo -Parent) 'launch_abaqus_with_env.ps1') -SkipVerify
$ErrorActionPreference = 'Continue'
Set-Location -LiteralPath $case
$sub = Join-Path $case 'vuamp_vmean5.f90'
$gateJob = $job + '_DATACHECK'
$gateOutput = & 'I:\SIMULIA\Commands\abq2025.bat' job=$gateJob input="$job.inp" user=$sub double=both datacheck ask_delete=ON interactive 2>&1
$gateExit = $LASTEXITCODE
$gateOutput | Out-File (Join-Path $out ($gateJob + '.log')) -Encoding ascii
$gatePassed = $gateExit -eq 0 -and ($gateOutput | Out-String) -notmatch '(?i)(Abaqus Error:|Problem during compilation|Analysis exited with errors)'
if (-not $gatePassed) { throw "$job datacheck failed; no dynamics launched" }
$started = Get-Date
$solverOutput = & 'I:\SIMULIA\Commands\abq2025.bat' job=$job input="$job.inp" user=$sub double=both cpus=1 ask_delete=ON interactive 2>&1
$solverExit = $LASTEXITCODE
$solverOutput | Out-File (Join-Path $out ($job + '_solver.log')) -Encoding ascii
if ($solverExit -ne 0) { throw "$job dynamics failed: $solverExit" }
$wallclock = ((Get-Date)-$started).TotalSeconds
$extract = Join-Path $repo 'calibration_analysis\FastPrecomputedNoReboundScreen\scripts\extract_fast_odb.py'
$extractOutput = & 'I:\SIMULIA\Commands\abq2025.bat' python $extract $job (Join-Path $case 'private') 2>&1
if ($LASTEXITCODE -ne 0) { throw "$job ODB extraction failed" }
$extractOutput | Out-File (Join-Path $out ($job + '_extract.log')) -Encoding ascii
$identity.status = 'SOLVED'
$identity.dynamics_run_count = 1
$identity | Add-Member wallclock_s $wallclock -Force
$identity | Add-Member cpus 1 -Force
$identity | Add-Member socket_calls 0 -Force
$identity | ConvertTo-Json -Depth 12 | Set-Content $identityPath -Encoding ascii
[ordered]@{job=$job; G_mT=$identity.gradient_mT; wallclock_s=$wallclock; datacheck_passed=$true} |
    ConvertTo-Json -Depth 5 | Set-Content (Join-Path $out 'refinement_run_status.json') -Encoding ascii
Write-Host ("COMPLETE {0}; wallclock={1:N1}s" -f $job,$wallclock)
