[CmdletBinding()]
param([Parameter(Mandatory=$true)][ValidateSet('B','C')][string]$Case)
$ErrorActionPreference='Stop'
$root='J:\abaqusfangzhen\abaqus_robot\calibration_analysis\TrueCELLongForwardTransit'
$job=if($Case -eq 'B'){'TRUECEL_A14P5_COARSETRAJ_CEL_COARSE_REPLAY_VVR'}else{'TRUECEL_A14P5_COARSETRAJ_CEL_HREFINE_REPLAY_VVR'}
$dir=Join-Path $root ('case\'+$job);$identityPath=Join-Path $dir 'case_identity.json'
$identity=Get-Content -Raw $identityPath|ConvertFrom-Json;$gate=Get-Content -Raw (Join-Path $root ($job+'_PreRun_Gate.json'))|ConvertFrom-Json
if($identity.status-ne'PREPARED'-or $identity.dynamics_run_count-ne 0){throw'Single-run gate failed'}
if((Get-FileHash (Join-Path $dir ($job+'.inp')) -Algorithm SHA256).Hash-ne$identity.input_sha256){throw'Input hash changed'}
if(-not $gate.driver_validated-or-not $gate.robot_part_byte_identical_to_coarse_original-or-not $gate.wall_parts_byte_identical_to_coarse_original){throw'Pre-run identity gate failed'}
if(Test-Path (Join-Path $dir ($job+'.sta'))){throw'Dynamics already attempted'}
. 'J:\abaqusfangzhen\launch_abaqus_with_env.ps1' -SkipVerify
$identity.status='RUNNING';$identity.dynamics_run_count=1;$identity|ConvertTo-Json -Depth 25|Set-Content $identityPath -Encoding ascii
$started=Get-Date
$p=Start-Process 'I:\SIMULIA\Commands\abq2025.bat' -ArgumentList @("job=$job","input=$job.inp",'user=vuamp_precomputed_truecel.f90','double=both','cpus=1','ask_delete=OFF','interactive') -WorkingDirectory $dir -WindowStyle Hidden -RedirectStandardOutput (Join-Path $root ($job+'_solver.log')) -RedirectStandardError (Join-Path $root ($job+'_solver_stderr.log')) -PassThru
$failure=$null
try{
 while(-not $p.HasExited){Start-Sleep 15;$p.Refresh();$msg=Join-Path $dir ($job+'.msg');$sta=Join-Path $dir ($job+'.sta');$recent=if(Test-Path $msg){(Get-Content $msg -Tail 250)-join"`n"}else{''};if($recent-match'(?i)InfoNodeDeepPenetFirst|negative volume|severe element distortion|\*\*\*ERROR'){$failure='Severe solver/contact event'};if(-not$failure-and(Test-Path $sta)){$dt=@(Get-Content $sta -Tail 100|%{if($_-match'^\s*\d+\s+[\d.E+-]+\s+[\d.E+-]+\s+\S+\s+([\d.E+-]+)\s+\d+'){[double]::Parse($Matches[1],[Globalization.CultureInfo]::InvariantCulture)}});if($dt.Count-ge 10-and@($dt|select -Last 10|?{$_-lt 2e-8}).Count-eq 10){$failure='Persistent timestep collapse'}};if($failure){&'I:\SIMULIA\Commands\abq2025.bat' "terminate job=$job"|Out-Null;break}}
 $p.WaitForExit();$completed=(Test-Path (Join-Path $dir ($job+'.sta')))-and((Get-Content -Raw (Join-Path $dir ($job+'.sta')))-match'THE ANALYSIS HAS COMPLETED SUCCESSFULLY');if($failure-or-not$completed-or$p.ExitCode-ne 0){throw"Dynamics invalid: $failure exit=$($p.ExitCode) completed=$completed"};$identity.status='SOLVED'
}catch{$identity.status='FAILED';$identity|Add-Member failure $_.Exception.Message -Force;throw}finally{$identity|Add-Member wallclock_s (((Get-Date)-$started).TotalSeconds) -Force;$identity|ConvertTo-Json -Depth 25|Set-Content $identityPath -Encoding ascii}
Write-Host "COMPLETED $job"
