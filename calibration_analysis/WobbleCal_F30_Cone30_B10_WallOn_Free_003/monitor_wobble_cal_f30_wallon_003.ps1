param([string]$JobName='WobbleCal_COMFixed_F30_Cone30_B10_WallOn_003',[switch]$Once)
$job=$JobName
$sta=Join-Path (Get-Location) ($job+'.sta')
Write-Host "=== $job monitor ==="
if(Test-Path $sta){
  $tail=Get-Content $sta -Tail 40
  $prog=$tail | Select-String -Pattern '^\s*\d+\s+[0-9.E+-]+\s+[0-9.E+-]+'
  if($prog){$prog | Select-Object -Last 1 | ForEach-Object {Write-Host $_.Line}}
  $frame=$tail | Select-String -Pattern 'Output Field Frame Number'
  if($frame){$frame | Select-Object -Last 1 | ForEach-Object {Write-Host $_.Line}}
}
$p=Get-CimInstance Win32_Process | Where-Object { $_.CommandLine -match [regex]::Escape($job) -and $_.Name -match 'explicit|SMAPython|python' }
$log=Join-Path (Get-Location) ($job+'_run.log')
$done=$false
foreach($f in @($sta,$log)){
  if(Test-Path $f){
    if(Select-String -Path $f -Pattern 'COMPLETED SUCCESSFULLY|COMPLETED|COMPLETED WITH ERRORS|Abaqus JOB .* COMPLETED' -Quiet){$done=$true}
  }
}
if($done){Write-Host 'Abaqus/Socket: COMPLETED'; exit 0}
if($p){Write-Host 'Abaqus/Socket: RUNNING'} else {Write-Host 'Abaqus/Socket: STOPPED_OR_NOT_DETECTED'; if(-not $Once){exit 2}}
if(-not $Once){
  while($true){
    Start-Sleep -Seconds 60
    & $PSCommandPath -JobName $job -Once
    if($LASTEXITCODE -ne 0){break}
    $finished=$false
    foreach($f in @($sta,$log)){
      if(Test-Path $f -and (Select-String -Path $f -Pattern 'COMPLETED SUCCESSFULLY|COMPLETED|COMPLETED WITH ERRORS|Abaqus JOB .* COMPLETED' -Quiet)){$finished=$true}
    }
    if($finished){break}
  }
}
