param(
  [string]$WorkDir='J:\abaqusfangzhen',
  [string]$JobName='WobbleCal_COMFixed_F30_Cone30_B10',
  [int]$Port=65492,
  [int]$IntervalSeconds=30,
  [switch]$Once
)
$ErrorActionPreference='SilentlyContinue'; Set-Location -LiteralPath $WorkDir
$sta=Join-Path $WorkDir "$JobName.sta"; $msg=Join-Path $WorkDir "$JobName.msg"; $lck=Join-Path $WorkDir "$JobName.lck"; $tel=Join-Path $WorkDir "${JobName}_telemetry.csv"
do {
  Clear-Host; $stamp=Get-Date -Format 'yyyy-MM-dd HH:mm:ss'; Write-Host "=== 30 Hz COM-fixed wobble monitor ===" -ForegroundColor Cyan; Write-Host "[$stamp] $JobName"
  $tail=@(); foreach($p in @($sta,$msg)){if(Test-Path $p){$tail+=Get-Content $p -Tail 120}}
  $frame=$tail|Select-String -Pattern 'Output Field Frame Number\s+([0-9]+).*step time\s+([0-9.Ee+-]+)'|Select-Object -Last 1
  $inc=$tail|Select-String -Pattern '^\s*([0-9]+)\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s+\S+\s+([0-9.Ee+-]+)'|Select-Object -Last 1
  $done=($tail -join "`n") -match 'ANALYSIS HAS COMPLETED SUCCESSFULLY|COMPLETED SUCCESSFULLY'
  $proc=Get-CimInstance Win32_Process|Where-Object{$_.CommandLine -and $_.CommandLine -like "*$JobName*" -and $_.Name -match 'explicit|SMALauncher|package|cmd|powershell'}
  if($proc){Write-Host 'Abaqus: RUNNING' -ForegroundColor Green}elseif($done){Write-Host 'Abaqus: COMPLETED SUCCESSFULLY' -ForegroundColor Green}elseif(Test-Path $lck){Write-Host 'Abaqus: LOCK PRESENT / PROCESS NOT VISIBLE' -ForegroundColor Yellow}else{Write-Host 'Abaqus: NOT FOUND or STOPPED' -ForegroundColor Yellow}
  if($frame){Write-Host ("Frame {0}; t={1} s" -f $frame.Matches[0].Groups[1].Value,$frame.Matches[0].Groups[2].Value)}
  if($inc){Write-Host ("Increment {0}; stable dt={1} s" -f $inc.Matches[0].Groups[1].Value,$inc.Matches[0].Groups[4].Value)}
  $tcp=Get-NetTCPConnection -LocalPort $Port -State Listen; if($tcp){Write-Host "Socket $Port`: LISTENING" -ForegroundColor Green}else{Write-Host "Socket $Port`: not listening" -ForegroundColor Yellow}
  if(Test-Path $tel){$n=(Get-Content $tel|Measure-Object -Line).Lines-1;$last=Import-Csv $tel|Select-Object -Last 1;Write-Host ("Telemetry rows={0}; last t={1:F6} s; UR=({2:F3},{3:F3},{4:F3})" -f $n,[double]$last.t_s,[double]$last.ur1,[double]$last.ur2,[double]$last.ur3)}
  if(-not$Once){Start-Sleep -Seconds $IntervalSeconds}
} while(-not$Once)
