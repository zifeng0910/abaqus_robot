[CmdletBinding()]
param([string]$WorkDir='J:\abaqusfangzhen',[int]$RefreshSeconds=30)
& (Join-Path $WorkDir 'monitor_wobble_f30_g6l45_wallon_0083.ps1') -WorkDir $WorkDir -RefreshSeconds $RefreshSeconds
