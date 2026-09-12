[CmdletBinding()]
param([string]$WorkDir='J:\abaqusfangzhen',[int]$Port=65501)
& (Join-Path $WorkDir 'run_wobble_f30_g6l45_wallon_0083.ps1') -WorkDir $WorkDir -Port $Port
