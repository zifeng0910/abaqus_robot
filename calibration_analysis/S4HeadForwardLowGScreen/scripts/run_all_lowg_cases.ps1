$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$jobs = @(
    @{Name='S4_HEADFORWARD_G0'; Port=65521},
    @{Name='S4_HEADFORWARD_G0P25'; Port=65522},
    @{Name='S4_HEADFORWARD_G0P5'; Port=65523},
    @{Name='S4_HEADFORWARD_G1P0'; Port=65524}
)
$processes = @()
foreach ($item in $jobs) {
    $log = Join-Path $here ("{0}_runner.log" -f $item.Name)
    $err = Join-Path $here ("{0}_runner.err" -f $item.Name)
    $arguments = @('-NoProfile','-ExecutionPolicy','Bypass','-File',(Join-Path $here 'run_lowg_case.ps1'),'-Job',$item.Name,'-Port',([string]$item.Port))
    $processes += Start-Process -FilePath 'powershell.exe' -ArgumentList $arguments -WindowStyle Hidden -RedirectStandardOutput $log -RedirectStandardError $err -PassThru
}
$processes | Select-Object Id,ProcessName
Write-Host 'Four head-forward low-gradient solves launched in parallel.'
