$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$jobs = @(
    @{Name='FAST_PLANAR_15HZ_A12_FORWARD'; Port=65511},
    @{Name='FAST_PLANAR_20HZ_A12_FORWARD'; Port=65512},
    @{Name='FAST_PLANAR_20HZ_A14P343_FORWARD'; Port=65513},
    @{Name='FAST_ELLIPTIC_20HZ_A14P343_X2P5_FORWARD'; Port=65514}
)
$processes = @()
foreach ($item in $jobs) {
    $log = Join-Path $here ("{0}_runner.log" -f $item.Name)
    $err = Join-Path $here ("{0}_runner.err" -f $item.Name)
    $arguments = @('-NoProfile','-ExecutionPolicy','Bypass','-File',(Join-Path $here 'run_fast_case.ps1'),
                   '-Job',$item.Name,'-Port',([string]$item.Port))
    $processes += Start-Process -FilePath 'powershell.exe' -ArgumentList $arguments -WindowStyle Hidden `
        -RedirectStandardOutput $log -RedirectStandardError $err -PassThru
}
$processes | Select-Object Id,ProcessName
Write-Host 'Four fast-screen solves launched in parallel.'
