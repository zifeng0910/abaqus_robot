$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$python = 'I:\SIMULIA\EstProducts\2025\win_b64\tools\SMApy\python3.10\python.exe'
$jobs = @('S4_HEADFORWARD_G0','S4_HEADFORWARD_G0P25','S4_HEADFORWARD_G0P5','S4_HEADFORWARD_G1P0')
while ($true) {
    $ready = $true
    foreach ($job in $jobs) {
        $case = Join-Path $root ('cases\' + $job)
        $identity = Get-Content -Raw (Join-Path $case 'case_identity.json') | ConvertFrom-Json
        if ($identity.status -ne 'SOLVED' -or -not (Test-Path (Join-Path $case ($job + '.odb'))) -or -not (Test-Path (Join-Path $case 'private\rp_history_private.npz'))) { $ready = $false }
    }
    if ($ready) { break }
    Start-Sleep -Seconds 30
}
Set-Location -LiteralPath $root
& $python (Join-Path $root 'scripts\analyze_and_render_lowg.py')
if ($LASTEXITCODE -ne 0) { throw "low-G analysis failed with exit code $LASTEXITCODE" }
Write-Host 'LOWG_SCREEN_FINISHED'
