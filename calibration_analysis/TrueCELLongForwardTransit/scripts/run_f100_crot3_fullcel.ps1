param([switch]$CheckOnly, [switch]$RunOnly,
      [ValidateSet('TRUECEL_B0P11_G2P20_A14P5_F100_CROT3_FULLCEL50',
                   'TRUECEL_B0P11_G2P20_A14P5_F100_CROT1_FULLCEL50',
                   'TRUECEL_B0P11_G2P20_A14P5_F100_CROT0P3_FULLCEL50')]
      [string]$Name = 'TRUECEL_B0P11_G2P20_A14P5_F100_CROT3_FULLCEL50')
$ErrorActionPreference = 'Stop'
$root = Split-Path $PSScriptRoot -Parent
$name = $Name
$case = Join-Path $root "case\$name"
$identityPath = Join-Path $case 'case_identity.json'
$identity = Get-Content -Raw -LiteralPath $identityPath | ConvertFrom-Json
if ($identity.status -ne 'PREPARED' -or $identity.dynamics_run_count -ne 0) {
    throw "One-run gate refused: $name"
}
. (Join-Path (Split-Path (Split-Path (Split-Path $root -Parent) -Parent) -Parent) 'launch_abaqus_with_env.ps1') -SkipVerify
$abaqus = 'I:\SIMULIA\Commands\abq2025.bat'
$subroutine = Join-Path $case 'vuamp_precomputed_truecel.f90'
$check = "${name}_CHECK"
Set-Location -LiteralPath $case
if (-not $RunOnly) {
    & $abaqus job=$check input="$name.inp" user=$subroutine double=both datacheck cpus=1 ask_delete=ON interactive *> "$check.log"
    $sta = Get-Content -Raw -LiteralPath "$check.sta" -ErrorAction SilentlyContinue
    $dat = Get-Content -Raw -LiteralPath "$check.dat" -ErrorAction SilentlyContinue
    if ($sta -notmatch 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY' -or $dat -match '\*\*\*ERROR') {
        throw "Datacheck failed: $check"
    }
    Write-Host "DATACHECK PASSED $name"
}
if ($CheckOnly) { return }
$checkSta = Get-Content -Raw -LiteralPath "$check.sta" -ErrorAction SilentlyContinue
if ($checkSta -notmatch 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY') {
    throw "A completed datacheck is required before dynamics"
}
$identity.status = 'RUNNING'
$identity.dynamics_run_count = 1
$identity | Add-Member -NotePropertyName run_started_at -NotePropertyValue (Get-Date).ToString('o') -Force
$identity | ConvertTo-Json -Depth 14 | Set-Content -LiteralPath $identityPath -Encoding utf8
try {
    & $abaqus job=$name input="$name.inp" user=$subroutine double=both cpus=1 ask_delete=ON interactive *> "$name.run.log"
    $sta = Get-Content -Raw -LiteralPath "$name.sta" -ErrorAction SilentlyContinue
    if ($sta -notmatch 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY') {
        throw "Dynamics did not complete normally; inspect $name.run.log and .sta"
    }
    $extract = Join-Path (Split-Path $root -Parent) 'ProductionLocalFrameValidation\scripts\extract_dynamic_odb.py'
    & $abaqus python $extract $name (Join-Path $case 'private') *> "$name.extract.log"
    if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath (Join-Path $case 'private\rp_history_private.npz'))) {
        throw 'ODB history extraction failed'
    }
    $identity.status = 'SOLVED'
    Write-Host "SOLVED $name"
} catch {
    $identity.status = 'FAILED'
    $identity | Add-Member -NotePropertyName failure -NotePropertyValue $_.Exception.Message -Force
    throw
} finally {
    $identity | Add-Member -NotePropertyName run_ended_at -NotePropertyValue (Get-Date).ToString('o') -Force
    $identity | ConvertTo-Json -Depth 14 | Set-Content -LiteralPath $identityPath -Encoding utf8
}
