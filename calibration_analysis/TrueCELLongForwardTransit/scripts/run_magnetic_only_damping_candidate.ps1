param([Parameter(Mandatory=$true)][string]$FactorTag)
$ErrorActionPreference = 'Stop'
$root = Split-Path $PSScriptRoot -Parent
$name = "TRUECEL_MAGNETIC_ONLY_DAMPING_${FactorTag}X"
$case = Join-Path $root "case\$name"
$identityPath = Join-Path $case 'case_identity.json'
$identity = Get-Content -Raw -LiteralPath $identityPath | ConvertFrom-Json
if ($identity.status -ne 'PREPARED' -or $identity.dynamics_run_count -ne 0) {
    throw "Candidate is not prepared for one fresh run: $name"
}
. (Join-Path (Split-Path (Split-Path (Split-Path $root -Parent) -Parent) -Parent) 'launch_abaqus_with_env.ps1') -SkipVerify
$abaqus = 'I:\SIMULIA\Commands\abq2025.bat'
$subroutine = Join-Path $case 'vuamp_precomputed_truecel.f90'
Set-Location -LiteralPath $case
$check = "${name}_CHECK"
& $abaqus job=$check input="$name.inp" user=$subroutine double=both datacheck ask_delete=ON interactive *> "$check.log"
$sta = Get-Content -Raw -LiteralPath "$check.sta" -ErrorAction SilentlyContinue
$dat = Get-Content -Raw -LiteralPath "$check.dat" -ErrorAction SilentlyContinue
if ($sta -notmatch 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY' -or $dat -match '\*\*\*ERROR') {
    throw "Datacheck failed: $name; inspect $check.log"
}
Write-Host "DATACHECK PASSED $name"
$identity.status = 'RUNNING'
$identity.dynamics_run_count = 1
$identity | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $identityPath -Encoding ascii
try {
    & $abaqus job=$name input="$name.inp" user=$subroutine double=both cpus=1 ask_delete=ON interactive *> "$name.log"
    $sta = Get-Content -Raw -LiteralPath "$name.sta" -ErrorAction SilentlyContinue
    if ($sta -notmatch 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY') {
        throw "Dynamics did not complete: $name; inspect $name.log"
    }
    $extract = Join-Path (Split-Path $root -Parent) 'ProductionLocalFrameValidation\scripts\extract_dynamic_odb.py'
    & $abaqus python $extract $name (Join-Path $case 'private') *> "$name.extract.log"
    if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath (Join-Path $case 'private\rp_history_private.npz'))) {
        throw "ODB extraction failed: $name"
    }
    $identity.status = 'SOLVED'
    Write-Host "SOLVED $name"
} catch {
    $identity.status = 'FAILED'
    $identity | Add-Member -NotePropertyName failure -NotePropertyValue $_.Exception.Message -Force
    throw
} finally {
    $identity | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $identityPath -Encoding ascii
}
