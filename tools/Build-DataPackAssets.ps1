# Seals every pack in data-packs/ with a release's own binary and writes the release assets:
#   DynSysLab-datapack-<name>-<release-date>-v<version>.zip and SHA256SUMS-datapacks.txt.
# Run it with the CLI from the published release ZIP, never a local development build.
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$Cli,
    [Parameter(Mandatory = $true)][string]$OutputDirectory,
    [string]$Workspace = (Join-Path ([IO.Path]::GetTempPath()) ('dynsyslab-datapacks-' + [Guid]::NewGuid().ToString('N'))),
    [double]$MinFreeGb = 6
)
$ErrorActionPreference = 'Stop'
$source = Join-Path (Split-Path -Parent $PSScriptRoot) 'data-packs'
New-Item -ItemType Directory -Force -Path $OutputDirectory, $Workspace | Out-Null
$lines = [Collections.Generic.List[string]]::new()
foreach ($pack in Get-ChildItem -LiteralPath $source -Directory) {
    $archives = Join-Path $Workspace ($pack.Name + '-run')
    $sealed = Join-Path $Workspace ($pack.Name + '-sealed')
    & $Cli datapack validate $pack.FullName
    if ($LASTEXITCODE -ne 0) { throw "Pack failed validation: $($pack.Name)" }
    & $Cli datapack run $pack.FullName --archive-root $archives --min-free-gb $MinFreeGb
    if ($LASTEXITCODE -ne 0) { throw "Pack run failed: $($pack.Name)" }
    & $Cli datapack seal $pack.FullName --archive-root $archives --output $sealed
    if ($LASTEXITCODE -ne 0) { throw "Sealing failed: $($pack.Name)" }
    & $Cli datapack verify $sealed --archive-root $archives
    if ($LASTEXITCODE -ne 0) { throw "Sealed pack does not verify: $($pack.Name)" }
    $export = (& $Cli datapack export $sealed --output $OutputDirectory | Out-String) | ConvertFrom-Json
    if ($LASTEXITCODE -ne 0) { throw "Export failed: $($pack.Name)" }
    $lines.Add($export.sha256 + '  ' + [IO.Path]::GetFileName($export.file))
    Remove-Item -LiteralPath $archives -Recurse -Force
}
$sums = Join-Path $OutputDirectory 'SHA256SUMS-datapacks.txt'
[IO.File]::WriteAllText($sums, (($lines | Sort-Object) -join "`n") + "`n", [Text.UTF8Encoding]::new($false))
Remove-Item -LiteralPath $Workspace -Recurse -Force
Write-Host "Data-pack assets written to $OutputDirectory"
