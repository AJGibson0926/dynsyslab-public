# Release gate: DynSysLab ships settings, parameters, schemas and instructions, never solver output.
#   -Tracked       checks the public export (tracked files minus export-ignore in .gitattributes)
#   -PackagePath   checks a release ZIP
# Fails on: archives/(solutions|derived|campaigns|hierarchies)/ paths; .npz .npy .mat .h5 .parquet;
# PNG/JPG outside native_ui/assets and resources/assets; .jsonl; and (tracked) text that names an
# archived run outside test code that creates its own temporary archives.
[CmdletBinding()]
param([switch]$Tracked, [string]$PackagePath)
$ErrorActionPreference = 'Stop'
if (-not $Tracked -and -not $PackagePath) { throw 'Give -Tracked and/or -PackagePath.' }
$sourceRoot = Split-Path -Parent $PSScriptRoot

function Get-Violation([string]$path) {
    $p = $path.Replace('\', '/')
    if ($p -match '(^|/)archives/(solutions|derived|campaigns|hierarchies)/') { return 'archive path' }
    if ($p -match '\.(npz|npy|mat|h5|parquet)$') { return 'numerical data file' }
    if ($p -match '\.(png|jpe?g)$' -and $p -notmatch '(^|/)(native_ui|resources)/assets/') { return 'image of results' }
    if ($p -match '\.jsonl$') { return 'JSON-lines record' }
    # Data packs are instructions only: job files, a manifest, checksums and a README.
    if ($p -match '(^|/)data-packs/' -and $p -notmatch '\.(json|md|sha256)$') { return 'non-text data-pack file' }
    return $null
}

$failures = [Collections.Generic.List[string]]::new()
if ($Tracked) {
    Push-Location $sourceRoot
    try {
        $files = @(git ls-files)
        # Exactly what `git archive` (the public export) contains.
        $published = @(git ls-files -- ':(attr:!export-ignore)')
        $ignored = @{}; foreach ($file in $files) { $ignored[$file] = $true }
        foreach ($file in $published) { $ignored.Remove($file) }
        foreach ($file in $published) {
            $why = Get-Violation $file
            if ($why) { $failures.Add("tracked ${why}: $file") }
        }
        # Text naming an archived run (outside test code that writes its own temporary archives).
        $text = $published | Where-Object { $_ -match '\.(json|md|html|txt|py|ps1|toml|ya?ml)$' -and
            $_ -notmatch '(^|/)(tests?|tester)/' }
        foreach ($file in $text) {
            if (Select-String -LiteralPath $file -Pattern 'archives/(solutions|derived)/[^" ]+\.npz' -Quiet) {
                $failures.Add("tracked file names an archived run: $file")
            }
        }
        Write-Host "Checked $($published.Count) published files ($($ignored.Count) private-only files left out)."
    } finally { Pop-Location }
}
if ($PackagePath) {
    Add-Type -AssemblyName System.IO.Compression, System.IO.Compression.FileSystem
    $zip = [IO.Compression.ZipFile]::OpenRead((Resolve-Path -LiteralPath $PackagePath).Path)
    try {
        foreach ($entry in $zip.Entries) {
            if (-not $entry.Name) { continue }
            $why = Get-Violation $entry.FullName
            if ($why) { $failures.Add("packaged ${why}: $($entry.FullName)") }
        }
        Write-Host "Checked $($zip.Entries.Count) packaged entries."
    } finally { $zip.Dispose() }
}
if ($failures.Count) {
    $failures | Select-Object -First 50 | ForEach-Object { Write-Host "  $_" }
    throw "Release gate: $($failures.Count) item(s) would ship solver output or results of past runs."
}
Write-Host 'Release gate passed: no solver output, archived runs or result images.'
