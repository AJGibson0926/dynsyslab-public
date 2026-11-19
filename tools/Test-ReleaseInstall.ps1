[CmdletBinding()]
param(
    [Parameter(Mandatory=$true)][string]$PackagePath,
    [Parameter(Mandatory=$true)][string]$ChecksumPath
)
$ErrorActionPreference = 'Stop'
$sourceRoot = Split-Path -Parent $PSScriptRoot
$installer = Join-Path $sourceRoot 'install.ps1'
$PackagePath = (Resolve-Path -LiteralPath $PackagePath).Path
$ChecksumPath = (Resolve-Path -LiteralPath $ChecksumPath).Path
$work = Join-Path $sourceRoot ('build\release-installer-test-' + [Guid]::NewGuid().ToString('N'))
$install = Join-Path $work 'app'
$data = Join-Path $work 'data'
New-Item -ItemType Directory -Force -Path $work, $data | Out-Null
$sentinel = Join-Path $data 'keep-my-scientific-data.txt'
[IO.File]::WriteAllText($sentinel, 'preserve this existing data')
$runtimePath = Join-Path $env:LOCALAPPDATA 'DynSysLab\runtime.json'
$oldRuntime = if (Test-Path -LiteralPath $runtimePath) { (Get-FileHash -LiteralPath $runtimePath).Hash } else { '' }

function Expect-Rejection([scriptblock]$Action, [string]$Expected) {
    $rejected = $false
    try { & $Action } catch {
        if ($_.Exception.Message -notlike "*$Expected*") { throw }
        $rejected = $true
    }
    if (-not $rejected) { throw "Installer did not reject: $Expected" }
}

$badChecksum = Join-Path $work 'bad-checksums.txt'
[IO.File]::WriteAllText($badChecksum, ('0' * 64) + '  ' + [IO.Path]::GetFileName($PackagePath))
Expect-Rejection {
    & $installer -PackagePath $PackagePath -ChecksumPath $badChecksum -InstallRoot $install -DataRoot $data -NoShortcuts -NoUserConfig
} 'checksum mismatch'
if (Test-Path -LiteralPath $install) { throw 'Rejected checksum wrote an installation.' }

Add-Type -AssemblyName System.IO.Compression.FileSystem
$unsafeZip = Join-Path $work 'unsafe.zip'
$zip = [IO.Compression.ZipFile]::Open($unsafeZip, [IO.Compression.ZipArchiveMode]::Create)
try {
    $entry = $zip.CreateEntry('../outside.txt')
    $writer = [IO.StreamWriter]::new($entry.Open())
    try { $writer.Write('must not escape') } finally { $writer.Dispose() }
} finally { $zip.Dispose() }
$unsafeChecksum = Join-Path $work 'unsafe-checksums.txt'
[IO.File]::WriteAllText($unsafeChecksum, (Get-FileHash -LiteralPath $unsafeZip).Hash + '  unsafe.zip')
Expect-Rejection {
    & $installer -PackagePath $unsafeZip -ChecksumPath $unsafeChecksum -InstallRoot $install -DataRoot $data -NoShortcuts -NoUserConfig
} 'unsafe path'

$installOutput = & $installer -PackagePath $PackagePath -ChecksumPath $ChecksumPath -InstallRoot $install -DataRoot $data -NoShortcuts -NoUserConfig *>&1 | Out-String
Write-Host $installOutput
$cli = Join-Path $install 'bin\dynsyslab-cli.exe'
$archives = Join-Path $data 'archives'
# The installer's own output is the first-run experience: welcome, example results, figures, no GUI.
foreach ($expected in @('# Welcome to DynSysLab', '## Quickstart results', '| Duffing double well |',
                        '## Figures', '## Things you can ask next', 'The applications were not started')) {
    if (-not $installOutput.Contains($expected)) { throw "Installer output is missing: $expected" }
}
if ($installOutput -match 'AI agents:') { throw 'Installer output addresses agents directly.' }
$installedPrefix = $install.TrimEnd('\') + '\'
$gui = @(Get-Process -ErrorAction SilentlyContinue | Where-Object {
    try { $_.Path -and $_.Path.StartsWith($installedPrefix, [StringComparison]::OrdinalIgnoreCase) } catch { $false } })
if ($gui.Count) { throw "Installer started: $($gui.Name -join ', ')" }
$record = Get-Content -LiteralPath (Join-Path $install 'release-install.json') -Raw | ConvertFrom-Json
foreach ($field in @('version', 'cli', 'welcome_document', 'quickstart_command', 'archive_root')) {
    if (-not $record.$field) { throw "release-install.json lacks $field" }
}
if (-not (Test-Path -LiteralPath $record.welcome_document)) { throw 'Installed welcome document is missing.' }
# quickstart: JSON form, default archive root from the install, distinct runs, long working directory.
$long = Join-Path $work ('long-working-directory-' + ('x' * 150))
New-Item -ItemType Directory -Force -Path $long | Out-Null
Push-Location -LiteralPath $long
try {
    $first = (& $cli quickstart --format json | Out-String) | ConvertFrom-Json
    if ($LASTEXITCODE -ne 0) { throw 'quickstart failed from a long working directory.' }
    $second = (& $cli quickstart --format json | Out-String) | ConvertFrom-Json
    if ($LASTEXITCODE -ne 0) { throw 'second quickstart failed.' }
} finally { Pop-Location }
if ($first.schema -ne 'dynsyslab.quickstart/v1' -or $first.runs.Count -lt 2) { throw 'quickstart JSON report is incomplete.' }
if ([IO.Path]::GetFullPath($first.archive_root) -ne [IO.Path]::GetFullPath($archives)) {
    throw "quickstart used $($first.archive_root) instead of the installed archive folder."
}
$paths = @($first.runs.archive) + @($second.runs.archive)
if (@($paths | Select-Object -Unique).Count -ne $paths.Count) { throw 'Repeated quickstart runs share an archive.' }
foreach ($path in $paths) { if (-not (Test-Path -LiteralPath $path)) { throw "Missing quickstart archive: $path" } }
foreach ($figure in @($first.runs.figures) + @($second.runs.figures)) {
    if (-not $figure.rendered -or -not (Test-Path -LiteralPath $figure.path)) { throw "Quickstart figure not rendered: $($figure.title)" }
}
if (Get-ChildItem -LiteralPath $long -Force) { throw 'quickstart wrote into the working directory.' }
$example = Join-Path $install 'resources\examples\duffing-free-solve.json'
$resultText = & $cli run $example --archive-root (Join-Path $data 'archives')
if ($LASTEXITCODE -ne 0) { throw 'Installed example solve failed.' }
$result = $resultText | ConvertFrom-Json
if ($result.schema -ne 'dynsyslab.workflow-result/v1' -or $result.archives.Count -lt 1) {
    throw 'Installed example returned no scientific archive.'
}
foreach ($app in @('DynSysLab.exe', 'DynSysLabViewer.exe')) {
    $exe = Join-Path $install ('bin\' + $app)
    & $exe --smoke-test --project-root $data --archive-root (Join-Path $data 'archives') --cache (Join-Path $work 'cache')
    if ($LASTEXITCODE -ne 0) { throw "Installed $app did not initialize." }
}
if ((Get-Content -LiteralPath $sentinel -Raw) -ne 'preserve this existing data') {
    throw 'Installer changed existing scientific data.'
}
$newRuntime = if (Test-Path -LiteralPath $runtimePath) { (Get-FileHash -LiteralPath $runtimePath).Hash } else { '' }
if ($oldRuntime -ne $newRuntime) { throw 'Isolated install changed the user runtime configuration.' }
Write-Host "Release installer checks passed: checksum, ZIP containment, first-run quickstart without GUI, distinct repeated runs, long working directory, native solve, both apps, data preservation. Evidence: $work"
