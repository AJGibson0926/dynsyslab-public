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

& $installer -PackagePath $PackagePath -ChecksumPath $ChecksumPath -InstallRoot $install -DataRoot $data -NoShortcuts -NoUserConfig
$cli = Join-Path $install 'bin\dynsyslab-cli.exe'
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
Write-Host "Release installer checks passed: checksum, ZIP containment, native solve, both apps, data preservation. Evidence: $work"
