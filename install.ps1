# Install a verified prebuilt Windows release. No compiler, Python, or admin rights.
[CmdletBinding()]
param(
    [string]$Repository = 'AJGibson0926/dynsyslab-public',
    [string]$Version = 'latest',
    [string]$InstallRoot = (Join-Path $env:LOCALAPPDATA 'Programs\DynSysLab'),
    [string]$DataRoot,
    [string]$PackagePath,
    [string]$ChecksumPath,
    [switch]$NoShortcuts,
    [switch]$NoUserConfig
)

$ErrorActionPreference = 'Stop'
if ($env:OS -ne 'Windows_NT' -or -not [Environment]::Is64BitOperatingSystem) {
    throw 'This release supports 64-bit Windows. See docs/BUILD_NATIVE.md for source builds.'
}
if ($Repository -notmatch '^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$') { throw 'Invalid GitHub repository name.' }
if ($Version -ne 'latest' -and $Version -notmatch '^v?[0-9]+\.[0-9]+\.[0-9]+(?:[-.][A-Za-z0-9.]+)?$') {
    throw 'Version must be latest or a release tag such as v0.6.0.'
}
if ([bool]$PackagePath -ne [bool]$ChecksumPath) {
    throw 'Offline installation requires both -PackagePath and -ChecksumPath.'
}
$InstallRoot = [IO.Path]::GetFullPath($InstallRoot)
if ($InstallRoot.TrimEnd('\') -eq [IO.Path]::GetPathRoot($InstallRoot).TrimEnd('\')) {
    throw 'InstallRoot must be an application directory, not a drive root.'
}
$runtimeFile = Join-Path $env:LOCALAPPDATA 'DynSysLab\runtime.json'
$explicitDataRoot = $PSBoundParameters.ContainsKey('DataRoot')
$previous = $null
if (-not $NoUserConfig -and (Test-Path -LiteralPath $runtimeFile -PathType Leaf)) {
    $previous = Get-Content -LiteralPath $runtimeFile -Raw | ConvertFrom-Json
}
if (-not $DataRoot) {
    $DataRoot = if ($previous -and $previous.project_root) { $previous.project_root } else {
        Join-Path ([Environment]::GetFolderPath('MyDocuments')) 'DynSysLab'
    }
}
$DataRoot = [IO.Path]::GetFullPath($DataRoot)
$archiveRoot = if (-not $explicitDataRoot -and $previous -and $previous.archive_root) {
    [IO.Path]::GetFullPath($previous.archive_root)
} else { Join-Path $DataRoot 'archives' }

$work = Join-Path ([IO.Path]::GetTempPath()) ('DynSysLab-install-' + [Guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $work | Out-Null
if (-not $PackagePath) {
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
    $endpoint = if ($Version -eq 'latest') { 'latest' } else { 'tags/' + $Version }
    $headers = @{ 'User-Agent' = 'DynSysLab-installer'; 'Accept' = 'application/vnd.github+json' }
    $release = Invoke-RestMethod -Uri "https://api.github.com/repos/$Repository/releases/$endpoint" -Headers $headers
    $assets = @($release.assets | Where-Object { $_.name -match '^DynSysLab-[0-9].*-windows-x64\.zip$' })
    $checksums = @($release.assets | Where-Object { $_.name -eq 'SHA256SUMS.txt' })
    if ($assets.Count -ne 1 -or $checksums.Count -ne 1) {
        throw 'The release must contain one Windows x64 ZIP and SHA256SUMS.txt.'
    }
    $PackagePath = Join-Path $work $assets[0].name
    $ChecksumPath = Join-Path $work 'SHA256SUMS.txt'
    Invoke-WebRequest -UseBasicParsing -Uri $assets[0].browser_download_url -OutFile $PackagePath
    Invoke-WebRequest -UseBasicParsing -Uri $checksums[0].browser_download_url -OutFile $ChecksumPath
}
$PackagePath = (Resolve-Path -LiteralPath $PackagePath).Path
$ChecksumPath = (Resolve-Path -LiteralPath $ChecksumPath).Path
$packageName = [IO.Path]::GetFileName($PackagePath)
$checksumLines = @(Get-Content -LiteralPath $ChecksumPath | Where-Object {
    $_ -match ('^[0-9a-fA-F]{64}\s+\*?' + [regex]::Escape($packageName) + '$')
})
if ($checksumLines.Count -ne 1) { throw "No unique checksum for $packageName." }
$expected = ($checksumLines[0] -split '\s+')[0]
$actual = (Get-FileHash -LiteralPath $PackagePath -Algorithm SHA256).Hash
if ($actual -ine $expected) { throw 'Release checksum mismatch. Nothing has been installed.' }

# Reject traversal and symbolic-link entries before extracting the verified ZIP.
Add-Type -AssemblyName System.IO.Compression.FileSystem
$expanded = Join-Path $work 'expanded'
$prefix = [IO.Path]::GetFullPath($expanded).TrimEnd('\') + '\'
$zip = [IO.Compression.ZipFile]::OpenRead($PackagePath)
try {
    foreach ($entry in $zip.Entries) {
        $destination = [IO.Path]::GetFullPath((Join-Path $expanded $entry.FullName))
        if (-not $destination.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) {
            throw 'Release ZIP contains an unsafe path.'
        }
        $mode = ($entry.ExternalAttributes -shr 16) -band 0xF000
        if ($mode -eq 0xA000) { throw 'Release ZIP contains a symbolic link.' }
    }
} finally { $zip.Dispose() }
Expand-Archive -LiteralPath $PackagePath -DestinationPath $expanded
$executables = @(Get-ChildItem -LiteralPath $expanded -Filter 'dynsyslab-cli.exe' -Recurse -File)
if ($executables.Count -ne 1) { throw 'Release must contain exactly one native CLI.' }
$packageRoot = Split-Path -Parent (Split-Path -Parent $executables[0].FullName)
foreach ($name in @('DynSysLab.exe', 'DynSysLabViewer.exe', 'DynSysLabViewerCli.exe')) {
    if (-not (Test-Path -LiteralPath (Join-Path $packageRoot "bin\$name") -PathType Leaf)) {
        throw "Release is missing $name."
    }
}
$installedPrefix = $InstallRoot.TrimEnd('\') + '\'
$running = @(Get-Process -ErrorAction SilentlyContinue | Where-Object {
    try { $_.Path -and $_.Path.StartsWith($installedPrefix, [StringComparison]::OrdinalIgnoreCase) } catch { $false }
})
if ($running.Count) { throw 'Close running DynSysLab applications before updating this installation.' }

New-Item -ItemType Directory -Force -Path $InstallRoot, $DataRoot, $archiveRoot | Out-Null
Get-ChildItem -LiteralPath $packageRoot -Force | ForEach-Object {
    Copy-Item -LiteralPath $_.FullName -Destination $InstallRoot -Recurse -Force
}
$config = [ordered]@{ schema='dynsyslab.runtime-paths/v1'; project_root=$DataRoot; archive_root=$archiveRoot } | ConvertTo-Json
$utf8 = [Text.UTF8Encoding]::new($false)
$installedRuntime = Join-Path $InstallRoot 'resources\runtime.json'
[IO.File]::WriteAllText($installedRuntime, $config, $utf8)
$cli = Join-Path $InstallRoot 'bin\dynsyslab-cli.exe'
& $cli doctor --archive-root $archiveRoot
if ($LASTEXITCODE -ne 0) { throw 'Installed native diagnostics failed.' }
if (-not $NoUserConfig) {
    New-Item -ItemType Directory -Force -Path (Split-Path -Parent $runtimeFile) | Out-Null
    [IO.File]::WriteAllText($runtimeFile, $config, $utf8)
}
if (-not $NoShortcuts) {
    $shell = New-Object -ComObject WScript.Shell
    $menu = Join-Path $env:APPDATA 'Microsoft\Windows\Start Menu\Programs\DynSysLab'
    New-Item -ItemType Directory -Force -Path $menu | Out-Null
    foreach ($app in @(@('DynSysLab','DynSysLab.exe','dynsyslab-app.ico'), @('Viewer','DynSysLabViewer.exe','dynsyslab-viewer.ico'))) {
        $shortcut = $shell.CreateShortcut((Join-Path $menu ($app[0] + '.lnk')))
        $shortcut.TargetPath = Join-Path $InstallRoot ('bin\' + $app[1])
        $shortcut.WorkingDirectory = Join-Path $InstallRoot 'bin'
        $shortcut.IconLocation = (Join-Path $InstallRoot ('resources\assets\' + $app[2])) + ',0'
        $shortcut.Save()
    }
}
[IO.File]::WriteAllText((Join-Path $InstallRoot 'release-install.json'),
    ([ordered]@{schema='dynsyslab.release-install/v1'; repository=$Repository; package=$packageName; sha256=$actual} | ConvertTo-Json), $utf8)
Write-Host "DynSysLab installed: $InstallRoot"
Write-Host "Scientific archives: $archiveRoot"
Write-Host "Launch: $(Join-Path $InstallRoot 'bin\DynSysLab.exe')"
Write-Host "Verified download retained at: $work"
