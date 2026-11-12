[CmdletBinding()]
param(
    [ValidateSet('Release','RelWithDebInfo','Debug')]
    [string]$Configuration = 'Release',
    [string]$InstallRoot = (Join-Path $env:LOCALAPPDATA 'Programs\DynSysLab'),
    [string]$ExistingBuildRoot,
    [string]$ArchiveRoot,
    [switch]$SkipTests,
    [switch]$PinToTaskbar
)

$ErrorActionPreference = 'Stop'
$sourceRoot = Split-Path -Parent $PSScriptRoot
$buildRoot = if ($ExistingBuildRoot) {
    (Resolve-Path -LiteralPath $ExistingBuildRoot).Path
} else {
    Join-Path $sourceRoot 'build\native-suite-fast'
}

if (-not $ExistingBuildRoot) {
    if (-not $env:VCPKG_ROOT) {
        throw 'VCPKG_ROOT must name a vcpkg installation before building DynSysLab.'
    }
    if ($Configuration -ne 'Release') {
        throw 'The standard native installer currently supports Release; use -ExistingBuildRoot for another configuration.'
    }
    cmake --preset windows-fast-release
    if ($LASTEXITCODE -ne 0) { throw 'DynSysLab CMake configuration failed.' }
    cmake --build --preset windows-fast-release --parallel
    if ($LASTEXITCODE -ne 0) { throw 'DynSysLab native build failed.' }
}
if (-not $SkipTests) {
    cmake --build $buildRoot --config $Configuration --target dynsyslab_test_binaries --parallel
    if ($LASTEXITCODE -ne 0) { throw 'DynSysLab native test build failed.' }
    ctest --test-dir $buildRoot -C $Configuration --output-on-failure
    if ($LASTEXITCODE -ne 0) { throw 'DynSysLab native release gates failed.' }
    cmake --build $buildRoot --config $Configuration --target dynsyslab_package_smoke
    if ($LASTEXITCODE -ne 0) { throw 'DynSysLab isolated package gate failed.' }
}
cmake --install $buildRoot --config $Configuration --prefix $InstallRoot
if ($LASTEXITCODE -ne 0) { throw 'DynSysLab native installation failed.' }

$resolvedArchiveRoot = if ($ArchiveRoot) {
    New-Item -ItemType Directory -Force -Path $ArchiveRoot | Out-Null
    (Resolve-Path -LiteralPath $ArchiveRoot).Path
} else {
    (Resolve-Path -LiteralPath (Join-Path $sourceRoot 'archives')).Path
}
$runtimeRoot = Join-Path $env:LOCALAPPDATA 'DynSysLab'
$runtimePath = Join-Path $runtimeRoot 'runtime.json'
New-Item -ItemType Directory -Force -Path $runtimeRoot | Out-Null
$runtimeConfiguration = [ordered]@{
    schema = 'dynsyslab.runtime-paths/v1'
    project_root = $sourceRoot
    archive_root = $resolvedArchiveRoot
} | ConvertTo-Json
[System.IO.File]::WriteAllText(
    $runtimePath, $runtimeConfiguration, [System.Text.UTF8Encoding]::new($false))
$installedRuntimePath = Join-Path $InstallRoot 'resources\runtime.json'
New-Item -ItemType Directory -Force -Path (Split-Path -Parent $installedRuntimePath) | Out-Null
[System.IO.File]::WriteAllText(
    $installedRuntimePath, $runtimeConfiguration, [System.Text.UTF8Encoding]::new($false))
# These launch-era JSON files predate the Archiver-owned SQLite catalogue.
# They are derived metadata, never scientific archives, and can make a fallback
# launch appear stale. Native apps now read the central index directly.
$obsoleteMetadata = @(
    (Join-Path $runtimeRoot 'Compute\catalog.json'),
    (Join-Path $runtimeRoot 'Compute\registry.json'),
    (Join-Path $runtimeRoot 'Viewer\catalog.json'),
    (Join-Path $InstallRoot 'resources\catalog.json')
)
foreach ($path in $obsoleteMetadata) {
    if (Test-Path -LiteralPath $path -PathType Leaf) {
        Remove-Item -LiteralPath $path -Force
    }
}
foreach ($cacheDirectory in @(
        (Join-Path $runtimeRoot 'Compute'),
        (Join-Path $runtimeRoot 'Viewer'))) {
    if (-not (Test-Path -LiteralPath $cacheDirectory -PathType Container)) { continue }
    Get-ChildItem -LiteralPath $cacheDirectory -File |
        Where-Object {
            $_.Name -like 'startup-diagnostics-*.json' -or
            $_.Name -eq 'startup-profile-debug.csv'
        } |
        ForEach-Object { Remove-Item -LiteralPath $_.FullName -Force }
}

$executable = Join-Path $InstallRoot 'bin\DynSysLab.exe'
$icon = Join-Path $InstallRoot 'resources\assets\dynsyslab-app.ico'
if (-not (Test-Path -LiteralPath $executable -PathType Leaf)) {
    throw "Installed executable is missing: $executable"
}

function Write-Shortcut([string]$Path) {
    $directory = Split-Path -Parent $Path
    New-Item -ItemType Directory -Force -Path $directory | Out-Null
    $shell = New-Object -ComObject WScript.Shell
    $shortcut = $shell.CreateShortcut($Path)
    $shortcut.TargetPath = $executable
    $shortcut.WorkingDirectory = Split-Path -Parent $executable
    $shortcut.IconLocation = "$icon,0"
    $shortcut.Description = 'DynSysLab Solver and Post-processor'
    $shortcut.Save()
}

$startMenu = Join-Path $env:APPDATA 'Microsoft\Windows\Start Menu\Programs\DynSysLab.lnk'
Write-Shortcut $startMenu
if ($PinToTaskbar) {
    $taskbar = Join-Path $env:APPDATA 'Microsoft\Internet Explorer\Quick Launch\User Pinned\TaskBar\DynSysLab.lnk'
    Write-Shortcut $taskbar
}

Write-Host "DynSysLab installed at $InstallRoot"
Write-Host "Launch: $executable"
Write-Host "Archives: $resolvedArchiveRoot"
