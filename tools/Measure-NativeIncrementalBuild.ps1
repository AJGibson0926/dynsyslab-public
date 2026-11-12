[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$Source,
    [string]$BuildDirectory = "build/native-suite-fast",
    [string]$Target = "dynsyslab_compute",
    [ValidateRange(1, 64)]
    [int]$Jobs = 4,
    [string]$JsonOutput
)

$ErrorActionPreference = "Stop"
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
function Resolve-ProjectPath([string]$Value) {
    if ([IO.Path]::IsPathRooted($Value)) { return $Value }
    return Join-Path $projectRoot $Value
}
$sourcePath = Resolve-ProjectPath $Source
$buildPath = Resolve-ProjectPath $BuildDirectory
if (-not (Test-Path -LiteralPath $sourcePath -PathType Leaf)) {
    throw "Incremental-build source does not exist: $sourcePath"
}
if (-not (Test-Path -LiteralPath $buildPath -PathType Container)) {
    throw "Native build directory does not exist: $buildPath"
}

$item = Get-Item -LiteralPath $sourcePath
$originalTime = $item.LastWriteTimeUtc
$stopwatch = [Diagnostics.Stopwatch]::StartNew()
$output = @()
$exitCode = 1
try {
    $item.LastWriteTimeUtc = [DateTime]::UtcNow.AddSeconds(1)
    $output = @(& cmake --build $buildPath --target $Target -j $Jobs 2>&1)
    $exitCode = $LASTEXITCODE
} finally {
    $stopwatch.Stop()
    (Get-Item -LiteralPath $sourcePath).LastWriteTimeUtc = $originalTime
}
$actions = @($output | Where-Object { $_ -match '^\[[0-9]+/[0-9]+\]' })
$compiles = @($actions | Where-Object { $_ -match 'Building (?:CXX|C) object' })
$ninja = Get-Command ninja -ErrorAction Stop
$graphCommands = @(& $ninja.Source -C $buildPath -t commands $Target 2>&1)
if ($LASTEXITCODE -ne 0) {
    throw "Ninja could not inspect $Target after the build"
}
$eigenTargets = [Collections.Generic.HashSet[string]]::new(
    [StringComparer]::OrdinalIgnoreCase)
foreach ($command in $graphCommands) {
    if ($command -match '(?i)[\\/]eigen3(?:[\\/]|\s|$)' -and
        $command -match 'CMakeFiles[\\/]([^\\/ ]+)\.dir[\\/]') {
        [void]$eigenTargets.Add($Matches[1])
    }
}
$eigenCompiles = @($compiles | Where-Object {
    if ($_ -notmatch 'CMakeFiles[\\/]([^\\/ ]+)\.dir[\\/]') { return $false }
    return $eigenTargets.Contains($Matches[1])
})
$result = [ordered]@{
    schema = "dynsyslab.native-incremental-build/v1"
    source = (Resolve-Path -LiteralPath $sourcePath).Path
    buildDirectory = (Resolve-Path -LiteralPath $buildPath).Path
    target = $Target
    jobs = $Jobs
    success = $exitCode -eq 0
    exitCode = $exitCode
    wallSeconds = [math]::Round($stopwatch.Elapsed.TotalSeconds, 3)
    actionCount = $actions.Count
    compileCount = $compiles.Count
    eigenCompileCount = $eigenCompiles.Count
    eigenCompileTargets = @($eigenTargets | Sort-Object)
    actions = $actions
}
$json = $result | ConvertTo-Json -Depth 6
if ($JsonOutput) {
    $outputPath = Resolve-ProjectPath $JsonOutput
    $parent = Split-Path -Parent $outputPath
    if ($parent) { New-Item -ItemType Directory -Force -Path $parent | Out-Null }
    Set-Content -LiteralPath $outputPath -Value $json -Encoding utf8
}
$json
if ($exitCode -ne 0) {
    Write-Error ($output -join [Environment]::NewLine)
    exit $exitCode
}
