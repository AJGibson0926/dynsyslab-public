[CmdletBinding()]
param(
    [string]$BuildDirectory = "build/native-suite-fast",
    [string[]]$Targets = @("dynsyslab_compute", "dynsyslab_viewer", "dynsyslab_cli"),
    [string]$JsonOutput
)

$ErrorActionPreference = "Stop"
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$buildPath = if ([IO.Path]::IsPathRooted($BuildDirectory)) {
    $BuildDirectory
} else {
    Join-Path $projectRoot $BuildDirectory
}
if (-not (Test-Path -LiteralPath $buildPath -PathType Container)) {
    throw "Native build directory does not exist: $buildPath"
}
$ninja = Get-Command ninja -ErrorAction Stop
$commands = @(& $ninja.Source -C $buildPath -t commands @Targets 2>&1)
if ($LASTEXITCODE -ne 0) {
    throw "Ninja could not inspect the requested targets: $($commands -join [Environment]::NewLine)"
}

$compileCommands = @($commands | Where-Object {
    $_ -match '(?:^|[\\/])(?:cl|clang-cl|clang\+\+|g\+\+)\.exe"?\s+' -and
    $_ -match '(?:\s/c\s|\s-c\s)'
})
$eigenCommands = @($compileCommands | Where-Object {
    $_ -match '(?i)[\\/]eigen3(?:[\\/]|\s|$)'
})

function Get-CompileItem([string]$Command) {
    $target = "unknown"
    $source = "unknown"
    if ($Command -match 'CMakeFiles[\\/]([^\\/ ]+)\.dir[\\/]') {
        $target = $Matches[1]
    }
    if ($Command -match '(?i)[\\/]([^\\/ ]+\.(?:cpp|cc|cxx))(?:"?\s|$)') {
        $source = $Matches[1]
    } elseif ($Command -match 'cmake_pch') {
        $source = "<precompiled-header>"
    }
    [pscustomobject]@{ target = $target; source = $source }
}

$eigenItems = @($eigenCommands | ForEach-Object { Get-CompileItem $_ } |
    Sort-Object target, source -Unique)
$currentObjects = [Collections.Generic.HashSet[string]]::new(
    [StringComparer]::OrdinalIgnoreCase)
foreach ($command in $compileCommands) {
    if ($command -match '(?:/Fo|-o\s+)"?([^"\s]+\.(?:obj|o))') {
        [void]$currentObjects.Add(($Matches[1] -replace '\\', '/'))
    }
}
$timings = @()
$logPath = Join-Path $buildPath ".ninja_log"
if (Test-Path -LiteralPath $logPath -PathType Leaf) {
    $latest = @{}
    foreach ($line in Get-Content -LiteralPath $logPath) {
        if (-not $line -or $line.StartsWith('#')) { continue }
        $fields = $line -split "`t"
        if ($fields.Count -lt 5) { continue }
        $start = 0L; $end = 0L
        if (-not [long]::TryParse($fields[0], [ref]$start) -or
            -not [long]::TryParse($fields[1], [ref]$end)) { continue }
        $latest[$fields[3]] = [pscustomobject]@{
            output = $fields[3]
            milliseconds = [math]::Max(0L, $end - $start)
        }
    }
    $timings = @($latest.Values | Where-Object {
            $currentObjects.Contains(($_.output -replace '\\', '/'))
        } |
        Sort-Object milliseconds -Descending | Select-Object -First 20)
}

$ratio = if ($compileCommands.Count) {
    [math]::Round(100.0 * $eigenCommands.Count / $compileCommands.Count, 2)
} else { 0.0 }
$result = [ordered]@{
    schema = "dynsyslab.native-build-graph/v1"
    projectRoot = $projectRoot
    buildDirectory = (Resolve-Path -LiteralPath $buildPath).Path
    targets = $Targets
    commandCount = $commands.Count
    compileCommandCount = $compileCommands.Count
    eigenExposedCompileCount = $eigenCommands.Count
    eigenExposedPercent = $ratio
    eigenCompileItems = $eigenItems
    longestRecentCompiles = $timings
}
$json = $result | ConvertTo-Json -Depth 8
if ($JsonOutput) {
    $outputPath = if ([IO.Path]::IsPathRooted($JsonOutput)) {
        $JsonOutput
    } else { Join-Path $projectRoot $JsonOutput }
    $parent = Split-Path -Parent $outputPath
    if ($parent) { New-Item -ItemType Directory -Force -Path $parent | Out-Null }
    Set-Content -LiteralPath $outputPath -Value $json -Encoding utf8
}
$json
