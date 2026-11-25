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

Add-Type -AssemblyName System.IO.Compression, System.IO.Compression.FileSystem
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

$installOutput = & $installer -PackagePath $PackagePath -ChecksumPath $ChecksumPath -InstallRoot $install -DataRoot $data -NoShortcuts -NoUserConfig -NoAgentSetup *>&1 | Out-String
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
# About 200 characters: long, yet below the Windows limit for a process working directory.
$long = Join-Path $work ('long-working-directory-' + ('x' * [Math]::Max(1, 200 - $work.Length - 24)))
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
# DynSysLab Hub and agent connection, sandboxed: --peek leaves the user's hub state alone and
# --home/--no-agent-cli keep the user's own Claude Code and Codex configuration untouched.
$hub = (& $cli hub --peek --format json | Out-String) | ConvertFrom-Json
if ($hub.schema -ne 'dynsyslab.hub/v1' -or $hub.options.Count -lt 5 -or -not $hub.greeting_markdown.Contains('What would you like to do?')) {
    throw 'hub greeting is incomplete.'
}
$sandbox = Join-Path $work 'agent-home'; $workspace = Join-Path $work 'workspace'
New-Item -ItemType Directory -Force -Path (Join-Path $sandbox '.claude'), (Join-Path $sandbox '.codex') | Out-Null
$setup = (& $cli agent-setup --home $sandbox --workspace $workspace --no-agent-cli --format json | Out-String) | ConvertFrom-Json
if ($LASTEXITCODE -ne 0 -or $setup.schema -ne 'dynsyslab.agent-setup/v1') { throw 'agent-setup failed.' }
if (-not (Select-String -Quiet -LiteralPath (Join-Path $sandbox '.claude\settings.json') -SimpleMatch 'hub --hook')) { throw 'Claude start-up hook missing.' }
if (-not (Select-String -Quiet -LiteralPath (Join-Path $sandbox '.codex\config.toml') -SimpleMatch '[mcp_servers.dynsyslab]')) { throw 'Codex MCP entry missing.' }
foreach ($file in @('AGENTS.md', 'CLAUDE.md', '.dynsyslab-workspace')) {
    if (-not (Test-Path -LiteralPath (Join-Path $workspace $file))) { throw "Workspace file missing: $file" }
}
$mcpInput = '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18"}}' + "`n" +
            '{"jsonrpc":"2.0","id":2,"method":"tools/list"}'
$mcp = @($mcpInput | & $cli mcp | ForEach-Object { $_ | ConvertFrom-Json })
if ($mcp.Count -ne 2 -or $mcp[0].result.serverInfo.name -ne 'dynsyslab' -or -not ($mcp[1].result.tools.name -contains 'hub')) {
    throw 'MCP server handshake failed.'
}
& $cli agent-setup --remove --home $sandbox --workspace $workspace --no-agent-cli | Out-Null
# Fresh install: every shipped example works with an empty archive folder (no shipped solver output).
$examplesRoot = Join-Path $install 'resources\examples'
$fresh = Join-Path $work 'fresh-archives'
New-Item -ItemType Directory -Force -Path $fresh | Out-Null
$longRunning = @((Get-Content -Raw -LiteralPath (Join-Path $examplesRoot 'long-running.json') | ConvertFrom-Json).long_running)
$queryCommands = @{ 'dynsyslab.analytic-query/v1' = 'analytic'; 'dynsyslab.koopman-spectral-query/v1' = 'koopman-spectral'
    'dynsyslab.observable-closure-query/v1' = 'observable-closure'; 'dynsyslab.observable-closure-query/v2' = 'observable-closure'
    'dynsyslab.observable-query/v1' = 'observable'; 'dynsyslab.phase-space-query/v1' = 'phase-space'
    'dynsyslab.special-function-query/v1' = 'special-function' }
function Invoke-Example([string[]]$arguments, [int]$seconds) {
    $info = [Diagnostics.ProcessStartInfo]::new($cli)
    $info.Arguments = ($arguments | ForEach-Object { '"' + $_ + '"' }) -join ' '
    $info.UseShellExecute = $false; $info.RedirectStandardOutput = $true; $info.RedirectStandardError = $true
    $process = [Diagnostics.Process]::Start($info)
    $stdout = $process.StandardOutput.ReadToEndAsync(); $stderr = $process.StandardError.ReadToEndAsync()
    if (-not $process.WaitForExit($seconds * 1000)) { $process.Kill(); return 'timeout' }
    if ($process.ExitCode -ne 0) { return "exit $($process.ExitCode): $($stderr.Result)" }
    return ''
}
$ran = 0; $planned = 0
foreach ($file in Get-ChildItem -LiteralPath $examplesRoot -Recurse -Filter *.json) {
    $relative = $file.FullName.Substring($examplesRoot.Length + 1).Replace('\', '/')
    $schema = (Get-Content -Raw -LiteralPath $file.FullName | ConvertFrom-Json).schema
    if ($schema -eq 'dynsyslab.compute-job/v1') {
        foreach ($command in @('validate', 'plan')) {
            $problem = Invoke-Example @($command, $file.FullName, '--archive-root', $fresh) 120
            if ($problem) { throw "Fresh install: $command $relative failed: $problem" }
        }
        $planned++
        if ($longRunning -contains $relative) { continue }
        $arguments = @('run', $file.FullName, '--archive-root', $fresh)
    } elseif ($queryCommands.ContainsKey($schema)) {
        $arguments = @($queryCommands[$schema], $file.FullName)
        if ($queryCommands[$schema] -eq 'phase-space') { $arguments += @('--archive-root', $fresh) }
    } else { continue }
    $problem = Invoke-Example $arguments 300
    if ($problem -eq 'timeout') { throw "Fresh install: $relative ran past 5 minutes; list it in examples/native/long-running.json or shorten it." }
    if ($problem) { throw "Fresh install: $relative failed: $problem" }
    $ran++
}
Write-Host "Fresh install: $planned compute jobs validated and planned; $ran examples run from an empty archive folder."
if (Select-String -Quiet -LiteralPath (Join-Path $sandbox '.claude\settings.json') -SimpleMatch 'hub --hook') { throw 'agent-setup --remove left the hook.' }
# Tracked data packs pinned to this release validate and plan under it (no pack is run here). A
# pack pinned to another release belongs to that release and is refused by this one, by design.
$packRoot = Join-Path $sourceRoot 'data-packs'
$releaseVersion = ([regex]::Match((Get-Content -Raw -LiteralPath (Join-Path $sourceRoot 'CMakeLists.txt')), 'VERSION (\d+\.\d+\.\d+)')).Groups[1].Value
if (Test-Path -LiteralPath $packRoot) {
    foreach ($pack in Get-ChildItem -LiteralPath $packRoot -Directory) {
        $pinned = (Get-Content -Raw -LiteralPath (Join-Path $pack.FullName 'pack.json') | ConvertFrom-Json).release.version
        if ($pinned -ne $releaseVersion) { Write-Host "Data pack $($pack.Name) is pinned to $pinned; skipped for $releaseVersion."; continue }
        foreach ($action in @('validate', 'plan')) {
            & $cli datapack $action $pack.FullName | Out-Null
            if ($LASTEXITCODE -ne 0) { throw "Data pack $action failed: $($pack.Name)" }
        }
    }
}
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
