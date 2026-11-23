<#
.SYNOPSIS
Plan and apply an explicit, metadata-checked cleanup of disposable files.
.DESCRIPTION
Policy roots must be reviewed by a person/agent before use. Age alone is not
evidence of disposability. No recursive deletion, reparse traversal, ACL changes,
or deletion of files added/changed after planning. Plans and results are private.
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory)][ValidateSet('Plan','Apply')][string]$Mode,
    [Parameter(Mandatory)][string]$PlanPath,
    [string]$PolicyPath,
    [string]$ResultPath
)
$ErrorActionPreference = 'Stop'

function Write-NewJson($Path, $Value) {
    $absolute = [IO.Path]::GetFullPath($Path)
    $parent = [IO.Path]::GetDirectoryName($absolute)
    [IO.Directory]::CreateDirectory($parent) | Out-Null
    $stream = [IO.File]::Open($absolute, 'CreateNew', 'Write', 'None')
    try {
        $writer = [IO.StreamWriter]::new($stream, [Text.UTF8Encoding]::new($false))
        try { $writer.Write(($Value | ConvertTo-Json -Depth 12)) }
        finally { $writer.Dispose() }
    } finally { $stream.Dispose() }
}

function Assert-PlainPath([string]$Path) {
    # Check every ancestor, including a junction that would hide outside the root.
    $item = Get-Item -LiteralPath $Path -Force
    $current = $item
    while ($null -ne $current) {
        if (($current.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) {
            throw "Reparse path refused: $($current.FullName)"
        }
        if ($current -is [IO.DirectoryInfo]) { $current = $current.Parent }
        else { $current = $current.Directory }
    }
    return $item
}

function Assert-Root($Rule) {
    $absolute = [IO.Path]::GetFullPath([string]$Rule.path).TrimEnd('\','/')
    if ($absolute.Length -le [IO.Path]::GetPathRoot($absolute).Length) {
        throw 'Volume-root cleanup is forbidden'
    }
    $item = Assert-PlainPath $absolute
    if (-not $item.PSIsContainer) { throw "Not a directory: $absolute" }
    if ($Rule.requireCMakeCache -and -not (Test-Path -LiteralPath (Join-Path $absolute 'CMakeCache.txt') -PathType Leaf)) {
        throw "Missing CMake cache: $absolute"
    }
    if (-not $Rule.reason) { throw 'Each root requires a reviewed reason' }
    return $item.FullName.TrimEnd('\','/')
}

function Assert-Child([string]$Root, [string]$Path) {
    $absolute = [IO.Path]::GetFullPath($Path)
    if (-not $absolute.StartsWith($Root + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Path outside reviewed root: $absolute"
    }
    return $absolute
}

if ($Mode -eq 'Plan') {
    if (-not $PolicyPath) { throw 'Plan requires PolicyPath' }
    $policy = Get-Content -LiteralPath $PolicyPath -Raw | ConvertFrom-Json
    $files = [Collections.Generic.List[object]]::new()
    $skips = [Collections.Generic.List[object]]::new()
    $seen = [Collections.Generic.HashSet[string]]::new([StringComparer]::OrdinalIgnoreCase)
    foreach ($rule in $policy.roots) {
        $rootPath = Assert-Root $rule
        $rule.path = $rootPath
        $cutoff = [DateTime]::UtcNow.AddDays(-[double]$rule.minAgeDays)
        $stack = [Collections.Generic.Stack[string]]::new()
        $stack.Push($rootPath)
        while ($stack.Count) {
            $directory = $stack.Pop()
            try { $entries = @(Get-ChildItem -LiteralPath $directory -Force) }
            catch { $skips.Add(@{path=$directory;reason=$_.Exception.Message}); continue }
            foreach ($entry in $entries) {
                if (($entry.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) {
                    $skips.Add(@{path=$entry.FullName;reason='reparse point'}); continue
                }
                if ($entry.PSIsContainer) { $stack.Push($entry.FullName); continue }
                if ($entry.LastWriteTimeUtc -gt $cutoff) { continue }
                if ($rule.extensions -and $entry.Extension.ToLowerInvariant() -notin $rule.extensions) { continue }
                $path = Assert-Child $rootPath $entry.FullName
                if ($seen.Add($path)) {
                    $files.Add(@{path=$path;root=$rootPath;length=$entry.Length;
                        modifiedTicks=$entry.LastWriteTimeUtc.Ticks;reason=$rule.reason})
                }
            }
        }
    }
    $total = ($files | Measure-Object -Property length -Sum).Sum
    Write-NewJson $PlanPath @{schema='dynsyslab.storage-cleanup/v1';createdUtc=[DateTime]::UtcNow.ToString('o');
        policy=$policy;files=@($files.ToArray());skipped=@($skips.ToArray());logicalBytes=$total}
    @{files=$files.Count;logicalBytes=$total;skipped=$skips.Count;plan=$PlanPath} | ConvertTo-Json
} else {
    if (-not $ResultPath) { throw 'Apply requires ResultPath' }
    if (Test-Path -LiteralPath $ResultPath) { throw 'Use a fresh result path' }
    $plan = Get-Content -LiteralPath $PlanPath -Raw | ConvertFrom-Json
    if ($plan.schema -ne 'dynsyslab.storage-cleanup/v1') { throw 'Unknown plan schema' }
    $roots = @{}
    foreach ($rule in $plan.policy.roots) { $roots[(Assert-Root $rule)] = $rule }
    # Validate the whole plan before the first mutation.
    foreach ($file in $plan.files) {
        if (-not $roots.ContainsKey([string]$file.root)) { throw 'Unreviewed root in plan' }
        Assert-Child $file.root $file.path | Out-Null
        $rule = $roots[[string]$file.root]
        if ($rule.extensions -and [IO.Path]::GetExtension($file.path).ToLowerInvariant() -notin $rule.extensions) {
            throw "Extension outside policy: $($file.path)"
        }
    }
    $results = [Collections.Generic.List[object]]::new()
    [long]$deletedBytes=0; [int]$deletedCount=0
    foreach ($file in $plan.files) {
        try {
            $item = Assert-PlainPath $file.path
            if ($item.PSIsContainer) { throw 'Expected a file' }
            if ($item.Length -ne $file.length -or $item.LastWriteTimeUtc.Ticks -ne $file.modifiedTicks) {
                throw 'File changed since planning'
            }
            Remove-Item -LiteralPath $item.FullName -Force
            $deletedBytes += $file.length; $deletedCount++
            $results.Add(@{path=$file.path;status='deleted';logicalBytes=$file.length})
        } catch { $results.Add(@{path=$file.path;status='skipped';reason=$_.Exception.Message}) }
    }
    Write-NewJson $ResultPath @{schema='dynsyslab.storage-cleanup-result/v1';createdUtc=[DateTime]::UtcNow.ToString('o');
        plan=[IO.Path]::GetFullPath($PlanPath);deletedFiles=$deletedCount;deletedLogicalBytes=$deletedBytes;
        skippedFiles=($results.Count-$deletedCount);files=@($results.ToArray())}
    @{deletedFiles=$deletedCount;deletedLogicalBytes=$deletedBytes;skippedFiles=($results.Count-$deletedCount);result=$ResultPath} | ConvertTo-Json
}
