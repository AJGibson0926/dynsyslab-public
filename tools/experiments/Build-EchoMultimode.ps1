$ErrorActionPreference = 'Stop'
$sourceRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $sourceRoot 'build/echo-multimode'
New-Item -ItemType Directory -Force -Path $buildDir | Out-Null
$vcvars = 'C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat'
$includes = "/I`"$sourceRoot/native/include`" /I`"$sourceRoot/native/tools/echo_input`""
$sources = "`"$PSScriptRoot/echo_multimode.cpp`" `"$PSScriptRoot/echo_multimode_control.cpp`" `"$PSScriptRoot/echo_multimode_audit.cpp`" `"$PSScriptRoot/echo_multimode_joint.cpp`" `"$sourceRoot/native/tools/echo_input/Experiment.cpp`""
$libs = 'dynsyslab_foundation.lib dynsyslab_koopman.lib dynsyslab_eigenfunctions.lib dynsyslab_history_koopman.lib dynsyslab_control.lib dynsyslab_control_maps.lib dynsyslab_numerics.lib'
Push-Location $buildDir
try {
    cmd /d /c "call `"$vcvars`" >nul && cl /nologo /EHsc /std:c++17 /O2 /MD $includes $sources /Fe:echo_multimode.exe /link /LIBPATH:`"$sourceRoot/build/native-suite-fast/lib`" $libs"
    if ($LASTEXITCODE -ne 0) { throw "Standalone experiment compile failed: $LASTEXITCODE" }
} finally { Pop-Location }
