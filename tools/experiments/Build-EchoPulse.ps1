$ErrorActionPreference = 'Stop'
$sourceRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$buildDir = Join-Path $sourceRoot 'build/echo-pulse'
New-Item -ItemType Directory -Force -Path $buildDir | Out-Null
$vcvars = 'C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat'
$includes = "/I`"$sourceRoot/native/include`" /I`"$sourceRoot/native/tools/echo_input`""
$files = @('echo_pulse.cpp','echo_pulse_acquire.cpp','echo_pulse_analyze.cpp','echo_pulse_fit.cpp','echo_pulse_dense.cpp','echo_multimode.cpp','echo_multimode_control.cpp','echo_multimode_joint.cpp')
$sources = ($files | ForEach-Object { "`"$PSScriptRoot/$_`"" }) -join ' '
$sources += " `"$sourceRoot/native/tools/echo_input/Experiment.cpp`""
$libs = 'dynsyslab_foundation.lib dynsyslab_koopman.lib dynsyslab_eigenfunctions.lib dynsyslab_history_koopman.lib dynsyslab_control.lib dynsyslab_control_maps.lib dynsyslab_numerics.lib'
Push-Location $buildDir
try {
    cmd /d /c "call `"$vcvars`" >nul && cl /nologo /EHsc /std:c++17 /O2 /MD /DDSL_ECHO_NO_MAIN $includes $sources /Fe:echo_pulse.exe /link /LIBPATH:`"$sourceRoot/build/native-suite-fast/lib`" $libs"
    if ($LASTEXITCODE -ne 0) { throw "Pulse experiment compile failed: $LASTEXITCODE" }
} finally { Pop-Location }
