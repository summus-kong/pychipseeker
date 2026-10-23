# Builds pyBigWig 0.3.x as a native Windows .pyd with a MinGW toolchain.
# See README.md in this folder for background.
#
# Parameters
#   -EnvRoot       conda env directory that contains python.exe, include/, libs/ ...
#   -PyBigWigSrc   extracted pyBigWig sdist (contains pyBigWig.c, libBigWig/)
#   -ZlibLib       directory containing a MinGW-static libz.a
#   -IODeviceSource optional patched libBigWig io.c to copy over the sdist copy
param(
    [Parameter(Mandatory = $true)][string]$EnvRoot,
    [Parameter(Mandatory = $true)][string]$PyBigWigSrc,
    [Parameter(Mandatory = $true)][string]$ZlibLib,
    [string]$IODeviceSource
)
$ErrorActionPreference = 'Stop'

$mingwBin = Join-Path $EnvRoot 'Library\mingw-w64\bin'
$gcc      = Join-Path $mingwBin 'gcc.exe'
$dlltool  = Join-Path $mingwBin 'dlltool.exe'
$objdump  = Join-Path $mingwBin 'objdump.exe'
$python   = Join-Path $EnvRoot 'python.exe'
if (-not (Test-Path $gcc))     { throw "gcc not found: $gcc" }
if (-not (Test-Path $dlltool)) { throw "dlltool not found: $dlltool" }

# cc1.exe fails to load unless mingw-w64 bin is on PATH.
$env:PATH = "$mingwBin;$env:PATH"

$ver = & $python -c "import sys; print('%d%d' % sys.version_info[:2])"
$dllName = "python$ver.dll"
$outName = "pyBigWig.cp$ver-win_amd64.pyd"

Write-Output "1/5 generating def file from python$ver.dll export table"
$peDump   = Join-Path $env:TEMP 'pyBigWig_build_pe.txt'
$defFile  = Join-Path $env:TEMP 'pyBigWig_build_python.def'
& $objdump -p (Join-Path $EnvRoot $dllName) > $peDump
$lines = Get-Content $peDump
$start = -1
for ($i = 0; $i -lt $lines.Count; $i++) {
    if ($lines[$i] -match 'Ordinal/Name Pointer') { $start = $i + 1; break }
}
if ($start -lt 0) { throw 'export table header not found in objdump output' }
$names = foreach ($l in $lines[$start..($lines.Count - 1)]) {
    if ($l -match '^\s*\[\s*\d+\]\s+(\S+)\s*$') { $Matches[1] } else { break }
}
if (-not $names) { throw 'no exported names parsed' }
Set-Content -Encoding ascii $defFile (@('LIBRARY ' + $dllName, 'EXPORTS') + $names)

Write-Output ("2/5 dlltool import library (" + $names.Count + " symbols)")
$importLib = Join-Path $env:TEMP 'pyBigWig_build_libpython.dll.a'
& $dlltool -d $defFile -D $dllName -l $importLib -m i386:x86-64
if ($LASTEXITCODE -ne 0) { throw 'dlltool failed' }

if ($IODeviceSource) {
    Write-Output '3/5 applying patched io.c'
    Copy-Item $IODeviceSource (Join-Path $PyBigWigSrc 'libBigWig\io.c') -Force
} else {
    Write-Output '3/5 using libBigWig/io.c as-is (non-ASCII paths will not work)'
}

Write-Output "4/5 compiling with -DNOCURL"
$buildDir = Join-Path $PyBigWigSrc 'build_before_pyd'
New-Item -ItemType Directory -Force $buildDir | Out-Null
$defines = @('-DMS_WIN64', '-DMS_WIN32', '-DNOCURL', '-D_CRT_SECURE_NO_WARNINGS')
$includes = @(
    '-I' + (Join-Path $PyBigWigSrc 'libBigWig'),
    '-I' + (Join-Path $EnvRoot 'include'),
    '-I' + (Join-Path $EnvRoot 'Library\include')
)
$sources = @(Get-ChildItem $PyBigWigSrc -Filter 'pyBigWig.c') +
           (Get-ChildItem (Join-Path $PyBigWigSrc 'libBigWig') -Filter '*.c')
$objs = foreach ($f in $sources) {
    $o = Join-Path $buildDir ($f.BaseName + '.o')
    & $gcc -O2 @includes @defines -c $f.FullName -o $o
    if ($LASTEXITCODE -ne 0) { throw ("compile failed: " + $f.Name) }
    $o
}

Write-Output "5/5 linking $outName"
$out = Join-Path $PyBigWigSrc $outName
& $gcc @('-shared', '-o', $out) $objs @($importLib, '-L' + $ZlibLib, '-lz')
if ($LASTEXITCODE -ne 0) { throw 'link failed' }

Write-Output "OK: $out"
& $python -c "
import sys, shutil, os
dst = os.path.join(r'$EnvRoot', 'Lib', 'site-packages')
shutil.copy(r'$out', dst)
sys.path.insert(0, dst)
import pyBigWig
print('import OK from site-packages, module file:', pyBigWig.__file__)
"
if ($LASTEXITCODE -ne 0) { throw 'post-build import check failed' }
