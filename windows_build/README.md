# Building pyBigWig on Windows from source

There are no prebuilt native Windows pyBigWig wheels on PyPI (manylinux only)
and bioconda only ships linux-64/osx-64/osx-arm64/linux-aarch64.  The
extension does build on Windows with a MinGW toolchain if you work around two
issues:

1. **Python data-symbol imports fail.**  MinGW `ld` (e.g. binutils 2.25 from
   `m2w64-gcc`) cannot resolve `__imp_*` data symbols (`PyExc_RuntimeError`,
   `_Py_NoneStruct`, ...) from the MSVC-format `python311.lib` import library,
   even though `nm` finds them.  Function imports link fine, data imports do
   not.  **Workaround:** generate a MinGW-native import library from
   `python311.dll` with `dlltool` (see `build_pybigwig_windows.ps1`; it parses
   the DLL export table with `objdump -p` and feeds a `*.def` file to
   `dlltool`).
2. **Non-ASCII local paths fail to open.**  libBigWig's `fopen()` interprets
   Python's UTF-8 file names through the ANSI codepage on Windows.  The patch
   in `io.c` (this folder) converts UTF-8 -> UTF-16 with
   `MultiByteToWideChar` and calls `_wfopen()`.  Copy the patched `io.c`
   over `libBigWig/io.c` before compiling.

## Reproduced environment (2026-10)

- conda env `pybw`: `conda create -y -n pybw -c conda-forge --override-channels python=3.11 m2w64-gcc zlib`
- gcc from `envs/pybw/Library/mingw-w64/bin/gcc.exe` (5.3.0, MSYS2 build);
  **its `bin` directory must be on `PATH` or internal `cc1.exe` fails to load**
- zlib static build from the upstream `zlib-1.3.1` tarball (`-O2`, `ar rcs`)
- pyBigWig 0.3.26 sdist, compiled with:
  `-DMS_WIN64 -DMS_WIN32 -DNOCURL -D_CRT_SECURE_NO_WARNINGS`
  (`NOCURL` = no remote http/ftp access, local files only; numpy support
  (`WITHNUMPY`) omitted by choice)

Result: `pyBigWig.cp311-win_amd64.pyd` imports and reads/writes bigWig files
correctly (validated with upstream `pyBigWigTest/test.bw`).

## Build script

`build_pybigwig_windows.ps1` automates def-file generation, import-library
creation, compilation and linking.  Example:

```powershell
powershell -ExecutionPolicy Bypass -File build_pybigwig_windows.ps1 `
    -EnvRoot <path-to-conda-env-with-m2w64-gcc> `
    -PyBigWigSrc <extracted-pybigwig-sdist> `
    -ZlibLib   <dir-containing-libz.a>
```

The produced `pyd` must match the Python ABI (name suffix
`cp311-win_amd64`); copy it into the env's `Lib\site-packages`.
