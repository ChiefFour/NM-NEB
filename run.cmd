@echo off
rem ============================================================
rem  NM-NEB launcher (Windows)
rem
rem  Requires:
rem    - a built nmneb.exe in the current directory (build.cmd)
rem    - Gaussian 16 installed by the USER (not bundled here)
rem
rem  The Gaussian executables are located via the environment
rem  variables GAUSS_EXE / FORMCHK_EXE; if they are not set,
rem  g16 / formchk are looked up in PATH.
rem ============================================================
setlocal

set "GAUSS="
if defined GAUSS_EXE (
  set "GAUSS=%GAUSS_EXE%"
) else (
  where g16 >nul 2>nul
  if not errorlevel 1 set "GAUSS=g16"
)
if not defined GAUSS (
  echo [NM-NEB] ERROR: Gaussian 16 executable not found.
  echo           Install Gaussian 16 yourself and either add g16 and formchk
  echo           to PATH, or set GAUSS_EXE / FORMCHK_EXE to their full paths.
  exit /b 1
)

set "FCHK="
if defined FORMCHK_EXE (
  set "FCHK=%FORMCHK_EXE%"
) else (
  where formchk >nul 2>nul
  if not errorlevel 1 set "FCHK=formchk"
)
if not defined FCHK (
  echo [NM-NEB] ERROR: formchk executable not found. Install Gaussian 16
  echo           or set FORMCHK_EXE to its full path.
  exit /b 1
)

if not exist "nmneb.exe" (
  echo [NM-NEB] ERROR: nmneb.exe not found in the current directory.
  echo           Build it first with build.cmd (see README, Installation).
  exit /b 1
)

echo [NM-NEB] Gaussian executable : %GAUSS%
echo [NM-NEB] formchk executable : %FCHK%
echo [NM-NEB] Running nmneb ...
nmneb.exe
exit /b %ERRORLEVEL%
