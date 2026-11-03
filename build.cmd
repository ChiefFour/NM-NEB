@echo off
rem ============================================================
rem  Build all NM-NEB executables (Windows, MinGW-w64 gfortran)
rem
rem  Requirements:
rem    - gfortran (MinGW-w64) on PATH
rem    - BLAS/LAPACK libraries (adjust -lblas -llapack below if
rem      you use OpenBLAS: replace with -lopenblas)
rem
rem  The DL-FIND kernel sources are compiled in directly, so no
rem  separate DL-FIND installation is needed. The result is a
rem  self-contained nmneb.exe.
rem ============================================================
setlocal enabledelayedexpansion
set "ROOT=%~dp0"
set "DLF=%ROOT%third_party\libdlfind\src"

echo [NM-NEB] Building nmneb ...
cd /d "%ROOT%src\nmneb"
set "DLFSRC="
for %%f in ("%DLF%\dlf_*.f90") do set "DLFSRC=!DLFSRC! "%%f""
gfortran -O2 -o nmneb.exe main.f90 get_input.f90 fortran_dl_find.f90 func.f90 subroutine.f90 molegeom_module.f90 "%DLF%\dl-find.f90" !DLFSRC! -lblas -llapack
if errorlevel 1 exit /b 1
del /q *.o *.mod 2>nul
echo     -^> src\nmneb\nmneb.exe

echo [NM-NEB] Building mode-align ...
cd /d "%ROOT%src\mode-align"
gfortran -O2 -o mode-align.exe main.f90 sub.f90
if errorlevel 1 exit /b 1
echo     -^> src\mode-align\mode-align.exe

echo [NM-NEB] Building mode-R / mode-P ...
cd /d "%ROOT%src\mode-displacement"
gfortran -O2 -o mode-R.exe mode-R.f90
if errorlevel 1 exit /b 1
gfortran -O2 -o mode-P.exe mode-P.f90
if errorlevel 1 exit /b 1
echo     -^> src\mode-displacement\mode-R.exe, mode-P.exe

echo [NM-NEB] Building z-interplot ...
cd /d "%ROOT%src\z-interplot"
gfortran -O2 -o z-interplot.exe z-matrix-interpolated.f90
if errorlevel 1 exit /b 1
echo     -^> src\z-interplot\z-interplot.exe

echo [NM-NEB] Building zmat_to_cart (ZtoC) ...
cd /d "%ROOT%src\z-interplot\ZtoC"
gfortran -O2 -o zmat_to_cart.exe main.f90 zmat_to_cart.f90
if errorlevel 1 exit /b 1
echo     -^> src\z-interplot\ZtoC\zmat_to_cart.exe

echo [NM-NEB] All executables built.
