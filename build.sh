#!/bin/sh
# ============================================================
#  Build all NM-NEB executables (Linux / macOS / WSL)
#
#  Requirements: gfortran, BLAS/LAPACK (development libraries).
#
#  The DL-FIND kernel sources are compiled in directly, so no
#  separate DL-FIND installation is needed. The result is a
#  self-contained nmneb binary.
# ============================================================
set -e
cd "$(dirname "$0")"
ROOT=$(pwd)
DLF="$ROOT/third_party/libdlfind/src"

echo "[NM-NEB] Building nmneb ..."
cd "$ROOT/src/nmneb"
gfortran -O2 -c main.f90 get_input.f90 fortran_dl_find.f90 func.f90 subroutine.f90 molegeom_module.f90
# DL-FIND kernel: dl-find.f90 plus all dlf_*.f90 modules.
# The files api.f90, main.f90, dlfind_crystal.f90, draw.f90, mod_globals.f90
# and mod_api.f90 are the Python-binding / standalone-entry files of the
# libdlfind repository and are intentionally NOT compiled in, because they
# provide their own dlf_get_params / dlf_get_gradient / dl_find symbols.
gfortran -O2 -c "$DLF/dl-find.f90" "$DLF"/dlf_*.f90
gfortran -O2 -o nmneb *.o -lblas -llapack
rm -f *.o *.mod
echo "    -> src/nmneb/nmneb"

echo "[NM-NEB] Building mode-align ..."
cd "$ROOT/src/mode-align"
gfortran -O2 -o mode-align main.f90 sub.f90
rm -f *.o *.mod
echo "    -> src/mode-align/mode-align"

echo "[NM-NEB] Building mode-R / mode-P ..."
cd "$ROOT/src/mode-displacement"
gfortran -O2 -o mode-R mode-R.f90
gfortran -O2 -o mode-P mode-P.f90
rm -f *.o *.mod
echo "    -> src/mode-displacement/mode-R, mode-P"

echo "[NM-NEB] Building z-interplot ..."
cd "$ROOT/src/z-interplot"
gfortran -O2 -o z-interplot z-matrix-interpolated.f90
rm -f *.o *.mod
echo "    -> src/z-interplot/z-interplot"

echo "[NM-NEB] Building zmat_to_cart (ZtoC) ..."
cd "$ROOT/src/z-interplot/ZtoC"
gfortran -O2 -o zmat_to_cart main.f90 zmat_to_cart.f90
rm -f *.o *.mod
echo "    -> src/z-interplot/ZtoC/zmat_to_cart"

echo "[NM-NEB] All executables built."
