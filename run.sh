#!/bin/sh
# ============================================================
#  NM-NEB launcher (Linux / macOS / WSL)
#
#  Requires:
#    - a built nmneb binary in the current directory (build.sh)
#    - Gaussian 16 installed by the USER (not bundled here)
#
#  The Gaussian executables are located via the environment
#  variables GAUSS_EXE / FORMCHK_EXE; if they are not set,
#  g16 / formchk are looked up in PATH.
# ============================================================
set -u

GAUSS=""
if [ -n "${GAUSS_EXE:-}" ]; then
  GAUSS="$GAUSS_EXE"
elif command -v g16 >/dev/null 2>&1; then
  GAUSS="g16"
fi
if [ -z "$GAUSS" ]; then
  echo "[NM-NEB] ERROR: Gaussian 16 executable not found." >&2
  echo "         Install Gaussian 16 yourself and either add g16 and formchk" >&2
  echo "         to PATH, or set GAUSS_EXE / FORMCHK_EXE to their full paths." >&2
  exit 1
fi

FCHK=""
if [ -n "${FORMCHK_EXE:-}" ]; then
  FCHK="$FORMCHK_EXE"
elif command -v formchk >/dev/null 2>&1; then
  FCHK="formchk"
fi
if [ -z "$FCHK" ]; then
  echo "[NM-NEB] ERROR: formchk executable not found. Install Gaussian 16" >&2
  echo "         or set FORMCHK_EXE to its full path." >&2
  exit 1
fi

if [ ! -x ./nmneb ]; then
  echo "[NM-NEB] ERROR: nmneb binary not found in the current directory." >&2
  echo "         Build it first with build.sh (see README, Installation)." >&2
  exit 1
fi

echo "[NM-NEB] Gaussian executable : $GAUSS"
echo "[NM-NEB] formchk executable : $FCHK"
echo "[NM-NEB] Running nmneb ..."
exec ./nmneb
