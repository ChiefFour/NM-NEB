# NM-NEB: Normal-Mode-Based Nudged Band Method

A normal-mode-guided nudged elastic band (NEB) implementation that retrieves
potential energy surface (PES) curvature information from the local Hessians of
the reactant and product minima to construct physically reasonable initial
paths and improve transition state (TS) search.

**Reference:** S. Liang, Z. Zhou, Z. Qu, *J. Phys. Chem. A*, 2026.
DOI: 10.1021/acs.jpca.6c02928.

## Overview

Conventional NEB methods interpolate linearly (in Cartesian or internal
coordinates) between the reactant and product minima, which often yields
unphysical initial paths when the endpoints are far apart. NM-NEB instead uses
the *local curvature information* contained in the endpoint Hessians:

1. **Frequency calculations** — Run Hessian/frequency calculations on the
   reactant (R) and product (P) minima with Gaussian 16. This produces
   `reactant.log` / `product.log` containing the normal-mode eigenvectors and
   vibrational frequencies.
2. **Mode matching** (`src/mode-align`) — Extract the modes, compute the
   per-mode displacement matrix **D** and a Franck–Condon-type overlap
   matrix **F**, and find the optimal one-to-one reactant–product mode
   assignment by a backtracking search (equivalent to the Hungarian
   assignment used in the paper).
3. **Endpoint displacement** (`src/mode-displacement`) — Project the full
   R→P Cartesian displacement onto the matched normal modes, scale the
   projections by empirical coefficients, and construct improved endpoints
   R′ and P′ that remain within the harmonic-approximation regime of the
   endpoint minima (Eqs. (16)–(18) of the paper).
4. **NEB optimization** (`src/nmneb`) — Build the initial path from the
   improved endpoints and optimize it with the NEB algorithm of
   [DL-FIND](https://www.chemshell.org/dl-find). Energies and gradients are
   supplied on the fly by Gaussian 16 single-point calculations.

An optional Z-matrix internal-coordinate interpolator
(`src/z-interplot`) is also provided for generating NEB initial paths in
internal coordinates.

## Repository structure

```
NM-NEB/
├── src/
│   ├── nmneb/                  # NEB driver: DL-FIND kernel + Gaussian 16 interface
│   │   ├── main.f90            #   entry point (get_input -> fortranb_dl_find)
│   │   ├── get_input.f90       #   reads the "input" file (geometry + DL-FIND params)
│   │   ├── fortran_dl_find.f90 #   assembles nvar/nvar2/nspec and calls dl_find
│   │   ├── func.f90            #   DL-FIND callbacks (dlf_get_params / dlf_get_gradient)
│   │   ├── subroutine.f90      #   single_point() -> Gaussian 16; atomic utilities
│   │   └── molegeom_module.f90 #   global molecular geometry data
│   ├── mode-align/             # reactant/product normal-mode matching
│   │   ├── main.f90            #   displacement D, FC overlap F, backtracking assignment
│   │   └── sub.f90             #   supporting routines
│   ├── mode-displacement/      # construction of improved endpoints R' and P'
│   │   ├── mode-R.f90          #   R-side variant (TS guess built from R)
│   │   └── mode-P.f90          #   P-side variant (TS guess built from P)
│   └── z-interplot/            # Z-matrix internal-coordinate interpolation
│       ├── z-matrix-interpolated.f90
│       ├── test.f90
│       └── ZtoC/               # zmat_to_cart: Z-matrix -> Cartesian converter
├── examples/                   # input templates for each component
│   ├── nmneb/                  #   input / input-z (NEB driver)
│   ├── lineb/                  #   input for linear-interpolation NEB runs
│   ├── mode-align/             #   reactant.com / product.com (Gaussian inputs)
│   ├── mode-displacement/      #   reactant.com / product.com
│   └── z-interplot/            #   info_connection, initio.dat, final.dat, ...
├── build.sh / build.cmd        # one-shot build scripts (all executables)
├── run.sh / run.cmd            # launchers with Gaussian environment check
├── docs/
│   └── input-format.md         # full description of the "input" file
└── third_party/
    └── libdlfind/              # DL-FIND as a git submodule (see Installation)
```

## Dependencies

| Dependency | Purpose | Notes |
|---|---|---|
| **DL-FIND** | NEB/optimization kernel (`dl_find`) | Sources ship as a git submodule ([libdlfind](https://github.com/digital-chemistry-laboratory/libdlfind)) and are compiled **directly into** the `nmneb` binary; **no separate DL-FIND installation is needed** |
| **BLAS / LAPACK** | linear-algebra backend required by the DL-FIND kernel | e.g. OpenBLAS; on Debian/Ubuntu: `sudo apt install libblas-dev liblapack-dev` |
| **Gaussian 16** | energies and gradients (`g16`, `formchk`) | **NOT bundled with this repository.** Each user must install and license Gaussian 16 themselves and make `g16` / `formchk` available (see [Environment requirements](#environment-requirements)) |
| **gfortran** | Fortran compiler | any recent gfortran (GCC) works |

The auxiliary components (`mode-align`, `mode-displacement`,
`z-interplot`, `zmat_to_cart`) require **only** a Fortran compiler.

## Environment requirements

NM-NEB is a method driver: it does not bundle, provide, or interface with any
quantum chemistry package.

**Gaussian 16 is strictly required at run time.** During the NEB optimization
the driver calls `g16` for every image single point and `formchk` to read the
gradient back from the checkpoint file. Without a working Gaussian 16
installation the program starts but cannot perform any calculation. Gaussian
16 is commercial software; it must be obtained, licensed and installed by
each user.

The driver locates the Gaussian executables in the following order:

1. Environment variable `GAUSS_EXE` (path to the `g16` executable).
2. Environment variable `FORMCHK_EXE` (path to the `formchk` executable).
3. `g16` / `formchk` found in `PATH` (the default Gaussian setup).

```bash
# Linux / macOS
export GAUSS_EXE=/opt/g16/g16
export FORMCHK_EXE=/opt/g16/formchk
```

```bat
rem Windows (cmd)
set GAUSS_EXE=C:\Gaussian\g16\g16.exe
set FORMCHK_EXE=C:\Gaussian\g16\formchk.exe
```

Paths containing spaces are handled automatically by the driver. The
`run.sh` / `run.cmd` launchers perform the same detection and refuse to start
when Gaussian is missing.

**Platform support:** the code is plain Fortran 90 and compiles with gfortran
on Linux, macOS and Windows (MinGW-w64). No pre-built binaries are provided;
see Installation below.

## Installation

### 1. Clone this repository (with the DL-FIND submodule)

```shell
git clone --recursive https://github.com/ChiefFour/NM-NEB.git
# or, if already cloned without --recursive:
git submodule update --init --recursive
```

The DL-FIND Fortran sources are then available under `third_party/libdlfind/`.

### 2. One-shot build

The DL-FIND kernel sources are compiled **directly into** the `nmneb` binary,
so no separate DL-FIND installation (cmake, shared library, or ChemShell
registration) is needed. Only gfortran and BLAS/LAPACK are required.

Linux / macOS / WSL:

```shell
./build.sh
```

Windows (MinGW-w64 gfortran on `PATH`):

```bat
build.cmd
```

Each script builds all executables:

| Executable | Location |
|---|---|
| `nmneb` | `src/nmneb/` |
| `mode-align` | `src/mode-align/` |
| `mode-R`, `mode-P` | `src/mode-displacement/` |
| `z-interplot` | `src/z-interplot/` |
| `zmat_to_cart` | `src/z-interplot/ZtoC/` |

> **What gets compiled in:** `build.sh` / `build.cmd` compile `dl-find.f90`
> plus the `dlf_*.f90` kernel modules from `third_party/libdlfind/src`. The
> files `api.f90`, `main.f90`, `dlfind_crystal.f90`, `draw.f90`,
> `mod_globals.f90` and `mod_api.f90` are the Python-binding /
> standalone-entry files of the libdlfind repository and are intentionally
> excluded, because they define their own `dlf_get_params` /
> `dlf_get_gradient` / `dl_find` symbols. The DL-FIND callbacks in
> `src/nmneb/func.f90` (including the `kiter` argument) are matched to the
> compiled kernel.
>
> **Manual compilation:** if you prefer to compile by hand, do the same
> steps directly: `cd src/nmneb`, compile the six NM-NEB sources together
> with `$DLF/dl-find.f90` and `$DLF/dlf_*.f90`, then link with
> `-lblas -llapack`. On Windows, add `-o nmneb.exe` and use MinGW-w64
> gfortran; adjust the BLAS/LAPACK flags if you use OpenBLAS
> (`-lopenblas`).

### 3. Alternative: link against an existing DL-FIND installation

If you already have DL-FIND installed (e.g. the official distribution from
the [ChemShell website](https://www.chemshell.org/dl-find)), you can compile
`nmneb` against it instead of compiling the kernel in. In that case make sure
the callback signatures in `src/nmneb/func.f90` (in particular the `kiter`
argument of `dlf_get_gradient`) match your DL-FIND version.

## Usage

The complete NM-NEB workflow is described below. Detailed input keywords are
documented in [docs/input-format.md](docs/input-format.md); an example `input`
file is provided in `examples/nmneb/`.

### Step 1 — Frequency calculations (Gaussian 16)

For the reactant and product minima, run a Gaussian frequency calculation
(e.g. `#force b3lyp/6-31G* em=GD3 nosymm`). Keep the `.log` files; they
contain the normal-mode vectors and frequencies consumed by the next two
steps.

### Step 2 — Mode matching (`src/mode-align`)

Place `reactant.log` and `product.log` in the working directory and run:

```shell
./mode-align
```

The program reads the equilibrium geometries and modes, computes the
displacement matrix **D** and the FC-type overlap matrix **F**, and prints the
optimal one-to-one mode assignment (row → column).

### Step 3 — Endpoint displacement (`src/mode-displacement`)

In the `P/` and `R/` working directories (containing `reactant.com`,
`product.com` and the mode data extracted from the `.log` files), run:

```shell
# in R/ : construct R' and the TS guess starting from the reactant
./mode-R
# in P/ : construct P' and the TS guess starting from the product
./mode-P
```

The programs grep the normal modes from `reactant.log`/`product.log`, project
the R→P displacement onto the matched modes, apply the empirical scaling
coefficients and write the displaced endpoint geometry used as NEB input.

> **Note:** the TS-guess displacement coefficients are currently hard-coded
> for the benchmark system (Natom = 8) used in the paper. Adapt
> `mode-R.f90` / `mode-P.f90` (or generalize them) for other systems.

### Step 4 — NEB optimization (`src/nmneb`)

Write an `input` file (see `examples/nmneb/input`) containing the two endpoint
geometries and the DL-FIND parameters, then run the launcher, which checks the
Gaussian environment first and then starts the driver:

```shell
./run.sh          # Linux / macOS / WSL
```

```bat
run.cmd           rem Windows
```

Equivalently, run the driver directly (`./nmneb` / `nmneb.exe`) once
`g16` / `formchk` are configured. During the optimization, each image geometry
is written to `single_point.gjf` and evaluated by Gaussian 16; energies and
gradients are read back from the Gaussian output. The NEB trajectory is
written to `neb_*.xyz`, and the optimized path is summarized in `nebinfo`.
The linear interpolation (LI-NEB) variant uses the same binary with the input
file shown in `examples/lineb/input`.

## Examples

| Directory | Contents |
|---|---|
| `examples/nmneb/input` | NEB driver input (Cartesian, icoord=100) with displaced endpoints |
| `examples/nmneb/input-z` | NEB driver input for internal-coordinate (Z-matrix) runs |
| `examples/lineb/input` | NEB input for the linear-interpolation variant (icoord=120) |
| `examples/mode-align/*.com` | Gaussian inputs that generate the `.log` files for mode matching |
| `examples/mode-displacement/*.com` | Gaussian inputs for the endpoint frequency calculations |
| `examples/z-interplot/` | Z-matrix connection table and parameter files for the interpolator |

## Citation

If you use NM-NEB in your work, please cite:

- **NM-NEB paper:** S. Liang, Z. Zhou, Z. Qu, *J. Phys. Chem. A*, 2026.
  DOI: 10.1021/acs.jpca.6c02928.
- **DL-FIND:** J. Kästner, J. M. Carr, T. W. Keal, W. Thiel, A. Wander,
  P. Sherwood, *J. Phys. Chem. A*, 2009, **113**, 11856–11865.
- **Gaussian 16:** M. J. Frisch et al., Gaussian 16, Revision C.01, Gaussian,
  Inc., Wallingford CT, 2016.

## License

The NM-NEB source code in this repository is distributed under the MIT
License (see [LICENSE](LICENSE)), except for the DL-FIND-derived parts.
DL-FIND, obtained via the
[libdlfind](https://github.com/digital-chemistry-laboratory/libdlfind)
submodule, is distributed under the GPL-3.0 / LGPL-3.0 licenses of its own
repository. The callback routines in `src/nmneb/func.f90` and
`src/nmneb/fortran_dl_find.f90` were adapted from the DL-FIND crystal
interface, and `build.sh` / `build.cmd` compile the DL-FIND kernel directly
into the `nmneb` binary. When you build and redistribute that combined
binary, the GPL-3.0 / LGPL-3.0 terms of DL-FIND apply to the combined work;
the complete source needed to rebuild it is this repository.
