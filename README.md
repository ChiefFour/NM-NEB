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
├── docs/
│   └── input-format.md         # full description of the "input" file
└── third_party/
    └── libdlfind/              # DL-FIND as a git submodule (see Installation)
```

## Dependencies

| Dependency | Purpose | Notes |
|---|---|---|
| **DL-FIND** | NEB/optimization kernel (`dl_find`) | Obtained as a git submodule from [libdlfind](https://github.com/digital-chemistry-laboratory/libdlfind) (contains the DL-FIND Fortran sources); see Installation |
| **BLAS / LAPACK** | linear-algebra backend required by the DL-FIND build | e.g. OpenBLAS; on Debian/Ubuntu: `sudo apt install libblas-dev liblapack-dev` |
| **Gaussian 16** | energies and gradients (`g16`, `formchk`) | single_point() shells out to `g16` and parses `.fchk`/`.log`; both executables must be in `PATH` |
| **gfortran** | Fortran compiler | any recent gfortran (GCC) works |

The three auxiliary components (`mode-align`, `mode-displacement`,
`z-interplot`) require **only** a Fortran compiler. Only `src/nmneb` links
against DL-FIND.

## Installation

### 1. Clone this repository with the DL-FIND submodule

```shell
git clone --recursive https://github.com/<your-account>/NM-NEB.git
# or, if already cloned without --recursive:
git submodule update --init --recursive
```

The DL-FIND Fortran sources are then available under `third_party/libdlfind/`.

### 2. Build and install DL-FIND

Build the DL-FIND shared library from the submodule and install it to a
prefix of your choice. This requires a Fortran compiler and BLAS/LAPACK:

```shell
cd third_party/libdlfind
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/.local
cmake --build build
cmake --install build
```

> If you prefer, you may instead obtain the original DL-FIND distribution
> from the [ChemShell website](https://www.chemshell.org/dl-find) (after
> registration) and install it in the usual way. Adjust the linker flags in
> the next step to match your installation.
>
> **Callback interface:** the DL-FIND callback routines in `src/nmneb/func.f90`
> (`dlf_get_params`, `dlf_get_gradient`) are matched to the API of the
> `libdlfind` submodule (including the `kiter` argument of
> `dlf_get_gradient`). If you link against a different DL-FIND distribution,
> check that its callback signatures match yours.

### 3. Compile the NM-NEB driver

```shell
cd src/nmneb
gfortran -O2 \
    main.f90 get_input.f90 fortran_dl_find.f90 \
    func.f90 subroutine.f90 molegeom_module.f90 \
    -o nmneb \
    -L$HOME/.local/lib -ldlfind     # adjust to your DL-FIND installation
```

### 4. Compile the auxiliary components

```shell
# mode matching
cd src/mode-align
gfortran -O2 main.f90 sub.f90 -o mode-align

# endpoint displacement (run in the P/ and R/ working directories)
cd src/mode-displacement
gfortran -O2 mode-R.f90 -o mode-R
gfortran -O2 mode-P.f90 -o mode-P

# Z-matrix interpolation
cd src/z-interplot
gfortran -O2 z-matrix-interpolated.f90 -o z-interplot

# Z-matrix -> Cartesian conversion (used for internal-coordinate NEB runs)
cd src/z-interplot/ZtoC
gfortran -O2 main.f90 zmat_to_cart.f90 -o zmat_to_cart
```

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
geometries and the DL-FIND parameters, then run:

```shell
./nmneb
```

During the optimization, each image geometry is written to
`single_point.gjf` and evaluated by Gaussian 16; energies and gradients are
read back from the Gaussian output. The NEB trajectory is written to
`neb_*.xyz`, and the optimized path is summarized in `nebinfo`. The linear
interpolation (LI-NEB) variant uses the same binary with the input file shown
in `examples/lineb/input`.

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
License (see [LICENSE](LICENSE)). DL-FIND, obtained via the
[libdlfind](https://github.com/digital-chemistry-laboratory/libdlfind)
submodule, is distributed under the GPL-3.0 / LGPL-3.0 licenses of its own
repository. Note the interface routines in `src/nmneb/func.f90` and
`src/nmneb/fortran_dl_find.f90` were adapted from the DL-FIND crystal
interface; if you redistribute these files together with DL-FIND, the
GPL-3.0 terms of DL-FIND apply to the combined work.
