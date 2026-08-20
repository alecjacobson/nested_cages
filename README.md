# Nested Cages 

This C++ project implements:

[Nested Cages] (http://www.cs.columbia.edu/cg/nested-cages/)  
ACM Transactions on Graphics, vol. 34, no. 6 (SIGGRAPH Asia 2015).  
Leonardo Sacht, Etienne Vouga and Alec Jacobson

> Get started with:
>
```bash
git clone https://github.com/alecjacobson/nested_cages.git
```

(No `--recursive` needed — all dependencies are fetched by CMake.)

[![Bunny teaser from "Nested Cages"](http://www.cs.columbia.edu/cg/nested-cages/bunny-shelf-teaser.jpg)](http://www.cs.columbia.edu/cg/nested-cages/)

## Compilation

This code has been tested on Linux and Mac OS X. In theory this should also
work on Windows.

To compile the command-line tool:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This will fetch and build all dependencies and produce the `nested_cages`
executable in `build/`.

### Dependencies

All third-party libraries are fetched automatically by CMake (via
`FetchContent`/`ExternalProject`) — there are no git submodules to initialize:

  * [libigl](https://github.com/libigl/libigl) `v2.6.0` (which in turn fetches
    CGAL and tetgen),
  * [Eigen](https://gitlab.com/libeigen/eigen) `3.4.0`,
  * [eltopo](https://github.com/leokollersacht/eltopo),
  * [collisiondetection](https://github.com/evouga/collisiondetection),
  * [MeshFix](https://github.com/alecjacobson/MeshFix-v2.1),
  * reference BLAS/LAPACK (built in isolation for eltopo), and
  * [nanobind](https://github.com/wjakob/nanobind) `v2.7.0` for the Python bindings.

The only system prerequisites are a C++17 compiler, a Fortran compiler
(for BLAS/LAPACK), CMake ≥ 3.16, and the Boost headers (used by CGAL). On
Debian/Ubuntu:

```bash
sudo apt-get install build-essential gfortran cmake libboost-dev
```

(GMP/MPFR are *not* required: libigl v2.6.0 builds CGAL in its GMP-free
`CGAL_DISABLE_GMP` configuration.)

By default the build fetches and compiles reference BLAS/LAPACK. To instead use
a system BLAS/LAPACK (e.g. OpenBLAS) for speed, configure with
`-DNESTED_CAGES_FETCH_BLAS=OFF`.

## Python bindings

The pipeline is also exposed as a Python module built with
[nanobind](https://github.com/wjakob/nanobind), mirroring the layout of the
[libigl Python bindings](https://github.com/libigl/libigl-python-bindings).
Installing with pip compiles the C++ code (and fetches all dependencies) via
scikit-build-core, so it needs the same system prerequisites listed above
(C++17 + Fortran compilers, CMake, Boost headers).

Install directly from GitHub — no need to clone first:

```bash
# latest on the default branch
pip install "git+https://github.com/alecjacobson/nested_cages.git"

# or pin a branch / tag / commit with the usual `@ref` suffix
pip install "git+https://github.com/alecjacobson/nested_cages.git@modernize-deps-and-python-bindings"
```

Or, from a local checkout of this repository:

```bash
pip install .            # regular install
pip install -e .         # editable/development install
```

The first build fetches and compiles all dependencies, so it can take several
minutes.

Then, from Python:

```python
import numpy as np
import nested_cages

# V: (#V,3) float64 vertices, F: (#F,3) int faces of the fine mesh
cages = nested_cages.nested_cages(
    V, F,
    num_faces=[1000, 500],      # one target face count per cage level
    quad_order=2,
    energy_expansion="None",
    energy_final="Volume",
    regular=[True, True],        # regular (vs adaptive) decimation per level
)
# cages is a list of (V_i, F_i) numpy arrays, one nested cage per level.
```

See `tests/regression.py` for a check that the Python bindings and the CLI
produce identical cages.

### Wheels / publishing to PyPI

`.github/workflows/wheels.yml` builds redistributable wheels with
[cibuildwheel](https://cibuildwheel.pypa.io) for Linux (`x86_64`, `aarch64`) and
macOS (`x86_64`, `arm64`), CPython 3.8–3.13, plus an sdist. It runs the full
matrix on tags/releases/manual dispatch and a single-version smoke build on
pushes and PRs. Each macOS arch is built natively (`macos-15-intel` for x86_64,
`macos-14` for arm64), because eltopo's Fortran BLAS/LAPACK cannot be
cross-compiled (Homebrew's `gfortran` only targets the host arch).

Note: Apple has retired the Intel architecture, and GitHub will drop macOS
`x86_64` runners after `macos-15-intel` is retired (~Fall 2027); at that point
Intel-Mac users would install from the sdist.

**Windows** is not built (eltopo links a Fortran BLAS/LAPACK, for which there is
no turn-key CI toolchain); Windows users install from the sdist.

Publishing to PyPI happens automatically when a GitHub Release is *published*,
via [trusted publishing](https://docs.pypi.org/trusted-publishers/) (OIDC, no
API token). To enable it once: on PyPI, add a trusted publisher for this project
pointing at `alecjacobson/nested_cages`, workflow `wheels.yml`, environment
`pypi`.


## Example usages

Help information

    ./nested_cages

Obtain 2 volume minimizing nested cages for `../gargo.off`: one regular with
1000 faces and the other regular with 500 faces. Output resulting cages to
`../test_1.off` and `../test_2.off`

    ./nested_cages ../gargo.off 2 1000r 500r None Volume ../test

The same as above, but outputs adaptive decimations (instead of regular)

    ./nested_cages ../gargo.off 2 1000 500 None Volume ../test

Obtain 2 nested cages for `../gargo.off` that minimize surface ARAP energy,
using as input decimations `../gargo_1000.off` and `gargo_500.off`

    ./nested_cages ../gargo.off 2 ../gargo_1000.off ../gargo_500.off SurfARAP None ../test


## Contact

If you have any comments or questions, please contact Leonardo Sacht by e-mail:
leo@mtm.ufsc.br
