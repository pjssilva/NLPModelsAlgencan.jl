# NLPModelsAlgencan.jl
[![Docs](https://img.shields.io/badge/docs-blue.svg)](https://pjssilva.github.io/NLPModelsAlgencan.jl/dev/)
![CI](https://github.com/pjssilva/NLPModelsAlgencan.jl/workflows/CI/badge.svg?branch=master)
[![codecov](https://codecov.io/gh/pjssilva/NLPModelsAlgencan.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/pjssilva/NLPModelsAlgencan.jl)

**NLPModelsAlgencan.jl** is a 
[NLPModels](https://github.com/JuliaSmoothOptimizers/NLPModels.jl)
interface to the [Algencan](https://www.ime.usp.br/~egbirgin/tango/codes.php)
nonlinear solver.

Algencan is a large scale high performance augmented Lagrangian solver written
by Ernesto Birgin and Mario Martínez. It has many special features like being
able to use the HSL library to speed up the sparse matrix linear algebra and
some smart acceleration strategies.

## How to cite

NLPModelsAlgencan.jl is based on Algencan that is a software from the [Tango
Project](https://www.ime.usp.br/~egbirgin/tango/). If you use this software in
your research, we kindly ask you to cite it according to [its
guidelines](https://www.ime.usp.br/~egbirgin/tango/license.php). In
particular, if you use Algencan we suggest citing:

1. R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt, "On
   Augmented Lagrangian methods with general lower-level constraints", SIAM
   Journal on Optimization 18, pp. 1286-1309, 2007.
1. R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt, "Augmented
   Lagrangian methods under the Constant Positive Linear Dependence constraint
   qualification", Mathematical Programming 111, pp. 5-32, 2008.

If your work uses Gencan, the suggested references are:

1. E. G. Birgin and J. M. Martínez, "Large-scale active-set box-constrained
   optimization method with spectral projected gradients", Computational
   Optimization and Applications 23, pp. 101-125, 2002.
1. M. Andretta, E. G. Birgin and J. M. Martínez, "Practical active-set
   Euclidian trust-region method with spectral projected gradients for
   bound-constrained minimization", Optimization 54, pp. 305-325, 2005.
1. E. G. Birgin and J. M. Martínez, "A box-constrained optimization algorithm
   with negative curvature directions and spectral projected gradients",
   Computing [Suppl] 15, pp. 49-60, 2001.

## Status

At this point this is beta software. It requires Julia 1.10 or later.

## Installation

```julia
(@v1.x) pkg> add NLPModelsAlgencan
julia> using NLPModelsAlgencan
```

That is all. The Algencan binary comes from
[`Algencan_jll`](https://github.com/JuliaBinaryWrappers/Algencan_jll.jl), so
nothing is compiled at installation time and you do not need a compiler, a
Fortran toolchain or a BLAS/Lapack development environment.

## Getting the most out of Algencan: using HSL

Algencan solves the sparse linear systems that arise in its subproblems much
faster when it can use an HSL linear solver. The binary in `Algencan_jll` is
built to do so: it looks for MA57 when it starts a solve and falls back to
Algencan's own truncated Newton solver when it does not find it. The choice is
made at run time, so one binary covers both cases and **nothing has to be
recompiled**.

HSL is proprietary and cannot be redistributed, so it is not included. To enable
MA57, obtain a licensed `libHSL` — it is free for academic use — from the [STFC
licences portal](https://licences.stfc.ac.uk/products/Software/HSL/LibHSL), and
install the `HSL_jll.jl` package that comes with it:

```julia
import Pkg
Pkg.develop(path = "/full/path/to/HSL_jll.jl")
```

Restart Julia and Algencan will use MA57 by itself. There is nothing to
configure and no compiler, Fortran toolchain or BLAS development environment is
involved. Algencan reports what it found in its own output:

```
 Available HSL subroutines = MA57
 lsslvr in TR           =            MA57/NONE
```

Only MA57 is supported at the moment. MA86 and MA97 are waiting on a fix to the
Fortran modules shipped in the public `HSL_jll`.

### Using your own Algencan build

If you would rather build Algencan yourself — for instance against a locally
patched MA57 — point the package at your shared library:

```julia
using NLPModelsAlgencan
set_algencan_library!("/path/to/libalgencan.so")
```

Then restart Julia. The path is stored as a preference of the active project, so
it applies to that project alone and survives restarts. Call
`set_algencan_library!(nothing)` to go back to the library from `Algencan_jll`.
The patches that build applies live in [`contrib/hsl`](contrib/hsl), and the
[wiki page on compiling HSL
libraries](https://github.com/pjssilva/NLPModelsAlgencan.jl/wiki/Compiling-HSL-Libraries-for-use-with-NLPModelsAlgencan.jl)
documents the process.

The `ALGENCAN_LIB_DIR` environment variable that earlier versions used still
works, so existing setups keep running, but it is deprecated and warns on load.
Prefer `set_algencan_library!`: environment variables are invisible to
precompilation, so changing one does not invalidate the cached module and can be
silently ignored.

## Contributing

See the [developer
notes](https://pjssilva.github.io/NLPModelsAlgencan.jl/dev/developer/) for how
the `Algencan_jll` binary is built and how to change its Yggdrasil recipe.


