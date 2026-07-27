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
faster when it is linked against an HSL linear solver. The binary in
`Algencan_jll` is *not*: the HSL solvers are proprietary and cannot be
redistributed, so the packaged build uses Algencan's own routines instead.

If you work on large or very sparse problems, it is worth compiling Algencan
yourself against HSL. We only support MA57 at this point. Grab your copy from
[HSL](http://www.hsl.rl.ac.uk/catalogue/hsl_ma57.html) — it has a free academic
license — and follow the [wiki page on compiling HSL
libraries](https://github.com/pjssilva/NLPModelsAlgencan.jl/wiki/Compiling-HSL-Libraries-for-use-with-NLPModelsAlgencan.jl),
which documents the whole process. The patches it refers to live in
[`contrib/hsl`](contrib/hsl).

Once you have your own shared library, point the package at it:

```julia
using NLPModelsAlgencan
set_algencan_library!("/path/to/libalgencan.so")
```

Then restart Julia. The path is stored as a preference of the active project, so
it applies to that project alone and survives restarts. Call
`set_algencan_library!(nothing)` to go back to the library from `Algencan_jll`.

The `ALGENCAN_LIB_DIR` environment variable that earlier versions used still
works, so existing setups keep running, but it is deprecated and warns on load.
Prefer `set_algencan_library!`: environment variables are invisible to
precompilation, so changing one does not invalidate the cached module and can be
silently ignored.

## Contributing

See the [developer
notes](https://pjssilva.github.io/NLPModelsAlgencan.jl/dev/developer/) for how
the `Algencan_jll` binary is built and how to change its Yggdrasil recipe.


