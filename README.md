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

At this point this is beta software. It will only work with Julia LTS or later. 

## Prerequisites

The package downloads and installs Algencan upon installation. Therefore, you
need to have a minimal development environment installed. You need at least
`gcc`, `gfortran`, `make` and a development version of a BLAS/Lapack libraries
(for example `libopenblas-dev`). The BLAS/Lapack implementation is important to
get good performance. Use a high quality one like Openblas or Intel MKL.

## Installation

There are three main modes of installation, depending on how you want to compile
Algencan.

### The preferred way: using HSL

Obs: We only give support for MA57 at this point.

The preferred way to use Algencan is to link it against an HSL function to solve
sparse linear systems. To do this you need to grab your copy of MA57 from
[HSL](http://www.hsl.rl.ac.uk/catalogue/hsl_ma57.html). It has a free academic
license. You should receive a file named `hsl_ma57-5.2.0.tar.gz`.

All you need to do is to create an environment variable named
`MA57_SOURCE` pointing to this file *before* installing NLPModelsAlgencan.jl. For
example, if the file is located at the `/tmp` folder, in bash you would do:
```bash
export MA57_SOURCE=/tmp/hsl_ma57-5.2.0.tar.gz
```

After that just install NLPModelsAlgencan.jl from Julia's REPL and import it to force
pre-compilation.

```julia
(@v1.x) pkg> add NLPModelsAlgencan
julia> using NLPModelsAlgencan
```

### The easy way

Just type
```julia
(@v1.x) pkg> add NLPModelsAlgencan
julia> using NLPModelsAlgencan
```
in package mode in Julia's REPL.

This will download Algencan's code, compile it and make it available to the
NLPModelsAlgencan.jl package. **However, there is a major caveat here. The
Algencan solver will be compiled without any HSL support. This will have a major
negative impact on its behavior and performance. You should use HSL whenever you
have access to it.**

### Precompiled `libalgencan.so`

If you have your own copy of `libalgencan.so` (note that you need a shared
library), you can use it. Just create an environment variable like below
pointing to the directory where the library resides *before* installing
NLPModelsAlgencan.jl.

```bash
export ALGENCAN_LIB_DIR=/path/where/algencan/libray/is
```

You can then proceed to install NLPModelsAlgencan.jl from the REPL
```julia
(@v1.x) pkg> add NLPModelsAlgencan
julia> using NLPModelsAlgencan
```

## Hints for self compiling Algencan with HSL libraries

This [wiki
page](https://github.com/pjssilva/NLPModelsAlgencan.jl/wiki/Compiling-HSL-Libraries-for-use-with-NLPModelsAlgencan.jl)
documents the steps I used to compile a version of `libalgencan.so` with HSL
support.


