# NLPModelsAlgencan.jl documentation

**NLPModelsAlgencan.jl** is a [NLPModels](https://github.com/JuliaSmoothOptimizers/NLPModels.jl)
interface to the [Algencan](https://www.ime.usp.br/~egbirgin/tango/codes.php)
nonlinear solver.

Algencan is a large scale high performance augmented Lagrangian solver written
by Ernesto Birgin and Mario Martínez. It has many special features like being
able to use the HSL library to speed up the sparse matrix linear algebra and
some smart acceleration strategies.

## Status

At this point this is beta software. It requires Julia 1.10 or later.

## Installation

```julia
(@v1.x) pkg> add NLPModelsAlgencan
julia> using NLPModelsAlgencan
```

That is all. The Algencan binary comes from `Algencan_jll`, so nothing is
compiled at installation time and you do not need a compiler, a Fortran
toolchain or a BLAS/LAPACK development environment.

### Getting the most out of Algencan: using HSL

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
`contrib/hsl` in this repository.

Once you have your own shared library, point the package at it:

```julia
using NLPModelsAlgencan
set_algencan_library!("/path/to/libalgencan.so")
```

Then restart Julia. The path is stored as a preference of the active project, so
it applies to that project alone and survives restarts. Call
`set_algencan_library!(nothing)` to go back to the library from `Algencan_jll`.

!!! note "The ALGENCAN\\_LIB\\_DIR environment variable"
    Earlier versions selected a custom library through the `ALGENCAN_LIB_DIR`
    environment variable. It still works, so existing setups keep running, but
    it is deprecated and warns on load. Prefer `set_algencan_library!`:
    environment variables are invisible to precompilation, so changing one does
    not invalidate the cached module and can be silently ignored.

## Usage

See [First steps](@ref) for the fundamental usage. In addition, if there is a
need to configure the solver algorithm, you may also check [Optional parameters](@ref).
