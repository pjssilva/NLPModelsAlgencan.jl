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

MA86 and MA97 are found the same way, and a licensed library reports all three:

```
 Available HSL subroutines = MA57 MA86 MA97
```

How the run-time switch works is described in the [developer notes](developer.md).

#### Choosing a linear solver

MA57 is the default and needs no configuration. MA86 and MA97 are alternatives
for the Newton line search and for the acceleration process, selected with a
keyword argument:

```julia
algencan(nlp; NEWTON_LINE_SEARCH_INNER_SOLVER = "MA86")
```

The value is the solver, optionally followed by a scaling, as in `"MA86 MC64"`.
The two systems are set independently, and the acceleration process does not
follow the Newton setting:

```julia
algencan(nlp; NEWTON_LINE_SEARCH_INNER_SOLVER = "MA86 MC64",
              LINEAR_SYSTEMS_SOLVER_IN_ACCELERATION_PROCESS = "MA97 MC64")
```

Algencan reports what it settled on, which is worth reading back: an
unrecognised keyword is passed through as a specification file line and quietly
ignored, so a misspelling looks like nothing happened.

```
 lsslvr in TR           =            MA57/NONE
 lsslvr in NW           =            MA86/MC64
 lsslvr in ACCPROC      =            MA97/MC64
```

The trust region accepts MA57 alone, whatever is asked for. Do not expect MA86
or MA97 to be faster in general, but they are worth trying on a problem where
MA57 struggles.

### Building Algencan yourself

Fully supported. Point the package at your own shared library:

```julia
using NLPModelsAlgencan
set_algencan_library!("/path/to/libalgencan.so")
```

Then restart Julia. The path is stored as a preference of the active project, so
it applies to that project alone and survives restarts. Call
`set_algencan_library!(nothing)` to go back to the library from `Algencan_jll`.

We must recall that such a library has to be built using a patched version of
MA57 from HSL as suggested in the original Algencan installation instructions.
This would avoid an extra (sparse) matrix times vector operation that the
patched version of libAlgencan_jll uses to avoid touching the HSL code.

The patches are in `contrib/hsl` in this repository, and the [wiki page on
compiling HSL
libraries](https://github.com/pjssilva/NLPModelsAlgencan.jl/wiki/Compiling-HSL-Libraries-for-use-with-NLPModelsAlgencan.jl)
documents the process.

!!! note "The ALGENCAN\\_LIB\\_DIR environment variable"
    Earlier versions selected a custom library through the `ALGENCAN_LIB_DIR`
    environment variable. It still works, so existing setups keep running, but
    it is deprecated and warns on load. Prefer `set_algencan_library!`:
    environment variables are invisible to precompilation, so changing one does
    not invalidate the cached module and can be silently ignored.

## Usage

See [First steps](@ref) for the fundamental usage. In addition, if there is a
need to configure the solver algorithm, you may also check [Optional parameters](@ref).
