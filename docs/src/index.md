# NLPModelsAlgencan.jl documentation

**NLPModelsAlgencan.jl** is a [NLPModels](https://github.com/JuliaSmoothOptimizers/NLPModels.jl)
interface to the [Algencan](https://www.ime.usp.br/~egbirgin/tango/codes.php)
nonlinear solver.

Algencan is a large scale high performance augmented Lagrangian solver written
by Ernesto Birgin and Mario Martínez. It has many special features like being
able to use the HSL library to speed up the sparse matrix linear algebra and
some smart acceleration strategies.

## Status

At this point this is alpha software. It will only work only with Julia 1.0 or later.

## Installation

We have not registered the NLPModelsAlgencan.jl package yet. Hence, to add it you
will need to give the full URL to the github directory. See examples below.

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
(@v1.x) pkg> add https://github.com/pjssilva/NLPModelsAlgencan.jl
julia> using NLPModelsAlgencan
```

### The easy way

Just type
```julia
(@v1.x) pkg> add https://github.com/pjssilva/NLPModelsAlgencan.jl
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
pointing to the directory where the library find resides *before* installing
NLPModelsAlgencan.jl.

```bash
export ALGENCAN_LIB_DIR=/path/where/algencan/libray/is
```

You can then proceed to install NLPModelsAlgencan.jl from the REPL
```julia
(@v1.x) pkg> add https://github.com/pjssilva/NLPModelsAlgencan.jl
julia> using NLPModelsAlgencan
```

## Hints to self compilation of Algencan with HSL libraries

This [wiki
page](https://github.com/pjssilva/NLPModelsAlgencan.jl/wiki/Compiling-HSL-Libraries-for-use-with-NLPModelsAlgencan.jl)
documents the steps I used to compile a version of `libalgencan.so` with HSL
support.

## Usage

See [First steps](@ref) for the fundamental usage. In addition, if there is a
need to configure the solver algorithm, you may also check [Optional parameters](@ref).
