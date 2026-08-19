# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A [NLPModels](https://github.com/JuliaSmoothOptimizers/NLPModels.jl) interface to
Algencan, an augmented Lagrangian nonlinear solver written in Fortran. The whole
Julia side is one file, `src/NLPModelsAlgencan.jl` (~800 lines); everything else
is tests, docs and maintainer material. The package compiles nothing: the binary
comes from `Algencan_jll`, built by a Yggdrasil recipe that lives outside this
repository.

## My preferences

After dealing with Claude Code I have noticed two biases that I consider
really annoying:

1. Claude likes to document the history of a change. Something along the
   lines: "we changed this to solve that bug" or "the new implementation based
   on X is faster than the older implementation that used Y". I don't like
   this. History belongs to `git`. I commit that comment may refer to the old
   implementation, to the old bug, but *documentation and code comments should
   always reflect the current, actual, code*. It should explain what, why, and
   how. But it should not refer to the past, to how it was. If, when making a
   change or fixing a bug, we learn about a pitfall that should be avoided in
   the future, probably the right way to put it is in `docs/src/developer.md`.

2. Claude likes to add comments to git commits that are way too long. They
   look like a prose that tells the story of the commit. I prefer shorter,
   direct comments that describe what changed. Avoid long comments. 

## Commands

```bash
# Full test suite
julia --project=. -e 'using Pkg; Pkg.test()'

# Iterating on tests: Pkg.test resolves test/Project.toml into a temp env, so it
# reprecompiles CUTEst/JuMP each run. For a faster loop, build the test env once
# and reuse it (Manifest.toml files are untracked, so this costs nothing):
julia --project=test -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=test test/runtests.jl

# There is no per-file test split and no test_args filtering: to run one
# testset, use the reusable env above and include just what you need, e.g.
julia --project=test -e '
    using NLPModelsAlgencan, ADNLPModels, Test
    include("test/problems/autodiff/hs12.jl")
    @test algencan(hs12()).status == :first_order'

# Format (yas style, 92 column margin, see .JuliaFormatter.toml)
julia -e 'using JuliaFormatter; format(".")'

# Build docs locally
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

CUTEst is a test dependency and compiles decoded SIF problems at run time, so
the test suite (not the package) needs `gfortran` on PATH.

CI runs on Linux x64, macOS aarch64 and Windows x64, Julia 1.10, with docs and
coverage deployed from the Linux job only.

## Architecture

### The single `ccall`

`algencan(nlp; kwargs...)` builds an `AlgencanSolver` (which translates the
NLPModel's problem data into Algencan's conventions) and hands it to
`SolverCore.solve!`, which makes exactly one `ccall` into `c_algencan`. Algencan
drives the optimization and calls back into Julia for every evaluation. The
result is unpacked into a `GenericExecutionStats`.

Four evaluation callbacks are wired up (`coded[7]`, `[8]`, `[10]`): `julia_fc`
(objective + constraints), `julia_gjac` (gradient + Jacobian), `julia_hl`
(Hessian of the Lagrangian) and `julia_hlp` (Hessian-vector product). The rest
are passed as `C_NULL`.

### Three things in `solve!` are load-bearing — read the comments before touching

1. **`_CURRENT_SOLVER` and the `_c_julia_*` trampolines.** `c_algencan` has no
   user-data pointer, so the solver cannot be passed to the callbacks: they are
   plain top-level functions reading a module-global `Ref`, set immediately
   before the `ccall` and cleared in a `finally`. They have to stay plain — the
   closure form of `@cfunction` needs a runtime trampoline Julia only has on
   x86, so it fails on Apple Silicon with "closures are not supported on this
   platform". Consequences: a solve is not reentrant or thread-safe, and the
   `Ref` must not outlive the solve (there is a regression test for both).

2. **dlopen/dlclose around every solve.** Algencan 3.1.1 leaks the linear system
   it hands to MA57, after which every later solve in the process silently runs
   without MA57. Loading and unloading the library each solve clears it. This
   only works because `Algencan_jll` declares the product with
   `dont_dlopen=true`; the `@assert`s on `Libdl.dllist()` guard that invariant.
   A fix is pending upstream — the unload and the JLL flag come out together,
   not separately.

3. **`ensure_lp64_blas!()` before every solve.** An HSL-backed Algencan reaches
   MA57 through `libhsl_subset`, which is LP64. Julia registers only an ILP64
   backend, and libblastrampoline answers an unmatched call by writing to stderr
   and returning the result *untouched* — no error, no crash, just a
   factorization on stale data and convergence to a wrong point. Runs per solve,
   not in `__init__`, because loading MKL.jl clears the registration.

### Constraint transformation

Algencan takes constraints as `c(x) <= 0` / `c(x) = 0` only, so
`treat_lower_bounds` sets up a remapping that the callbacks then apply
consistently:

- lower-bound-only constraints (`meta.jlow`) are negated — tracked in `g_sense`;
- range constraints (`meta.jrng`) are **duplicated**: the extra rows are
  appended after the original `m`, mapped by `g_two_smap` / `g_two_sinvmap`.

So the `m` Algencan sees is `solver.m + length(g_two_smap)`, and `jcnnzmax` is
`2 * nnzj`. `solve!` folds the duplicate multipliers back into `solver.mult`
before reporting. Any change to a callback has to keep the `g_has_lb` branch and
the non-`g_has_lb` fast path in agreement. Indices are stored 0-based (`.- 1` in
the constructor) because Algencan's C interface expects them that way.

### Options

Everything not in the fixed tolerance set (`epsfeas`, `epsopt`, `efstain`,
`eostain`, `efacc`, `eoacc`, `outputfnm`, `specfnm`) is passed through
`option2vparam`, which turns `key => value` into the string `"KEY-WITH-DASHES
value"` that Algencan's specification-file parser reads. That is why keyword
arguments use underscores for Algencan's dashes
(`outer_iterations_limit=50` → `OUTER-ITERATIONS-LIMIT 50`). A `specfnm` file is
additionally parsed on the Julia side by
`read_options_from_specification_file`, but only for the six tolerances, since
those are `ccall` arguments rather than vparams.

### Library selection

`algencan_lib_path()` resolves, in order: the `libalgencan_path` preference (set
by `set_algencan_library!`, read at precompile time so Preferences invalidates
the cache), the deprecated `ALGENCAN_LIB_DIR` environment variable (warns in
`__init__`), then `Algencan_jll`. Prefer the preference in anything new —
environment variables are invisible to precompilation.

## Maintainer material

- `docs/src/developer.md` is the real reference for the `Algencan_jll` build:
  the Yggdrasil recipe, why the source tarball is mirrored as a release asset of
  this repo (do not delete the `algencan-3.1.1` release), the run-time HSL
  detection patch, and the recipe details that are easy to break.
- `contrib/hsl/` holds patches for compiling Algencan against HSL by hand; the
  wiki page walks through it. Nothing applies them automatically.
- `contrib/hsl-check/` is a CUTEst sweep harness for comparing an HSL-backed
  build against the fallback. Small problems do not exercise MA57's blocked BLAS
  path, so only a real sweep catches the LP64 failure mode above.
- `deps/` and `debug/` are untracked leftovers from the pre-`Algencan_jll` days,
  when the package built Algencan at install time. They are not part of the
  build.
