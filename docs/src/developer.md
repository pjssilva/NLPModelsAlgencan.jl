# Developer notes

Notes for maintaining NLPModelsAlgencan.jl itself. Users of the package do not
need anything here.

## Where the Algencan binary comes from

NLPModelsAlgencan.jl does not compile Algencan. It depends on `Algencan_jll`, a
binary package built by [Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil),
the community build system that produces every `*_jll` package in the General
registry. The build recipe lives there:

```
Yggdrasil/A/Algencan/build_tarballs.jl
```

That file is the only authoritative copy. It is deliberately **not** duplicated
in this repository: Yggdrasil maintainers periodically sweep every recipe in the
repository when BinaryBuilder APIs change or new platforms are added, so a local
copy would silently fall out of date, and edits made to it would never reach
Yggdrasil anyway.

The recipe links the *public* `HSL_jll`, whose solvers are stubs, and patches
Algencan so that it decides at run time which of MA57, MA86 and MA97 are really
present. One binary therefore serves both cases: users without an HSL licence
get the truncated Newton inner solver, and users who override the `HSL_jll`
artifact with their licensed build get the solvers that build provides. No
licensed code enters the build or the tarballs. See [HSL support](@ref) for how
that works and what it requires.

### The source tarball mirror

The recipe does not download Algencan from
`https://www.ime.usp.br/~egbirgin/tango/sources/`. It fetches a byte-identical
copy (same SHA256) published as an asset of the `algencan-3.1.1` release of this
repository, because Yggdrasil rebuilds recipes years after they are merged and a
personal academic URL is not a dependable long-term build input. Redistributing
the tarball this way is permitted: Algencan is GPL-2.0-or-later.

!!! warning "Do not delete the `algencan-3.1.1` release"
    That release is a permanent build input of `Algencan_jll`. Every Yggdrasil
    rebuild refetches the asset and verifies its SHA256. Deleting the release,
    removing the asset, or renaming the tag makes `Algencan_jll` unbuildable.

## Changing the recipe

All changes — a new Algencan version, an added platform, a fix requested by a
reviewer — follow the same route:

1. Fork and clone <https://github.com/JuliaPackaging/Yggdrasil>.
2. Edit `A/Algencan/build_tarballs.jl` in your clone. If you are bumping the
   Algencan version, publish the new upstream tarball as a release asset here
   first (see above) and update both the URL and the SHA256.
3. Test locally before opening the PR (next section). This matters: Yggdrasil CI
   builds every platform, and a recipe that fails there wastes reviewer time.
4. Push and open a pull request against `master`, titled `[Algencan] ...`.
5. Yggdrasil CI builds all platforms on the PR. Fix whatever fails and push again.
6. On merge, Yggdrasil builds the tarballs, uploads them to a release of
   `JuliaBinaryWrappers/Algencan_jll.jl`, and registers the new `Algencan_jll`
   version in the General registry, usually within an hour.
7. Once the new `Algencan_jll` is registered, bump its `[compat]` entry in this
   package's `Project.toml` if the version bound needs to change.

## Building and testing the recipe locally

BinaryBuilder runs the build in a sandbox. On Linux it uses unprivileged user
namespaces directly; Docker is only needed on macOS. It needs Julia 1.12 or
newer.

Build from the Yggdrasil clone's own environment, which pins the BinaryBuilder
version the recipe builds under. Instantiate it after cloning, and again after
updating to a newer upstream master; otherwise the build fails on missing
dependencies, which looks like a broken recipe and is not:

```bash
cd /path/to/Yggdrasil
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then build a single platform, running from the directory holding the recipe so
that `build/` and `products/` land beside it:

```bash
cd /path/to/Yggdrasil/A/Algencan
julia --project=/path/to/Yggdrasil build_tarballs.jl --verbose x86_64-linux-gnu-libgfortran5
```

On Ubuntu and derivatives the sandbox is blocked by default, and BinaryBuilder
falls back to asking for `sudo`:

```bash
sudo sysctl -w kernel.apparmor_restrict_unprivileged_userns=0
```

That reverts on reboot.

Omit the triplet to build every platform in the recipe, which takes considerably
longer and downloads a compiler toolchain per target.

Add `--debug` to drop into a shell inside the sandbox when a step fails. Do not
use it in a non-interactive or backgrounded run — it will wait forever for input.

Building the macOS targets fails locally with `macOS SDK not installable` unless
you let BinaryBuilder download Apple's SDK, which means accepting Apple's
license:

```bash
export BINARYBUILDER_AUTOMATIC_APPLE=true
```

Yggdrasil's CI already has the SDK, so macOS is built on the pull request either
way; this only matters if you want to test those targets on your own machine.

### Checking the result

A successful build is not the same as a usable one. BinaryBuilder reports audit
problems as *warnings*, so the process still exits 0 when it produces a library
that cannot be loaded. Read the audit section of the log, and check the artifact
directly:

```bash
tar xzf products/Algencan.v3.1.1.x86_64-linux-gnu-libgfortran5.tar.gz -C /tmp/check
nm -D /tmp/check/lib/libalgencan.so | grep c_algencan   # the entry point we ccall
objdump -p /tmp/check/lib/libalgencan.so | grep NEEDED  # expect libgfortran.so.5
```

### Testing against this package

`--deploy=local` builds the JLL and writes it to `~/.julia/dev/Algencan_jll`
without uploading anything:

```bash
julia --project=~/tmp/bb build_tarballs.jl --verbose --deploy=local x86_64-linux-gnu-libgfortran5
```

Then point this package at it and run the test suite:

```julia
using Pkg
Pkg.activate("/path/to/NLPModelsAlgencan.jl")
Pkg.develop(path=expanduser("~/.julia/dev/Algencan_jll"))
Pkg.test()
```

## Recipe details worth knowing

These are easy to get wrong or to "clean up" by mistake.

### `dont_dlopen=true` is load-bearing

The library product is declared as:

```julia
products = [
    LibraryProduct("libalgencan", :libalgencan; dont_dlopen=true),
]
```

Algencan 3.1.1 leaks the linear system it hands to MA57, and once that happens
every later solve in the same process quietly runs without MA57. This package
therefore loads the library and unloads it around *every* solve, which clears the
leak. A JLL that opened the library in its own `__init__` would keep it resident
and defeat that, and the `@assert`s guarding the load/unload cycle in
`src/NLPModelsAlgencan.jl` would fail immediately. Do not remove this flag while
the unload is there; both go together, once Algencan itself is fixed.

### Command-line triplets bypass the `platforms` variable

Passing a bare triplet such as `x86_64-linux-gnu` overrides the recipe's
`platforms` list entirely, including the `expand_gfortran_versions` call. The
build then falls back to BinaryBuilder's default toolchain, which links
`libgfortran.so.3` — an ABI that `CompilerSupportLibraries_jll` no longer ships,
producing a library that cannot be `dlopen`ed. Always name the tagged triplet,
`x86_64-linux-gnu-libgfortran5`, when testing a single platform.

### `lssma97.f90` has CRLF line endings

The recipe normalises them before applying the patch:

```bash
sed -i 's/\r$//' sources/algencan/lssma97.f90
```

Without it the patch fails with "different line endings", because Yggdrasil's
root `.gitattributes` (`* text=auto eol=lf`) rewrites the CRLF context lines out
of the stored patch on checkout. Keep the patch plain LF and let the recipe do
the normalising.

### `HSL_jll` carries no version bound

```julia
Dependency("HSL_jll"),
```

A JLL's `[compat]` is frozen into the registered version at build time and
cannot differ between build numbers of one version, so a bound here would lock
users out of future HSL_jll releases and force a new Algencan version to lift.
If a bound is ever genuinely needed, it belongs in this package's
`Project.toml`, which can be revised in an ordinary release.

### libgfortran ABI tagging

```julia
platforms = expand_gfortran_versions(platforms)
```

`libalgencan` links `libgfortran`, whose ABI breaks across GCC major versions, so
each tarball has to be tagged with the version it binds to. Without this line,
Pkg has no way to hand Julia a matching build.

## HSL support

Algencan solves its subproblems much faster with an HSL linear solver, but HSL
is licensed and cannot be redistributed. `Algencan_jll` gets around that by
deciding at run time, so a single binary covers both cases.

### How the run-time switch works

`HSL_jll`'s `libhsl_subset` exports `ma57_available`, `ma86_available` and
`ma97_available`, public variables of the `hsl_maNN_double` modules. Each is
`.false.` in the freely distributed stub and `.true.` in a licensed build that
implements that solver. A patch, carried in the recipe's `bundled/patches`,
changes Algencan's `lss_ma57()`, `lss_ma86()` and `lss_ma97()` to return the
matching variable rather than a compile time constant, which is enough to turn
Algencan's own solver selection into a run time decision.

MA57 stays the default and is the only solver the trust region accepts. MA86 and
MA97 are reachable for the Newton line search and the acceleration process
through the specification file, which from Julia means a keyword argument; the
[main page](index.md) shows the form. The recipe symlinks all three
`hsl_maNN_double.mod` files into the directory `sources/algencan/Makefile`
searches, which is what makes it compile the real `lssmaNN.o` for each rather
than the stub.

The same patch removes Algencan's use of `finfo%pivot`. That field is not part
of `ma57_finfo` in a standard HSL distribution, it is added by the local MA57
patch in `contrib/hsl`, so a build against stock HSL does not compile without
this change. Its only consumer is the More-Sorensen safeguard in `moresor.f90`,

```fortran
ls = max( l + (d / ueucn2), ls )
```

which raises the lower bound on the trust region multiplier so that the
iteration can leave the region where the matrix is indefinite. The value can be
computed rather than queried: with the vector `u` that `scalcu` already returns,
the pivot the factorization rejected is the quadratic form `u'(B + lI)u` over
the leading block.

This is a reconstruction rather than a reproduction. It matches the real
`finfo%pivot` to eight or more significant digits while the pivot is well
scaled, and loses it once the pivot falls to where summing the quadratic form
cancels — at which point the safeguard it feeds is inactive anyway. It also
costs a pass over the matrix entries each time a pivot is rejected. A build
against a locally patched MA57 avoids both, which is one reason
`set_algencan_library!` remains supported.

A licensed user enables MA57 by pointing a `Pkg` artifact override at their own
`HSL_jll` build, the same way `HSL.jl` and other JuliaSmoothOptimizers packages
do. Nothing has to be reinstalled or recompiled.

### An LP64 BLAS has to be forwarded

!!! warning "Without this, MA57 silently returns wrong answers"
    `libhsl_subset` is LP64: it calls `dgemm_` and friends with 32 bit integers.
    Julia registers only an ILP64 backend, and libblastrampoline answers a call
    it cannot match by writing a line to stderr and returning with the result
    untouched. There is no error and no crash. The factorization simply works on
    stale data, MA57 reports as indefinite matrices that are not, and the solver
    converges somewhere wrong.

`ensure_lp64_blas!` registers an LP64 backend from `OpenBLAS32_jll` when none is
present, and runs before every solve. Once at load time is not enough: loading
MKL.jl, or anything else that reconfigures the trampoline, drops the
registration afterwards. A backend that is already there, MKL's for instance, is
left alone. Any LP64 BLAS will do, the interface is what matters and not the
implementation.

This belongs in this package rather than in the JLL, as it does in `HSL.jl` and
`Ipopt.jl`: forwarding is a run time decision about the current session, which a
binary artifact cannot make.

## Using a different Algencan library

If you compiled Algencan yourself, for instance against a patched MA57 in the
old way, you can use it instead of the one from `Algencan_jll`:

```julia
using NLPModelsAlgencan
set_algencan_library!("/path/to/libalgencan.so")
```

The path is stored as a preference of the active project, so it applies to that
project alone and survives restarts. Julia must be restarted for the change to
take effect. Pass `nothing` to go back to the library from `Algencan_jll`.

The `ALGENCAN_LIB_DIR` environment variable is honoured too, but it is deprecated
and warns on load: environment variables are invisible to precompilation, so a
change to one does not invalidate the cached module and can be silently ignored.

### The HSL patches in `contrib/hsl`

These are needed to compile Algencan against HSL by hand, for the library
`set_algencan_library!` then points at. Nothing applies them automatically. The
[wiki page on compiling HSL libraries](https://github.com/pjssilva/NLPModelsAlgencan.jl/wiki/Compiling-HSL-Libraries-for-use-with-NLPModelsAlgencan.jl)
walks through the whole process.

`patch_algencan.txt` applies to the Algencan 3.1.1 source tree with `patch -p1`
and edits its top-level `Makefile` to:

- use `gcc`/`g++` rather than the hardcoded `gcc-4.9`/`g++-4.9`;
- add `-fPIC` to `CFLAGS` and `FFLAGS`, needed to link a shared library;
- set `BLAS_LAPACK := 0`, so the system BLAS/LAPACK are used instead of the
  reference copies bundled with Algencan;
- read the METIS and MA57 locations from the `METISPATH` and `MA57PATH`
  environment variables;
- comment out the CUTEst paths, which otherwise point at directories that do not
  exist.

`patch_ma57.txt` applies to the HSL MA57 5.2.0 source tree with `patch -p1` and
changes `src/ddeps.f` and `src/hsl_ma57d.f90` so that `RINFO(20)` is set to the
pivot value on the error and warning paths where MA57 would otherwise leave it
undefined. Algencan reads that value to decide how to react to a failed
factorization.
