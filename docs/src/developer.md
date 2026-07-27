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

The recipe builds Algencan *without* HSL, matching what a plain source build of
Algencan produces. Distributing an HSL-linked binary is not possible, since the
HSL linear solvers are proprietary. See
[Using a different Algencan library](@ref) below for how to use your own build.

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
namespaces directly; Docker is only needed on macOS. These notes were written
against Julia 1.10 and BinaryBuilder 0.6.6 on `x86_64-linux`.

Create a scratch environment with BinaryBuilder in it — keep it separate from
this package's environment:

```bash
mkdir -p ~/tmp/bb && cd ~/tmp/bb
julia --project=. -e 'using Pkg; Pkg.add("BinaryBuilder")'
```

Then build a single platform, running from the directory holding the recipe so
that `build/` and `products/` land beside it:

```bash
cd /path/to/Yggdrasil/A/Algencan
julia --project=~/tmp/bb build_tarballs.jl --verbose x86_64-linux-gnu-libgfortran5
```

Omit the triplet to build every platform in the recipe, which takes considerably
longer and downloads a compiler toolchain per target.

Add `--debug` to drop into a shell inside the sandbox when a step fails. Do not
use it in a non-interactive or backgrounded run — it will wait forever for input.

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

Three things in the recipe are easy to get wrong or to "clean up" by mistake.

### `dont_dlopen=true` is load-bearing

The library product is declared as:

```julia
products = [
    LibraryProduct("libalgencan", :libalgencan; dont_dlopen=true),
]
```

Algencan is Fortran code that keeps state in common blocks, so this package
loads the library and unloads it around *every* solve to start from a clean
state. A JLL that opened the library in its own `__init__` would keep it
resident and defeat that, and the `@assert`s guarding the load/unload cycle in
`src/NLPModelsAlgencan.jl` would fail immediately. Do not remove this flag.

### Command-line triplets bypass the `platforms` variable

Passing a bare triplet such as `x86_64-linux-gnu` overrides the recipe's
`platforms` list entirely, including the `expand_gfortran_versions` call. The
build then falls back to BinaryBuilder's default toolchain, which links
`libgfortran.so.3` — an ABI that `CompilerSupportLibraries_jll` no longer ships,
producing a library that cannot be `dlopen`ed. Always name the tagged triplet,
`x86_64-linux-gnu-libgfortran5`, when testing a single platform.

### libgfortran ABI tagging

```julia
platforms = expand_gfortran_versions(platforms)
```

`libalgencan` links `libgfortran`, whose ABI breaks across GCC major versions, so
each tarball has to be tagged with the version it binds to. Without this line,
Pkg has no way to hand Julia a matching build.

## Using a different Algencan library

The `Algencan_jll` build has no HSL support, which costs a lot of performance on
larger problems. To use your own build — typically one linked against MA57:

```julia
using NLPModelsAlgencan
set_algencan_library!("/path/to/libalgencan.so")
```

The path is stored as a preference of the active project, so it applies to that
project alone and survives restarts. Julia must be restarted for the change to
take effect. Pass `nothing` to go back to the library from `Algencan_jll`.

The `ALGENCAN_LIB_DIR` environment variable is still honoured, so setups
predating `Algencan_jll` keep working, but it is deprecated and warns on load:
environment variables are invisible to precompilation, so a change to one does
not invalidate the cached module and can be silently ignored.

### The HSL patches in `contrib/hsl`

These are the patches the old, pre-`Algencan_jll` build script applied when
compiling Algencan against HSL. They are kept because they are still needed for
a manual HSL build; nothing applies them automatically any more. The
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
