# Shipping ALGENCAN as a JLL for NLPModelsAlgencan.jl

Goal: prebuilt ALGENCAN binaries for Julia, **no compilation at install time**,
that use HSL's MA57 when the user has a licensed HSL and fall back to truncated
Newton (TN) otherwise.

~~Betra (trust-region inner solver) is **disabled**, because it needs a locally
patched MA57 that cannot be redistributed.~~ Superseded, see point 1 below:
Betra is kept and the pivot it needs is computed instead.

> **Status: the design works and is in the Yggdrasil pull request, but three
> conclusions below were superseded by later testing.** Read this section first.
>
> 1. **Betra is no longer disabled.** Part A.1 gives up the trust region inner
>    solver because it needs `finfo%pivot`. That turned out to be unnecessary
>    *and* harmful: with MA57 present, `TR` is the solver Algencan normally
>    picks, so disabling it moved licensed users onto `NW`, which on CUTEst
>    hs104 was worse than having no HSL at all. The pivot is instead computed
>    in `moresor.f90` as `u'(B + lI)u`, using the vector `scalcu` already
>    returns. Verified against a patched MA57 that reports the real pivot:
>    agreement to 7 to 11 significant digits over 55 samples. The patch is now
>    2 files, `lssma57.f90` and `moresor.f90`, and `algencan.f90` is untouched.
>
> 2. **The claim that A.3 is safe because nothing reads `pval` is wrong.**
>    Part A.3 lists three callers of `lssfac`; there are four, and it misses
>    `newtd.f90:187`, which is on the `NW` path that HSL enables. The
>    conclusion survives, none of those four read the value, but the evidence
>    as written did not cover it. `moresor` does read it, which is why it is
>    computed there rather than zeroed.
>
> 3. **The verification in C.6 was too weak.** A 2x2 solve never reaches
>    MA57's blocked BLAS path, so it could not expose the problem that actually
>    broke the first JLL: `libhsl_subset` is LP64, Julia registers only an
>    ILP64 BLAS backend, and libblastrampoline answers an unmatched call by
>    writing to stderr and returning without computing. MA57 then factorizes
>    stale data and reports healthy matrices as indefinite. On CUTEst SWOPF
>    that meant 91629 factorizations and an infeasible answer, against 659 and
>    the correct answer once an LP64 backend was forwarded. See
>    "Part H" at the end.
>
> Everything else, in particular the `ma57_available` discovery in C.4, the
> `HSL_jll` layout in C.1 to C.5, and the build recipe in Part B, held up.
> Part E (two JLLs plus a package extension) remains unnecessary.

---

## Part A — ALGENCAN source changes, as a patch

Upstream (Birgin's ALGENCAN 3.1.1) is not forked. The changes live in a single
patch, `algencan-3.1.1-runtime-hsl.patch`, applied by the Yggdrasil recipe to
the pristine upstream tarball. **4 hunks, 2 files, no code removed.**

As finally shipped, after the corrections at the top of this document:

```
sources/algencan/lssma57.f90    3 hunks  runtime MA57 detection + drop finfo%pivot
sources/algencan/moresor.f90    2 hunks  compute the pivot, bail out if lambda stalls
```

The version originally planned here, kept for the record:

```
sources/algencan/algencan.f90   1 hunk   disable the TR inner solver
sources/algencan/lssma57.f90    3 hunks  runtime MA57 detection + drop finfo%pivot
```

Verified: `patch -p1` applies cleanly to a pristine tree and reproduces the
edited tree exactly. In the recipe:

```julia
atomic_patch -p1 ${WORKSPACE}/srcdir/patches/algencan-3.1.1-runtime-hsl.patch
```

All line numbers below refer to the stock 3.1.1 tree.

### A.1 Disable Betra / trust-region inner solver — **NOT DONE, superseded**

> This was the original plan and it is wrong. Disabling the trust region solver
> takes MA57 users off the solver Algencan would normally choose for them, and
> on hs104 that was worse than having no HSL at all. `algencan.f90` is not
> patched; the pivot Betra needs is computed in `moresor.f90` instead. The rest
> of this subsection describes the abandoned approach.

`sources/algencan/algencan.f90:141-147` enables the TR path whenever MA57 is
linked. The patch replaces that `if/else` with an unconditional disable — the
operative lines being:

```fortran
  call setalgparam(val_lsssubTR = 'NONE')
  call setalgparam(val_sclsubTR = 'NONE')
```

(in the patch these are preceded by an eight-line comment explaining why, for
the benefit of anyone reading the diff upstream).

That alone suffices: the selection cascade at `algencan.f90:173-182` guards the
`'TR'` branch with `lsssubTR .ne. 'NONE'`, and `betra` is called from exactly
one place, `gencan.f90:551`, under `if ( inniter .eq. 'TR' )`.

### A.2 The spec file cannot re-enable TR — no patch needed

Investigated and confirmed unnecessary:

- `ikey .eq. 4` (`LINEAR-SYSTEMS-SOLVER-IN-TRUST-REGIONS`) **has no branch at
  all** in `fparam.f90`'s if-chain — it is already a silent no-op upstream.
- `ikey .eq. 3` (`TRUST-REGIONS-INNER-SOLVER`) is guarded at `fparam.f90:487`
  by `seconde .and. lsssubTR .ne. 'NONE'`, and `procpar` runs at
  `algencan.f90:224` — *after* the defaults set at 141-182.

So A.1 alone makes key 3 fall into its existing "option not available"
warning branch. `fparam.f90` is untouched by the patch.

### A.3 Remove the dependence on the patched MA57 — **confirmed necessary**

The patch in `sources/hsl/fix-hsl_ma57.sh` injects a `pivot` field into
`ma57_finfo`. **Verified against the real HSL_jll module file: stock
`hsl_ma57_double` has no `pivot` field** (zero occurrences in
`modules/hsl_ma57_double.mod`). So `lssma57.f90` will not compile against stock
HSL until this is changed.

Both reads are in `lssfac_ma57` (`lssma57.f90:301` and `:309`):

```fortran
       pind    = finfo%more
       pval    = abs( finfo%pivot )      ! <-- does not exist in stock HSL
```

The patch changes both to a bare:

```fortran
       pind    = finfo%more
       pval    = 0.0d0
```

with a six-line comment above the first occurrence recording why the value is
not needed (the second refers back to it).

~~**Safe because** `lssfac` has exactly three callers: `moresor.f90` (TR path,
dead after A.1), `newtonkkt.f90:690` and `accpro.f90:968` — the latter two
receive `pind`/`pval` and never read them, branching only on
`lssinfo`/`nneigv`/`nrank`.~~

**Corrected.** `lssfac` has **four** callers, not three: the missing one is
`newtd.f90:187`, which computes the Newton direction and so sits on the `NW`
path that HSL enables — the worst place to overlook one. Re-traced, all four of
`newtd`, `newtonkkt`, `accpro` and `moresor` declare and pass `pind`/`pval`;
the first three never read them. `moresor` does read `pval`, in the
More-Sorensen safeguard, which is why the shipped patch computes the value
there rather than leaving it zero.

### A.4 Nothing is deleted

`betra.o`/`moresor.o` and every TR-only entry point in `lssma57.f90` stay in
the build, simply unreached. Keeping the patch small and local matters more
than shaving a few kB off the archive, and it keeps the diff reviewable by
upstream.

---

## Part B — Building it (verified, no build-system patch)

The `sources/hsl/` machinery — which compiles raw HSL source and applies the
local MA57 patch — must be bypassed. It turns out **no Makefile patch is
needed**: invoke the core library's own makefile directly and point `HSLSRC`
at a directory holding only `hsl_ma57_double.mod`. That single file drives the
solver selection in `sources/algencan/Makefile`, so MA57 gets the real
`lssma57.o` while MA86/MA97 fall back to their stubs:

```bash
HSL_MODS=<HSL_jll artifact>/modules

mkdir -p hsldetect
ln -sf $HSL_MODS/hsl_ma57_double.mod hsldetect/

make -C sources/algencan lib \
     FC=gfortran-14 AR=ar \
     FFLAGS="-O3 -fPIC -I$HSL_MODS" \
     HSLSRC=$PWD/hsldetect
```

Verified to produce exactly:

```
ar rcs libalgencan.a … lssma57.o dummy_lssma86.o dummy_lssma97.o
```

This also sidesteps the root `Makefile`'s `hsl` target entirely (the `algencan`
target depends on it), so neither Makefile is touched by the patch.

A JLL then needs a shared library rather than the static `libalgencan.a`:

```
gfortran -shared -o ${libdir}/libalgencan.${dlext} \
    -Wl,--whole-archive lib/libalgencan.a -Wl,--no-whole-archive
```

(`--whole-archive` keeps `c_algencan` from being dropped; macOS `-Wl,-all_load`,
Windows needs `--out-implib`.) `c_algencan` (`calgencan.f90`) is already
`bind(C)`, so no extra shim is needed for `ccall`.

---

## Part C — What `HSL_jll` actually contains (verified)

Inspected: `HSL_jll` v4.0.6+0, artifact
`~/.julia/artifacts/1a92e5b95830b425db1ad8f555e08f6479e0cc7f`.

### C.1 Layout

```
lib/libhsl.so           18 KB    23 symbols   <- dummy/stub
lib/libhsl_subset.so   810 KB  1467 symbols   <- real F90 infrastructure
lib/libhsl_subset_64.so 858 KB               <- ILP64 variant
modules/                208 *.mod files
include/                hsl_ma57.h, ..., libhsl.h
```

### C.2 Fortran `.mod` files **are** shipped

`modules/hsl_ma57_double.mod` exists and declares exactly the derived types
`lssma57.f90` uses — `ma57_ainfo`, `ma57_control`, `ma57_factors`,
`ma57_finfo`. **So `use hsl_ma57_double` compiles, and the feared rewrite
against the F77 entry points is not needed.** This was the biggest risk in the
plan; it is retired.

**Compiler constraint:** those `.mod` files are `GFORTRAN module version '15'`.
The system gfortran 15.2 emits version `'16'` and **cannot read them**.
**gfortran-14 emits version `'15'` — verified working.** Build in Yggdrasil
with `preferred_gcc_version = v"14"` (or any GCC in the libgfortran5 range that
still emits `'15'`); use `gfortran-14` locally.

**The licensed package ships no `.mod` files at all** (0 found) — only C
headers. That is fine and is the intended split: compile against the *public*
`HSL_jll` modules, run against whichever library is on the path. Verified in
C.6, including across the version skew (public v4.0.6 modules vs licensed
v2025.7.21 binary).

### C.3 The licensed cores are no-op stubs, and that is the whole trick

- `libhsl.so` (dummy) defines `LIBHSL_isfunctional` plus flat F77 entry points
  `ma27ad_`, `ma57ad_`, `ma57bd_`, `ma57cd_`, `ma57dd_`, `ma57ed_`, `ma57id_`,
  and a few `ma86_*_d`. Every one disassembles to a bare `ret`.
- `libhsl_subset.so` carries the **real F90 module layer** —
  `__hsl_ma57_double_MOD_ma57_analyse`, `_ma57_factorize`, `_ma57_solve1/2`,
  `_ma57_finalize`, `_ma57_get_factors`, `_ma57_enquire`, … — but the
  computational cores it calls (`ma57ad_` … `ma57ed_`) are **also bare `ret`
  stubs** in the unlicensed build.

So **linking succeeds with no licence**, which is what makes the whole
artifact-override design work. Confirmed at runtime:

```
LIBHSL_isfunctional() = 0
```

### C.4 The key finding: `ma57_available` is a public module variable

`hsl_ma57_double` exports a public variable `ma57_available`
(`__hsl_ma57_double_MOD_ma57_available`, in `.bss`). On the unlicensed artifact
it reads **`false`**.

That is precisely the runtime probe you originally asked about — and far
cleaner than `dlopen`/`dlsym`, because it needs no ABI gymnastics. **No `use`
statement has to be added**: `lssma57.f90:15` already carries a bare
`use hsl_ma57_double`, so `ma57_available` is in scope. The patch changes one
line in `lss_ma57()`:

```fortran
-    lss_ma57 = .true.
+    lss_ma57 = ma57_available
```

`lss_ma57()` stops being a compile-time constant and becomes a genuine runtime
query. `algencan.f90:141-182` already branches on it, so **one binary** then
selects NW/MA57 when the user has a licensed HSL and TN when they do not.

### C.5 The licensed build — verified

Inspected `HSL_jll.jl.v2025.7.21`, `override/lib/x86_64-linux-gnu-libgfortran5`:

| | public v4.0.6 | licensed v2025.7.21 |
|---|---|---|
| `libhsl_subset.so` | 810 KB, `ma57ad_` = bare `ret` | **8.7 MB, real `ma57ad_`…`ma57zd_`** |
| → `ma57_available` | **`false`** (`.bss`) | **`true`** (`.data`, statically initialised) |
| `libhsl.so` | 18 KB stub, `LIBHSL_isfunctional()=0` | 15 MB, **also needs `libmetis.so`** |
| `.mod` files | 208 shipped | none |

Conclusions:

- **Link `libhsl_subset`, not `libhsl`.** It is the only library present in
  *both* builds that carries `__hsl_ma57_double_MOD_*` and `ma57_available`,
  and it does not drag in METIS. (`libhsl.so` failed to load here for exactly
  that reason.) Its only extra dependency is `libblastrampoline.so.5`, which
  Julia itself provides.
- `ma57_available` lives in `.data`, statically initialised — so it is readable
  the moment the library is loaded, with no init call first. This resolves the
  timing concern about reading it early at `algencan.f90:141`.
- `LIBHSL_isfunctional` exists **only in `libhsl.so`**, so it is *not* the flag
  to use from ALGENCAN. `ma57_available` is.

### C.6 End-to-end proof: one binary, both libraries

A 25-line Fortran program (`use hsl_ma57_double`, factorise
`[[2,1],[1,3]]`, solve for `x = (1,1)`), compiled **once** with
`gfortran-14 -I<public HSL_jll modules> -lhsl_subset`, then run against each
library via `LD_LIBRARY_PATH`:

```
=== LICENSED libhsl_subset ===
 ma57_available = T
 analyse   flag =  0
 factorize flag =  0   neig= 0   rank= 2
 solve     flag =  0
 x =  0.99999999999999978  1.0000000000000002

=== PUBLIC (unlicensed) libhsl_subset, same binary ===
 ma57_available = F
 analyse   flag = -29
 factorize flag = -29  neig= -1  rank= -1
 solve     flag = -29
 x =  3.0000000000000000   4.0000000000000000
```

This is the whole design working: same binary, no relink, correct solve with a
licence, and **a clean error code (`-29`) rather than a crash or silent garbage**
without one. So `ma57_available` is the primary switch and the `-29` return is a
free second line of defence.

(Note for the recipe: modern `ld` writes `DT_RUNPATH`, which is *not* inherited
by dependencies — an `-Wl,-rpath` on the executable did not let
`libhsl_subset.so` find `libblastrampoline`. BinaryBuilder handles this, but
keep it in mind when testing by hand.)

---

## Part D — The plan: one JLL, real runtime detection

Confirmed viable by C.5/C.6.

1. Apply `algencan-3.1.1-runtime-hsl.patch` (Part A — it contains A.1, A.3 and
   C.4 together).
   - `lss_ma86()`/`lss_ma97()` stay on their dummy stubs — support **MA57
     only**. (`libgomp` comes in via `libhsl_subset` regardless, so that is not
     the reason; the reason is simply fewer moving parts, and MA57 is what the
     NW path prefers first at `algencan.f90:149`.)
   - Use the **LP64** `libhsl_subset` (not `libhsl_subset_64`): ALGENCAN is
     32-bit `integer` throughout.
   - Compile against the **public** `HSL_jll` `.mod` files with GCC 14.
2. Build **one** `Algencan_jll` with `[deps] HSL_jll`,
   `CompilerSupportLibraries_jll`. ALGENCAN itself needs no BLAS — only HSL
   does, via `libblastrampoline`, which `HSL_jll` already pulls in.
3. `NLPModelsAlgencan.jl` depends on that single JLL. Nothing conditional:
   ALGENCAN decides internally, per solve. Licensed users get MA57 by pointing
   Julia's artifact override at the licensed package's `override/` directory —
   no reinstall, no recompile.

No package extension, no "or" dependency, no duplicate libraries. Users without
a licence get TN automatically; users who override the `HSL_jll` artifact get
MA57 with no reinstall of anything.

### D.1 Belt and braces: force TN from Julia

Independently of C.4, the existing `vparam` mechanism can force TN at runtime —
`c_algencan` already takes `nvparam`/`vparam`. Key 7,
`TRUNCATED-NEWTON-LINE-SEARCH-INNER-SOLVER` (`fparam.f90:269`), sets
`innslvr = 'TN'` (`fparam.f90:609`).

```
"TRUNCATED-NEWTON-LINE-SEARCH-INNER-SOLVER"
"SKIP-ACCELERATION-PROCESS"
```

The second matters: the acceleration process (`accpro.f90`/`newtonkkt.f90`)
uses the linear solver independently of the inner solver, so forcing TN alone
does not eliminate every HSL path.

Worth exposing as a user-facing option (forcing TN is sometimes desirable on
its own), but **not needed as a safety net** — C.6 showed the in-library
detection works, and the stubs return a clean `-29` even if it ever did not.

⚠️ **Do not** reach for `LIBHSL_isfunctional()` as the Julia-side check.
It lives only in `libhsl.so`, and the *licensed* `libhsl.so` is the 15 MB
build that additionally needs `libmetis.so` — it **failed to `dlopen` on this
machine** for exactly that reason. Probing it would make licensed users worse
off than unlicensed ones. `ma57_available`, read inside ALGENCAN against
`libhsl_subset`, is the right and only check.

---

## Part E — Two JLLs with a package extension (NOT NEEDED — kept for reference)

**Superseded by Part D.** C.5/C.6 confirmed the single-JLL design works, so
none of this is required. It is retained only to answer the original question —
*can `NLPModelsAlgencan.jl` depend on one JLL or the other?* — and in case the
single-JLL route ever has to be abandoned.

**Julia has no "or" dependency** — every `[deps]` entry is mandatory. The
equivalent is a package extension (Julia ≥ 1.9; this machine runs 1.10.11):

```toml
[deps]
Algencan_jll = "<uuid>"          # non-HSL build, always installed

[weakdeps]
AlgencanHSL_jll = "<uuid>"

[extensions]
NLPModelsAlgencanHSLExt = "AlgencanHSL_jll"

[compat]
julia = "1.9"
```

Plain `add NLPModelsAlgencan` → TN-only build. `add AlgencanHSL_jll` → the
extension loads and switches backend. Cost: both libraries on disk.

`ccall((:sym, lib), …)` needs `lib` constant-foldable, which makes swapping
libraries from an extension awkward. Resolve the symbol to a `Ptr{Cvoid}` once
and `ccall` the pointer, so the long `c_algencan` signature is written once:

```julia
module NLPModelsAlgencan
using Algencan_jll, Libdl

const _ALGENCAN_FPTR = Ref{Ptr{Cvoid}}(C_NULL)

__init__() = (_ALGENCAN_FPTR[] = dlsym(dlopen(Algencan_jll.libalgencan), :c_algencan))

_call_algencan(args...) = ccall(_ALGENCAN_FPTR[], Cvoid, (#= … =#), args...)
end
```

```julia
module NLPModelsAlgencanHSLExt
using NLPModelsAlgencan, AlgencanHSL_jll, HSL_jll, Libdl

function __init__()
    if ccall((:LIBHSL_isfunctional, HSL_jll.libhsl), Cint, ()) != 0
        NLPModelsAlgencan._ALGENCAN_FPTR[] =
            dlsym(dlopen(AlgencanHSL_jll.libalgencan), :c_algencan)
    end
end
end
```

Extension `__init__` runs after the parent's, so the override is well-ordered.

---

## Summary of effort

| Task | Status |
|---|---|
| A. `algencan-3.1.1-runtime-hsl.patch` (4 hunks, 2 files) | **done** — applies with `-p1`, builds, runs |
| B. build without touching either Makefile | **done** — command in Part B, verified |
| ~~`fparam.f90` change~~ | **not needed** — A.2 |
| ~~C.5 licensed-artifact behaviour~~ | **answered** — C.5/C.6 |
| ~~E. two-JLL fallback~~ | **not needed** — single JLL proven |
| ~~F77 rewrite of `lssma57.f90`~~ | **not needed** — `.mod` files ship |
| ~~deleting dead TR code~~ | **not doing** — A.4, patch stays local |
| D. `build_tarballs.jl` recipe | **remaining work** — pin `preferred_gcc_version = v"14"` |
| `NLPModelsAlgencan.jl` wiring | **remaining work** — plain dep, nothing conditional |

### Whole design validated end-to-end

`libalgencan.a` built **once** (gfortran-14, public `HSL_jll` modules), linked
into `toyprob` against `-lhsl_subset`, then run against each library by
swapping `LD_LIBRARY_PATH` only:

| | public (unlicensed) | licensed v2025.7.21 |
|---|---|---|
| `innslvr` | **`TN`** | **`NW`** |
| result | `Flag of ALGENCAN: Solution was found.` | `Flag of ALGENCAN: Solution was found.` |
| `fsub` calls | 38 | 62 |

Same binary, no relink, no recompile — exactly the behaviour the JLL needs.

(`toyprob.f90` ships with `checkder = .true.`, whose derivative checker reads
from stdin and aborts when run non-interactively; set it to `.false.` for
scripted testing.)

**Next step:** the BinaryBuilder recipe — `build_tarballs.jl` with the upstream
tarball plus `bundled/patches/algencan-3.1.1-runtime-hsl.patch`,
`preferred_gcc_version = v"14"`, `HSL_jll` as a dependency, and the shared-lib
link step from Part B.

---

## Part F — Reproduction

Everything below was actually run on `x86_64-linux-gnu`, Julia 1.10.11,
gfortran-14, Ubuntu. Reproduce it before trusting any of the claims above.

### F.0 Prerequisites and paths

```bash
# gfortran-14 is required: the system gfortran 15.x emits GFORTRAN module
# version '16' and CANNOT read HSL_jll's version '15' .mod files.
sudo apt install gfortran-14

# Public (unlicensed) HSL_jll — provides the .mod files used at compile time.
julia -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add("HSL_jll");
          using HSL_jll; println(HSL_jll.artifact_dir)'
```

Set these for the commands that follow:

```bash
PUB=<HSL_jll artifact_dir printed above>          # has lib/ and modules/
MODS=$PUB/modules
LIC=<licensed HSL_jll.jl>/override/lib/x86_64-linux-gnu-libgfortran5
JLIB=$(julia -e 'print(joinpath(Sys.BINDIR,"..","lib","julia"))')
```

On the machine this was developed on:

```
PUB  = ~/.julia/artifacts/1a92e5b95830b425db1ad8f555e08f6479e0cc7f
LIC  = ~/documentos/programas/HSL/HSL_jll.jl/HSL_jll.jl.v2025.7.21/override/lib/x86_64-linux-gnu-libgfortran5
JLIB = ~/.julia/juliaup/julia-1.10.11+0.x64.linux.gnu/lib/julia
```

`libhsl_subset.so` needs `libblastrampoline.so.5`, which Julia ships in
`$JLIB` — hence its presence on every link and run line below.

### F.1 Sanity check the two HSL builds (no ALGENCAN involved)

```bash
julia -e '
using Libdl
dlopen(joinpath(Sys.BINDIR,"..","lib","julia","libblastrampoline.so.5"), RTLD_GLOBAL|RTLD_LAZY)
for p in ARGS
    h = dlopen(p, RTLD_GLOBAL|RTLD_LAZY)
    v = unsafe_load(Ptr{Int32}(dlsym(h, :__hsl_ma57_double_MOD_ma57_available)))
    println("ma57_available=", rpad(Bool(v), 6), "  ", p)
end' $PUB/lib/libhsl_subset.so $LIC/libhsl_subset.so
```

Expected: `false` for the public library, `true` for the licensed one.

### F.2 Apply the patch and build ALGENCAN

From a pristine ALGENCAN 3.1.1 tree:

```bash
patch -p1 < algencan-3.1.1-runtime-hsl.patch      # 2 files, 4 hunks

# Only hsl_ma57_double.mod goes in the detection dir, so that the makefile
# picks the real lssma57.o but leaves MA86/MA97 on their dummy stubs.
mkdir -p hsldetect && ln -sf $MODS/hsl_ma57_double.mod hsldetect/

make -C sources/algencan lib \
     FC=gfortran-14 AR=ar \
     FFLAGS="-O3 -fPIC -I$MODS" \
     HSLSRC=$PWD/hsldetect
```

The final `ar` line must end with `lssma57.o dummy_lssma86.o dummy_lssma97.o`.

### F.3 The decisive test

```bash
# toyprob ships with checkder = .true.; its derivative checker reads stdin and
# aborts when run non-interactively.
sed 's/checkder = \.true\./checkder = .false./' \
    sources/examples/f90/toyprob.f90 > toyprob_nc.f90

gfortran-14 -O3 toyprob_nc.f90 \
    -L$PWD/sources/algencan -lalgencan \
    -L$LIC -lhsl_subset -L$JLIB -lblastrampoline -o toyprob_nc

for P in $PUB/lib $LIC; do
  echo "--- $P ---"
  LD_LIBRARY_PATH=$P:$JLIB ./toyprob_nc | grep -E "innslvr|Flag of ALGENCAN"
done
```

Expected — **same binary, only `LD_LIBRARY_PATH` differs**:

```
--- <public>/lib ---
 innslvr                =                   TN
 Flag of ALGENCAN: Solution was found.
--- <licensed> ---
 innslvr                =                   NW
 Flag of ALGENCAN: Solution was found.
```

Note the link uses `-L$LIC` but that only resolves the soname; the library
actually used is whichever `LD_LIBRARY_PATH` finds at run time. Linking against
the public one instead works identically.

(Modern `ld` writes `DT_RUNPATH`, which dependencies do **not** inherit, so
`-Wl,-rpath` on the executable will not let `libhsl_subset.so` find
`libblastrampoline` — use `LD_LIBRARY_PATH` when testing by hand.)

### F.4 Optional: standalone MA57 probe

Useful to isolate HSL behaviour from ALGENCAN entirely.

```fortran
program ma57test
  use hsl_zd11_double
  use hsl_ma57_double
  implicit none
  type(zd11_type)    :: m
  type(ma57_control) :: cntl
  type(ma57_factors) :: fact
  type(ma57_ainfo)   :: ainfo
  type(ma57_finfo)   :: finfo
  type(ma57_sinfo)   :: sinfo
  real(kind=8) :: rhs(2)

  write(*,*) 'ma57_available =', ma57_available
  m%n = 2; m%ne = 3
  allocate(m%row(3), m%col(3), m%val(3))
  m%row = (/1,1,2/); m%col = (/1,2,2/)
  m%val = (/2.0d0,1.0d0,3.0d0/)          ! [[2,1],[1,3]]
  call ma57_initialize(control=cntl)
  call ma57_initialize(factors=fact)
  call ma57_analyse(m,fact,cntl,ainfo)
  write(*,*) 'analyse   flag =', ainfo%flag
  call ma57_factorize(m,fact,cntl,finfo)
  write(*,*) 'factorize flag =', finfo%flag
  rhs = (/3.0d0, 4.0d0/)                 ! exact solution x = (1,1)
  call ma57_solve(m,fact,rhs,cntl,sinfo)
  write(*,*) 'solve     flag =', sinfo%flag
  write(*,*) 'x =', rhs
end program ma57test
```

```bash
gfortran-14 -I$MODS ma57test.f90 -L$LIC -lhsl_subset -L$JLIB -lblastrampoline -o ma57test
LD_LIBRARY_PATH=$LIC:$JLIB     ./ma57test     # flags 0, x = (1,1)
LD_LIBRARY_PATH=$PUB/lib:$JLIB ./ma57test     # flags -29, rhs untouched
```

### F.5 What was NOT tested — verify these before merging

In rough order of risk:

1. **Only `x86_64-linux-gnu`.** Yggdrasil also builds musl, macOS
   (`.dylib`), FreeBSD and Windows. On Windows the HSL libraries live under
   `override/bin/`, not `lib/`. The `.mod` availability, the `libhsl_subset`
   naming and the `ma57_available` symbol all need re-checking per platform.
2. **GCC version pin may be off by one.** `preferred_gcc_version = v"14"` is
   recommended above because gfortran-14 emits module version `'15'`. But a
   runtime backtrace leaked `/workspace/srcdir/gcc-13.2.0/libgfortran/...`,
   suggesting HSL_jll's own toolchain is **GCC 13**. Both emit `'15'`, so
   either should work — confirm rather than inherit the guess.
3. **`pval = 0.0d0` (A.3)** is the one genuine semantic change. The claim that
   `newtonkkt.f90:690` and `accpro.f90:968` never read `pind`/`pval` was
   established by reading those two call sites; re-trace it independently,
   because an error there fails silently.
4. **`toyprob` is a 2-variable problem.** It does select NW, but exercises MA57
   trivially. Run the `chap13-*` examples or a CUTEst problem to stress
   factorization and the acceleration process.
5. **`HSL_jll` `[compat]`** in the recipe must admit both the public release
   line (v4.0.6 here) and whatever version a licensed override reports
   (v2025.7.21 here).
6. **Dependency-footprint change to the Yggdrasil PR.** The existing PR builds
   a non-HSL ALGENCAN; this makes the JLL depend on `HSL_jll` unconditionally.
   See Part G.

---

## Part G — Notes for the Yggdrasil PR reviewer

The existing PR builds a plain, non-HSL ALGENCAN. This revision adds a patch
and a dependency. The points below are the ones a reviewer is most likely to
raise; stating them up front should save a round trip.

### G.1 No licensed code is redistributed

This is the point to lead with.

- The build links **only** against the public `HSL_jll`, whose MA57 cores are
  no-op stubs (`ma57ad_` … `ma57ed_` each disassemble to a bare `ret`). It is
  the ordinary, freely redistributable JuliaBinaryWrappers artifact.
- **No HSL source or licensed binary enters the build, the recipe, or the
  produced tarballs.** Nothing from an HSL licence is needed to build or ship
  `Algencan_jll`.
- Users who hold an HSL licence supply it the standard way — a `Pkg` artifact
  override pointing `HSL_jll` at their own build — exactly as `HSL.jl` and
  other JuliaSmoothOptimizers packages already do. Nothing about that flow is
  invented here.
- ALGENCAN's own `sources/hsl/` machinery, which compiles user-supplied HSL
  source and applies a local patch to MA57, is **bypassed entirely**. See G.5.

### G.2 Why `HSL_jll` is a hard dependency

A reviewer may reasonably ask why a build that works fine without HSL now
depends on it.

- Julia has **no optional or "or" dependencies** — every `[deps]` entry is
  mandatory. The alternative is two separate JLLs plus a package extension in
  `NLPModelsAlgencan.jl` (documented in Part E), which means two binaries to
  build, publish and keep in lockstep, and both libraries on disk for HSL
  users.
- Because the public `HSL_jll` artifact is a linkable stub, one binary covers
  both cases. ALGENCAN queries `ma57_available` (a public module variable of
  `hsl_ma57_double`) at run time and selects MA57 or truncated Newton
  accordingly.
- Cost to users without a licence: the public `HSL_jll` artifact, ~3 MB
  unpacked (1.7 MB libraries, 1.3 MB Fortran modules), compressed for
  download.
- **Behaviour for those users is unchanged from the current PR**: they get the
  truncated Newton inner solver, the same as a build with no HSL at all. This
  revision is a strict superset — it only adds a capability that activates
  when a licensed HSL is present.

### G.3 The link line does not imply a licensed dependency

In the manual reproduction (Part F.3) the link uses `-L<licensed>`, which can
look alarming out of context. It is not meaningful:

- `-lhsl_subset` resolves a **soname** at link time. Which library is actually
  used is decided at run time by the loader (`LD_LIBRARY_PATH` when testing by
  hand; the artifact or its override in a real Julia session).
- Linking against the **public** stub instead produces an identical binary that
  behaves identically. It was only convenient to point at the licensed copy
  while testing both paths.
- In the recipe proper, BinaryBuilder links against `HSL_jll` from the build
  prefix — i.e. the public stub. There is no path by which a licensed artifact
  participates in the build.

### G.4 Functional differences from a conventional ALGENCAN build

Worth stating explicitly so nobody is surprised later:

- **The trust-region inner solver (BETRA) is disabled.** It reads the offending
  pivot via `finfo%pivot`, a field that exists only in a locally patched MA57
  and is absent from every standard HSL distribution. ALGENCAN's own solver
  cascade already falls back from `TR` to `NW`/`TN`, so this is a
  configuration change, not a code removal — `betra.f90` and `moresor.f90` are
  still compiled, just never selected.
- **MA57 only.** MA86 and MA97 keep ALGENCAN's dummy stubs. MA57 is what the
  Newton path prefers first anyway (`algencan.f90:149`).
- **LP64**: links `libhsl_subset`, not `libhsl_subset_64`. ALGENCAN uses
  default 32-bit `integer` throughout.
- Links `libhsl_subset` rather than `libhsl`. `libhsl` additionally requires
  `libmetis`, and in the licensed build it failed to load here for that reason;
  `libhsl_subset` is the library that carries the `hsl_ma57_double` module
  symbols in **both** the public and licensed builds.

### G.5 Relationship to upstream

- ALGENCAN 3.1.1 is **not forked**. The recipe fetches the pristine upstream
  tarball and applies one patch, `algencan-3.1.1-runtime-hsl.patch`: **4 hunks
  across 2 files, nothing deleted**, each hunk carrying a comment explaining
  itself. It is deliberately small enough to be read and, if he wishes,
  adopted or rejected by the upstream author.
- The root `Makefile` and `sources/algencan/Makefile` are **not patched**. The
  recipe drives `sources/algencan`'s makefile directly and points `HSLSRC` at a
  directory containing only `hsl_ma57_double.mod`, which is what the existing
  makefile logic uses to choose the real `lssma57.o` over the stub (Part B).

### G.6 Failure modes

- If HSL is absent or unlicensed, `ma57_available` is `.false.` and ALGENCAN
  never calls a stub. Even if it did, the stubs return a clean error flag
  (`-29`) rather than crashing or returning silent garbage — verified in F.4.
- If a future `HSL_jll` ever renamed or dropped `ma57_available`, the build
  **fails loudly at compile time** (unresolved name in the module), not
  silently at run time.
- `preferred_gcc_version` matters: gfortran reads `.mod` files only from a
  compatible module-format version (`'15'` for the current `HSL_jll`; GCC 15
  emits `'16'` and cannot read them). A mismatch is a hard compile error, again
  not a silent failure. See F.5 item 2 for the exact pin to confirm.

---

## Part H — The LP64 BLAS trap (added after the first JLL failed)

The first HSL aware `Algencan_jll` produced wrong answers, and the cause was
neither Algencan, nor the patch, nor MA57. It cost several days to find, so it
is written up in full here.

### H.1 Symptom

Against the licensed HSL, CUTEst SWOPF returned `infeasible` at `-6.8984` after
91629 factorizations. A conventional Algencan built with a locally patched MA57
returned `first_order` at `0.06786018` in 659. Across 60 CUTEst problems the JLL
build differed on 13, almost all of them worse.

### H.2 What it was not

Each of these was tested and ruled out, in this order:

| Suspect | How it was excluded |
|---|---|
| The patch | Applied to the MA57 5.2.0 build, SWOPF still solved, 12 digit match |
| BinaryBuilder toolchain | A system gfortran-14 build reproduced the failure bit for bit |
| JLL packaging | Same |
| NLPModelsAlgencan version | Algencan parameter dumps identical between runs |
| Derived type ABI skew, 4.0.6 modules against a 2025.7.21 library | `ma57_control` and `ma57_finfo` are byte identical between the two |
| MA57 definite mode | A probe with `pivoting = 2` on an indefinite matrix returns `flag = -6` correctly |
| METIS version | Forcing `ordering = 0`, which bypasses METIS, changed nothing in either build |
| The MA57 source | Compiling the 2025.7.21 subset sources locally and linking them gave the correct answer |

That last row is what finally isolated it: identical objects, relinked, differing
only in which BLAS they bound to.

### H.3 The cause

`libhsl_subset` is LP64. It calls `dgemm_`, `dgemv_` and `dtpsv_` with 32 bit
integer arguments. A stock Julia session reports

```
LBTConfig([ILP64] libopenblas64_.so)
```

so there is no LP64 backend registered. libblastrampoline answers an unmatched
call by writing `Error: no BLAS/LAPACK library loaded for dgemm_()` to stderr
and returning with the output untouched. The frontal matrix updates therefore
never happen, MA57 factorizes stale data, healthy matrices come back indefinite,
More-Sorensen never converges, and Algencan stops somewhere wrong.

The evidence had been in the logs all along: the failing run contains 3766989 of
those lines, buried under Algencan's own output. The working run contains none.

Two things made it hard to see. A 2x2 test matrix, as used in C.6, never reaches
MA57's blocked BLAS path, so small probes all passed. And `LD_PRELOAD` of a real
BLAS does not help, because `libhsl_subset`'s `DT_NEEDED` binds those symbols to
libblastrampoline directly.

### H.4 The fix

Register an LP64 backend, which is what `HSL.jl` and `Ipopt.jl` do:

```julia
using LinearAlgebra, OpenBLAS32_jll

if !any(lib -> lib.interface === :lp64, BLAS.get_config().loaded_libs)
    BLAS.lbt_forward(OpenBLAS32_jll.libopenblas)
end
```

In `NLPModelsAlgencan.jl` this is `ensure_lp64_blas!`, called before every solve
rather than once at load, because loading MKL.jl or anything else that
reconfigures the trampoline drops the registration. Any LP64 BLAS serves; MKL
was tested and gives the same results as OpenBLAS32. The two interfaces coexist,
so Julia keeps using ILP64 OpenBLAS for its own linear algebra.

With that in place, 55 CUTEst problems match the conventional HSL build, with
one problem the JLL solves that the source build does not.

### H.5 Minimal reproducer

No HSL and no Algencan needed:

```julia
using LinearAlgebra
A = [1.0 2.0; 3.0 4.0]; B = [5.0 6.0; 7.0 8.0]; C = zeros(2,2)
n = Int32(2); alpha, beta = 1.0, 0.0
ccall((:dgemm_, "libblastrampoline"), Cvoid,
      (Ptr{UInt8},Ptr{UInt8},Ref{Int32},Ref{Int32},Ref{Int32},Ref{Float64},
       Ptr{Float64},Ref{Int32},Ptr{Float64},Ref{Int32},Ref{Float64},
       Ptr{Float64},Ref{Int32}),
      "N","N",n,n,n,alpha,A,n,B,n,beta,C,n)
@show C      # zeros, silently
@show A*B    # the answer it should have produced
```
