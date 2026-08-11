# Reaching HSL's MA57 from Algencan_jll

How `NLPModelsAlgencan` gets MA57 without anyone compiling anything at install
time, and without a licensed HSL ever entering a build.

This describes the solution as it ships. The history of how it was arrived at,
including several approaches that were tried and abandoned, is in the git log.

## The problem

Algencan chooses its linear solver at compile time. `lss_ma57()` in
`sources/algencan/lssma57.f90` returns `.true.` in the real file and `.false.`
in `dummy_lssma57.f90`, and the build system picks one. A binary is therefore
either an MA57 build or a truncated Newton build, decided when it is compiled.

That is incompatible with shipping a JLL. HSL is distributed to licensees only,
so the binary has to be built and published without it, and pick it up later if
the user turns out to have one.

## The solution

One binary, one decision made at run time.

`libhsl_subset` exports `ma57_available`, a public module variable of
`hsl_ma57_double` that is `.true.` only when the loaded library really
implements MA57. The public HSL_jll artifact is a linkable stub whose
computational cores are bare `ret` instructions and whose `ma57_available` is
`.false.`; a licensed HSL has the real routines and `.true.`. Returning that
variable from `lss_ma57()` turns Algencan's existing selection cascade into a
run-time query, with no other change to the algorithm.

So:

- `Algencan_jll` is built against the **public** HSL_jll and depends on it
  unconditionally. Nothing licensed enters the build, the recipe or the
  tarballs.
- Users without a licence get truncated Newton, exactly as a build with no HSL
  at all.
- Users with a licence point a `Pkg` artifact override at their own HSL_jll,
  and get MA57 with no reinstall and no recompilation. This is the flow HSL.jll
  and the other JuliaSmoothOptimizers packages already use.

`ma57_available` lives in `.data`, statically initialised, so it is readable the
moment the library loads — there is no init call to order against.

## The patch

`algencan-3.1.1-runtime-hsl.patch` in the Yggdrasil recipe. Upstream 3.1.1 is
not forked: the recipe fetches the pristine tarball and applies one patch,
**2 files, 8 hunks, nothing deleted**.

`sources/algencan/lssma57.f90` (3 hunks)

- `lss_ma57 = .true.` becomes `lss_ma57 = ma57_available`. No `use` statement is
  needed; the file already has a bare `use hsl_ma57_double`.
- The two reads of `finfo%pivot` become `pval = 0.0d0`. That field is not part
  of `ma57_finfo` in any standard HSL distribution — it is added by a local
  patch to MA57 that cannot be redistributed — so `lssma57.f90` will not
  compile against stock HSL until this changes.

`sources/algencan/moresor.f90` (5 hunks)

`pval` has exactly one consumer, the Moré-Sorensen safeguard

```fortran
ls = max( l + (d / ueucn2), ls )
```

which raises the lower bound on the trust-region multiplier so the iteration can
leave the region where the matrix is indefinite. Rather than give that up, the
value is reconstructed in `scalcu`, where the factorization has just failed:
with `u(idx) = 1` and the leading block solved, the pivot the factorization
rejected is the quadratic form `u'(A + lI)u` over the first `idx` unknowns, and
`alin`/`acol` are still in the ordering `u` refers to.

Two things this patch deliberately does **not** do:

- **It does not disable the trust-region inner solver.** An earlier version did,
  on the grounds that it needed `finfo%pivot`. That was wrong and harmful: with
  MA57 present, `TR` is the solver Algencan normally picks, so disabling it
  moved licensed users onto `NW`, which on CUTEst hs104 was worse than having no
  HSL at all.
- **It does not add a fallback when `ls` fails to advance.** An earlier version
  returned `msinfo = 5` in that case so `betra` would take a dogleg direction
  rather than a null step. That was written when the pivot was left at zero and
  `ls` therefore could not advance at all; once the pivot is reconstructed it
  advances normally and the test is redundant. It was also actively harmful: it
  pre-empts the `iter > maxit` bound that already terminates the loop, and
  `betra` re-enters on the same condition. Instrumented on CUTEst `OPTPRLOC` it
  fired 71,939,942 times in 240 seconds while the iteration limit was never
  reached. Removing it took `LAUNCH` from a 30 minute timeout to 3.4 seconds and
  `OPTPRLOC` to 2.7. `moresor.f90` keeps Algencan's own two `msinfo = 5` exits
  and nothing else.

## Building

No Makefile is patched. `sources/algencan/Makefile` selects the real
`lssmaNN.o` over the stub for each solver whose module it finds in `HSLSRC`, so
the recipe simply creates a directory holding only the modules it wants found:

```bash
mkdir -p hsldetect
ln -s ${prefix}/modules/hsl_ma57_double.mod hsldetect/

make -C sources/algencan lib \
     FC=gfortran AR=ar \
     FFLAGS="-O3 -ffree-form -fPIC -I${prefix}/modules" \
     HSLSRC=${PWD}/hsldetect
```

The resulting `ar` line ends with `lssma57.o dummy_lssma86.o dummy_lssma97.o`.
This also sidesteps the root `Makefile`'s `hsl` target, which would otherwise
compile user-supplied HSL source and apply the local MA57 patch.

A JLL needs a shared library, and `--whole-archive` is required so `c_algencan`
survives:

```bash
${FC} -shared -o "${libdir}/libalgencan.${dlext}" \
    -Wl,--whole-archive sources/algencan/libalgencan.a -Wl,--no-whole-archive \
    -L${libdir} -L${prefix}/lib -lhsl_subset -lgfortran
```

macOS uses `-Wl,-all_load`, Windows adds `--export-all-symbols`. `c_algencan`
is already `bind(C)`, so no shim is needed for `ccall`.

## Why `libhsl_subset`

HSL ships several libraries. `libhsl_subset` is the right one and the others
are not:

- It is the only library present in **both** the public and licensed builds
  that carries the `__hsl_ma57_double_MOD_*` symbols and `ma57_available`.
- `libhsl` additionally requires `libmetis`, and the licensed `libhsl` failed to
  load here for that reason. Reaching it would also mean calling the F77 entry
  points directly and abandoning the Fortran modules — a refactor of
  `lssma57.f90` and friends. This was put to @amontoison and confirmed:
  `libhsl_subset` is what a Fortran consumer should link.
- `libhsl_subset_64` is ILP64. Algencan uses default 32-bit `integer`
  throughout, so the LP64 build is the correct one.

Note that Ipopt and Uno both take the other route — they declare the structures
themselves in C, hold the factorization handle as an opaque pointer, dlopen the
user's `libhsl` and resolve the C API. Neither links HSL at build time and
neither uses a Fortran module. That is worth knowing, because it means Algencan
is the only consumer exercising HSL_jll's Fortran modules, which is why the
defects described below went unnoticed.

The public HSL_jll ships 208 `.mod` files; the licensed package ships none, only
libraries and C headers. That split is intentional and is what makes the design
work: compile against the public modules, run against whichever library is on
the path.

## MA86 and MA97

Not shipped yet, waiting on a release of `HSL_jll`.

A JLL is compiled against the public artifact's Fortran modules and only meets
the licensed library at run time, which makes those modules an ABI contract.
For MA86 and MA97 they did not honour it, and a build against them corrupted
memory rather than failing: over 308 CUTEst problems MA86 segfaulted on 107 and
MA97 on 293. The derived types were wrong first, and once those were corrected
the routine interfaces still differed — five MA86 routines missing a trailing
`scale`, and three MA97 mismatches in both directions. Both rounds were fixed
upstream in `ralna/libHSL`.

The lesson worth keeping is how to check it. Compare the argument lists and
intents of every routine in the dummy `hsl_subset` sources against the ones in
the licensed libHSL tarball; they must agree exactly, and MA57 agreeing while
MA86 did not is what localised the fault. Then run something large: a 500x500
solve fails where a 3x3 succeeds, because a shifted argument only crashes when
the address it lands on happens to matter.

Once the modules are fixed, the work here is small and has been tested:

- `lssma86.f90`: `lss_ma86 = .true.` becomes `lss_ma86 = ma86_available`.
- `lssma97.f90`: the same with `ma97_available`.
- `build_tarballs.jl`: symlink `hsl_ma86_double.mod` and `hsl_ma97_double.mod`
  into `hsldetect` beside MA57.
- Do **not** add `-fopenmp`. Both carry `!$omp threadprivate` directives, which
  are comments unless the compiler is told otherwise, and Algencan keeps state
  in common blocks. A serial build is the conservative choice and is what was
  tested.
- `lssma97.f90` ships with CRLF line endings, so a patch touching it carries
  CRLF context lines, and Yggdrasil's root `.gitattributes` (`* text=auto
  eol=lf`) will strip them and silently break it. Add a
  `bundled/patches/.gitattributes` containing `<patch name> -text diff`, as
  `A/algoim/bundled/patches` does.

How to reach them, which is not obvious. The solver is chosen through the
specification file or `vparam`, as `SOLVER [SCALING]` — two words, the scaling
optional and defaulting to `MC64`:

| slot | keyword | solvers |
|---|---|---|
| trust regions | `LINEAR-SYSTEMS-SOLVER-IN-TRUST-REGIONS` | **MA57 only** |
| Newton line search | `LINEAR-SYSTEMS-SOLVER-IN-NEWTON-LINE-SEARCH` | MA57, MA86, MA97 |
| acceleration | `LINEAR-SYSTEMS-SOLVER-IN-ACCELERATION-PROCESS` | MA57, MA86, MA97 |

The `...-INNER-SOLVER` variants set the inner solver as well as the linear one.
From Julia any unrecognised keyword argument becomes a `vparam` line with
underscores turned into dashes, so
`algencan(nlp; NEWTON_LINE_SEARCH_INNER_SOLVER = "MA86 MC77")` is the way to ask.

MA86 and MA97 are never defaults and never reach the trust region: the cascade
in `algencan.f90` gives `lsssubTR` only `MA57` or `NONE`, and `fparam.f90`
contains exactly one `setalgparam(val_lsssubTR = ...)`, guarded by `MA57`. That
is consistent with `lssfac_ma86` and `lssfac_ma97` never assigning `pval`, so
they cannot feed the Moré-Sorensen safeguard. Note also that acceleration does
not follow the Newton setting — `algencan.f90` sets `lsssubACC = lsssubNW`
before `procpar` reads the specification file, so overriding only the Newton
keyword leaves acceleration on MA57.

With faithful modules, a build against the licensed `libhsl_subset` solves
CUTEst `HS106` with MA86 to the true optimum, 7049.24802053, where MA57 stops at
an infeasible 7239.49565125.

Depending on a fixed HSL_jll will also need a version bound, which is awkward
because the two distributions use different schemes: the licensed packages are
versioned by date (`2025.7.21`) and the registered one is `4.0.6`. Ipopt.jl
expresses this as `HSL_jll = "3, 4, 2023, 2024, 2025"`. Our floor is higher than
theirs — `2023.11.7` ships no `libhsl_subset` at all and exports no
`*_available` symbols, so a licensee on it gets a hard
`libhsl_subset.so => not found` rather than a fallback to truncated Newton.

## Pitfalls

**The LP64 BLAS trap.** This is the one that matters most. `libhsl_subset` is
LP64 and calls `dgemm_`, `dgemv_` and `dtpsv_` with 32-bit integer arguments.
Julia registers only an ILP64 backend, and libblastrampoline answers an
unmatched call by writing `Error: no BLAS/LAPACK library loaded for dgemm_()` to
stderr and **returning with the output untouched**. No error, no crash: the
frontal matrix updates never happen, MA57 factorizes stale data, healthy
matrices come back indefinite, and the solver stops somewhere wrong. On CUTEst
SWOPF that was 91629 factorizations and an infeasible answer, against 659 and
the correct one with an LP64 backend forwarded. The failing run contained
3766989 of those stderr lines, buried under Algencan's own output.

`ensure_lp64_blas!` forwards `OpenBLAS32_jll` before **every** solve, not once
at load, because loading MKL.jl or anything else that reconfigures the
trampoline drops the registration. Any LP64 BLAS serves. The two interfaces
coexist, so Julia keeps using ILP64 OpenBLAS for its own linear algebra. If HSL
results ever look wrong again, check `BLAS.get_config()` for an `[LP64]` entry
and grep stderr for `no BLAS/LAPACK library loaded` before suspecting anything
else.

**Small probes prove nothing.** A 2x2 solve never reaches MA57's blocked BLAS
path, so it cannot expose the trap above; every small test passed while the
build was badly broken. Use CUTEst, or at least the `chap13-*` examples.

**`LD_PRELOAD` does not override HSL's BLAS.** `libhsl_subset` has
libblastrampoline in its `DT_NEEDED`, so those symbols bind there directly.

**`DT_RUNPATH` is not inherited by dependencies.** Modern `ld` writes `RUNPATH`
rather than `RPATH`, so `-Wl,-rpath` on an executable will not let
`libhsl_subset.so` find `libblastrampoline`. Use `LD_LIBRARY_PATH` when testing
by hand; BinaryBuilder handles it in the recipe.

**Outside Julia there is no BLAS backend at all.** A bare Fortran program
linked against `libhsl_subset` jumps through an unregistered trampoline and
segfaults for reasons that have nothing to do with the solver. Preload a real
BLAS or put a Julia `lib/julia` on `LD_LIBRARY_PATH`.

**Do not probe `LIBHSL_isfunctional()` from Julia.** It exists only in
`libhsl.so`, and the licensed `libhsl.so` needs `libmetis` and failed to load
here. Probing it would make licensed users worse off than unlicensed ones.
`ma57_available`, read inside Algencan against `libhsl_subset`, is the check.

**Everything HSL related is verified only on `x86_64-linux-gnu`.** CI proves the
other platforms build, not that MA57 works there. Licensed HSL ships no Windows,
i686, armv6l/v7l or riscv64 build, so those users get truncated Newton
regardless. On Windows the HSL libraries live under `override/bin/`, not `lib/`.

**A missing `ma57_available` fails loudly.** If a future HSL_jll renamed or
dropped it, the build breaks at compile time with an unresolved name, not
silently at run time. Likewise a module-format mismatch is a hard compile error.
The recipe pins no GCC version and does not need to: HSL_jll's modules are
`GFORTRAN module version '15'` and gfortran 15.2 reads them without complaint.
