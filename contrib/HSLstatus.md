# Reaching HSL from NLPModelsAlgencan: where things stand

Status of the work to reach HSL's linear solvers through `Algencan_jll`: what
is done, what is left, and what MA86 and MA97 are still waiting on. Last
updated 7 August 2026. Keep until that work is finished, then fold whatever is
still true into the manual and delete the rest.

For how the pieces actually work, read the "Developer notes" page in the manual
instead; for the full HSL investigation,
`contrib/patched_algencan_for_jll.md`.

## What was done

**Released.** v0.8.1, carrying @joaquimg's fix for `@cfunction` closures being
unsupported on aarch64 — the package could not run at all on Apple Silicon or
ARM Linux — plus a follow up so the callback reference stops keeping the first
solver alive for the whole session.

**Branch `jll-migration`, pushed.** Moves from compiling Algencan at install
time to `Algencan_jll`:

- `Algencan_jll` and a Preferences based override replace `deps/build.jl`;
  `set_algencan_library!` selects a custom library, `ALGENCAN_LIB_DIR` still
  works but is deprecated.
- `ensure_lp64_blas!` forwards an LP64 BLAS before each solve. Not optional,
  see below.
- Docs, CI and the developer notes updated. `deps/build.jl` deleted, the old
  HSL patches moved to `contrib/hsl`.
- `contrib/hsl-check/` runs a CUTEst sweep, for comparing an HSL build against
  one without, serially or several problems at a time.

**Yggdrasil PR #14301**, branch `algencan` in the fork clone at
`~/documentos/programas/Yggdrasil_Algencan`, commit `4ab57f90`, pushed and
**accepted**. Builds Algencan with run time HSL detection, MA57 only.

The open question about linking `libhsl` instead of `libhsl_subset` is closed.
@amontoison: `libhsl_subset` is the right one precisely because Algencan uses
the Fortran modules; reaching `libhsl` would mean calling the F77 symbols
directly, which is a refactor of `lssma57.f90` and friends. Ipopt and Uno both
do exactly that — they declare the structures themselves in C, hold the
factorization handle as an opaque pointer, dlopen the user's `libhsl` at run
time and resolve the C API. Neither links HSL at build time and neither uses a
Fortran module. That is worth knowing for two reasons: it is the reason nobody
had noticed the module problems described below, and it is the shape of any
future rewrite, should the module route ever become untenable.

## What is left

1. **Done.** #14301 is accepted and `Algencan_jll` 3.1.1+0 is in General.
2. **Done, 7 August 2026.** `Pkg.free("Algencan_jll")` here and in the sweep
   environment, re-resolved, tests pass against the registered build rather
   than a locally deployed one.

   Worth having checked, since the registry carries a single build and it was
   not obvious whether it predated the fix: it does not. With the registered
   JLL and a licensed HSL, `OPTPRLOC` solves in 1.6 s and `LAUNCH` in 2.3 s,
   to the same objectives as the local build, so the merged recipe includes
   `4ab57f90` and the dogleg fallback is gone from what users get. Had it
   timed out instead, the recipe would have needed rebuilding.
3. **Done.** `contrib/yggdrasil/` is gone. The recipe's home is the Yggdrasil
   fork at `~/documentos/programas/Yggdrasil_Algencan`, and a second copy here
   could silently drift from the one that actually ships: during the August
   2026 work both copies were edited and were briefly out of step, with a
   broken patch in one of them. The investigation writeup it contained was
   worth keeping and is now `contrib/patched_algencan_for_jll.md`. Old
   revisions of the recipe copy remain in git history if ever needed.
4. Push `jll-migration`, merge to master, register 0.9.0.
5. Optionally enable the macOS and Windows CI matrix entries, now that the JLL
   covers them and no compiler is needed.
6. Run the wider CUTEst sweep, on a machine with cores to spare. **Done, on
   quorra, 5 August 2026**: two runs, 6 jobs, all 308 problems, on a locally
   deployed JLL, once with MA57 reachable and once without. The result files
   are scratch output of one machine and are not kept here; see "Numbers worth
   keeping" below for what came out of them.

   ```bash
   cd contrib/hsl-check
   tmux new-session -d -s cutest \
     'CUTEST_CHECK_JOBS=8 CUTEST_CHECK_LIMIT=308 ./setup_and_run.sh 2>&1 | tee tmux.log'
   ```

   One problem per process, dynamically scheduled, with a 30 minute limit each.
   Results are appended as they finish and reruns resume, so it can be left
   alone and looked at later. Of the five problems that never finished before
   — `LAUNCH`, `LHAIFAM`, `NASH`, `OPTPRLOC`, `PALMER5ANE` — three now finish
   in seconds, so their earlier non-completion was the old serial run being
   fragile rather than anything real. `LAUNCH` and `OPTPRLOC` timed out on the
   first sweep and no longer do; that was the patch, see below. What still
   times out is the `DMN*` family, `DMN37143` in every configuration tried and
   `DMN15103` in the last one.

   One trap, which cost a wasted hour-long run here. Leaving the HSL argument
   off does not disable HSL: `Pkg.develop` on this repo also carries over the
   `HSL_jll` developed in its own `Manifest.toml`, so the licensed one arrives
   unasked and the sweep silently repeats the HSL run while printing a warning
   saying the opposite. `setup_and_run.sh` now calls `Pkg.free("HSL_jll")`, and
   the results header carries an `ma57` line saying which solver actually ran.
   Read that line, not the `HSL_jll` version: the public stub and the licensed
   release can share a version string and both export the MA57 symbols.

7. Add MA86 and MA97, once @amontoison has fixed the Fortran modules in the
   public `HSL_jll`. Agreed with him on 7 August 2026: ship #14301 with MA57
   only, he fixes the modules over the following weeks, then we try again. What
   that will take is in `contrib/patched_algencan_for_jll.md`.

## MA86 and MA97

Blocked on the public `HSL_jll`, not on Algencan and not on HSL. The licensed
`libhsl_subset` implements all three solvers correctly; three of the Fortran
modules shipped in the public artifact do not match it, so a binary built
against them misbehaves against the licensed library — `hsl_mc69_double` will
not compile, `ma86_factor` returns a corrupted `info%flag`, and `ma97_akeep` is
truncated to the point of corrupting the heap.

Reported to @amontoison on 7 August 2026 and acknowledged; he expects to fix
them over the following weeks. The reproduction sent to him is in
`~/documentos/programas/hsl-jll-module-report`.

The numbers that identify each defect, what to change here once they are fixed,
and how the solvers are selected at all, are in
`contrib/patched_algencan_for_jll.md` under "MA86 and MA97". That document is
permanent; this one is not, so the details live there.

## The one thing not to forget

`libhsl_subset` is LP64 and Julia registers only an ILP64 BLAS backend. An
unmatched call returns **with the output untouched**, no error and no crash, so
MA57 factorizes stale data and the solver stops somewhere wrong. If HSL results
ever look wrong again, check `BLAS.get_config()` for an `[LP64]` entry and grep
stderr for `no BLAS/LAPACK library loaded` before suspecting anything else.

The full account, and the rest of the traps, are under "Pitfalls" in
`contrib/patched_algencan_for_jll.md`.

## Numbers worth keeping

- 55 CUTEst problems, HSL build against a conventional Algencan with locally
  patched MA57: no status regressions, one problem the JLL solves that the
  source build does not. Objective differences around 1e-5 are expected,
  Algencan stops on a gradient criterion.

- 308 CUTEst problems, the same JLL with MA57 reachable and with it not, on
  quorra, 6 jobs. Nothing crashed in either: no `CRASHED` rows, no Julia level
  exceptions, 308 of 308 recorded both times. MA57 solves 201 to first order
  against 196, 25 problems change status, and of the 190 both solve, 183 agree
  in the objective to better than 1e-5. Of the seven that do not, three are
  marginal and four are plainly different local minima, with MA57 finding the
  better point on `OET6`, `OET7` and `ROBOT`. `SWOPF`, which is what exposed
  the ILP64 disaster, now agrees to eight significant figures in 2.6 seconds.

  Read this for what it is. Two sound methods, a trust region over MA57 and
  truncated Newton, have no reason to agree problem by problem or for either
  to dominate over a selection this broad, so the 25 changes of status are the
  expected outcome and not a finding. Eleven problems go the way of MA57 and
  eight the other way; that split says little on its own. What the sweep is
  evidence for is narrower and worth more: the HSL path does not crash, does
  not hang the sweep and does not compute nonsense. It is not evidence of "no
  status regressions", which was measured against a conventional build with a
  patched MA57 and would need a source build to reproduce at this size.

  `LAUNCH` and `OPTPRLOC` time out at 30 minutes here and finish without MA57.
  That turned out to be a defect in the patch rather than anything about MA57,
  and is fixed; the entry below has it. `DMN37143` also times out, but it does
  so in every configuration tried, including conventional builds, so it is
  simply a slow problem.
- 308 CUTEst problems, the JLL against the conventional Algencan with the
  locally patched MA57 that the 55 problem comparison used, same selection, 6
  jobs, the library selected through `ALGENCAN_LIB_DIR`. Eight problems change
  status. The JLL solves `HIMMELBJ` and `HS106`, which the source build does
  not; the source build solves `ACOPR30` and `NGONE`, and `LAUNCH` and
  `OPTPRLOC`, which the JLL runs into the 30 minute limit. `AGG` and `CRESC50`
  go from `exception` to `infeasible`. Of the 199 both solve, 197 agree in the
  objective to better than 1e-5, the exceptions being `KISSING` at 2e-4 and
  `POLYGON` at 5e-4.

  This is the closest thing to an apples to apples comparison here, and the
  strongest evidence for the Moré-Sorensen substitution, that being the only
  algorithmic difference between the two builds. It is not spotless: "no status
  regressions" does not survive at this size, four problems go the wrong way.

  **Superseded by the entry on the removed fallback below.** Four of those
  eight, including both timeouts, were an unnecessary safeguard in the patch
  rather than the substitution, and are gone.

- The MA57 version was ruled out as the cause of those four. Rebuilding the
  conventional Algencan against libHSL 2025.7.21, so that both sides use the
  same MA57 and the same METIS 5 and only the patch differs, gives **zero**
  status changes against the 5.2.0 build over all 308, and 6 objectives of 203
  differing at all, one of them above 1e-5. The eight changes above are
  therefore the patch, not the library version.

  With that noise removed the remaining disagreement is 199 common solves, 35
  differing at some level and 2 above 1e-5, `KISSING` at 2e-4 and `POLYGON` at
  1e-4, the largest gap having come down from 5e-4. `POLYGON` is simply a
  sensitive problem: it is also the one that moves between the two conventional
  builds, so most of its earlier disagreement was MA57 version noise rather
  than the substitution. Objectives differing somewhere in the fifth digit on a
  third of the problems is what computing the pivot a different way should look
  like; nothing grows beyond that.

  The recipe, should it be needed again: the pivot patch ports to 2025.7.21
  unchanged in substance, three `RINFO(20) = PIVOT` in `ma57/ma57d.f` at the
  `-5`, `-6` and rank deficient sites, and the `pivot` field plus
  `finfo%pivot = rinfo(20)` in `hsl_ma57/hsl_ma57d.f90`. MA57 needs `mc21`,
  `mc22`, `mc34`, `mc47`, `mc59`, `mc64`, `mc71`, `hsl_zd11` and one of the
  `metis/` wrappers, all of which 5.2.0 had bundled in `ddeps.f`. Compile them
  into `<path>/src`, `ar` them into `<path>/lib/libhsl_ma57.a`, and point
  `MA57PATH` there; leave `METISPATH` unset and link METIS at the final shared
  library step instead, which avoids needing a static METIS.

- **The dogleg fallback has been removed from the patch, and it is what
  `LAUNCH` and `OPTPRLOC` were.** The patch used to return `msinfo = 5` when
  `ls` failed to advance, so that `betra` would take a dogleg direction rather
  than a null step. That was written when the pivot was left at zero and `ls`
  therefore could not advance at all; once the pivot is reconstructed it
  advances normally and the test is redundant. It also pre-empts the
  `iter > maxit` bound that already terminates the loop: instrumented on
  `OPTPRLOC` it fired 71,939,942 times in 240 seconds while the iteration limit
  was reached zero times, a livelock where ALGENCAN had a terminating loop.
  Removing it leaves the patch doing only run-time MA57 detection and the pivot
  reconstruction, which is all it was ever meant to do.

  Rebuilt and swept again, against the same conventional build: seven status
  changes instead of eight, and 203 solved to first order, the same as the
  conventional build. `LAUNCH` goes from a 30 minute timeout to 3.4 seconds and
  `OPTPRLOC` to 2.7 seconds, both matching the conventional objective, and
  `NGONE` returns to `first_order`. `HS106` and `HIMMELBJ`, which the JLL
  solves and the source build does not, are unaffected. Of the 201 both solve,
  200 agree to better than 1e-5, `POLYGON` alone above it at 1e-4, down from
  two problems and 2e-4.

  Two differences appear that were not there before, both mild. `HS99EXP` has
  the same objective to nine digits in every run and only flips status, like
  `ACOPR30`, so it is a feasibility tolerance boundary rather than a different
  answer. `DMN15103` crosses the 30 minute limit, having taken 248 seconds
  conventionally and 920 with the old patch; without the early exit `moresor`
  runs to its natural bound, and the `DMN*` family was already the slow corner
  of this selection.

  What remains, and is genuinely the substitution: `ACOPR30` and `HS99EXP` flip
  status at an unchanged objective, `AGG` and `CRESC50` move between two
  failure statuses, and `HS106` and `HIMMELBJ` are solved that were not.

- The computed Moré-Sorensen pivot agrees with a patched MA57's real
  `finfo%pivot` to 7–11 significant digits over 55 samples. Measured again at
  308 problem scale, by building a conventional Algencan that computes both and
  prints them, that holds for the early iterations and in the well conditioned
  range: `NGONE`, `SWOPF` and `EXPFITA` agree to 1e-13 or better, several
  samples bit identical. It does not hold everywhere. When the true pivot falls
  to the level where summing the quadratic form cancels, the reconstruction
  loses it: on `HS106` MA57 reports 5.68e-14 and the formula returns exactly
  zero. `LAUNCH` sits in between, agreeing to 1e-6 rather than 1e-13. The
  earlier 55 samples evidently missed both tails. This matters less than it
  sounds, since a pivot that small leaves the safeguard it feeds inactive, but
  it is not a formula that reproduces `finfo%pivot` unconditionally.
- Everything HSL related is verified only on `x86_64-linux-gnu`. CI proves the
  other platforms build, not that MA57 works there. Licensed HSL ships no
  Windows, i686, armv6l/v7l or riscv64 build, so those users get truncated
  Newton regardless.

## Picking up on another machine

`~/documentos` syncs, so the repo, the Yggdrasil clone and the licensed HSL are
there. `/tmp` and `~/.julia` are not: the deployed `Algencan_jll`, the CUTEst
setup and every scratch build stay behind. Simplest is to wait for the JLL to be
registered, then `contrib/hsl-check/setup_and_run.sh` builds what it needs.

To test the unregistered build instead, deploy it first:

```bash
cd ~/documentos/programas/Yggdrasil_Algencan/A/Algencan
julia --project=<a project with BinaryBuilder> build_tarballs.jl \
      --verbose --deploy=local x86_64-linux-gnu-libgfortran5
```

then pass `ALGENCAN_JLL_PATH=~/.julia/dev/Algencan_jll` to the sweep. Needs
Julia 1.12 or newer for BinaryBuilder.

Two traps, both of which wasted time here:

- A bare triplet on the BinaryBuilder command line overrides the recipe's
  `platforms`, including `expand_gfortran_versions`, and silently produces a
  libgfortran3 build that cannot be dlopened. Always write
  `x86_64-linux-gnu-libgfortran5`.
- `ALGENCAN_LIB_DIR` is read into a `const` at precompile time in 0.8.1, so
  setting it does nothing if the package is already precompiled. Delete
  `~/.julia/compiled/v1.10/NLPModelsAlgencan` after changing it. Fixed in 0.9.0,
  which resolves the path at call time.
