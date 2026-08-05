# Where things stand — session notes, 2026-07-29

Working notes for picking this up on another machine. Temporary: delete once
`Algencan_jll` is registered and 0.9.0 is out. For how the pieces actually work,
read the "Developer notes" page in the manual instead; for the full HSL
investigation, `contrib/yggdrasil/patched_algencan_for_jll.md`.

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
`~/documentos/programas/Yggdrasil_Algencan`, commit `0a45f6fd`, pushed. Builds
Algencan with run time HSL detection. CI green on all platforms, mergeable.

Quiet since 29 July. The description has been rewritten to describe the HSL
aware recipe, since it still described the original one, and a comment points
@imciner2 and @odow at what changed. One question is open, asked of
@amontoison by email rather than in the thread: whether to link `libhsl`
instead of `libhsl_subset`. Today that does not work, because the public
`libhsl` exports no `__hsl_ma57_double_MOD_*` symbols and neither build
exports `ma57_available`, but he could add them to the dummy as he did for
MA86. Switching would change the link line and one bullet of the description,
nothing else.

## What is left

1. Wait for #14301 to merge, then for `Algencan_jll` to appear in General.
2. `Pkg.free("Algencan_jll")` in this repo, re-resolve, re-run the tests. Every
   test so far used a locally deployed JLL, not the registered one.
3. Delete `contrib/yggdrasil/` — the recipe's home is Yggdrasil, and a copy here
   would go stale. Keep `contrib/hsl` and `contrib/hsl-check`.
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
   fragile rather than anything real. `LAUNCH` and `OPTPRLOC` still time out,
   as does `DMN37143`, and all three finish without HSL.

   One trap, which cost a wasted hour-long run here. Leaving the HSL argument
   off does not disable HSL: `Pkg.develop` on this repo also carries over the
   `HSL_jll` developed in its own `Manifest.toml`, so the licensed one arrives
   unasked and the sweep silently repeats the HSL run while printing a warning
   saying the opposite. `setup_and_run.sh` now calls `Pkg.free("HSL_jll")`, and
   the results header carries an `ma57` line saying which solver actually ran.
   Read that line, not the `HSL_jll` version: the public stub and the licensed
   release can share a version string and both export the MA57 symbols.

## The one thing not to forget

`libhsl_subset` is LP64. Julia registers only an ILP64 BLAS backend, and
libblastrampoline answers a call it cannot match by writing a line to stderr and
returning **with the output untouched**. No error, no crash: MA57 factorizes
stale data, healthy matrices come back indefinite, and the solver stops
somewhere wrong. On CUTEst SWOPF that was 91629 factorizations and an infeasible
answer, against 659 and the right one with an LP64 backend forwarded.

This cost most of a day and looked, in turn, like a bug in the patch, in MA57,
in the ABI, and in METIS. If HSL results ever look wrong again, check
`BLAS.get_config()` for an `[LP64]` entry and grep stderr for
`no BLAS/LAPACK library loaded` before suspecting anything else.

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

  The one thing worth returning to is that `LAUNCH`, `OPTPRLOC` and `DMN37143`
  time out at 30 minutes under MA57 and all three finish without it. That is a
  cost question rather than a correctness one, and `DMN37143` is the least
  suspicious of the three, the whole `DMN*` family being slow here.
- The computed Moré-Sorensen pivot agrees with a patched MA57's real
  `finfo%pivot` to 7–11 significant digits over 55 samples.
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
