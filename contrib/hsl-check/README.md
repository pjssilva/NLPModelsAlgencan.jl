# Checking Algencan against CUTEst, with and without HSL

Runs NLPModelsAlgencan over a list of CUTEst problems and records status,
objective and time for each, so that a build reaching MA57 through a licensed
HSL can be compared against one that falls back to truncated Newton, or against
a conventional Algencan compiled the old way with a locally patched MA57.

This exists because the HSL path is easy to get subtly wrong: a missing LP64
BLAS backend makes MA57 return wrong answers rather than errors, and small test
problems do not notice, because they never reach MA57's blocked BLAS path. Only
a real sweep shows it. See the "HSL support" section of the developer notes.

## Running it

```bash
cd contrib/hsl-check
./setup_and_run.sh
```

The environment is built in `~/.julia/environments/cutest-check` and reused.
Results go to `results-<hostname>.tsv` beside the script, appended and flushed
after every problem, so an interrupted run is still useful and rerunning skips
whatever is already recorded.

To leave it going after disconnecting:

```bash
tmux new-session -d -s cutest './setup_and_run.sh 2>&1 | tee tmux.log'
tmux attach -t cutest          # look in
tmux kill-session -t cutest    # stop it
```

The machine has to stay awake. A laptop that sleeps suspends the run with it.

## Running it in parallel

```bash
CUTEST_CHECK_JOBS=8 CUTEST_CHECK_LIMIT=308 ./setup_and_run.sh
```

Whole problems are handed to separate processes, several at a time. Processes
and not threads, because Algencan keeps state in Fortran common blocks and the
library is loaded and unloaded around every solve; and one process per problem
so that the driver can impose a hard time limit, since a solve stuck inside
Algencan cannot be interrupted from Julia but the process can be killed.
`CUTEST_CHECK_TIMEOUT`, 1800 seconds by default, sets that limit, and problems
that hit it are recorded as `TIMEOUT`.

The parallel run decodes serially first. CUTEst.jl decodes and compiles inside a
single directory of its own and changes into it to do so, so several processes
decoding at once trample each other's Fortran scratch files and fail with
`File cannot be deleted`. Decoding is cached, so this happens once and reruns
skip it.

Julia startup costs a few seconds per problem, which is negligible against
solves that take minutes, and buys robustness: a crash or a hang costs one
problem rather than the run.

## Settings

Both scripts read their configuration from the environment:

| Variable | Default | Meaning |
|---|---|---|
| `CUTEST_CHECK_LIST` | `example/cutest_selection.txt` | problems to run, one per line |
| `CUTEST_CHECK_LIMIT` | `60` | how many of them |
| `CUTEST_CHECK_OUT` | `results-<hostname>.tsv` | where results go |
| `CUTEST_CHECK_TIME_LIMIT` | `600` | warn about problems slower than this |
| `CUTEST_CHECK_ENV` | `~/.julia/environments/cutest-check` | environment to build |
| `ALGENCAN_JLL_PATH` | unset | dev an unregistered `Algencan_jll` instead of taking it from the registry |

The licensed HSL is the first argument, defaulting to
`~/documentos/programas/HSL/HSL_jll.jl/HSL_jll.jl.v2025.7.21`. Without it the
run still works, but `ma57_available` is false and every problem takes the
truncated Newton path, which is a useful comparison in its own right. Pass a
path that does not exist, not an empty string: the default fills in for an
empty argument and you get the licensed HSL back without noticing.

Leaving the argument off is not by itself enough to disable HSL, and used not
to be handled at all. `Pkg.develop` on this repo also carries over whatever
`HSL_jll` is developed in its own `Manifest.toml`, so a licensed HSL can enter
the environment without ever being asked for; a whole sweep then silently
repeats the HSL run while the script prints a warning saying the opposite.
`setup_and_run.sh` now calls `Pkg.free("HSL_jll")` in that branch to put it
back on the public stub. Confirm it worked by reading the `ma57` line of the
results header rather than the warning.

The whole selection is 308 problems; `CUTEST_CHECK_LIMIT=308` runs all of them.

## Reading the output

The header records host, Julia version, `Algencan_jll` and `HSL_jll` versions,
whether MA57 is really reachable, and the BLAS configuration. Check it first.

The `ma57` line is the one that says which solver actually ran. It reads the
`ma57_available` module variable out of the `libhsl_subset` that got loaded,
which is what Algencan itself tests at run time, so `true` means MA57 and
`false` means truncated Newton. Do not infer this from the `HSL_jll` version:
the public stub and the licensed release can carry the same version string,
and both export the MA57 symbols, so neither the version nor `nm` tells them
apart.

Then check the BLAS line: if it shows only an `[ILP64]` backend and the run is
supposed to be using MA57, the numbers are not trustworthy, and `stderr` will
be full of

```
Error: no BLAS/LAPACK library loaded for dgemm_()
```

`NLPModelsAlgencan` forwards an LP64 backend before each solve, so this should
not happen, but it is the first thing to rule out if results look wrong.

To compare two runs:

```bash
join -t $'\t' -j1 \
  <(grep -v '^#' results-A.tsv | sort -k1,1) \
  <(grep -v '^#' results-B.tsv | sort -k1,1) \
  | awk -F'\t' '$2 != $5 { print }'
```

Differing objectives are not automatically a problem. Algencan stops on a
gradient criterion, not on function values, so agreement to about 1e-5 in the
objective is as much as the stopping test guarantees. Changes in `status` are
what matter.

## Known slow problems

Three problems hit the 30 minute limit on an HSL enabled build and finish on
one without it: `LAUNCH`, `OPTPRLOC` and `DMN37143`. Whether they are slow or
stuck is still unknown, but the limit now bounds them either way and they are
recorded as `TIMEOUT` instead of holding up the sweep.

`LAUNCH` is the long standing case, over 40 minutes on HSL against a
conventional build that solves it. `DMN37143` is the least suspicious of the
three: the whole `DMN*` family is slow here, several of them taking 900 to
1500 seconds, so it may simply be over the line rather than stuck.
