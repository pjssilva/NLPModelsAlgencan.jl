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
truncated Newton path, which is a useful comparison in its own right.

The whole selection is 308 problems; `CUTEST_CHECK_LIMIT=308` runs all of them.

## Reading the output

The header records host, Julia version, `Algencan_jll` and `HSL_jll` versions,
and the BLAS configuration. Check that first: if it shows only an `[ILP64]`
backend and the run is supposed to be using MA57, the numbers are not
trustworthy, and `stderr` will be full of

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

## Known slow problem

`LAUNCH` did not finish within 40 minutes on an HSL enabled build while a
conventional one solved it. It was never run to completion, so whether it is
slow or stuck is unknown. It is worth watching.
