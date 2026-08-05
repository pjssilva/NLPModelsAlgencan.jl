#!/usr/bin/env bash
# Runs the CUTEst sweep across several processes, one problem per process.
#
# Processes rather than threads because Algencan keeps state in Fortran common
# blocks, and one process per problem because a solve stuck inside Algencan
# cannot be interrupted from Julia -- the process can be killed, so a hung
# problem costs a slot rather than the run.
#
# Expects the environment built by setup_and_run.sh. Usually invoked through it
# with CUTEST_CHECK_JOBS set, but can be run directly.

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENVDIR="${CUTEST_CHECK_ENV:-$HOME/.julia/environments/cutest-check}"
LIST="${CUTEST_CHECK_LIST:-$HERE/../../example/cutest_selection.txt}"
LIMIT="${CUTEST_CHECK_LIMIT:-60}"
OUT="${CUTEST_CHECK_OUT:-$HERE/results-$(hostname).tsv}"
JOBS="${CUTEST_CHECK_JOBS:-4}"
# Hard limit per problem. LAUNCH ran over 40 minutes without finishing once.
TIMEOUT="${CUTEST_CHECK_TIMEOUT:-1800}"

if [ ! -f "$OUT" ]; then
    {
        echo "# host          $(hostname)"
        echo "# jobs          $JOBS"
        echo "# timeout       ${TIMEOUT}s per problem"
        echo "# started       $(date -Iseconds)"
        julia --project="$ENVDIR" -e '
            using NLPModelsAlgencan, LinearAlgebra, Libdl
            NLPModelsAlgencan.ensure_lp64_blas!()
            # Deliberately not a loop: assigning to an outer variable inside a
            # top level for loop lands in soft scope, so hsl would stay nothing
            # and the ma57 line below would always read unknown.
            jll(n, u) = try Base.require(Base.PkgId(Base.UUID(u), n)) catch; nothing end
            ver(m) = m === nothing ? "none" : string(pkgversion(m))
            algencan = jll("Algencan_jll", "07ede149-d6eb-53b6-8e3c-1a25465d123c")
            hsl      = jll("HSL_jll", "017b0a0e-03f4-516a-9b91-836bbd1904dd")
            println("# ", rpad("Algencan_jll", 13), " ", ver(algencan))
            println("# ", rpad("HSL_jll", 13), " ", ver(hsl))
            # The version alone does not say whether MA57 is really there: the
            # public stub and the licensed release can carry the same version,
            # and both export the MA57 symbols. Only this module variable is
            # .true. when the loaded library actually implements MA57, and it
            # is exactly what Algencan tests at run time to choose between
            # MA57 and the truncated Newton inner solver.
            println("# ", rpad("ma57", 13), " ", try
                s = dlsym(dlopen(hsl.libhsl_subset), :__hsl_ma57_double_MOD_ma57_available)
                unsafe_load(Ptr{Int32}(s)) != 0
            catch; "unknown" end)
            println("# BLAS          ", BLAS.get_config())'
        echo -e "# problem\tstatus\tobjective\tseconds"
    } > "$OUT"
fi

# Skip anything already recorded, so a rerun resumes.
todo=$(grep -v '^#' "$OUT" | cut -f1 | sort -u > /tmp/.cutest-done.$$;
       head -n "$LIMIT" "$LIST" | tr -d ' \r' | grep -v '^$' | sort -u |
       comm -23 - /tmp/.cutest-done.$$)
rm -f /tmp/.cutest-done.$$

n=$(printf '%s\n' "$todo" | grep -c . || true)
echo "$n problems to run, $JOBS at a time, ${TIMEOUT}s limit each"
[ "$n" -eq 0 ] && { echo "nothing to do"; exit 0; }

# Decode first, serially. CUTEst.jl decodes inside one shared directory and
# changes into it, so parallel decoding fails; solving afterwards with
# decode=false only loads the cached libraries. Skips what is already built.
echo
echo "decoding (serial, cached between runs)..."
printf '%s\n' "$todo" | xargs julia --project="$ENVDIR" "$HERE/decode_all.jl"
echo

printf '%s\n' "$todo" | xargs -P "$JOBS" -I@ bash -c '
    out="$1"; env="$2"; here="$3"; t="$4"; p="$5"
    # A directory of its own: the SIF decoder writes fort.NN scratch files into
    # the working directory, and parallel runs sharing one collide.
    wd=$(mktemp -d "${TMPDIR:-/tmp}/cutest-$p-XXXXXX")
    timeout "$t" julia --project="$env" "$here/run_one.jl" "$p" "$wd/line" \
        >/dev/null 2>&1
    rc=$?
    # Each process writes its own file; only the shell appends to the shared
    # one, in a single small write. Letting several Julia processes append
    # directly shredded lines, because buffered IO splits the write.
    if [ -s "$wd/line" ]; then
        cat "$wd/line" >> "$out"
    elif [ $rc -eq 124 ]; then
        printf "%s\tTIMEOUT\tNaN\t%s\n" "$p" "$t" >> "$out"
    else
        printf "%s\tCRASHED(rc=%s)\tNaN\t0\n" "$p" "$rc" >> "$out"
    fi
    rm -rf "$wd"
    tail -n 1 "$out"
' _ "$OUT" "$ENVDIR" "$HERE" "$TIMEOUT" @

echo "# finished       $(date -Iseconds)" >> "$OUT"
echo
echo "Done. Results in $OUT"
