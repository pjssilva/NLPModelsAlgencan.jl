#!/usr/bin/env bash
# Builds the environment the CUTEst check needs and runs it. Safe to rerun: the
# environment is reused and problems already recorded are skipped.
#
#   ./setup_and_run.sh                      # licensed HSL at the default path
#   ./setup_and_run.sh /path/to/HSL_jll.jl  # or point it somewhere else
#
# To leave it running after you disconnect:
#
#   tmux new-session -d -s cutest './setup_and_run.sh 2>&1 | tee tmux.log'
#   tmux attach -t cutest      # to look in
#   tmux kill-session -t cutest
#
# See README.md.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
HSL="${1:-$HOME/documentos/programas/HSL/HSL_jll.jl/HSL_jll.jl.v2025.7.21}"
ENVDIR="${CUTEST_CHECK_ENV:-$HOME/.julia/environments/cutest-check}"

# A dev'd Algencan_jll, for testing a build that is not registered yet. Leave
# unset to take Algencan_jll from the registry, which is what users will get.
ALGENCAN_JLL="${ALGENCAN_JLL_PATH:-}"

echo "repo:        $REPO"
echo "environment: $ENVDIR"
echo "licensed HSL: ${HSL:-none}"
[ -n "$ALGENCAN_JLL" ] && echo "Algencan_jll: $ALGENCAN_JLL (dev)"

mkdir -p "$ENVDIR"

julia --project="$ENVDIR" -e '
using Pkg
repo, hsl, ajll = ARGS[1], ARGS[2], ARGS[3]
Pkg.develop(path=repo)
isempty(ajll) || Pkg.develop(path=ajll)
if !isempty(hsl) && isdir(hsl)
    Pkg.develop(path=hsl)        # licensed HSL, overriding the public stub
else
    @warn "no licensed HSL at $hsl; MA57 will be unavailable and Algencan will use truncated Newton"
end
Pkg.add("CUTEst")
Pkg.instantiate()
Pkg.status()
' "$REPO" "$HSL" "$ALGENCAN_JLL"

echo
echo "starting the sweep; results are appended after every problem"
echo

# CUTEST_CHECK_JOBS > 1 runs one process per problem, several at a time, with a
# hard time limit each. See run_parallel.sh.
if [ "${CUTEST_CHECK_JOBS:-1}" -gt 1 ]; then
    exec "$HERE/run_parallel.sh"
else
    exec julia --project="$ENVDIR" "$HERE/run_cutest_check.jl"
fi
