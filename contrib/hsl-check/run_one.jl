# Solves a single CUTEst problem and appends one TSV line to a results file.
#
#   julia --project=<env> run_one.jl PROBLEM RESULTS.tsv
#
# Kept as its own process so that the driver can impose a hard time limit: a
# solve stuck inside Algencan cannot be interrupted from Julia, but the process
# can be killed. A crash then costs one problem rather than the run.
#
# The result goes to the file rather than to stdout because Algencan prints a
# great deal to stdout, which the driver discards.
#
# Expects the problem to have been decoded already, by decode_all.jl.

using Printf
using NLPModelsAlgencan, CUTEst

function main(name, out)
    t0 = time()
    status, obj = "EXCEPTION", NaN
    nlp = nothing
    try
        # decode=false: decode_all.jl has already built and cached the problem
        # library. Decoding here would put every process into CUTEst.jl's single
        # shared build directory at once, which fails.
        nlp = CUTEstModel(name; decode=false)
        # CUTEST_CHECK_SOLVER picks the linear solver for the Newton line
        # search, as "MA86" or "MA86 MC64". Unset leaves Algencan's own choice,
        # which is MA57 whenever it is available.
        solver = get(ENV, "CUTEST_CHECK_SOLVER", "")
        stats = isempty(solver) ? algencan(nlp) :
                algencan(nlp; NEWTON_LINE_SEARCH_INNER_SOLVER=solver)
        status, obj = string(stats.status), stats.objective
    catch e
        status = "EXCEPTION:" * replace(first(sprint(showerror, e), 40), r"[\t\n]" => " ")
    finally
        nlp === nothing || (try finalize(nlp) catch end)
    end
    # A private file, which the driver appends to the shared results. Writing
    # to the shared file from here looked atomic and was not: buffered IO
    # splits the write, and concurrent processes shredded each other's lines.
    open(out, "w") do io
        @printf(io, "%s\t%s\t%.12g\t%.1f\n", name, status, obj, time() - t0)
    end
end

main(ARGS[1], ARGS[2])
