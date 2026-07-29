# Decodes CUTEst problems, one after another, so that the parallel solve phase
# does not have to.
#
#   julia --project=<env> decode_all.jl PROBLEM...
#
# CUTEst.jl decodes and compiles inside a single directory of its own,
# `deps/files` in the installed package, and changes into it to do so. That
# directory is not configurable, so several processes decoding at once trample
# each other's Fortran scratch files and fail with
#
#   Fortran runtime error: File cannot be deleted
#
# Decoding serially first, and then solving with `decode=false`, avoids it: the
# compiled problem libraries are cached, so the parallel phase only loads them.
#
# Already decoded problems are skipped, so this is cheap on a rerun.

using CUTEst

failed = String[]
for (i, name) in enumerate(ARGS)
    print("[$i/$(length(ARGS))] $name ... ")
    try
        nlp = CUTEstModel(name)
        finalize(nlp)
        println("ok")
    catch e
        push!(failed, name)
        println("FAILED: ", first(sprint(showerror, e), 60))
    end
    flush(stdout)
end

if !isempty(failed)
    println("\n$(length(failed)) problem(s) could not be decoded and will be ",
            "reported as exceptions by the sweep:")
    foreach(p -> println("  ", p), failed)
end
