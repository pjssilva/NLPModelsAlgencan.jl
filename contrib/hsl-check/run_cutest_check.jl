# Runs NLPModelsAlgencan over a list of CUTEst problems and records what came
# out, so that a build using HSL's MA57 can be compared against one that does
# not, or against a conventional Algencan compiled with a locally patched MA57.
#
# It is meant to survive the session that started it: results are appended and
# flushed after every problem, so killing it half way still leaves usable
# output, and rerunning it skips whatever is already recorded.
#
# See contrib/hsl-check/README.md for how to run it.

using Printf

const OUT = get(ENV, "CUTEST_CHECK_OUT",
                joinpath(@__DIR__, "results-$(gethostname()).tsv"))
const LIST = get(ENV, "CUTEST_CHECK_LIST",
                 joinpath(@__DIR__, "..", "..", "example", "cutest_selection.txt"))
const LIMIT = parse(Int, get(ENV, "CUTEST_CHECK_LIMIT", "60"))
const TIME_LIMIT = parse(Float64, get(ENV, "CUTEST_CHECK_TIME_LIMIT", "600"))

using NLPModelsAlgencan, CUTEst, LinearAlgebra, Libdl

# Which HSL is in play, and whether MA57 can actually be reached. Without an
# LP64 backend registered libblastrampoline answers MA57's calls by writing to
# stderr and returning without computing, which produces wrong answers rather
# than errors, so the count of those messages is worth watching.
function environment_banner()
    # NLPModelsAlgencan does this before each solve; do it here too so that the
    # BLAS line below reports what the solves will actually see rather than the
    # bare ILP64 configuration Julia starts with.
    NLPModelsAlgencan.ensure_lp64_blas!()
    hsl, ma57 = "none", "unknown"
    try
        m = Base.require(Base.PkgId(Base.UUID("017b0a0e-03f4-516a-9b91-836bbd1904dd"), "HSL_jll"))
        hsl = string(pkgversion(m))
        # The version alone does not say whether MA57 is really there: the
        # public stub and the licensed release can carry the same version, and
        # both export the MA57 symbols. Only this module variable is .true.
        # when the loaded library actually implements MA57, and it is exactly
        # what Algencan tests at run time to choose between MA57 and the
        # truncated Newton inner solver.
        s = dlsym(dlopen(m.libhsl_subset), :__hsl_ma57_double_MOD_ma57_available)
        ma57 = string(unsafe_load(Ptr{Int32}(s)) != 0)
    catch
    end
    algencan = "unknown"
    try
        m = Base.require(Base.PkgId(Base.UUID("07ede149-d6eb-53b6-8e3c-1a25465d123c"), "Algencan_jll"))
        algencan = string(pkgversion(m))
    catch
    end
    return """
    # host          $(gethostname())
    # julia         $(VERSION)
    # Algencan_jll  $algencan
    # HSL_jll       $hsl
    # ma57          $ma57
    # BLAS          $(BLAS.get_config())
    # started       $(round(Int, time()))
    """
end

already_done() = isfile(OUT) ?
    Set(split(l, '\t')[1] for l in eachline(OUT) if !startswith(l, "#")) : Set{String}()

function main()
    problems = [strip(l) for l in eachline(LIST) if !isempty(strip(l))]
    length(problems) > LIMIT && (problems = problems[1:LIMIT])
    done = already_done()

    if !isfile(OUT)
        open(OUT, "w") do io
            print(io, environment_banner())
            println(io, "# problem\tstatus\tobjective\tseconds")
        end
    end

    for p in problems
        p in done && continue
        t0 = time()
        status, obj = "EXCEPTION", NaN
        nlp = nothing
        try
            nlp = CUTEstModel(p)
            stats = algencan(nlp)
            status, obj = string(stats.status), stats.objective
        catch e
            status = "EXCEPTION:" * first(sprint(showerror, e), 40)
        finally
            nlp === nothing || (try finalize(nlp) catch end)
        end
        el = time() - t0
        open(OUT, "a") do io
            @printf(io, "%s\t%s\t%.12g\t%.1f\n", p, status, obj, el)
        end
        @printf("%-14s %-14s %18.12g  %6.1fs\n", p, status, obj, el)
        flush(stdout)
        el > TIME_LIMIT && @warn "$p took $(round(el, digits=1))s, over the soft limit"
    end

    open(OUT, "a") do io
        println(io, "# finished\t$(round(Int, time()))")
    end
    println("\nDone. Results in $OUT")
end

main()
