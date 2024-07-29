using Printf
using NLPModels
using CUTEst
using NLPModelsAlgencan

"""
Solve a single CUTEst problem given as a parameter in the command line.
"""
function solve_cutest(pname)
    # solver = AlgencanSolver(epsfeas=1.0e-5, epsopt=1.0e-5,
    #     efstain=3.162278e-03, eostain=3.162278e-08, efacc=3.162278e-3,
    #     eoacc=3.162278e-3,
    #     ITERATIONS_OUTPUT_DETAIL=10, NUMBER_OF_ARRAYS_COMPONENTS_IN_OUTPUT=0)
    nlp = CUTEstModel(pname)
    bench_data = @timed status = algencan(nlp)
    finalize(nlp)
    println("Solver status = ", status)
    n_fc, n_ggrad, n_hl, n_hlp = (
        status.solver_specific[:nfc], status.solver_specific[:ngjac], 
        status.solver_specific[:nhl], status.solver_specific[:nhlp]
    )
    println("Perfomance data (time, n_fc, n_ggrad, n_hl, n_nlp)")
    @printf("%.4e\t%10d\t%10d\t%10d\t%10d\n", bench_data[2], n_fc, n_ggrad, n_hl, n_hlp)
end

# Calls a simple problem to compile the Julia code in order to get a reasonable
# timing information at the end.
set_mastsif()
solve_cutest("HS110")
println("\n\n", '*'^40, " Solving ", ARGS[1], "\n")
solve_cutest(ARGS[1])
