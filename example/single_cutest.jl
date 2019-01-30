using Printf
using MathProgBase
MPB = MathProgBase
using NLPModels
using NLPModelsJuMP
using CUTEst
using Algencan

"""
Solve a single CUTEst problem given as a parameter in the command line.
"""
function solve_cutest(pname)
    # solver = AlgencanSolver(epsfeas=1.0e-5, epsopt=1.0e-5,
    #     efstain=3.162278e-03, eostain=3.162278e-08, efacc=3.162278e-3,
    #     eoacc=3.162278e-3,
    #     ITERATIONS_OUTPUT_DETAIL=10, NUMBER_OF_ARRAYS_COMPONENTS_IN_OUTPUT=0)
    solver = AlgencanSolver(epsfeas=1.0e-5, epsopt=1.0e-5, specfnm="algencan.dat")
    nlp = CUTEstModel(pname)
    model = NLPtoMPB(nlp, solver)
    bench_data = @timed status = MPB.optimize!(model)
    finalize(nlp)
    println("Solver status = ", status)
    println("Perfomance data (time, n_fc, n_ggrad, n_hl, n_nlp)")
    n_fc, n_ggrad, n_hl, n_hlp = getnfevals(model)
    @printf("%.4e\t%10d\t%10d\t%10d\t%10d\n", bench_data[2], n_fc, n_ggrad, n_hl, n_hlp)
end

# Calls a simple problem to compile the Julia code in order to get a reasonable
# timing information at the end.
solve_cutest("HS110")
println("\n\n", '*'^40, " Solving ", ARGS[1], "\n")
solve_cutest(ARGS[1])
