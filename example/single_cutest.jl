using MathProgBase
using NLPModels
using CUTEst
using Algencan

"""
Solve a single CUTEst problem given as a parameter in the command line.
"""

solver = AlgencanSolver(epsfeas=1.0e-5, epsopt=1.0e-5,
    efstain=3.162278e-03, eostain=3.162278e-08, efacc=3.162278e-3,
    eoacc=3.162278e-3,
    ITERATIONS_OUTPUT_DETAIL=10, NUMBER_OF_ARRAYS_COMPONENTS_IN_OUTPUT=0)
nlp = CUTEstModel(ARGS[1])
model = NLPtoMPB(nlp, solver)
MathProgBase.optimize!(model)
println("Solver status = ", status(model))
println("Number of function evaluations = ", getnfevals(model))
finalize(nlp)
