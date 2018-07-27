using MathProgBase
using NLPModels
using CUTEst
using Algencan

solver = AlgencanSolver()
nlp = CUTEstModel("HS6")
model = NLPtoMPB(nlp, solver)
MathProgBase.optimize!(model)
finalize(nlp)

solver = AlgencanSolver(epsfeas=1.0e-5, epsopt=1.0e-5,
    efstain=3.162278e-03, eostain=3.162278e-08, efacc=3.162278e-3,
    eoacc=3.162278e-3,
    ITERATIONS_OUTPUT_DETAIL=10, NUMBER_OF_ARRAYS_COMPONENTS_IN_OUTPUT=0)
nlp = CUTEstModel("HS6")
model = NLPtoMPB(nlp, solver)
MathProgBase.optimize!(model)
println("Solver status = ", status(model))
finalize(nlp)
