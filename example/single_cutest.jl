push!(LOAD_PATH, "../src/")
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
    efstain=sqrt(1.0e-5), eostain=1.0e-5^1.5, efacc=sqrt(1.0e-5),
    eoacc=sqrt(1.0e-5),
    ITERATIONS_OUTPUT_DETAIL=15, NUMBER_OF_ARRAYS_COMPONENTS_IN_OUTPUT=0)
nlp = CUTEstModel("MPC7")
model = NLPtoMPB(nlp, solver)
MathProgBase.optimize!(model)
println("Solver status = ", status(model))
finalize(nlp)
