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

solver = AlgencanSolver(ITERATIONS_OUTPUT_DETAIL=10,epsfeas=1.0e-5, epsopt=1.0e-5)
nlp = CUTEstModel("MPC7")
model = NLPtoMPB(nlp, solver)
MathProgBase.optimize!(model)
println("Solver status = ", status(model))
finalize(nlp)
