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

solver = AlgencanSolver()
nlp = CUTEstModel("A4X12")
model = NLPtoMPB(nlp, solver)
MathProgBase.optimize!(model)
finalize(nlp)
