push!(LOAD_PATH, "../src/")
using MathProgBase
using NLPModels
using CUTEst
using Algencan

solver = AlgencanSolver(epsfeas=1.0e-5, epsopt=1.0e-5)
nlp = CUTEstModel("HS6")
model = NLPtoMPB(nlp, solver)
MathProgBase.optimize!(model)
