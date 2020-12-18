using JuMP
using NLPModels, NLPModelsJuMP, NLPModelsAlgencan

println("Building the first model...\n")
m = Model()
@variable(m, 3.0 <= var[1:2] <= 5.0)
@NLobjective(m, Min, 123.45*(var[1] - π)^2 + sin(var[1]) + (var[2] - ℯ)^2 + sin(π - var[2]))
println(m)

println("Create the bridge to NLPModels..\n")
nlp = MathOptNLPModel(m)

println("Solving...")
status = algencan(nlp, epsfeas=1.0e-5, epsopt=1.0e-5)
println("Solution status: $status.")
print("(Primal) Solution to first problem: ")
println(status.solution, "\n\n")

println("Building second model...\n")
m2 = Model()
@variable(m2, 2.0 <= var2[1:3] <= 5.0)
@NLobjective(m2, Min, (var2[1] - 4)^2 + (var2[2] - 4)^2 + var2[3])
@constraint(m2, -5 <= -(var2[1] + var2[2]))
println(m2)

println("Create the bridge to NLPModels..\n")
cnlp = MathOptNLPModel(m2)

println("Solving...")
status = algencan(cnlp, epsfeas=1.0e-8, epsopt=1.0e-6)

println("Solution status: $status.")
print("Primal solution to second problem: ")
println(status.solution)
print("Multipliers: ")
println(status.multipliers)
println("Full solve time = ", status.elapsed_time, "s")
