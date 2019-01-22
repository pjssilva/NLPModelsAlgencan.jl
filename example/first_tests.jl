using JuMP
using Algencan
import MathProgBase
MPB = MathProgBase

println("Building the first model...\n")
m = Model(solver=AlgencanSolver(epsopt=1.0e-5, epsfeas=1.0e-5))
@variable(m, 3.0 <= var[1:2] <= 5.0)
@NLobjective(m, :Min, 123.45*(var[1] - π)^2 + sin(var[1]) + (var[2] - ℯ)^2 + sin(π - var[2]))
println(m)

println("Solving...")
status = solve(m)
println("Solution status: $status.")
println("(Primal) Solution to first problem: ")
println(getvalue(var), "\n\n")

println("Building second model...\n")
m2 = Model(solver=Algencan.AlgencanSolver(epsfeas=1.0e-8, epsopt=1.0e-5))
@variable(m2, 2.0 <= var2[1:3] <= 5.0)
@NLobjective(m2, :Min, (var2[1] - 4)^2 + (var2[2] - 4)^2 + var2[3])
@constraint(m2, -5 <= -(var2[1] + var2[2]))
println(m2)

println("Solving...")
status = solve(m2)
println("Solution status: $status.")
println("Primal-dual solution to second problem: ")
println(getvalue(var2))
println(MPB.getconstrduals(internalmodel(m2)))
println("Full solve time = ", MPB.getsolvetime(internalmodel(m2)), "s")
