
push!(LOAD_PATH, "../src")
using JuMP
using Ipopt
using Algencan

m = Model(solver=AlgencanSolver(epsopt=1.0e-5, epsfeas=1.0e-5))
# m = Model(solver=IpoptSolver())
@variable(m, 3.0 <= var[1:2] <= 5.0)
@NLobjective(m, :Min, 123.45*(var[1] - pi)^2 + sin(var[1]) + (var[2] - e)^2 + sin(pi - var[2]))

m2 = Model(solver=Algencan.AlgencanSolver(epsfeas=1.0e-8, epsopt=1.0e-5))
# m2 = Model(solver=IpoptSolver())
@variable(m2, 2.0 <= var2[1:3] <= 5.0)
@NLobjective(m2, :Max, (var2[1] - 4)^2 + (var2[2] - 4)^2 + var2[3])
@constraint(m2, var2[1]  + var2[2] == 5)

status = solve(m)
println("Solution status: $status.")
println("Solution to first problem: ")
println(getvalue(var))

status = solve(m2)
println("Solution status: $status.")
println("Primal-dual pair to second problem: ")
println(getvalue(var2))
println(getconstrduals(internalmodel(m2)))
