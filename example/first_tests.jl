
push!(LOAD_PATH, "../src")
using JuMP
using Ipopt
import Algencan

m = Model(solver=Algencan.AlgencanSolver())
@variable(m, 3.0 <= var[1:2] <= 5.0)
@NLobjective(m, :Min, 123.45*(var[1] - pi)^2 + sin(var[1]) + (var[2] - e)^2 + sin(pi - var[2]))

m2 = Model(solver=Algencan.AlgencanSolver())
# m2 = Model(solver=IpoptSolver())
@variable(m2, 3.0 <= var2[1:3] <= 5.0)
@NLobjective(m2, :Max, (var2[1] - 4)^2 + (var2[2] - 4)^2 + var2[3])
@constraint(m2, var2[1]  + var2[2] == 7)


solve(m)
println("Solução do primeiro problema")
println(getvalue(var))

solve(m2)
println("Solução do segundo problema ")
println(getvalue(var2))
