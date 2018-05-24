
push!(LOAD_PATH, "../src")
using JuMP
import Algencan

m = Model(solver=Algencan.AlgencanSolver())
@variable(m, 3.0 <= var[1:2] <= 5.0)
@NLobjective(m, Min, 123.45*(var[1] - pi)^2 + sin(var[1]) + (var[2] - e)^2 + sin(pi - var[2]))

m2 = Model(solver=Algencan.AlgencanSolver())
@variable(m2, 3.0 <= var[1:3] <= 5.0)
@NLobjective(m2, Min, 123.45*(var[1] - 1)^2 + (var[2] - 3)^2 + var[3])

solve(m)

solve(m2)
