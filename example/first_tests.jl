
push!(LOAD_PATH, "../src")
using JuMP
using Ipopt
import Algencan

m = Model(solver=Algencan.AlgencanSolver(OPTIMALITY_TOLERANCE=1.0e-5, FEASIBILITY_TOLERANCE=1.0e-5))
@variable(m, 3.0 <= var[1:2] <= 5.0)
@NLobjective(m, :Min, 123.45*(var[1] - pi)^2 + sin(var[1]) + (var[2] - e)^2 + sin(pi - var[2]))

m2 = Model(solver=Algencan.AlgencanSolver(OPTIMALITY_TOLERANCE=1.0e-5, FEASIBILITY_TOLERANCE=1.0e-5))
# m2 = Model(solver=IpoptSolver())
@variable(m2, 2 <= var2[1:3] <= 5.0)
@NLobjective(m2, :Max, (var2[1] - 4)^2 + (var2[2] - 4)^2 + var2[3])
@constraint(m2, var2[1]  + var2[2] == 5)

return_code = solve(m)
println("Solução do primeiro problema")
println("Código de retorno = ", return_code)
println(getvalue(var))

return_code = solve(m2)
println("Solução do segundo problema ")
println("Código de retorno = ", return_code)
println(getvalue(var2))
println(m2.internalModel.mult)
