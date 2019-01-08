using Algencan, JuMP, Test

@testset "Basic model test" begin
  model = Model(solver = AlgencanSolver())
  @variable(model, x[1:2])
  @objective(model, Min, 4 * x[1]^2 + 9 * x[2]^2)
  @constraint(model, c, x[1] + x[2] == 1.0)

  status = solve(model)
  @test status == :Optimal

  x = getvalue(x)
  @test x ≈ [9 / 13; 4 / 13]

  λ = getdual(c)
  @test λ ≈ -72 / 13
end
