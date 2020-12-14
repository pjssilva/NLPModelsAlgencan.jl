using NLPModels, NLPModelsAlgencan, Test

@testset "Basic model test" begin
  f(x) = 4 * x[1]^2 + 9 * x[2]^2
  x0   = zeros(2)
  c(x) = [x[1] + x[2]]
  lcon = [1.0]
  ucon = [1.0]

  nlp = ADNLPModel(f, x0, c, lcon, ucon)

  output = algencan(nlp)
  @test output.status == :first_order

  @test output.solution ≈ [9 / 13; 4 / 13]

  @test output.multipliers ≈ [-72 / 13]
end
