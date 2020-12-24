using NLPModels, NLPModelsAlgencan, Test, LinearAlgebra, CUTEst, NLPModelsIpopt

@testset "Autodiff tests" begin
  include("problems/autodiff/hs12.jl")
  include("problems/autodiff/hs52.jl")
  include("problems/autodiff/hs104.jl")

  stats = algencan(hs12())
  @test stats.status     == :first_order
  @test stats.objective   ≈ -30.0 rtol=1.0e-6
  @test stats.solution    ≈ [2.0, 3.0] rtol=1.0e-6
  @test stats.multipliers ≈ [-0.5] rtol=1.0e-6

  stats = algencan(hs52())
  @test stats.status     == :first_order
  @test stats.objective   ≈ 1859/349 rtol=1.0e-6
  @test stats.solution    ≈ [-0.0945558739311759, 0.03151862464668532,
                            0.5157593121810534, -0.4527220628824326,
                            0.031518624649457136] rtol=1.0e-6
  @test stats.multipliers ≈ [3.277936962780388, 2.905444126051696,
                             -7.747851002665478] rtol=1.0e-6

  stats = algencan(hs104())
  @test stats.status     == :first_order
  @test stats.objective   ≈ 3.951163439 rtol=1.0e-6
  @test stats.solution    ≈ [6.465114030593916, 2.232708645245993,
                             0.6673974911544844, 0.5957564229231442,
                             5.932675677556293, 5.527234565979492,
                             1.01332200852064, 0.4006682289818261] rtol=1.0e-6
  @test stats.multipliers ≈ [-2.3596898286423866, -6.205501994097753,
                             -0.9276198475360146, -0.8472008503830812,
                              0.0] rtol=1.0e-6

end
