using ADNLPModels, NLPModelsAlgencan, Test, CUTEst

include("problems/autodiff/hs12.jl")
include("problems/autodiff/hs52.jl")
include("problems/autodiff/hs104.jl")

@testset "Autodiff basic tests" begin
  stats = algencan(hs104())
  @test stats.status == :first_order
  @test stats.objective ≈ 3.951163439 rtol = 1.0e-6
  @test stats.solution ≈ [6.46511403059392, 2.232708645245993,
    0.66739749115448, 0.595756422923144,
    5.93267567755629, 5.527234565979492,
    1.01332200852064, 0.4006682289818261] rtol = 1.0e-6
  @test stats.multipliers ≈ [-2.3596898286423866, -6.205501994097753,
    -0.9276198475360146, -0.8472008503830812,
    0.0] rtol = 1.0e-6
end

@testset "Autodiff params tests" begin
  # Improves precision
  stats = algencan(hs12(); epsfeas=1.0e-13, epsopt=1.0e-13)
  @test stats.status == :first_order
  @test stats.objective ≈ -30.0 rtol = 1.0e-12
  @test stats.solution ≈ [2.0, 3.0] rtol = 1.0e-12
  @test stats.multipliers ≈ [-0.5] rtol = 1.0e-12
end


@testset "Autodiff specification file tests" begin
  stats = algencan(hs52(); specfnm=string(@__DIR__) * "/spec.dat")
  @test stats.status == :first_order
  @test stats.objective ≈ 1859 / 349 rtol = 1.0e-4
  @test stats.solution ≈ [-0.0945558739311759, 0.03151862464668532,
    0.5157593121810534, -0.4527220628824326,
    0.031518624649457136] rtol = 1.0e-4
  @test stats.multipliers ≈ [3.277936962780388, 2.905444126051696,
    -7.747851002665478] rtol = 1.0e-4
end


@testset "CUTEst unconstrained test" begin
  nlp = CUTEstModel("BROWNDEN")
  stats = algencan(nlp)

  @test stats.status == :first_order
  @test stats.objective ≈ 85822.2013524326 rtol = 1.0e-6
  @test stats.solution ≈ [-11.594439894369563, 13.203630043985436,
    -0.40343948846314803, 0.23677877472091655] rtol = 1.0e-6

  finalize(nlp)
end

@testset "CUTEst linear test" begin
  nlp = CUTEstModel("S277-280")
  stats = algencan(nlp)

  @test stats.status == :first_order
  @test stats.objective ≈ 5.076190464 rtol = 1.0e-5
  @test stats.solution ≈ [1.0, 1.0, 1.0, 1.0] rtol = 1.0e-5
  @test stats.multipliers ≈ [-1.0, -1.0, -1.0, -1.0] rtol = 1.0e-5

  finalize(nlp)
end

@testset "CUTEst equality constraints test" begin
  nlp = CUTEstModel("BT9")
  stats = algencan(nlp)

  @test stats.status == :first_order
  @test stats.objective ≈ -1.0 rtol = 1.0e-6
  @test stats.solution ≈ [1.0, 0.0, 0.0, 1.0] rtol = 1.0e-6
  @test stats.multipliers ≈ [-1.0, -1.0] rtol = 1.0e-6

  finalize(nlp)
end

@testset "CUTEst double bounded constraints test" begin
  nlp = CUTEstModel("HS83")
  stats = algencan(nlp;)

  @test stats.status == :first_order
  @test stats.objective ≈ -30665.53867 rtol = 1.0e-6
  @test stats.solution ≈ [78.0, 36.7758129, 29.995256,
    45.0, 33.0] rtol = 1.0e-6
  @test stats.multipliers ≈ [403.26887955845586, 0.0,
    -809.4250334347122] rtol = 1.0e-6

  finalize(nlp)
end

@testset "CUTEst larger problem test" begin
  nlp = CUTEstModel("SWOPF")
  stats = algencan(nlp)

  @test stats.status == :first_order
  @test stats.objective ≈ 0.06786018 rtol = 1.0e-6
  @test stats.solution ≈ [1.10000000, 0.00000000, -0.51719452, 1.06418307,
    -0.02550540, -1.16246134, 1.09999993, 0.00038566, 0.45449012, 1.08906193,
    -0.03016740, 1.18456331, 1.09995970, 0.00941532, 0.64553505, 1.07848290,
    -0.03588621, -1.06848617, 1.01126996, -0.12692320, -0.45745730, 1.87934821,
    -2.23746789, -1.87934821, 2.23746789, 2.06728303, 2.46121468, -2.05703806,
    -2.33314191, 0.05986148, 0.25775443, -0.05151223, 0.03364814, 0.05712945,
    -0.27582469, -0.05676256, -0.03444036, 0.00871582, 0.36724047, -0.00370526,
    -0.17344842, -0.00009139, -0.39103339, 0.00119723, 0.18900785, 2.60397993,
    -0.69449373, -2.60397994, 0.69449373, 2.86410992, 0.76494730, -2.85684648,
    -0.67779138, 0.77773350, -0.28671072, -0.75636917, 0.57235585, 0.85564926,
    0.28878353, -0.83753866, -0.48280548, 1.94858398, -0.73248917, -1.94858398,
    0.73248917, 2.13646724, 0.82405512, -2.12780075, -0.72004976, 1.08313902,
    1.03877642, 1.16441318, 1.21000000, 1.18696596, 1.21000000, 1.13313613,
    1.21000000, 2.06728303, 3.46410992, 2.63646724, 2.46121468, 0.84494730,
    0.87405512] rtol = 1.0e-6

  finalize(nlp)
end
