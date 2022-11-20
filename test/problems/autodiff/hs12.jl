#Problem 12 in the Hock-Schittkowski suite
function hs12()

  x0 = [0.0, 0.0]
  f(x) = (x[1]^2) / 2 + x[2]^2 - x[1] * x[2] - 7 * x[1] - 7 * x[2]
  c(x) = [-4 * x[1]^2 - x[2]^2]
  lcon = [-25.0]
  ucon = [Inf]

  return ADNLPModel(f, x0, c, lcon, ucon, name="hs12_autodiff")

end
