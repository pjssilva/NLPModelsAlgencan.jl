#Problem 52 in the Hock-Schittkowski suite
function hs52()

  x0   = [2.0; 2.0; 2.0; 2.0; 2.0]
  f(x) = (4*x[1]-x[2])^2 + (x[2]+x[3]-2)^2 + (x[4]-1)^2 + (x[5]-1)^2
  c(x) = [x[1] + 3*x[2]; x[3] +   x[4] - 2*x[5]; x[2] -   x[5]]
  lcon = [0.0; 0.0; 0.0]
  ucon = [0.0; 0.0; 0.0]

  return ADNLPModel(f, x0, c, lcon, ucon, name="hs52_autodiff")

end
