#Problem 104 in the Hock-Schittkowski suite
function hs104()

  x0   = [6.0, 3.0, 0.4, 0.2, 6.0, 6.0, 1.0, 0.5]
  f(x) = .4*x[1]^.67*x[7]^(-.67)+.4*x[2]^.67*x[8]^(-.67)+10-x[1]-x[2]
  lvar =  0.1 * ones(8)
  uvar = 10.0 * ones(8)
  c(x) = [1-.0588*x[5]*x[7]-.1*x[1];
          1-.0588*x[6]*x[8]-.1*x[1]-.1*x[2];
          1-4*x[3]/x[5]-2/(x[3]^.71*x[5])-.0588*x[7]/x[3]^1.3;
          1-4*x[4]/x[6]-2/(x[4]^.71*x[6])-.0588*x[8]/x[4]^1.3;
          .4*x[1]^.67*x[7]^(-.67)+.4*x[2]^.67*x[8]^(-.67)+10-x[1]-x[2]
          ]
  lcon = [0.0; 0.0; 0.0; 0.0; .1]
  ucon = [Inf; Inf; Inf; Inf; 4.2]

  return ADNLPModel(f, x0, lvar, uvar, c, lcon, ucon, name="hs104_autodiff")

end
