# Examples of usage

## Using NLPModels

```@example
using ADNLPModels
using NLPModelsAlgencan

# Define the HS52 problem
x0 = [2.0; 2.0; 2.0; 2.0; 2.0]
f(x) = (4 * x[1] - x[2])^2 + (x[2] + x[3] - 2)^2 + (x[4] - 1)^2 + (x[5] - 1)^2
c(x) = [x[1] + 3 * x[2]; x[3] + x[4] - 2 * x[5]; x[2] - x[5]]
lcon = [0.0; 0.0; 0.0]
ucon = [0.0; 0.0; 0.0]
model = ADNLPModel(f, x0, c, lcon, ucon; name="hs52_autodiff")
stats = algencan(model)
@show stats.solution
``` 

## Using JuMP

```@example
using JuMP
using NLPModelsJuMP
using NLPModelsAlgencan

# Define a problem usng JuMP
model = Model(NLPModelsJuMP.Optimizer)
set_attribute(model, "solver", NLPModelsAlgencan.AlgencanSolver)
@variable(model, x[i=1:2], start = [-1.2; 1.0][i])
@objective(model, Min, (x[1] - 1)^2 + 100 * (x[2] - x[1]^2)^2)
@constraint(model, x[1]^2 + x[2]^2 == 1)

# Solve problem
optimize!(model)
@show values.(x)
```

## Using [CUTEst](https://github.com/JuliaSmoothOptimizers/CUTEst.jl) Models

```@example
using NLPModels
using CUTEst
using NLPModelsAlgencan

set_mastsif()

nlp = CUTEstModel("HS6")

stats = algencan(nlp)

print(stats)

finalize(nlp)
```

## Utilizing a specification file

**spec.dat**
```
FEASIBILITY-TOLERANCE 1.0e-5
OPTIMALITY-TOLERANCE 1.0e-5
STAINF-FEASIBILITY-TOLERANCE 3.162278e-03
STAINF-OPTIMALITY-TOLERANCE 3.162278e-08
ACC-FEASIBILITY-THRESHOLD 3.162278e-3
ACC-OPTIMALITY-THRESHOLD 3.162278e-3
ITERATIONS-OUTPUT-DETAIL 10
NUMBER-OF-ARRAYS-COMPONENTS-IN-OUTPUT 0
SKIP-ACCELERATION-PROCESS 1
```

**Code**
```@example
io = open("specfnm.dat", "w") # hide
println(io, "FEASIBILITY-TOLERANCE 1.0e-5") # hide
println(io, "OPTIMALITY-TOLERANCE 1.0e-5") # hide
println(io, "STAINF-FEASIBILITY-TOLERANCE 3.162278e-03") # hide
println(io, "STAINF-OPTIMALITY-TOLERANCE 3.162278e-08") # hide
println(io, "ACC-FEASIBILITY-THRESHOLD 3.162278e-3") # hide
println(io, "ACC-OPTIMALITY-THRESHOLD 3.162278e-3") # hide
println(io, "ITERATIONS-OUTPUT-DETAIL 10") # hide
println(io, "NUMBER-OF-ARRAYS-COMPONENTS-IN-OUTPUT 0") # hide
println(io, "SKIP-ACCELERATION-PROCESS 1") # hide
close(io) # hide
using NLPModels
using CUTEst
using NLPModelsAlgencan

set_mastsif()

nlp = CUTEstModel("HS6")

stats = algencan(nlp, specfnm="specfnm.dat")

print(stats)

finalize(nlp)
```
