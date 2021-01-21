# Examples of usage

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

**specfnm.dat**
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
