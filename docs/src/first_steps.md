# First steps

## Using [NLPModels](https://github.com/JuliaSmoothOptimizers/NLPModels.jl) standard modeling interface

NLPModelsAlgencan is naturally built as an NLPModels interface. Hence, one needs
an implementation of an AbstractNLPModel to call the Algencan solver.

```@docs
algencan
```

For instance, the NLPModels package provides an automatic differentiation NLPModel,
called ADNLPModel, using ForwardDiff to compute the derivatives. For further details,
it is recommended to access [this tutorial]
(https://juliasmoothoptimizers.github.io/NLPModels.jl/stable/tutorial/#ADNLPModel-Tutorial).
As an example, the following optimization problem is defined as:

```math
\begin{aligned}
\min \quad & f(x_1, x_2) = x_1x_2 + 5 \\
\textrm{s.t.} \quad & x_1+x_2 \leq 5, \\
& x_1^2+x_2^2=10, \\
& 0 \leq x_1, x_2 \leq 5,
\end{aligned}
```

And could be modeled as an ADNLModel as follows:

```julia
using ADNLPModels, NLPModelsAlgencan

x0   = [1.0; 1.0]
f(x) = x[1]*x[2] + 5
lvar = [0.0; 0.0]
uvar = [5.0; 5.0]
c(x) = [x[1] + x[2]; x[1]^2 + x[2]^2]
lcon = [-Inf; 10.0]
ucon = [5.0;  10.0]

nlp  = ADNLPModel(f, x0, lvar, uvar, c, lcon, ucon)
```

Moreover, it is solved with Algencan and displayed the problem statistics, which is
provided by a return instance of a [GenericExecutionStats]
(https://juliasmoothoptimizers.github.io/SolverTools.jl/dev/api/#SolverTools.GenericExecutionStats):

```@example
using ADNLPModels, NLPModelsAlgencan # hide
x0   = [1.0;1.0] # hide
f(x) = x[1]*x[2] + 5 # hide
lvar = [0.0; 0.0] # hide
uvar = [5.0; 5.0] # hide
c(x) = [x[1]+x[2]; x[1]^2+x[2]^2] # hide
lcon = [-Inf; 10.0] # hide
ucon = [5.0; 10.0] # hide
nlp  = ADNLPModel(f, x0, lvar, uvar, c, lcon, ucon) # hide
stats = algencan(nlp)
print(stats)
```

## Using [JuMP](https://jump.dev/JuMP.jl/stable/) modeling interface

Alternatively, one can use the JuMP interface to model the problem and then
solve it with Algencan. This solution is provided by the [NLPModelsJuMP]
(https://github.com/JuliaSmoothOptimizers/NLPModelsJuMP.jl) package. To
model the same problem as before, it could be done as follows:

```@example
using JuMP, NLPModelsJuMP, NLPModelsAlgencan

# Create a new JuMP Model and set Algencan as solver
model = Model(NLPModelsJuMP.Optimizer)
set_attribute(model, "solver", NLPModelsAlgencan.AlgencanSolver)

# Define the model
@variable(model, 0 ≤ x[1:2] ≤ 5)
set_start_value(x[1], 1.0)
set_start_value(x[2], 1.0)
@objective(model, Min, x[1]*x[2] + 5)
@constraint(model, x[1] + x[2] ≤ 5)
@constraint(model, x[1]^2 + x[2]^2 == 10)

# Solve the model and show the solution
optimize!(model)
@show value.(x)
```
