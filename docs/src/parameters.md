
# Optional parameters

Algencan has a set of parameters that can adjust the solver algorithm.
* `epsfeas`: a Float, ``\epsilon_{feas}``, to specify the feasibility tolerance of the Algencan main stopping criterion of convergence.
* `epsopt`: a Float, ``\epsilon_{opt}``, to specify the optimality tolerance of the Algencan main stopping criterion of convergence.
* `efstain`: a Float, ``\epsilon_{fstain}``, to specify the feasibility tolerance of the Algencan main stopping criterion to an infeasible point. A possible value for this parameter is ``\sqrt{\epsilon_{feas}}``.
* `eostain`: a Float, ``\epsilon_{ostain}``, to specify the optimality tolerance of the Algencan main stopping criterion to an infeasible point. A possible value for this parameter is ``\epsilon_{opt}^{1.5}``.
* `efacc`: a Float, ``\epsilon_{facc}``, to specify the feasibility tolerance of the Algencan threshold to launch the acceleration process.
* `eoacc`: a Float, ``\epsilon_{oacc}``, to specify the optimality tolerance of the Algencan threshold to launch the acceleration process.
* `outputfnm`: a String to specify the output filename.
* `specfnm`: a String to specify the specification filename, which stores the additional parameters as shown below.

These parameters can be defined as optional keywords in the solver function call, e.g:

```julia
algencan(nlp; epsfeas=1.0e-10, epsopt=1.0e-10, outputfnm='output.dat')
```

## Additional optional parameters

Algencan also provides a set of additional parameters (case-insensitive) that
can be defined in the specification file provided by `specfnm`. The table below
lists these parameters.


Parameter                                    | Additional value       |
:--------------------------------------------|:-----------------------|
SKIP-ACCELERATION-PROCESS                    |                        |
LINEAR-SYSTEMS-SOLVER-IN-ACCELERATION-PROCESS| `String`               |
TRUST-REGIONS-INNER-SOLVER                   | `String` (not required)|
LINEAR-SYSTEMS-SOLVER-IN-TRUST-REGIONS       | `String`               |
NEWTON-LINE-SEARCH-INNER-SOLVER              | `String` (not required)|
LINEAR-SYSTEMS-SOLVER-IN-NEWTON-LINE-SEARCH  | `String`               |
TRUNCATED-NEWTON-LINE-SEARCH-INNER-SOLVER    | `String` (not required)|
MATRIX-VECTOR-PRODUCT-IN-TRUNCATED-NEWTON-LS | `String`               |
FIXED-VARIABLES-REMOVAL-AVOIDED              |                        |
ADD-SLACKS                                   |                        |
OBJECTIVE-AND-CONSTRAINTS-SCALING-AVOIDED    |                        |
IGNORE-OBJECTIVE-FUNCTION                    |                        |
ITERATIONS-OUTPUT-DETAIL                     | `Int`                  |
NUMBER-OF-ARRAYS-COMPONENTS-IN-OUTPUT        | `Int`                  |
SOLUTION-FILENAME                            | `String`               |
ACCELERATION-PROCESS-ITERATIONS-LIMIT        | `Int`                  |
INNER-ITERATIONS-LIMIT                       | `Int`                  |
OUTER-ITERATIONS-LIMIT                       | `Int`                  |
PENALTY-PARAMETER-INITIAL-VALUE              | `Float`                |
LARGEST-PENALTY-PARAMETER-ALLOWED            | `Float`                |

As an illustration, the following expressions may form part of a valid specification
file:

```
# single line comments with '#' are allowed
SKIP-ACCELERATION-PROCESS
ITERATIONS-OUTPUT-DETAIL 10
# the syntax is case-insensitive
inner-iterations-limit 30
```

Furthermore, there is the possibility of these additional parameters being passed by
optional keywords in the solver call. In that case, the parameter's name has the
character `-` substituted by `_`. For instance:

```julia
# Assign 50 to the parameter `OUTER-ITERATIONS-LIMIT`.
algencan(nlp; outer_iterations_limit=50)
```

When there is no additional value associated with the parameter, you can assign an
empty string value to the respective keyword. For example:

```julia
algencan(nlp; skip_acceleration_process='')
```
