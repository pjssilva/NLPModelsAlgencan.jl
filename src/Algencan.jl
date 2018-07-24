__precompile__()

"""
Algencan interface to MathProgBase and JuMP.

See its [GitHub page](https://github.com/pjssilva/Algencan.jl)
"""
module Algencan

# TODO: This looks like things to allow for automatic download and
#       compilation of dependencies. Deal with it later.
# if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
#     include("../deps/deps.jl")
# else
#     error("Ipopt not properly installed. Please run Pkg.build(\"Ipopt\")")
# end
# amplexe = joinpath(dirname(libipopt), "..", "bin", "ipopt")

# Compiles to the *static* path of the algencan library
if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
    const algencan_lib_path = libalgencan
else
    const algencan_lib_path = string(joinpath(ENV["ALGENCAN_LIB_DIR"],
        "libalgencan.so"))
end


# Global to avoid closures as Algencan does not allow to send user information
# back to call backs. Long names to avoid conflicts.
global current_algencan_problem

# Standard LP interface
import MathProgBase
importall MathProgBase.SolverInterface
import Base.copy

###############################################################################
# Solver objects
export AlgencanSolver

"Algencan solver that stores options"
struct AlgencanSolver <: AbstractMathProgSolver
    options

    function AlgencanSolver(kwargs)
        options = Dict(
            :epsfeas=>1.0e-08,
            :epsopt=>1.0e-08,
            :efstain=>sqrt.(1.0e-8),
            :eostain=>(1.0e-8)^1.5,
            :efacc=>sqrt(1.0e-8),
            :eoacc=>sqrt(1.0e-8),
            :outputfnm=>"",
            :specfnm=>""
        )
        for i in kwargs
            key, value = i
            options[key] = value
        end
        solver = new(options)
    end
end
"Store keyworded argments as options"
AlgencanSolver(;kwargs...) = AlgencanSolver(kwargs)

"Algencan model, that storing solution data"
mutable struct AlgencanMathProgModel <: AbstractNonlinearModel
    # Problem data
    n::Int                          # Num variables
    lb::Vector{Float64}             # Lower bounds on vars
    ub::Vector{Float64}             # Upper bounds on var
    sense::Float64                  # 1.0 for :Min and -1 for :Max
    m::Int                          # Num constraints
    g_ub::Vector{Float64}           # Upper bounds on constraints
    g_lb::Vector{Float64}           # Lower bound on constrains
    g_sense::Vector{Int}
    g_two_sides::Vector{Bool}
    g_two_smap::Vector{Int}
    g_has_lb::Bool                  # true if at least one constraint has lower
                                    # bound
    evaluator::AbstractNLPEvaluator # Evaluator for functions
    j_row_inds::Vector{Int}         # NNZ row indexes of sparse const. Jacobian.
    j_col_inds::Vector{Int}         # NNZ col indexes of sparse const. Jacobian.
    h_row_inds::Vector{Int}         # NNZ row indexes of sparse Lag. Hessian.
    h_col_inds::Vector{Int}         # NNZ col indexes of sparse Lag. Hessian.
    is_equality::Vector{UInt8}      # 1 if equality, 0 if inequality
    is_g_linear::Vector{UInt8}      # 1 if constraint is linear, 0 otherwise

    # Decision variables and constraints values.
    x::Vector{Float64}              # Starting and final solution
    # TODO: See how to allow for initial multipliers
    g::Vector{Float64}              # Final constraint values
    mult::Vector{Float64}           # Final Lagrange multipliers on constraints
    obj_val::Float64                # Final objective
    status::Symbol                  # Final status
    solve_time::Float64             # Total solution time

    # Options to be passed to the solver
    options

    function AlgencanMathProgModel(options)
        model = new()

        # Set up model with dummy values to mark that is not initialized
        model.n, model.m = 0, 0
        model.lb, model.ub = Float64[], Float64[]
        model.g_lb, model.g_ub = Float64[], Float64[]
        model.g_sense, model.g_two_sides = Float64[], Bool[]
        model.g_has_lb = false
        model.g_two_smap = Int[]
        model.sense = 1.0
        model.j_row_inds, model.j_col_inds = Int[], Int[]
        model.h_row_inds, model.h_col_inds = Int[], Int[]
        model.x, model.g, model.mult = Float64[], Float64[], Float64[]
        model.obj_val, model.status = 0.0, :Undefined
        model.options = options
        model
    end
end
NonlinearModel(s::AlgencanSolver) = AlgencanMathProgModel(s.options)
LinearQuadraticModel(s::AlgencanSolver) = NonlinearToLPQPBridge(NonlinearModel(s))

###############################################################################
# Begin interface implementation
# TODO: Verify if I need to explicitly export something or not

# Simple functions
"Return the current primal point, the solution after calling optimize"
getsolution(model::AlgencanMathProgModel) = model.x
export getsolution

"Get objective value at current primal point"
getobjval(model::AlgencanMathProgModel) = model.sense*model.obj_val
export getobjval

"Get current model status"
status(model::AlgencanMathProgModel) = model.status
export status

"Get best bound on (local) optimal value"
getobjbound(m::AlgencanMathProgModel) = model.status == :optimal ? getobjval(m) : m.sense*Inf
export getobjbound

"Get gap to (local) optimality"
getobjgap(m::AlgencanMathProgModel) = model.status == :Optimal ? 0.0 : Inf
export getobjgap

"There is no inner solver, all functionality is exposed by the default interface"
getrawsolver(m::AlgencanMathProgModel) = nothing
export getrawsolver

"Get the solution time"
getsolvetime(m::AlgencanMathProgModel) = m.solve_time
export getsolvetime

"Change optimization sense, either :Min or :Max"
function setsense!(m::AlgencanMathProgModel, sense)
    m.sense = sense == :Max ? -1 : 1
end
export setsense!

"Return problem sense"
getsense(model::AlgencanMathProgModel) = (model.sense == 1 ? :Min : :Max)
export getsense

"Return the number of decision variables"
numvar(model::AlgencanMathProgModel) = model.n
export numvar

"Return the number of constraints"
numconstr(model::AlgencanMathProgModel) = model.m
export numconstr

"Return a copy of the current model"
copy(m::AlgencanMathProgModel) = deepcopy(m)
export copy

"Algencan only deals with continuous variable, so inform if this is not the case"
 function setvartype!(m::AlgencanMathProgModel, v::Vector{Symbol})
     if !all(v .== :Cont)
         throw("Algencan ony deals with continuous variables")
     end
 end
export setvartype!

"Return variable types: they are all continuous"
function getvartype(model::AlgencanMathProgModel)
    types = Vector{Symbol}(model.m)
    types .= :Cont
    return types
end
export getvartype

"""
Set parameter for a solver, that will be default for all next models

You can use any Algencan parameter that can be set in a specification file plus
:epsfeas, :epsopt, :efstain, :eostain, efacc, :eoac, :outputfnm, :specfnm.

See more details in Algencan documentation.
"""
function setparameters!(m::Union{AlgencanSolver, AlgencanMathProgModel};
    kwargs...)

    for (key, value) in kwargs
        m.options[key] = value
    end
end
export setparameters

"Set an initial value for the decision variables (warmstart)"
function setwarmstart!(model::AlgencanMathProgModel, x)
    model.x = copy(x)
end
export setwarmstart!

function getconstrduals(model::AlgencanMathProgModel)
    v = model.mult
    scale!(v, model.sense)
    return v
end
export getconstrduals

"Return mutipliers associated with bound cosntraints"
function getreducedcosts(model::AlgencanMathProgModel)
    if model.status != :Optimal
        throw("Cannot compute reduced costs: problem not solved to optimality")
    end

    # Compute the gradient of the Lagranfian
    # TODO: It seems like the eval_jac_prod_t is not implemented for
    # the JuMP evaluator or something. I need to implement this
    # using eval_jac_ǵ
    return zeros(model.n)

    # # Objective function
    # grad_lag = Vector{Float64}(model.n)
    # MathProgBase.eval_grad_f(model.evaluator, grad_lag, model.x)
    # grad_lag *= model.sense
    #
    # # Add the multiplier combination of the contraints gradients
    # grad_g = Vector{Float64}(model.n)
    # MathProgBase.eval_jac_prod_t(model.evaluator, grad_g, model.x, model.mult)
    #
    # grad_lag .+= grad_g
    #
    # reduced_costs = zeros(model.n)
    # in_lower = model.x .== model.lb
    # reduced_costs[in_lower] = max.(0.0, grad_lag[in_lower])
    # in_upper = model.x .== model.ub
    # reduced_costs[in_upper] = min.(0.0, grad_lag[in_upper])
    # reduced_costs .*= model.sense
    #
    # return reduced_costs
end
export getreducedcosts

# Simple function only defined for Algencan models

"Set an inital value for the constaint mutipliers (warmstart)"
function setmultwarmstart!(model::AlgencanMathProgModel, mult)
    model.mult = copy(mult)
end
export setmultwarmstart!

# More complex funcitons

"Loads the problem with its basic data and functions in a NLPEvaluator"
function loadproblem!(model::AlgencanMathProgModel, numVar::Integer,
    numConstr::Integer, x_l, x_u, g_lb, g_ub, sense::Symbol,
    d::AbstractNLPEvaluator)

    @assert sense == :Min || sense == :Max

    # Link the model with the problem representation
    # Initialize the evaluator with the right features
    features = features_available(d)
    has_hessian = (:Hess in features)
    init_feat = [:Grad]
    has_hessian && push!(init_feat, :Hess)
    numConstr > 0 && push!(init_feat, :Jac)
    initialize(d, init_feat)

    # Copy data
    model.n = numVar
    model.m = numConstr
    model.lb, model.ub = float(x_l), float(x_u)
    g_lb, g_ub = float(g_lb), float(g_ub)
    model.g_sense, model.g_two_sides, model.g_two_smap = treat_lower_bounds(
        g_lb, g_ub)
    model.g_has_lb = (model.m > 0 && (minimum(model.g_sense) == -1 ||
        maximum(model.g_two_sides)))
    g_only_low = (model.g_sense .== -1.0)
    model.g_lb, model.g_ub = g_lb, g_ub
    # println("lb = ", g_lb[model.g_two_sides])
    # println("ub = ", g_ub[model.g_two_sides])

    model.g_ub[g_only_low] = -g_lb[g_only_low]
    model.g_lb[g_only_low] = -Inf

    model.sense = sense == :Min ? 1.0 : -1.0
    model.evaluator = d

    # Constraints types
    model.is_equality = zeros(UInt8, numConstr)
    model.is_equality[model.g_lb .== model.g_ub] = 1
    model.is_g_linear = zeros(UInt8, numConstr)
    for i in 1:numConstr
        if MathProgBase.isconstrlinear(d, i)
            model.is_g_linear[i] = 1
        end
    end

    # Get strutural indices of Jacobian and Hessian.
    j_row_inds, j_col_inds = MathProgBase.jac_structure(d)
    h_row_inds, h_col_inds = MathProgBase.hesslag_structure(d)

    # C indices start in 0
    model.j_row_inds, model.j_col_inds = j_row_inds - 1, j_col_inds - 1
    model.h_row_inds, model.h_col_inds = h_row_inds - 1, h_col_inds - 1

    # Initial values
    model.x = zeros(numVar)
    model.g = zeros(numConstr)
    model.mult = zeros(numConstr)
    model.obj_val, model.status = 0.0, :Undefined
end
export loadproblem!

"Analyse the lower and upper bounds on the constraints and prepare the
 data structure to treat lower bounds."
function treat_lower_bounds(lb, ub)
    m = length(lb)
    sense = ones(m)
    only_lower = (-Inf .< lb) .& (ub .== Inf)
    sense[only_lower] = -1.0

    two_sides = -Inf .< lb .< ub .< Inf

    new_ind = 1
    two_smap = zeros(m)
    for i = 1:m
        if two_sides[i]
            two_smap[i] = new_ind
            new_ind += 1
        end
    end

    return sense, two_sides, two_smap
end

###########################################################################
# Algencan callbacks
###########################################################################

"Compute objective and constraints as required by Algencan"
function julia_fc(n::Cint, x_ptr::Ptr{Float64}, obj_ptr::Ptr{Float64},
    m::Cint, g_ptr::Ptr{Float64}, flag_ptr::Ptr{Cint})
    model::AlgencanMathProgModel = current_algencan_model

    # Evaluate objective and constraints
    x = unsafe_wrap(Array, x_ptr, Int(n))
    obj_val = MathProgBase.eval_f(model.evaluator, x)
    unsafe_store!(obj_ptr, model.sense*obj_val)
    g = unsafe_wrap(Array, g_ptr, Int(m))
    MathProgBase.eval_g(model.evaluator, g, x)

    # Treat lower bounds and two-sided constraints
    if model.g_has_lb
        first_g = view(g, 1:model.m)
        first_g .*= model.g_sense
        g[model.m + 1:m] = -first_g[model.g_two_sides] +
            model.g_lb[model.g_two_sides]
        first_g .-= model.g_ub
    else
        g .-= model.g_ub
    end

    # Report that evaluation of successful
    unsafe_store!(flag_ptr, Cint(0))
    nothing
end

"Compute objective gradient and constraints Jacobian as required by Algencan"
function julia_gjac(n::Cint, x_ptr::Ptr{Float64}, f_grad_ptr::Ptr{Float64},
    m::Cint, jrow_ptr::Ptr{Cint}, jcol_ptr::Ptr{Cint},
    jval_ptr::Ptr{Float64}, jnnz_ptr::Ptr{Cint}, lim::Cint,
    lmem_ptr::Ptr{UInt8}, flag_ptr::Ptr{Cint})
    model::AlgencanMathProgModel = current_algencan_model

    # Compute gradient of the objective
    x = unsafe_wrap(Array, x_ptr, Int(n))
    f_grad = unsafe_wrap(Array, f_grad_ptr, Int(n))
    MathProgBase.eval_grad_f(model.evaluator, f_grad, x)
    scale!(f_grad, model.sense)

    # Find structure of the constraints Jacobian
    nnz = length(model.j_row_inds)
    if nnz > Int(lim)
        unsafe_store!(lmem_ptr, Cint(1))
        unsafe_store!(flag_ptr, Cint(1))
        return nothing
    else
        unsafe_store!(lmem_ptr, Cint(0))
    end
    jcol_ind = unsafe_wrap(Array, jcol_ptr, Int(lim))
    jrow_ind = unsafe_wrap(Array, jrow_ptr, Int(lim))
    jrow_ind[1:nnz] .= model.j_row_inds
    jcol_ind[1:nnz] .= model.j_col_inds

    # Compute the constraints Jacobian
    J = unsafe_wrap(Array, jval_ptr, Int(lim))
    MathProgBase.eval_jac_g(model.evaluator, J, x)

    # Treat the presence of lower bound in the constraints
    if model.g_has_lb
        @inbounds for i = 1:length(model.j_row_inds)
            # +1, -1 to translate from C indexing to Julia indexing
            rind, cind = model.j_row_inds[i] + 1, model.j_col_inds[i]
            if model.g_two_sides[rind]
                nnz += 1
                jrow_ind[nnz] = model.m + model.g_two_smap[rind] - 1
                jcol_ind[nnz] = cind
                J[nnz] = -J[i]
            else
                J[i] *= model.g_sense[rind]
            end
        end
    end
    unsafe_store!(jnnz_ptr, nnz)

    # Declare success
    unsafe_store!(flag_ptr, Cint(0))
    nothing
end

"Compute the Hessian of the Lagrangian as required by Algencan"
function julia_hl(n::Cint, x_ptr::Ptr{Float64}, m::Cint,
    mult_ptr::Ptr{Float64}, scale_f::Float64, scale_g_ptr::Ptr{Float64},
    hrow_ptr::Ptr{Cint}, hcol_ptr::Ptr{Cint},
    hval_ptr::Ptr{Float64}, hnnz_ptr::Ptr{Cint}, lim::Cint,
    lmem_ptr::Ptr{UInt8}, flag_ptr::Ptr{Cint})
    model::AlgencanMathProgModel = current_algencan_model

    # Get nonzero indexes.
    nnz = length(model.h_row_inds)
    if nnz > Int(lim)
        unsafe_store!(lmem_ptr, 1)
        unsafe_store!(flag_ptr, 1)
        return nothing
    else
        unsafe_store!(lmem_ptr, 0)
    end
    unsafe_store!(hnnz_ptr, Cint(nnz))
    hcol_ind = unsafe_wrap(Array, hcol_ptr, Int(lim))
    hrow_ind = unsafe_wrap(Array, hrow_ptr, Int(lim))
    hrow_ind[1:nnz] .= model.h_row_inds
    hcol_ind[1:nnz] .= model.h_col_inds

    # Compute scaled multipliers
    σ = scale_f*model.sense
    alg_mult = unsafe_wrap(Array, mult_ptr, Int(m))
    scale_g = unsafe_wrap(Array, scale_g_ptr, Int(m))
    if !model.g_has_lb
        μ = alg_mult .* scale_g
    else
        μ = alg_mult[1:model.m] .* scale_g[1:model.m]
        μ[model.g_two_sides] -= scale_g[model.m + 1:m] .* alg_mult[model.m + 1:m]
    end

    # Evaluate the Hessian
    x = unsafe_wrap(Array, x_ptr, Int(n))
    H = unsafe_wrap(Array, hval_ptr, Int(lim))
    MathProgBase.eval_hesslag(model.evaluator, H, x, σ, μ)
    unsafe_store!(flag_ptr, Cint(0))
    nothing
end

"Compute the Hessian of the Lagrangian times p as required by Algencan"
function julia_hlp(n::Cint, x_ptr::Ptr{Float64}, m::Cint,
    mult_ptr::Ptr{Float64}, scale_f::Float64, scale_g_ptr::Ptr{Float64},
    p_ptr::Ptr{Float64}, hp_ptr::Ptr{Float64}, goth_ptr::Ptr{UInt8},
    flag_ptr::Ptr{Cint})
    model::AlgencanMathProgModel = current_algencan_model

    # Compute scaled multipliers
    σ = scale_f*model.sense
    alg_mult = unsafe_wrap(Array, mult_ptr, Int(m))
    scale_g = unsafe_wrap(Array, scale_g_ptr, Int(m))
    if !model.g_has_lb
        μ = alg_mult .* scale_g
    else
        μ = alg_mult[1:model.m] .* scale_g[1:model.m]
        μ[model.g_two_sides] -= scale_g[model.m + 1:m] .* alg_mult[model.m + 1:m]
    end

    # Evaluate Hessian times p
    x = unsafe_wrap(Array, x_ptr, Int(n))
    p = unsafe_wrap(Array, p_ptr, Int(n))
    hp = unsafe_wrap(Array, hp_ptr, Int(n))
    MathProgBase.eval_hesslag_prod(model.evaluator, hp, x, p,
        σ, μ)
    unsafe_store!(flag_ptr, Cint(0))
    nothing
end

function optimize!(model::AlgencanMathProgModel)
    # TODO: Allow warm start primal and specially dual
    #copy!(model.inner.x, model.warmstart) # set warmstart
    # TODO: No options for now
    # for (name,value) in model.options
    #     sname = string(name)
    #     if match(r"(^resto_)", sname) != nothing
    #         sname = replace(sname, r"(^resto_)", "resto.")
    #     end
    #     addOption(model.inner, sname, value)
    # end
    global current_algencan_model = model

    tic()
    ###########################################################################
    # Algencan callback function wrappers
    ###########################################################################
    const c_julia_fc = cfunction(julia_fc, Void, (Cint, Ptr{Float64},
        Ptr{Float64}, Cint, Ptr{Float64}, Ptr{Cint}))

    const c_julia_gjac = cfunction(julia_gjac, Void, (Cint, Ptr{Float64},
        Ptr{Float64}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Cint,
        Ptr{UInt8}, Ptr{Cint}))

    const c_julia_hl = cfunction(julia_hl, Void, (Cint, Ptr{Float64}, Cint,
        Ptr{Float64}, Float64, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, Ptr{Float64},
        Ptr{Cint}, Cint, Ptr{UInt8}, Ptr{Cint}))

    const c_julia_hlp = cfunction(julia_hlp, Void, (Cint, Ptr{Float64}, Cint,
            Ptr{Float64}, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
            Ptr{UInt8}, Ptr{Cint}))

    # Call Algencan. I will do it slowly, first I create variables for all
    # Algencan's parameters, that are a lot, and then I call it.

    # Basic callbacks
    myevalf = C_NULL
    myevalg = C_NULL
    myevalh = C_NULL
    myevalc = C_NULL
    myevaljac = C_NULL
    myevalhc = C_NULL
    myevalfc = c_julia_fc
    myevalgjac = c_julia_gjac
    myevalgjacp = C_NULL
    myevalhl = c_julia_hl
    myevalhlp = c_julia_hlp
    jcnnzmax = 2*length(model.j_row_inds)
    # TODO: Compute the right value here
    hnnzmax = 10*(length(model.h_row_inds) + 10*jcnnzmax)
    coded = zeros(UInt8, 11)
    coded[7] = UInt8(1)
    coded[8] = UInt8(1)
    coded[10] = UInt8(1)
    coded[11] = UInt8(0)
    checkder = UInt8(0)

    # Deal with lower bounds
    m = model.m + sum(model.g_two_sides)
    mult = zeros(m)
    is_equality = zeros(UInt8, m)
    is_equality[1:model.m] .= model.is_equality
    is_g_linear = zeros(UInt8, m)
    is_g_linear[1:model.m] .= model.is_g_linear
    is_g_linear[model.m + 1:m] .= model.is_g_linear[model.g_two_sides]

    # Parameters controling precision
    epsfeas = [model.options[:epsfeas]]
    epsopt = [model.options[:epsopt]]
    efstain = [model.options[:efstain]]
    eostain = [model.options[:eostain]]
    efacc  = [model.options[:efacc]]
    eoacc  = [model.options[:eoacc]]

    # Extra parameters
    outputfnm = model.options[:outputfnm]
    specfnm   = model.options[:specfnm]
    vparam = option2vparam(model)
    nvparam = length(vparam)

    # Return information
    f = [0.0]
    cnorm = [0.0]
    snorm = [0.0]
    nlpsupn = [0.0]
    inform = Vector{Cint}([0])

    ccall(
        (:c_algencan, algencan_lib_path),                # library
        Void,                                            # Return type
        (                                                # Parameters types
            Ptr{Void},                                   # *myevalf,
            Ptr{Void},                                   # *myevalg,
            Ptr{Void},                                   # *myevalh,
            Ptr{Void},                                   # *myevalc,
            Ptr{Void},                                   # *myevaljac,
            Ptr{Void},                                   # *myevalhc,
            Ptr{Void},                                   # *myevalfc,
            Ptr{Void},                                   # *myevalgjac,
            Ptr{Void},                                   # *myevalgjacp,
            Ptr{Void},                                   # *myevalhl,
            Ptr{Void},                                   # *myevalhlp,
            Cint,                                        # jcnnzmax,
            Cint,                                        # hnnzmax,
            Ref{Cdouble},                                # *epsfeas,
            Ref{Cdouble},                                # *epsopt,
            Ref{Cdouble},                                # *efstain,
            Ref{Cdouble},                                # *eostain,
            Ref{Cdouble},                                # *efacc,
            Ref{Cdouble},                                # *eoacc,
            Cstring,                                     # *outputfnm,
            Cstring,                                     # *specfnm,
            Cint,                                        # nvparam,
            Ptr{Ptr{UInt8}},                             # **vparam,
            Cint,                                        # int n,
            Ref{Cdouble},                                # *x,
            Ref{Cdouble},                                # *l,
            Ref{Cdouble},                                # *u,
            Cint,                                        #  m,
            Ref{Cdouble},                                # double *lambda,
            Ref{UInt8},                                  # *equatn,
            Ref{UInt8},                                  # _Bool *linear,
            Ref{UInt8},                                  # _Bool *coded,
            UInt8,                                       # _Bool checkder,
            Ref{Cdouble},                                # double *f,
            Ref{Cdouble},                                # double *cnorm,
            Ref{Cdouble},                                # double *snorm,
            Ref{Cdouble},                                # double *nlpsupn,
            Ref{Cint}                                    # int *inform
        ),
        myevalf, myevalg, myevalh, myevalc, myevaljac, myevalhc, myevalfc,
        myevalgjac, myevalgjacp, myevalhl, myevalhlp, jcnnzmax, hnnzmax,
        epsfeas, epsopt, efstain, eostain, efacc, eoacc, outputfnm, specfnm,
        nvparam, vparam, model.n, model.x, model.lb, model.ub, m, mult,
        is_equality, is_g_linear, coded, checkder, f, cnorm, snorm,
        nlpsupn, inform
    )

    # Fix sign of objetive function
    model.obj_val = model.sense*f[1]

    # Deal with lower bound and two-sided contraints
    model.mult = model.g_sense .* mult[1:model.m]
    model.mult[model.g_two_sides] -= mult[model.m + 1:m]

    # Recover status information
    model.status = find_status(model, cnorm[1], snorm[1], nlpsupn[1],
        Int(inform[1]))

    model.solve_time = toc()
    return Int(inform[1])

end
export optimize!

# Local, auxiliary functions

"Transform the option dictionary into a vparam string array"
function option2vparam(model::AlgencanMathProgModel)
    parameters = [:epsfeas, :epsopt, :efstain, :eostain, :efacc, :eoacc,
        :outputfnm, :specfnm]
    vparam = Vector{String}(0)
    for option in model.options
        key, value = option
        if key in parameters
            continue
        end
        key = replace(string(key), "_", "-")
        push!(vparam, "$key $value")
    end
    return vparam
end

"Find final status of Algencan"
function find_status(model::AlgencanMathProgModel, cnorm::Float64, snorm::Float64,
    nlpsupn::Float64, inform::Int)

    if inform != 0
        return :Error
    end

    # These constants come from Algencan code
    max_multiplier, fmin = 1.0e+20, -1.0e+20
    bounded_obj = model.sense*model.obj_val > fmin

    # Optimality thresholds
    epsopt, epsfeas = model.options[:epsopt], model.options[:epsfeas]

    # Conditions for constrained problems
    if model.m > 0
        bounded_mult = maximum(abs.(model.mult)) < max_multiplier
        feasible = cnorm <= epsfeas
        if feasible && (!bounded_mult || !bounded_obj)
            return :Unbounded
        elseif feasible && nlpsupn <= epsopt && snorm <= epsopt
            return :Optimal
        elseif !feasible
            return :Infeasible
        else
            return :Error
        end
    else
        if nlpsupn <= epsopt && bounded_obj
            return :Optimal
        elseif !bounded_obj
            return :Unbounded
        else
            return :Error
        end
    end
end

# Hints to precompile
# TODO: Automate this with proper packages
precompile(loadproblem!, (AlgencanMathProgModel, Integer, Integer,
    Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Symbol,
    AbstractNLPEvaluator))
precompile(option2vparam, (AlgencanMathProgModel,))
precompile(julia_fc, (Cint, Ptr{Float64}, Ptr{Float64}, Cint, Ptr{Float64},
    Ptr{Cint}))
precompile(julia_gjac, (Cint, Ptr{Float64}, Ptr{Float64}, Cint, Ptr{Cint},
    Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Cint, Ptr{UInt8}, Ptr{Cint}))
precompile(julia_hl, (Cint, Ptr{Float64}, Cint, Ptr{Float64}, Float64,
    Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Cint,
    Ptr{UInt8}, Ptr{Cint}))
precompile(julia_hlp, (Cint, Ptr{Float64}, Cint, Ptr{Float64}, Float64,
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{UInt8}, Ptr{Cint}))
precompile(optimize!, (AlgencanMathProgModel,))
precompile(find_status, (AlgencanMathProgModel, Float64, Float64, Float64, Int))

end # module
