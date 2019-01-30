"""
Algencan interface to MathProgBase and JuMP.

See its [GitHub page](https://github.com/pjssilva/Algencan.jl)
"""
module Algencan

using LinearAlgebra
import Libdl

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
    const algencan_lib_path = libalgencan
elseif "ALGENCAN_LIB_DIR" in keys(ENV)
    const algencan_lib_path = string(joinpath(ENV["ALGENCAN_LIB_DIR"],
        "libalgencan.so"))
else
    error("Algencan not properly installed. Please run Pkg.build(\"Algenacn\")")
    error("or set the LAGENCAL_LIB_DIR enviroment variable.")
end

# Compiles to the *static* path of the algencan library
if "ALGENCAN_LIB_DIR" in keys(ENV)
    const algencan_lib_path = string(joinpath(ENV["ALGENCAN_LIB_DIR"],
        "libalgencan.so"))
elseif isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
    const algencan_lib_path = libalgencan
end

# Standard LP interface
import MathProgBase
const MPB = MathProgBase
import Base.copy

###############################################################################
# Solver objects
export AlgencanSolver

"Algencan solver that stores options"
struct AlgencanSolver <: MPB.AbstractMathProgSolver
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
mutable struct AlgencanMathProgModel <: MPB.AbstractNonlinearModel
    # Problem data
    n::Int                          # Num variables
    lb::Vector{Float64}             # Lower bounds on vars
    ub::Vector{Float64}             # Upper bounds on var
    sense::Float64                  # 1.0 for :Min and -1 for :Max
    m::Int                          # Num constraints
    g_ub::Vector{Float64}           # Upper bounds on constraints
    g_lb::Vector{Float64}           # Lower bound on constrains
    g_sense::Vector{Int}
    g_two_sinvmap::Vector{Int}
    g_two_smap::Vector{Int}
    g_has_lb::Bool                  # true if at least one constraint has lower
                                    # bound
    evaluator::MPB.AbstractNLPEvaluator # Evaluator for functions
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

    # Performance information
    solve_time::Float64             # Total solution time
    n_fc::Int                       # Number of funciton evaluations
    n_gjac::Int                     # Number of gradients and Jacobian evaluations
    n_hl::Int                       # Number of Hessian evaluations
    n_hlp::Int                      # Number of Hessian times vector

    # Options to be passed to the solver
    options

    function AlgencanMathProgModel(options)
        model = new()

        # Set up model with dummy values to mark that is not initialized
        model.n, model.m = 0, 0
        model.lb, model.ub = Float64[], Float64[]
        model.g_lb, model.g_ub = Float64[], Float64[]
        model.g_sense, model.g_two_sinvmap, model.g_two_smap = Float64[], Int[], Int[]
        model.g_has_lb = false
        model.sense = 1.0
        model.j_row_inds, model.j_col_inds = Int[], Int[]
        model.h_row_inds, model.h_col_inds = Int[], Int[]
        model.x, model.g, model.mult = Float64[], Float64[], Float64[]
        model.obj_val, model.status = 0.0, :Undefined
        model.solve_time = 0.0
        model.n_fc, model.n_gjac, model.n_hl, model.n_hlp = 0, 0, 0, 0
        model.options = options
        model
    end
end
MPB.NonlinearModel(s::AlgencanSolver) = AlgencanMathProgModel(s.options)
MPB.LinearQuadraticModel(s::AlgencanSolver) = MPB.NonlinearToLPQPBridge(MPB.NonlinearModel(s))

###############################################################################
# Begin interface implementation
# TODO: Verify if really need not to explicitly export the MPB functions

# Simple functions
"Return the current primal point, the solution after calling optimize"
MPB.getsolution(model::AlgencanMathProgModel) = model.x

"Get objective value at current primal point"
MPB.getobjval(model::AlgencanMathProgModel) = model.sense*model.obj_val

"Get current model status"
MPB.status(model::AlgencanMathProgModel) = model.status

"Get best bound on (local) optimal value"
MPB.getobjbound(m::AlgencanMathProgModel) = model.status == :optimal ? getobjval(m) : m.sense*Inf

"Get gap to (local) optimality"
MPB.getobjgap(m::AlgencanMathProgModel) = model.status == :Optimal ? 0.0 : Inf

"There is no inner solver, all functionality is exposed by the default interface"
MPB.getrawsolver(m::AlgencanMathProgModel) = nothing

"Get the solution time"
MPB.getsolvetime(m::AlgencanMathProgModel) = m.solve_time

"Change optimization sense, either :Min or :Max"
function MPB.setsense!(m::AlgencanMathProgModel, sense)
    m.sense = sense == :Max ? -1 : 1
end

"Return problem sense"
MPB.getsense(model::AlgencanMathProgModel) = (model.sense == 1 ? :Min : :Max)

"Return the number of decision variables"
MPB.numvar(model::AlgencanMathProgModel) = model.n

"Return the number of constraints"
MPB.numconstr(model::AlgencanMathProgModel) = model.m

"Return a copy of the current model"
MPB.copy(m::AlgencanMathProgModel) = deepcopy(m)

"Algencan only deals with continuous variable, so inform if this is not the case"
function MPB.setvartype!(m::AlgencanMathProgModel, v::Vector{Symbol})
    if !all(v .== :Cont)
        throw("Algencan ony deals with continuous variables")
    end
end

"Return variable types: they are all continuous"
function MPB.getvartype(model::AlgencanMathProgModel)
    types = Vector{Symbol}(model.m)
    types .= :Cont
    return types
end

"""
Set parameter for a solver, that will be default for all next models

You can use any Algencan parameter that can be set in a specification file plus
:epsfeas, :epsopt, :efstain, :eostain, efacc, :eoac, :outputfnm, :specfnm.

See more details in Algencan documentation.
"""
function MPB.setparameters!(m::Union{AlgencanSolver, AlgencanMathProgModel};
    kwargs...)

    for (key, value) in kwargs
        m.options[key] = value
    end
end

"Set an initial value for the decision variables (warmstart)"
function MPB.setwarmstart!(model::AlgencanMathProgModel, x)
    model.x = copy(x)
end

function MPB.getconstrduals(model::AlgencanMathProgModel)
    v = model.mult
    v .*= model.sense
    return v
end

"Return mutipliers associated with bound cosntraints"
function MPB.getreducedcosts(model::AlgencanMathProgModel)
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
    # MPB.eval_grad_f(model.evaluator, grad_lag, model.x)
    # grad_lag *= model.sense
    #
    # # Add the multiplier combination of the contraints gradients
    # grad_g = Vector{Float64}(model.n)
    # MPB.eval_jac_prod_t(model.evaluator, grad_g, model.x, model.mult)
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

# Simple function only defined for Algencan models

"Set an inital value for the constaint mutipliers (warmstart)"
function setmultwarmstart!(model::AlgencanMathProgModel, mult)
    model.mult = copy(mult)
end
export setmultwarmstart!

function getnfevals(model::AlgencanMathProgModel)
    return model.n_fc, model.n_gjac, model.n_hl, model.n_hlp
end
export getnfevals

function resetnfevals(model::AlgencanMathProgModel)
   model.n_fc, model.n_gjac, model.n_hl, model.n_hlp = 0, 0, 0, 0
   nothing
end
export resetnfevals

# More complex funcitons from MPB interface

"Loads the problem with its basic data and functions in a NLPEvaluator"
function MPB.loadproblem!(model::AlgencanMathProgModel, numVar::Integer,
    numConstr::Integer, x_l, x_u, g_lb, g_ub, sense::Symbol,
    d::MPB.AbstractNLPEvaluator)

    @assert sense == :Min || sense == :Max

    # Link the model with the problem representation
    # Initialize the evaluator with the right features
    features = MPB.features_available(d)
    has_hessian = (:Hess in features)
    init_feat = [:Grad]
    has_hessian && push!(init_feat, :Hess)
    numConstr > 0 && push!(init_feat, :Jac)
    MPB.initialize(d, init_feat)

    # Copy data
    model.n = numVar
    model.m = numConstr
    model.lb, model.ub = float(x_l), float(x_u)
    g_lb, g_ub = float(g_lb), float(g_ub)
    model.g_sense, model.g_two_sinvmap, model.g_two_smap = treat_lower_bounds(
        g_lb, g_ub)
    model.g_has_lb = (model.m > 0 && (minimum(model.g_sense) == -1 ||
        length(model.g_two_smap) > 0))
    model.g_lb, model.g_ub = g_lb, g_ub

    # Contraints with only lower bound will be multiplied by -1.0, hence
    # their lower bounds becomes minus the upper bound
    g_only_low = (model.g_sense .== -1.0)
    model.g_ub[g_only_low] .= -g_lb[g_only_low]
    model.g_lb[g_only_low] .= -Inf

    model.sense = (sense == :Min ? 1.0 : -1.0)
    model.evaluator = d

    # Constraints types
    model.is_equality = zeros(UInt8, numConstr)
    model.is_equality[model.g_lb .== model.g_ub] .= 1
    model.is_g_linear = zeros(UInt8, numConstr)
    for i in 1:numConstr
        if MPB.isconstrlinear(d, i)
            model.is_g_linear[i] = 1
        end
    end

    # Get strutural indices of Jacobian and Hessian.
    j_row_inds, j_col_inds = MPB.jac_structure(d)
    h_row_inds, h_col_inds = MPB.hesslag_structure(d)

    # C indices start in 0
    model.j_row_inds, model.j_col_inds = j_row_inds .- 1, j_col_inds .- 1
    model.h_row_inds, model.h_col_inds = h_row_inds .- 1, h_col_inds .- 1

    # Initial values
    model.x = zeros(numVar)
    model.g = zeros(numConstr)
    model.mult = zeros(numConstr)
    model.obj_val, model.status = 0.0, :Undefined
    model.solve_time = 0.0
    model.n_fc, model.n_gjac, model.n_hl, model.n_hlp = 0, 0, 0, 0
end

"Analyse the lower and upper bounds on the constraints and prepare the
 data structure to treat lower bounds."
function treat_lower_bounds(lb, ub)
    m = length(lb)

    # Allow to treat constraints that have only lower bounds
    sense = ones(m)
    only_lower = (-Inf .< lb) .& (ub .== Inf)
    sense[only_lower] .= -1.0

    # Treat two side constraints
    two_sides = -Inf .< lb .!= ub .< Inf
    two_smap = (1:m)[two_sides]
    new_ind = 1
    two_sinvmap = zeros(m)
    for i = 1:m
        if two_sides[i]
            two_sinvmap[i] = new_ind
            new_ind += 1
        end
    end

    return sense, two_sinvmap, two_smap
end

###########################################################################
# Algencan callbacks
###########################################################################

"Compute objective and constraints as required by Algencan"
function julia_fc(model::AlgencanMathProgModel, n::Cint, x_ptr::Ptr{Float64}, obj_ptr::Ptr{Float64},
    m::Cint, g_ptr::Ptr{Float64}, flag_ptr::Ptr{Cint})

    # Evaluate objective and constraints
    model.n_fc += 1
    x = unsafe_wrap(Array, x_ptr, Int(n))
    obj_val = MPB.eval_f(model.evaluator, x)
    unsafe_store!(obj_ptr, model.sense*obj_val)
    g = unsafe_wrap(Array, g_ptr, Int(m))
    MPB.eval_g(model.evaluator, g, x)

    # Treat lower bounds and two-sided constraints
    if model.g_has_lb
        first_g = view(g, 1:model.m)
        first_g .*= model.g_sense
        for i = model.m + 1:m
            mapped_i = model.g_two_smap[i - model.m]
            g[i] = -first_g[mapped_i] + model.g_lb[mapped_i]
        end
        first_g .-= model.g_ub
    else
        g .-= model.g_ub
    end

    # Report that evaluation of successful
    unsafe_store!(flag_ptr, Cint(0))
    nothing
end

"Compute objective gradient and constraints Jacobian as required by Algencan"
function julia_gjac(model::AlgencanMathProgModel, n::Cint, x_ptr::Ptr{Float64}, f_grad_ptr::Ptr{Float64},
    m::Cint, jrow_ptr::Ptr{Cint}, jcol_ptr::Ptr{Cint},
    jval_ptr::Ptr{Float64}, jnnz_ptr::Ptr{Cint}, lim::Cint,
    lmem_ptr::Ptr{UInt8}, flag_ptr::Ptr{Cint})

    # Compute gradient of the objective
    model.n_gjac += 1
    x = unsafe_wrap(Array, x_ptr, Int(n))
    f_grad = unsafe_wrap(Array, f_grad_ptr, Int(n))
    MPB.eval_grad_f(model.evaluator, f_grad, x)
    f_grad .*= model.sense

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
    MPB.eval_jac_g(model.evaluator, J, x)

    # Treat the presence of lower bound in the constraints
    if model.g_has_lb
        @inbounds for i = 1:length(model.j_row_inds)
            # +1, -1 to translate from C indexing to Julia indexing
            rind, cind = model.j_row_inds[i] + 1, model.j_col_inds[i]
            if model.g_two_sinvmap[rind] > 0
                nnz += 1
                jrow_ind[nnz] = model.m + model.g_two_sinvmap[rind] - 1
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
function julia_hl(model::AlgencanMathProgModel, n::Cint, x_ptr::Ptr{Float64}, m::Cint,
    mult_ptr::Ptr{Float64}, scale_f::Float64, scale_g_ptr::Ptr{Float64},
    hrow_ptr::Ptr{Cint}, hcol_ptr::Ptr{Cint},
    hval_ptr::Ptr{Float64}, hnnz_ptr::Ptr{Cint}, lim::Cint,
    lmem_ptr::Ptr{UInt8}, flag_ptr::Ptr{Cint})

    model.n_hl += 1
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
    σ = model.sense*scale_f
    alg_mult = unsafe_wrap(Array, mult_ptr, Int(m))
    scale_g = unsafe_wrap(Array, scale_g_ptr, Int(m))
    if !model.g_has_lb
        μ = scale_g .* alg_mult
    else
        μ = model.g_sense .* scale_g[1:model.m] .* alg_mult[1:model.m]
        for i = 1:length(model.g_two_smap)
            μ[model.g_two_smap[i]] -= scale_g[model.m + i] * alg_mult[model.m + i]
        end
    end

    # Evaluate the Hessian
    x = unsafe_wrap(Array, x_ptr, Int(n))
    H = unsafe_wrap(Array, hval_ptr, Int(lim))
    MPB.eval_hesslag(model.evaluator, H, x, σ, μ)

    # Declare success
    unsafe_store!(flag_ptr, Cint(0))
    nothing
end

"Compute the Hessian of the Lagrangian times p as required by Algencan"
function julia_hlp(model::AlgencanMathProgModel, n::Cint, x_ptr::Ptr{Float64}, m::Cint,
    mult_ptr::Ptr{Float64}, scale_f::Float64, scale_g_ptr::Ptr{Float64},
    p_ptr::Ptr{Float64}, hp_ptr::Ptr{Float64}, goth_ptr::Ptr{UInt8},
    flag_ptr::Ptr{Cint})

    model.n_hlp += 1
    # Compute scaled multipliers
    σ = scale_f*model.sense
    alg_mult = unsafe_wrap(Array, mult_ptr, Int(m))
    scale_g = unsafe_wrap(Array, scale_g_ptr, Int(m))
    if !model.g_has_lb
        μ = model.g_sense .* alg_mult .* scale_g
    else
        μ = model.g_sense .* alg_mult[1:model.m] .* scale_g[1:model.m]
        for i = 1:length(model.g_two_smap)
            μ[model.g_two_smap[i]] -= scale_g[model.m + i] * alg_mult[model.m + i]
        end
    end

    # Evaluate Hessian times p
    x = unsafe_wrap(Array, x_ptr, Int(n))
    p = unsafe_wrap(Array, p_ptr, Int(n))
    hp = unsafe_wrap(Array, hp_ptr, Int(n))
    MPB.eval_hesslag_prod(model.evaluator, hp, x, p, σ, μ)

    # Declare success
    unsafe_store!(flag_ptr, Cint(0))
    nothing
end

function MPB.optimize!(model::AlgencanMathProgModel)
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

    start_time = time_ns()
    ###########################################################################
    # Algencan callback function wrappers
    ###########################################################################
    local_julia_fc = (n, x_ptr, obj_ptr, m, g_ptr, flag_ptr) -> julia_fc(model, n, x_ptr, obj_ptr, m, g_ptr, flag_ptr)
    c_julia_fc = @cfunction($local_julia_fc,
        Nothing, (Cint, Ptr{Float64}, Ptr{Float64}, Cint, Ptr{Float64}, Ptr{Cint}))

    local_julia_gjac = (n, x_ptr, f_grad_ptr, m, jrow_ptr, jcol_ptr, jval_ptr, jnnz_ptr, lim, lmem_ptr, flag_ptr) -> julia_gjac(model, n, x_ptr, f_grad_ptr, m, jrow_ptr, jcol_ptr, jval_ptr, jnnz_ptr, lim, lmem_ptr, flag_ptr)
    c_julia_gjac = @cfunction($local_julia_gjac,
        Nothing, (Cint, Ptr{Float64}, Ptr{Float64}, Cint, Ptr{Cint}, Ptr{Cint},
        Ptr{Float64}, Ptr{Cint}, Cint, Ptr{UInt8}, Ptr{Cint}))

    local_julia_hl = (n, x_ptr, m, mult_ptr, scale_f, scale_g_ptr, hrow_ptr, hcol_ptr, hval_ptr, hnnz_ptr, lim, lmem_ptr, flag_ptr) -> julia_hl(model, n, x_ptr, m, mult_ptr, scale_f, scale_g_ptr, hrow_ptr, hcol_ptr, hval_ptr, hnnz_ptr, lim, lmem_ptr, flag_ptr)
    c_julia_hl = @cfunction($local_julia_hl, Nothing, (Cint, Ptr{Float64}, Cint,
        Ptr{Float64}, Float64, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, Ptr{Float64},
        Ptr{Cint}, Cint, Ptr{UInt8}, Ptr{Cint}))

    local_julia_hlp = (n, x_ptr, m, mult_ptr, scale_f, scale_g_ptr, p_ptr, hp_ptr, goth_ptr, flag_ptr) -> julia_hlp(model, n, x_ptr, m, mult_ptr, scale_f, scale_g_ptr, p_ptr, hp_ptr, goth_ptr, flag_ptr)
    c_julia_hlp = @cfunction($local_julia_hlp, Nothing, (Cint, Ptr{Float64}, Cint,
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
    m = model.m + length(model.g_two_smap)
    mult = zeros(m)
    is_equality = zeros(UInt8, m)
    is_equality[1:model.m] .= model.is_equality
    is_g_linear = zeros(UInt8, m)
    is_g_linear[1:model.m] .= model.is_g_linear
    for i = 1:length(model.g_two_smap)
        is_g_linear[model.m + i] = model.is_g_linear[model.g_two_smap[i]]
    end

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

    @assert !(algencan_lib_path in Libdl.dllist())
    algencandl = Libdl.dlopen(algencan_lib_path)
    @assert algencan_lib_path in Libdl.dllist()
    algencansym = Libdl.dlsym(algencandl, :c_algencan)
    ccall(
        algencansym,                                     # function
        Nothing,                                         # Return type
        (                                                # Parameters types
            Ptr{Nothing},                                # *myevalf,
            Ptr{Nothing},                                # *myevalg,
            Ptr{Nothing},                                # *myevalh,
            Ptr{Nothing},                                # *myevalc,
            Ptr{Nothing},                                # *myevaljac,
            Ptr{Nothing},                                # *myevalhc,
            Ptr{Nothing},                                # *myevalfc,
            Ptr{Nothing},                                # *myevalgjac,
            Ptr{Nothing},                                # *myevalgjacp,
            Ptr{Nothing},                                # *myevalhl,
            Ptr{Nothing},                                # *myevalhlp,
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
    Libdl.dlclose(algencandl)
    @assert !(algencan_lib_path in Libdl.dllist())

    # Fix sign of objetive function
    model.obj_val = model.sense*f[1]

    # Deal with lower bound and two-sided contraints
    model.mult = model.g_sense .* mult[1:model.m]
    for i = 1:length(model.g_two_smap)
        model.mult[model.g_two_smap[i]] -= mult[model.m + i]
    end

    # Recover status information
    model.status = find_status(model, cnorm[1], snorm[1], nlpsupn[1],
        Int(inform[1]))

    model.solve_time = (time_ns() - start_time) / 1.0e+9
    return Int(inform[1])
end

# Local, auxiliary functions

"Transform the option dictionary into a vparam string array"
function option2vparam(model::AlgencanMathProgModel)
    parameters = [:epsfeas, :epsopt, :efstain, :eostain, :efacc, :eoacc,
        :outputfnm, :specfnm]
    vparam = Array{String}(undef, 0)
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
    # TODO: Deal with the possibility that specfnm may overrid these constans
    epsopt, epsfeas = model.options[:epsopt], model.options[:epsfeas]

    # Conditions for constrained problems
    if model.m > 0
        bounded_mult = maximum(abs.(model.mult)) < max_multiplier
        feasible = cnorm <= epsfeas
        if !feasible
            return :Infeasible
        elseif !(bounded_obj && bounded_mult)
            return :Unbounded
        elseif nlpsupn <= epsopt && snorm <= epsopt
            return :Optimal
        else
            return :Error
        end
    else
        if !bounded_obj
            return :Unbounded
        elseif nlpsupn <= epsopt
            return :Optimal
        else
            return :Error
        end
    end
end

# Hints to precompile
# TODO: Automate this with proper packages
# precompile(MPB.loadproblem!, (AlgencanMathProgModel, Integer, Integer,
#     Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Symbol,
#     MPB.AbstractNLPEvaluator))
# precompile(option2vparam, (AlgencanMathProgModel,))
# precompile(julia_fc, (Cint, Ptr{Float64}, Ptr{Float64}, Cint, Ptr{Float64},
#     Ptr{Cint}))
# precompile(julia_gjac, (Cint, Ptr{Float64}, Ptr{Float64}, Cint, Ptr{Cint},
#     Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Cint, Ptr{UInt8}, Ptr{Cint}))
# precompile(julia_hl, (Cint, Ptr{Float64}, Cint, Ptr{Float64}, Float64,
#     Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Cint,
#     Ptr{UInt8}, Ptr{Cint}))
# precompile(julia_hlp, (Cint, Ptr{Float64}, Cint, Ptr{Float64}, Float64,
#     Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{UInt8}, Ptr{Cint}))
# precompile(MPB.optimize!, (AlgencanMathProgModel,))
# precompile(find_status, (AlgencanMathProgModel, Float64, Float64, Float64, Int))
# include("precompile_Algencan.jl")
# _precompile_()
end # module
