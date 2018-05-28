__precompile__()

"""
Algencan interface to MathProgBase and JuMP.

See its [Git hub page](https://github.com/pjssilva/Algencan.jl)
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
const algencan_lib_path = string(joinpath(ENV["ALGENCAN_LIB_DIR"],
    "libalgencan.so"))

# Global to avoid closures as Algencan does not allow to send user information
# back to call backs. Long names to avoid conflicts.
global current_algencan_problem

# Standard LP interface
import MathProgBase
importall MathProgBase.SolverInterface

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
            :efstin=>sqrt.(1.0e-8),
            :eostin=>(1.0e-8)^1.5,
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

###############################################################################
# Model objects objects

export getconstrduals

"Algencan model, that storing solution data"
mutable struct AlgencanMathProgModel <: AbstractNonlinearModel
    # Problem data
    n::Int                          # Num variables
    m::Int                          # Num constraints
    lb::Vector{Float64}             # Lower bounds on vars
    ub::Vector{Float64}             # Upper bounds on var
    g_ub::Vector{Float64}           # Upper bounds on constraints
    g_lb::Vector{Float64}           # Lower bound on constrains
    sense::Float64                  # 1.0 for :Min and -1 for :Max
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

    # Options to be passed to the solver
    options

    function AlgencanMathProgModel(options)
        model = new()

        # Set up model with dummy values to mark that is not initialized
        model.n, model.m = 0, 0
        model.lb, model.ub = Float64[], Float64[]
        model.g_lb, model.g_ub = Float64[], Float64[]
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

"Loads the problem with its basica data and functions in a NLPEvaluator"
function loadproblem!(model::AlgencanMathProgModel, numVar::Integer,
    numConstr::Integer, x_l::Vector{Float64}, x_u::Vector{Float64},
    g_lb::Vector{Float64}, g_ub::Vector{Float64}, sense::Symbol,
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
    model.g_lb, model.g_ub = float(g_lb), float(g_ub)
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
    model.g, model.mult = zeros(numConstr), zeros(numConstr)
    model.obj_val, model.status = 0.0, :Undefined
end

# Simple functions
getsense(model::AlgencanMathProgModel) = (model.sense == 1 ? :Min : :Max)

numvar(model::AlgencanMathProgModel) = model.n

numconstr(model::AlgencanMathProgModel) = model.m

status(model::AlgencanMathProgModel) = model.status

getobjval(model::AlgencanMathProgModel) = model.obj_val * model.sense

getsolution(model::AlgencanMathProgModel) = model.x

function getreducedcosts(model::AlgencanMathProgModel)
    # TODO: Verify, I am not thinking about constraints yet.
    # sense = model.inner.sense
    # redcost = model.inner.mult_x_U - model.inner.mult_x_L
    # return sense == :Max ? redcost : -redcost
    return zeros(model.n)
end

function getconstrduals(model::AlgencanMathProgModel)
    v = model.mult
    scale!(v, model.sense)
    return v
end

getrawsolver(model::AlgencanMathProgModel) = nothing

setwarmstart!(model::AlgencanMathProgModel, x) = (model.x = x)

"Transform the option dictionary to a vparam string array"
function option2vparam(model::AlgencanMathProgModel)
    parameters = [:epsfeas, :epsopt, :efstin, :eostin, :efacc, :eoacc,
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

###########################################################################
# Algencan callbacks
###########################################################################

"Compute objective and constraints as required by Algencan"
function julia_fc(n::Cint, x_ptr::Ptr{Float64}, obj_ptr::Ptr{Float64},
    m::Cint, g_ptr::Ptr{Float64}, flag_ptr::Ptr{Cint})
    x = unsafe_wrap(Array, x_ptr, Int(n))
    obj_val = MathProgBase.eval_f(current_algencan_model.evaluator, x)
    unsafe_store!(obj_ptr, current_algencan_model.sense*obj_val)
    g = unsafe_wrap(Array, g_ptr, Int(m))
    MathProgBase.eval_g(current_algencan_model.evaluator, g, x)
    g .-= current_algencan_model.g_ub
    unsafe_store!(flag_ptr, Cint(0))
    nothing
end

"Compute objective gradient and constraints Jacobian as required by Algencan"
function julia_gjac(n::Cint, x_ptr::Ptr{Float64}, f_grad_ptr::Ptr{Float64},
    m::Cint, jrow_ptr::Ptr{Cint}, jcol_ptr::Ptr{Cint},
    jval_ptr::Ptr{Float64}, jnnz_ptr::Ptr{Cint}, lim::Cint,
    lmem_ptr::Ptr{UInt8}, flag_ptr::Ptr{Cint})

    # Compute gradient of the objective
    x = unsafe_wrap(Array, x_ptr, Int(n))
    f_grad = unsafe_wrap(Array, f_grad_ptr, Int(n))
    MathProgBase.eval_grad_f(current_algencan_model.evaluator, f_grad, x)
    scale!(f_grad, current_algencan_model.sense)

    # Find structure of the constraints Jacobian
    nnz = length(current_algencan_model.j_row_inds)
    if nnz > Int(lim)
        unsafe_store!(lmem_ptr, Cint(1))
        unsafe_store!(flag_ptr, Cint(1))
        return nothing
    else
        unsafe_store!(lmem_ptr, Cint(0))
    end
    unsafe_store!(jnnz_ptr, nnz)
    jcol_ind = unsafe_wrap(Array, jcol_ptr, Int(lim))
    jrow_ind = unsafe_wrap(Array, jrow_ptr, Int(lim))
    jrow_ind[1:nnz] = current_algencan_model.j_row_inds
    jcol_ind[1:nnz] = current_algencan_model.j_col_inds

    # Compute the constraints Jacobian
    J = unsafe_wrap(Array, jval_ptr, Int(lim))
    MathProgBase.eval_jac_g(current_algencan_model.evaluator, J, x)

    unsafe_store!(flag_ptr, Cint(0))
    nothing
end

"Compute the Hessian of the Lagrangian as required by Algencan"
function julia_hl(n::Cint, x_ptr::Ptr{Float64}, m::Cint,
    mult_ptr::Ptr{Float64}, scale_f::Float64, scale_g_ptr::Ptr{Float64},
    hrow_ptr::Ptr{Cint}, hcol_ptr::Ptr{Cint},
    hval_ptr::Ptr{Float64}, hnnz_ptr::Ptr{Cint}, lim::Cint,
    lmem_ptr::Ptr{UInt8}, flag_ptr::Ptr{Cint})
    # Get nonzero indexes.
    nnz = length(current_algencan_model.h_row_inds)
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
    hrow_ind[1:nnz] = current_algencan_model.h_row_inds
    hcol_ind[1:nnz] = current_algencan_model.h_col_inds

    # Compute the Hessian (for now objective function only)
    σ = scale_f*current_algencan_model.sense
    μ = unsafe_wrap(Array, mult_ptr, Int(m))
    scale_g = unsafe_wrap(Array, scale_g_ptr, Int(m))
    μ .*= scale_g
    x = unsafe_wrap(Array, x_ptr, Int(n))
    H = unsafe_wrap(Array, hval_ptr, Int(lim))
    MathProgBase.eval_hesslag(current_algencan_model.evaluator, H, x, σ, μ)
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
    myevalhlp = C_NULL
    jcnnzmax = length(model.j_row_inds)
    hnnzmax = length(model.h_row_inds) + 10*jcnnzmax
    coded = zeros(UInt8, 11)
    coded[7] = UInt8(1)
    coded[8] = UInt8(1)
    coded[10] = UInt8(1)
    checkder = UInt8(0)

    # Parameters controling precision
    epsfeas = [model.options[:epsfeas]]
    epsopt = [model.options[:epsopt]]
    efstin = [model.options[:efstin]]
    eostin = [model.options[:eostin]]
    efacc  = [model.options[:efacc]]
    eoacc  = [model.options[:eoacc]]

    # Extra parameters
    outputfnm = model.options[:outputfnm]
    specfnm   = model.options[:specfnm]
    vparam = Vector{String}(0) # option2vparam(model)
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
            Ref{Cdouble},                                # *efstin,
            Ref{Cdouble},                                # *eostin,
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
        epsfeas, epsopt, efstin, eostin, efacc, eoacc, outputfnm, specfnm,
        nvparam, vparam, model.n, model.x, model.lb, model.ub, model.m,
        model.mult, model.is_equality, model.is_g_linear,
        coded, checkder, f, cnorm, snorm, nlpsupn, inform
    )

    model.obj_val = model.sense*f[1]
    model.status = find_status(model, cnorm[1], snorm[1], nlpsupn[1],
        Int(inform[1]))
    return Int(inform[1])

end

function find_status(model::AlgencanMathProgModel, cnorm::Float64, snorm::Float64,
    nlpsupn::Float64, inform::Int)

    if inform != 0
        return :Error
    end

    # Constant comes from Algencan code
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

end # module
