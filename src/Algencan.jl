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
struct AlgencanSolver <: AbstractMathProgSolver
    options
end
AlgencanSolver(;kwargs...) = AlgencanSolver(kwargs)

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

    # Decision variables and constraints values.
    x::Vector{Float64}              # Starting and final solution
    # TODO: See how to allow for initial multiplies
    g::Vector{Float64}              # Final constraint values
    mult_g::Vector{Float64}         # Final Lagrange multipliers on constraints
    obj_val::Float64                # Final objective
    status::Int                     # Final status

    # Options to be passed to the solver
    options

    function AlgencanMathProgModel(;options...)
        model = new()

        # Set up model with dummy values to mark that is not initialized
        model.n, model.m = 0, 0
        model.lb, model.ub = Float64[], Float64[]
        model.g_lb, model.g_ub = Float64[], Float64[]
        model.sense = 1.0
        model.j_row_inds, model.j_col_inds = Int[], Int[]
        model.h_row_inds, model.h_col_inds = Int[], Int[]
        model.x, model.g, model.mult_g = Float64[], Float64[], Float64[]
        model.obj_val, model.status = 0.0, 0
        # Pass on the options
        model.options = options
        model
    end
end
NonlinearModel(s::AlgencanSolver) = AlgencanMathProgModel(;s.options...)
LinearQuadraticModel(s::AlgencanSolver) = NonlinearToLPQPBridge(NonlinearModel(s))

###############################################################################
# Begin interface implementation

# generic nonlinear interface
function loadproblem!(model::AlgencanMathProgModel, numVar::Integer,
    numConstr::Integer, x_l, x_u, g_lb, g_ub, sense::Symbol,
    d::AbstractNLPEvaluator)

    # Initialize the evaluator with the right features
    features = features_available(d)
    has_hessian = (:Hess in features)
    init_feat = [:Grad]
    has_hessian && push!(init_feat, :Hess)
    numConstr > 0 && push!(init_feat, :Jac)
    initialize(d, init_feat)

    # Link the model with the problem representation
    @assert sense == :Min || sense == :Max

    # Copy data
    model.n = numVar
    model.m = numConstr
    model.lb, model.ub = float(x_l), float(x_u)
    model.g_lb, model.g_ub = float(g_lb), float(g_ub)
    model.sense = sense == :Min ? 1.0 : -1.0
    model.evaluator = d

    # Get strutural indices of Jacobian and Hessian.
    j_row_inds, j_col_inds = MathProgBase.jac_structure(d)
    h_row_inds, h_col_inds = MathProgBase.hesslag_structure(d)

    model.j_row_inds, model.j_col_inds = j_row_inds, j_col_inds
    model.h_row_inds, model.h_col_inds = h_row_inds, h_col_inds
    model.x, model.g, model.mult_g = zeros(numVar), zeros(numConstr),
        zeros(numConstr)
    model.obj_val, model.status = 0.0, 0

end

# Simple functions
getsense(model::AlgencanMathProgModel) = model.sense

numvar(model::AlgencanMathProgModel) = model.n

numconstr(model::AlgencanMathProgModel) = model.m

# TODO: Implement this
numlinconstr(model::AlgencanMathProgModel) = 0

numquadconstr(model::AlgencanMathProgModel) = 0

function status(model::AlgencanMathProgModel)
    # TODO: I need to actually verify the KKT conditions to define the final
    # status as Algencan does not inform it (it only appears in the screen)
    return :Optimal
end

getobjval(model::AlgencanMathProgModel) = model.obj_val * (model.sense == :Max ? -1 : +1)

getsolution(model::AlgencanMathProgModel) = model.x

function getreducedcosts(model::AlgencanMathProgModel)
    # TODO: Verify, I am not thinking about constraints yet.
    # sense = model.inner.sense
    # redcost = model.inner.mult_x_U - model.inner.mult_x_L
    # return sense == :Max ? redcost : -redcost
    return zeros(model.m)
end

function getconstrduals(model::AlgencanMathProgModel)
    # TODO: Verify, I am not thinking about constraints yet.

    v = model.mult_g # return multipliers for all constraints
    return model.sense == :Max ? copy(v) : -v
end

getrawsolver(model::AlgencanMathProgModel) = nothing

setwarmstart!(model::AlgencanMathProgModel, x) = (model.x = x)

"Transform the option dictionariy to a vparam string array"
function option2vparam(model::AlgencanMathProgModel)
    vparam = Vector{String}(0)
    for option in model.options
        key, value = option
        key = replace(string(key), "_", "-")
        push!(vparam, "$key $value")
    end
    return vparam
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
    # TODO: Use model.sense to define "direction" of optimization
    function julia_fc(n::Cint, x_ptr::Ptr{Float64}, obj_ptr::Ptr{Float64},
        m::Cint, g_ptr::Ptr{Float64}, flag_ptr::Ptr{Cint})
        x = unsafe_wrap(Array, x_ptr, Int(n))
        obj_val = MathProgBase.eval_f(current_algencan_model.evaluator, x)
        unsafe_store!(obj_ptr, current_algencan_model.sense*obj_val)
        g = unsafe_wrap(Array, g_ptr, Int(m))
        MathProgBase.eval_g(current_algencan_model.evaluator, g, x)
        g .-= current_algencan_model.g_ub
        unsafe_store!(flag_ptr, 0)
        nothing
    end

    const c_julia_fc = cfunction(julia_fc, Void, (Cint, Ptr{Float64},
        Ptr{Float64}, Cint, Ptr{Float64}, Ptr{Cint}))

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
            unsafe_store!(lmem_ptr, 1)
            unsafe_store!(flag_ptr, 1)
            return nothing
        else
            unsafe_store!(lmem_ptr, 0)
        end
        unsafe_store!(jnnz_ptr, nnz)
        jcol_ind = unsafe_wrap(Array, jcol_ptr, Int(lim))
        jrow_ind = unsafe_wrap(Array, jrow_ptr, Int(lim))

        # TODO: This is extra work tha only happens because I am calling C instead
        # of Fortran.
        # TODO: Try to move to direct fortran calling.
        jrow_ind[1:nnz] = current_algencan_model.j_row_inds - 1
        jcol_ind[1:nnz] = current_algencan_model.j_col_inds - 1

        # Compute the constraints Jacobian
        J = unsafe_wrap(Array, jval_ptr, Int(lim))
        MathProgBase.eval_jac_g(current_algencan_model.evaluator, J, x)

        unsafe_store!(flag_ptr, 0)
        nothing
    end

    const c_julia_gjac = cfunction(julia_gjac, Void, (Cint, Ptr{Float64},
        Ptr{Float64}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Cint,
        Ptr{UInt8}, Ptr{Cint}))

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
        unsafe_store!(hnnz_ptr, nnz)
        hcol_ind = unsafe_wrap(Array, hcol_ptr, Int(lim))
        hrow_ind = unsafe_wrap(Array, hrow_ptr, Int(lim))
        # TODO: This is extra work tha only happens because I am calling C instead
        # of Fortran.
        # TODO: Try to move to direct fortran calling.
        hrow_ind[1:nnz] = current_algencan_model.h_row_inds - 1
        hcol_ind[1:nnz] = current_algencan_model.h_col_inds - 1

        # Compute the Hessian (for now objective function only)
        σ = scale_f*current_algencan_model.sense
        μ = unsafe_wrap(Array, mult_ptr, Int(m))
        scale_g = unsafe_wrap(Array, scale_g_ptr, Int(m))
        μ .*= scale_g
        x = unsafe_wrap(Array, x_ptr, Int(n))
        H = unsafe_wrap(Array, hval_ptr, Int(lim))
        MathProgBase.eval_hesslag(current_algencan_model.evaluator, H, x, σ, μ)
        unsafe_store!(flag_ptr, 0)
        nothing
    end

    const c_julia_hl = cfunction(julia_hl, Void, (Cint, Ptr{Float64}, Cint,
        Ptr{Float64}, Float64, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, Ptr{Float64},
        Ptr{Cint}, Cint, Ptr{UInt8}, Ptr{Cint}))

    # Now CallAlgencan. I will do it slowly, first I create variables for all
    # Algencan's parameters, that are a lot, and then I call it.
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
    hnnzmax = length(model.h_row_inds)
    coded = zeros(UInt8, 11)
    coded[7] = 1
    coded[8] = 1
    coded[10] = 1

    # Parameters setting
    epsfeas = [1.0e-08]
    epsopt  = [1.0e-08]
    efstin  = sqrt.(epsfeas)
    eostin  = epsopt.^1.5
    efacc   = sqrt.(epsfeas)
    eoacc   = sqrt.(epsopt)

    outputfnm = ""
    specfnm   = ""

    vparam = option2vparam(model)
    nvparam = length(vparam)

    n = model.n
    l = model.lb
    u = model.ub

    # Information on the constraints that do not exist
    m = model.m
    lambda = zeros(Float64, m)
    equatn = zeros(UInt8, m)
    equatn[model.g_lb .== model.g_ub] = 1
    linear = zeros(UInt8, m)

    checkder = UInt8(0)
    f = [0.0]
    cnorm = [0.0]
    snorm = [0.0]
    nlpsupn = [0.0]

    inform = Vector{Cint}([0])

    b = ccall(
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
        myevalf,
        myevalg,
        myevalh,
        myevalc,
        myevaljac,
        myevalhc,
        myevalfc,
        myevalgjac,
        myevalgjacp,
        myevalhl,
        myevalhlp,
        jcnnzmax,
        hnnzmax,
        epsfeas,
        epsopt,
        efstin,
        eostin,
        efacc,
        eoacc,
        outputfnm,
        specfnm,
        nvparam,
        vparam,
        n,
        model.x,
        model.lb,
        model.ub,
        m,
        lambda,
        equatn,
        linear,
        coded,
        checkder,
        f,
        cnorm,
        snorm,
        nlpsupn,
        inform
    )

    model.obj_val = f[1]
    model.status = inform[1]

    println("Direct inform = ", inform)

    return inform[1]
end

end # module
