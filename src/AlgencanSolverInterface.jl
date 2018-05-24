__precompile__()

"""
Algencan interface to MathProgBase and JuMP.

WARNING: This module is not threadsafe, you can not solve two models at
the same time. BAD THINGS WILL HAPPEN IF YOU TRY TO DO IT.

This modules creates an interface to allow you to use
[Algencan](https://www.ime.usp.br/~egbirgin/tango/codes.php) to solve nonlinear
optimization methods defined in JuMP.

Algencan is a high performance and large scale Augmented Lagrangian solver
written by Ernesto Birgin and Mario Martínez. It has many special features like
being able to use HSL librarry functions to speed up linear algebra with sparse
matrices and some smart accelaration strategies. To use it, at least for now,
you need to compile Algencan yourself. Go to Tango's project
[website](https://www.ime.usp.br/~egbirgin/tango/codes.php) grab Algencan and
compile it. However, you need to make the small changes below to be able to get
a dynamic library from it as Algencan creates a static library by default.

1. Add the option `-f PIC` to  `CFLAGS` in the top of the main Makefile.

2. After compiling, you'll get a libalgencan.a file in the `lib` subdir. In
order to create a dynamic library try the following:

```
gcc -shared -o libalgencan.so -Wl,--whole-archive libalgencan.a \\
    -Wl,--no-whole-archive -lgfortran
```

It should create a file named `libalgencan.so` in the `lib` dir.

3. Create a environmental library named `ALGENCAN_LIB_DIR` pointing to the
`lib` dir path.
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


# Globals to avoid closures as Algencan does not allow to send user information
# back to call backs. Long names to avoid conflicts.
global current_algencan_problem
const algencan_lib_path = string(joinpath(ENV["ALGENCAN_LIB_DIR"],
    "libalgencan.so"))

# Standard LP interface
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
    sense::Symbol                   # max or min
    evaluator::AbstractNLPEvaluator # Evaluator for functions
    g_row_inds::Vector{Int}         # NNZ row indexes of sparse const. Jacobian.
    g_col_inds::Vector{Int}         # NNZ col indexes of sparse const. Jacobian.
    h_row_inds::Vector{Int}         # NNZ row indexes of sparse Lag. Hessian.
    h_col_inds::Vector{Int}         # NNZ col indexes of sparse Lag. Hessian.

    # Decision variables and constraints values.
    x::Vector{Float64}              # Starting and final solution
    # TODO: No constraints for now and see how to allow for initial multiplies
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
        model.sense = :Min
        model.g_row_inds, model.g_col_inds = Int[], Int[]
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
function loadproblem!(m::AlgencanMathProgModel, numVar::Integer,
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
    model.g_lb, mode.g_ub = float(g_lb), float(g_ub)
    model.sense = sense
    model.evaluator = d

    # Get strutural indices of Jacobian and Hessian.
    g_row_inds, g_col_inds = MathProgBase.jac_structure(evaluator)
    h_row_inds, h_col_inds = MathProgBase.hesslag_structure(evaluator)

    model.g_row_inds, model.g_col_inds = g_row_inds, g_col_inds
    model.h_row_inds, model.h_col_inds = h_row_inds, h_col_inds
    model.x, model.g, model.mult_g = zeros(numVar), zeros(numConstr),
        zeros(numConstr)
    model.obj_val, model.status = 0.0, 0

end

# Simple functions
getsense(m::AlgencanMathProgModel) = m.sense

numvar(m::AlgencanMathProgModel) = m.n

numconstr(m::AlgencanMathProgModel) = m.m

# TODO: Implement this
numlinconstr(m::AlgencanMathProgModel) = 0

numquadconstr(m::AlgencanMathProgModel) = 0

function status(m::AlgencanMathProgModel)
    # TODO: I still need to map status, return that everyhing was OK
    return :Optimal

    # # Map all the possible return codes, as enumerated in
    # # Ipopt.ApplicationReturnStatus, to the MPB statuses:
    # # :Optimal, :Infeasible, :Unbounded, :UserLimit, and :Error
    # stat_sym = ApplicationReturnStatus[m.inner.status]
    # if  stat_sym == :Solve_Succeeded ||
    #     stat_sym == :Solved_To_Acceptable_Level
    #     return :Optimal
    # elseif stat_sym == :Infeasible_Problem_Detected
    #     return :Infeasible
    # elseif stat_sym == :Diverging_Iterates
    #     return :Unbounded
    #     # Things that are more likely to be fixable by changing
    #     # a parameter will be treated as UserLimit, although
    #     # some are error-like too.
    # elseif stat_sym == :User_Requested_Stop ||
    #     stat_sym == :Maximum_Iterations_Exceeded ||
    #     stat_sym == :Maximum_CpuTime_Exceeded
    #     return :UserLimit
    # else
    #     # Default is to not mislead user that it worked
    #     # Includes:
    #     #   :Search_Direction_Becomes_Too_Small
    #     #   :Feasible_Point_Found
    #     #   :Restoration_Failed
    #     #   :Error_In_Step_Computation
    #     #   :Not_Enough_Degrees_Of_Freedom
    #     #   :Invalid_Problem_Definition
    #     #   :Invalid_Option
    #     #   :Invalid_Number_Detected
    #     #   :Unrecoverable_Exception
    #     #   :NonIpopt_Exception_Thrown
    #     #   :Insufficient_Memory
    #     #   :Internal_Error
    #     warn("Ipopt finished with status $stat_sym")
    #     return :Error
    # end
    #
end

getobjval(m::AlgencanMathProgModel) = m.obj_val * (m.inner.sense == :Max ? -1 : +1)

getsolution(m::AlgencanMathProgModel) = m.x

function getreducedcosts(m::AlgencanMathProgModel)
    # TODO: Verify, I am not thinking about constraints yet.
    # sense = m.inner.sense
    # redcost = m.inner.mult_x_U - m.inner.mult_x_L
    # return sense == :Max ? redcost : -redcost
    return zeros(m.inner.m)
end

function getconstrduals(m::AlgencanMathProgModel)
    # TODO: Verify, I am not thinking about constraints yet.

    v = m.mult_g # return multipliers for all constraints
    return m.inner.sense == :Max ? copy(v) : -v
end

getrawsolver(m::AlgencanMathProgModel) = nothing

setwarmstart!(m::AlgencanMathProgModel, x) = (m.x = x)

function optimize!(m::AlgencanMathProgModel)
    # TODO: Allow warm start primal and specially dual
    #copy!(m.inner.x, m.warmstart) # set warmstart
    # TODO: No options for now
    # for (name,value) in m.options
    #     sname = string(name)
    #     if match(r"(^resto_)", sname) != nothing
    #         sname = replace(sname, r"(^resto_)", "resto.")
    #     end
    #     addOption(m.inner, sname, value)
    # end
    global current_algencan_model = m

    ###########################################################################
    # Algencan callback function wrappers
    ###########################################################################
    # TODO: Use model.sense to define "direction" of optimization
    function my_f(n::Cint, x_ptr::Ptr{Float64}, obj_ptr::Ptr{Float64},
        flag_ptr::Ptr{Cint})
        x = unsafe_wrap(Array, x_ptr, Int(n))
        obj_val = MathProgBase.eval_f(current_algencan_model.evaluator, x)
        unsafe_store!(obj_ptr, obj_val)
        unsafe_store!(flag_ptr, 0)
        nothing
    end

    const my_f_c = cfunction(my_f, Void, (Cint, Ptr{Float64}, Ptr{Float64},
        Ptr{Cint}))

    function my_evalg(n::Cint, x_ptr::Ptr{Float64}, g_ptr::Ptr{Float64},
        flag_ptr::Ptr{Cint})
        x = unsafe_wrap(Array, x_ptr, Int(n))
        g = unsafe_wrap(Array, g_ptr, Int(n))
        MathProgBase.eval_grad_f(current_algencan_model.evaluator, g, x)
        unsafe_store!(flag_ptr, 0)
        nothing
    end

    const my_evalg_c = cfunction(my_evalg, Void, (Cint, Ptr{Float64}, Ptr{Float64},
        Ptr{Cint}))

    function my_evalh(n::Cint, x_ptr::Ptr{Float64}, hrow_ptr::Ptr{Cint},
        hcol_ptr::Ptr{Cint}, hval_ptr::Ptr{Float64}, hnnz_ptr::Ptr{Cint}, lim::Cint,
        lmem_ptr::Ptr{UInt8}, flag_ptr::Ptr{Cint})
        # Get nonzero indexes.
        nnz = length(current_algencan_model.h_row_inds)
        if nnz > Int(lim)
            unsafe_store!(lmem_ptr, 1)
            unsafe_store!(flag_ptr, 1)
            return nothing
        end
        unsafe_store!(hnnz_ptr, nnz)
        hcol_ind = unsafe_wrap(Array, hcol_ptr, Int(lim))
        hrow_ind = unsafe_wrap(Array, hrow_ptr, Int(lim))
        # TODO: This is extra work tha only happens because I am calling C instead
        # of Fortran.
        # TODO: Try to move to direct fortran calling.
        hrow_ind[1:nnz] = current_algencan_problem.h_row_inds - 1
        hcol_ind[1:nnz] = current_algencan_problem.h_col_inds - 1

        # Compute the Hessian (for now objective function only)
        σ = 1.0
        μ = Vector{Float64}(0)
        x = unsafe_wrap(Array, x_ptr, Int(n))
        H = unsafe_wrap(Array, hval_ptr, Int(lim))
        MathProgBase.eval_hesslag(current_algencan_problem.evaluator, H, x, σ, μ)
        unsafe_store!(flag_ptr, 0)
        nothing
    end

    const my_evalh_c = cfunction(my_evalh, Void, (Cint, Ptr{Float64}, Ptr{Cint},
        Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Cint, Ptr{UInt8}, Ptr{Cint}))

    # Now CallAlgencan. I will do it slowly, first I create variables for all
    # Algencan's parameters, that are a lot, and then I call it.
    myevalf = my_f_c
    myevalg = my_evalg_c
    myevalh = my_evalh_c
    myevalc = C_NULL
    myevaljac = C_NULL
    myevalhc = C_NULL
    myevalfc = C_NULL
    myevalgjac = C_NULL
    myevalgjacp = C_NULL
    myevalhl = C_NULL
    myevalhlp = C_NULL
    jcnnzmax = length(prob.g_row_inds)
    hnnzmax = length(prob.h_row_inds)

    # Parameters setting
    epsfeas = [1.0e-08]
    epsopt  = [1.0e-08]
    efstin  = sqrt.(epsfeas)
    eostin  = epsopt.^1.5
    efacc   = sqrt.(epsfeas)
    eoacc   = sqrt.(epsopt)

    outputfnm = "algencan.out"
    specfnm   = ""

    vparam = ["ITERATIONS-OUTPUT-DETAIL 1"] #Vector{String}(nvparam)
    nvparam = length(vparam)

    n = prob.n
    l = prob.lb
    u = prob.ub

    # Information on the constraints that do not exist
    m = prob.m
    lambda = Vector{Float64}(m)
    equatn = Vector{UInt8}(m)
    linear = Vector{UInt8}(m)

    coded = zeros(UInt8, 11)
    coded[1:3] = 1
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
        prob.x,
        prob.lb,
        prob.ub,
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
    println("x final ", prob.x)

    prob.obj_val = f[1]
    prob.status = inform[1]

    return inform[1]
end
