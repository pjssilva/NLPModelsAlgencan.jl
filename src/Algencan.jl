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

# Globals to avoid closures as Algencan does not allow to send user information
# back to call backs. Long names to avoid conflicts.
global current_algencan_problem
const algencan_lib_path = string(joinpath(ENV["ALGENCAN_LIB_DIR"],
    "libalgencan.so"))

# Imports
import MathProgBase
import MathProgBase.AbstractNLPEvaluator

# TODO: This looks like things to allow for automatic download and
#       compilation of dependencies. Deal with it later.
# if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
#     include("../deps/deps.jl")
# else
#     error("Ipopt not properly installed. Please run Pkg.build(\"Ipopt\")")
# end
# amplexe = joinpath(dirname(libipopt), "..", "bin", "ipopt")

# Exports
export createProblem, addOption
export openOutputFile, setProblemScaling
export solveProblem
export AlgencanProblem

"""
Represents an optimization problem that is solvable by Algencan.
"""
mutable struct AlgencanProblem
    # Problem data
    n::Int                          # Num variables
    m::Int                          # Num constraints
    lb::Vector{Float64}             # Lower bound on vars
    ub::Vector{Float64}             # Upper bound on var
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

    function AlgencanProblem(n, m, lb, ub, sense, evaluator)
        # Get strutural indices of Jacobian and Hessian.
        g_row_inds, g_col_inds = MathProgBase.jac_structure(evaluator)
        h_row_inds, h_col_inds = MathProgBase.hesslag_structure(evaluator)

        # TODO: Verify if I should allow initial points here.
        prob = new(n, m, lb, ub, sense, evaluator, g_row_inds, g_col_inds,
            h_row_inds, h_col_inds, zeros(n), zeros(m), zeros(m), 0.0, 0)
        # Return the object we just made
        prob
    end
end

# TODO: I have to create this dictionary
# # From Ipopt/src/Interfaces/IpReturnCodes_inc.h
# ApplicationReturnStatus = Dict(
# 0=>:Solve_Succeeded,
# 1=>:Solved_To_Acceptable_Level,
# 2=>:Infeasible_Problem_Detected,
# 3=>:Search_Direction_Becomes_Too_Small,
# 4=>:Diverging_Iterates,
# 5=>:User_Requested_Stop,
# 6=>:Feasible_Point_Found,
# -1=>:Maximum_Iterations_Exceeded,
# -2=>:Restoration_Failed,
# -3=>:Error_In_Step_Computation,
# -4=>:Maximum_CpuTime_Exceeded,
# -10=>:Not_Enough_Degrees_Of_Freedom,
# -11=>:Invalid_Problem_Definition,
# -12=>:Invalid_Option,
# -13=>:Invalid_Number_Detected,
# -100=>:Unrecoverable_Exception,
# -101=>:NonIpopt_Exception_Thrown,
# -102=>:Insufficient_Memory,
# -199=>:Internal_Error)

function createProblem(n::Int, x_L::Vector{Float64}, x_U::Vector{Float64},
    m::Int, g_L::Vector{Float64}, g_U::Vector{Float64}, sense::Symbol,
    d::AbstractNLPEvaluator)
    @assert n == length(x_L) == length(x_U)
    @assert m == length(g_L) == length(g_U)

    return(AlgencanProblem(n, m, x_L, x_U, sense, d))
end

# TODO: Implement add option functions
function addOption(prob::AlgencanProblem, keyword::String, value::String)
    # #/** Function for adding a string option.  Returns FALSE the option
    # # *  could not be set (e.g., if keyword is unknown) */
    # if !(isascii(keyword) && isascii(value))
    #     error("IPOPT: Non ASCII parameters not supported")
    # end
    # ret = ccall((:AddIpoptStrOption, libipopt),
    # Cint, (Ptr{Void}, Ptr{UInt8}, Ptr{UInt8}),
    # prob.ref, keyword, value)
    # if ret == 0
    #     error("IPOPT: Couldn't set option '$keyword' to value '$value'.")
    # end
    nothing
end

# TODO: Do I need this?
function openOutputFile(prob::AlgencanProblem, file_name::String, print_level::Int)
    # #/** Function for opening an output file for a given name with given
    # # *  printlevel.  Returns false, if there was a problem opening the
    # # *  file. */
    # if !isascii(file_name)
    #     error("IPOPT: Non ASCII parameters not supported")
    # end
    # ret = ccall((:OpenIpoptOutputFile, libipopt),
    # Cint, (Ptr{Void}, Ptr{UInt8}, Cint),
    # prob.ref, file_name, print_level)
    # if ret == 0
    #     error("IPOPT: Couldn't open output file.")
    # end
    nothing
end

# TODO: Do I need this?
function setProblemScaling(prob::AlgencanProblem, obj_scaling::Float64,
    x_scaling = nothing,
    g_scaling = nothing)
    # #/** Optional function for setting scaling parameter for the NLP.
    # # *  This corresponds to the get_scaling_parameters method in TNLP.
    # # *  If the pointers x_scaling or g_scaling are NULL, then no scaling
    # # *  for x resp. g is done. */
    # x_scale_arg = (x_scaling == nothing) ? C_NULL : x_scaling
    # g_scale_arg = (g_scaling == nothing) ? C_NULL : g_scaling
    # ret = ccall((:SetIpoptProblemScaling, libipopt),
    # Cint, (Ptr{Void}, Float64, Ptr{Float64}, Ptr{Float64}),
    # prob.ref, obj_scaling, x_scale_arg, g_scale_arg)
    # if ret == 0
    #     error("IPOPT: Error setting problem scaling.")
    # end
    nothing
end

function solveProblem(prob::AlgencanProblem)
    global current_algencan_problem = prob

    ###########################################################################
    # Algencan callback function wrappers
    ###########################################################################
    # TODO: Use prob.sense to define "direction" of optimization
    function my_f(n::Cint, x_ptr::Ptr{Float64}, obj_ptr::Ptr{Float64},
        flag_ptr::Ptr{Cint})
        x = unsafe_wrap(Array, x_ptr, Int(n))
        obj_val = MathProgBase.eval_f(current_algencan_problem.evaluator, x)
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
        MathProgBase.eval_grad_f(current_algencan_problem.evaluator, g, x)
        unsafe_store!(flag_ptr, 0)
        nothing
    end

    const my_evalg_c = cfunction(my_evalg, Void, (Cint, Ptr{Float64}, Ptr{Float64},
        Ptr{Cint}))

    function my_evalh(n::Cint, x_ptr::Ptr{Float64}, hrow_ptr::Ptr{Cint},
        hcol_ptr::Ptr{Cint}, hval_ptr::Ptr{Float64}, hnnz_ptr::Ptr{Cint}, lim::Cint,
        lmem_ptr::Ptr{UInt8}, flag_ptr::Ptr{Cint})
        # Get nonzero indexes.
        nnz = length(current_algencan_problem.h_row_inds)
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
    coded[1:1] = 1
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

include("AlgencanSolverInterface.jl")

end # module
