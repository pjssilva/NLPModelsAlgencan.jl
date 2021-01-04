
"""
Algencan interface to NLPModels.
See its [GitHub page](https://github.com/pjssilva/NLPModelsAlgencan.jl)
"""
module NLPModelsAlgencan

using LinearAlgebra, SparseArrays, NLPModels, SolverTools
import Libdl

# Gets the path of the Algencan library
if "ALGENCAN_LIB_DIR" in keys(ENV)
    const algencan_lib_path = string(joinpath(ENV["ALGENCAN_LIB_DIR"],
        "libalgencan.so"))
elseif isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
    const algencan_lib_path = libalgencan
else
    error("Algencan not properly installed. Please run Pkg.build(\"Algencan\")")
    error("or set the ALGENCAN_LIB_DIR enviroment variable.")
end

export algencan

"Processes and stores data from an AbstractNLPModel to be used by Algencan
subroutines"
mutable struct AlgencanModelData
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
    j_row_inds::Vector{Int}         # NNZ row indexes of sparse const. Jacobian.
    j_col_inds::Vector{Int}         # NNZ col indexes of sparse const. Jacobian.
    h_row_inds::Vector{Int}         # NNZ row indexes of sparse Lag. Hessian.
    h_col_inds::Vector{Int}         # NNZ col indexes of sparse Lag. Hessian.
    is_equality::Vector{UInt8}      # 1 if equality, 0 if inequality
    is_g_linear::Vector{UInt8}      # 1 if constraint is linear, 0 otherwise

    # Decision variables and constraints values.
    x::Vector{Float64}              # Starting and final solution
    g::Vector{Float64}              # Final constraint values
    mult::Vector{Float64}           # Final Lagrange multipliers on constraints
    obj_val::Float64                # Final objective
    status::Symbol                  # Final status

    nlp::AbstractNLPModel

    # Options to be passed to the solver
    options

    function AlgencanModelData(nlp::AbstractNLPModel)
        model = new()

        model.nlp = nlp;

        model.sense = nlp.meta.minimize ? 1.0 : -1.0
        model.n = nlp.meta.nvar
        model.m = nlp.meta.ncon
        model.x = copy(nlp.meta.x0)
        model.lb = copy(nlp.meta.lvar)
        model.ub = copy(nlp.meta.uvar)
        model.mult = copy(nlp.meta.y0)

        jrow_inds, jcol_inds = jac_structure(nlp)
        model.j_row_inds, model.j_col_inds = jrow_inds .- 1, jcol_inds .- 1
        hrow_inds, hcol_inds = hess_structure(nlp)
        model.h_row_inds, model.h_col_inds = hrow_inds .- 1, hcol_inds .- 1

        model.options = Dict(
            :epsfeas=>1.0e-08,
            :epsopt=>1.0e-08,
            :efstain=>sqrt.(1.0e-8),
            :eostain=>(1.0e-8)^1.5,
            :efacc=>sqrt(1.0e-8),
            :eoacc=>sqrt(1.0e-8),
            :outputfnm=>"",
            :specfnm=>""
        )

        g_lb, g_ub = float(copy(nlp.meta.lcon)), float(copy(nlp.meta.ucon))
        model.g_sense, model.g_two_sinvmap,
                        model.g_two_smap = treat_lower_bounds(nlp, g_lb, g_ub)
        model.g_has_lb = length(model.nlp.meta.jrng) +
                            length(model.nlp.meta.jlow) > 0
        model.g_lb, model.g_ub = g_lb, g_ub

        # Contraints with only lower bound will be multiplied by -1.0, hence
        # their lower bounds becomes minus the upper bound
        model.g_ub[nlp.meta.jlow] .= -g_lb[nlp.meta.jlow]
        model.g_lb[nlp.meta.jlow] .= -Inf

        # Constraints types
        model.is_equality = zeros(UInt8, model.m)
        model.is_equality[nlp.meta.jfix] .= 1
        model.is_g_linear = zeros(UInt8, model.m)
        model.is_g_linear[nlp.meta.lin] .= 1

        return model
    end
end

# TODO: describe the keywords better
"""`output = algencan(nlp; kwargs...)`
Solves the `NLPModel` problem `nlp` using `Algencan`.

# Optional keyword arguments
* `epsfeas`: feasibility tolerance
* `epsopt`: optimality tolerance
* `efstain`: stainf feasibility tolerance
* `eostain`: stainf output tolerance
* `efacc`: acc feasibility tolerance
* `eoacc`: acc optimality tolerance
* `outputfnm`: output filename
* `specfnm`: specification filename

All other keyword arguments will be passed to Algencan as an option.
See Birgin and Martínez [1] on page 120 for the list of options accepted.

[1] Birgin, Ernesto G., and José Mario Martínez. Practical augmented Lagrangian
    methods for constrained optimization. Society for Industrial and Applied
    Mathematics, 2014. https://books.google.com.br/books?id=og1_AwAAQBAJ.
"""
function algencan(nlp::AbstractNLPModel; kwargs...)
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

    model = AlgencanModelData(nlp)

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
    # Since I can duplicate the two sided constraints, I need
    # to take that into account here
    jcnnzmax = 2*model.nlp.meta.nnzj
    # Using the same workaround as Birgin in the CUTEst interface.
    # See comment below:
    # hnnzmax must be an upper bound on the number of elements of the
    # Hessian of the Lagrangian plus the new elements that appear when
    # adding (to the Hessian of the Lagrangian) matrix rho \sum_j \nabla
    # cj(x) \nabla cj(x)^t to obtain the Hessian of the Augmented
    # Lagrangian. But this additional space is only need when Algencan
    # uses an Euclidian trust-region approach that is recommended only
    # for problems with no more than 500 variables. Therefore, since we
    # are not able to estimate this quantity here (when using CUTEst),
    # we arbitrarily add to hnnzmax the quantity 1,000,000 >= 500^2.
    hnnzmax = model.nlp.meta.nnzh + 1000000
    coded = zeros(UInt8, 11)
    coded[7] = UInt8(1)
    coded[8] = UInt8(1)
    coded[10] = UInt8(1)
    coded[11] = UInt8(0)
    checkder = UInt8(0)

    # Deal with lower bounds
    m = model.m + length(model.g_two_smap)
    mult = zeros(m)
    mult[1:model.m] .= model.mult[1:model.m]
    is_equality = zeros(UInt8, m)
    is_equality[1:model.m] .= model.is_equality
    is_g_linear = zeros(UInt8, m)
    is_g_linear[1:model.m] .= model.is_g_linear
    for i = 1:length(model.g_two_smap)
        is_g_linear[model.m + i] = model.is_g_linear[model.g_two_smap[i]]
        mult[model.m + i] = -mult[model.g_two_smap[i]]
    end

    # Optional keyword arguments
    for (key, value) in kwargs
        model.options[key] = value
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

    if !isempty(specfnm)
        read_options_from_specification_file(model)
    end

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
    try
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
    finally
        Libdl.dlclose(algencandl)
        @assert !(algencan_lib_path in Libdl.dllist())
    end

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

    Δt = (time_ns() - start_time) / 1.0e+9

    return GenericExecutionStats(model.status, model.nlp, solution=model.x,
                                 objective=model.obj_val,
                                 dual_feas=max(nlpsupn[1], snorm[1]),
                                 primal_feas=cnorm[1],
                                 elapsed_time=Δt,
                                 multipliers=model.mult[1:model.m]
                                )
end

"Read additional parameters present in the specification file"
function read_options_from_specification_file(model::AlgencanModelData)
    # Options that are present in our data model
    spec_params = Dict(
        "FEASIBILITY-TOLERANCE"=>:epsfeas,
        "OPTIMALITY-TOLERANCE"=>:epsopt,
        "STAINF-FEASIBILITY-TOLERANCE"=>:efstain,
        "STAINF-OPTIMALITY-TOLERANCE"=>:eostain,
        "ACC-FEASIBILITY-THRESHOLD"=>:efacc,
        "ACC-OPTIMALITY-THRESHOLD"=>:eoacc
    )

    specfnm = model.options[:specfnm]
    open(specfnm, "r") do io
        for line in split(read(io, String), "\n")
            # Skip comments
            line  = split(line, ['#', '*'])[1]

            # Only key and value
            words = split(line, ' ')
            if length(words) == 2
                key = uppercase(words[1])

                # Set the options that are present
                if key in keys(spec_params)
                    value = parse(Float64, words[2])
                    opt_key = spec_params[key]
                    model.options[opt_key] = value
                end
            end
        end
    end
end

"Analyse the lower and upper bounds on the constraints and prepare the
 data structure to treat lower bounds."
function treat_lower_bounds(nlp::AbstractNLPModel, lb, ub)
    m = length(lb)

    # Allow to treat constraints that have only lower bounds
    sense = ones(m)
    sense[nlp.meta.jlow] .= -1.0

    # Treat two side constraints
    two_smap = nlp.meta.jrng
    new_ind = 1
    two_sinvmap = zeros(m)
    for i = two_smap
        two_sinvmap[i] = new_ind
        new_ind += 1
    end

    return sense, two_sinvmap, two_smap
end

"Transform the option dictionary into a vparam string array"
function option2vparam(model::AlgencanModelData)
    parameters = [:epsfeas, :epsopt, :efstain, :eostain, :efacc, :eoacc,
        :outputfnm, :specfnm]
    vparam = Array{String}(undef, 0)
    for option in model.options
        key, value = option
        if key in parameters
            continue
        end
        key = replace(string(key), "_" => "-")
        push!(vparam, "$key $value")
    end
    return vparam
end

"Find final status of Algencan"
function find_status(model::AlgencanModelData, cnorm::Float64, snorm::Float64,
    nlpsupn::Float64, inform::Int)

    if inform != 0
        return :exception
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
        if !feasible
            return :infeasible
        elseif !(bounded_obj && bounded_mult)
            return :unbounded
        elseif nlpsupn <= epsopt && snorm <= epsopt
            return :first_order
        else
            return :exception # :exception or :unknown ?
        end
    else
        if !bounded_obj
            return :unbounded
        elseif nlpsupn <= epsopt
            return :first_order
        else
            return :exception
        end
    end
end

###########################################################################
# Algencan callbacks
###########################################################################

"Compute objective and constraints as required by Algencan"
function julia_fc(model::AlgencanModelData, n::Cint, x_ptr::Ptr{Float64}, obj_ptr::Ptr{Float64},
    m::Cint, g_ptr::Ptr{Float64}, flag_ptr::Ptr{Cint})

    # Evaluate objective and constraints
    x = unsafe_wrap(Array, x_ptr, Int(n))
    obj_val = obj(model.nlp, x)
    unsafe_store!(obj_ptr, model.sense * obj_val)
    g = unsafe_wrap(Array, g_ptr, Int(m))
    g[1:model.m] .= cons(model.nlp, x)

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
function julia_gjac(model::AlgencanModelData, n::Cint, x_ptr::Ptr{Float64}, f_grad_ptr::Ptr{Float64},
    m::Cint, jrow_ptr::Ptr{Cint}, jcol_ptr::Ptr{Cint},
    jval_ptr::Ptr{Float64}, jnnz_ptr::Ptr{Cint}, lim::Cint,
    lmem_ptr::Ptr{UInt8}, flag_ptr::Ptr{Cint})

    # Compute gradient of the objective
    x = unsafe_wrap(Array, x_ptr, Int(n))
    f_grad = unsafe_wrap(Array, f_grad_ptr, Int(n))

    grad!(model.nlp, x, f_grad)
    f_grad[1:Int(n)] .*= model.sense

    # Find structure of the constraints Jacobian
    nnz = model.nlp.meta.nnzj
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
    J[1:nnz] .= jac_coord(model.nlp, x)

    # Treat the presence of lower bound in the constraints
    if model.g_has_lb
        @inbounds for i = 1:model.nlp.meta.nnzj
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
function julia_hl(model::AlgencanModelData, n::Cint, x_ptr::Ptr{Float64}, m::Cint,
    mult_ptr::Ptr{Float64}, scale_f::Float64, scale_g_ptr::Ptr{Float64},
    hrow_ptr::Ptr{Cint}, hcol_ptr::Ptr{Cint},
    hval_ptr::Ptr{Float64}, hnnz_ptr::Ptr{Cint}, lim::Cint,
    lmem_ptr::Ptr{UInt8}, flag_ptr::Ptr{Cint})

    # Get nonzero indexes.
    nnz = model.nlp.meta.nnzh
    if nnz > Int(lim)
        unsafe_store!(lmem_ptr, Cint(1))
        unsafe_store!(flag_ptr, Cint(1))
        return nothing
    else
        unsafe_store!(lmem_ptr, Cint(0))
    end
    unsafe_store!(hnnz_ptr, Cint(nnz))
    hcol_ind = unsafe_wrap(Array, hcol_ptr, Int(lim))
    hrow_ind = unsafe_wrap(Array, hrow_ptr, Int(lim))
    hrow_ind[1:nnz] .= model.h_row_inds
    hcol_ind[1:nnz] .= model.h_col_inds

    # Compute scaled multipliers
    σ = model.sense * scale_f
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
    H[1:nnz] = hess_coord(model.nlp, x, μ; obj_weight = σ)

    # Declare success
    unsafe_store!(flag_ptr, Cint(0))
    nothing
end


"Compute the Hessian of the Lagrangian times p as required by Algencan"
function julia_hlp(nlp::AlgencanModelData, n::Cint, x_ptr::Ptr{Float64}, m::Cint,
    mult_ptr::Ptr{Float64}, scale_f::Float64, scale_g_ptr::Ptr{Float64},
    p_ptr::Ptr{Float64}, hp_ptr::Ptr{Float64}, goth_ptr::Ptr{UInt8},
    flag_ptr::Ptr{Cint})

    # Compute scaled multipliers
    σ =  model.sense * scale_f
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
    hprod!(model.nlp, x, μ, p, hp; obj_weight=σ)

    # Declare success
    unsafe_store!(flag_ptr, Cint(0))
    nothing
end

end
