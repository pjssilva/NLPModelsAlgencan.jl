__precompile__()

"""
Dirty and quick interface to use Algencan to solve NLPModels
"""
module NLPMAlgencan

export optimize

using NLPModels


# Compiles to the *static* path of the algencan library
const algencan_lib_path = string(joinpath(ENV["ALGENCAN_LIB_DIR"],
    "libalgencan.so"))

# Global to avoid closures as Algencan does not allow to send user information
# back to call backs. Long names to avoid conflicts.
global current_algencan_nlpmodel

"Compute objective and constraints as required by Algencan"
function nlpm_fc(n::Cint, x_ptr::Ptr{Float64}, obj_ptr::Ptr{Float64},
    m::Cint, g_ptr::Ptr{Float64}, flag_ptr::Ptr{Cint})
    x = unsafe_wrap(Array, x_ptr, Int(n))
    obj_val = obj(current_algencan_nlpmodel, x)
    unsafe_store!(obj_ptr, obj_val)
    g = unsafe_wrap(Array, g_ptr, Int(m))
    cons!(current_algencan_nlpmodel, x, g)
    g .-= current_algencan_nlpmodel.meta.ucon
    unsafe_store!(flag_ptr, 0)
    nothing
end

"Compute objective gradient and constraints Jacobian as required by Algencan"
function nlpm_gjac(n::Cint, x_ptr::Ptr{Float64}, f_grad_ptr::Ptr{Float64},
    m::Cint, jrow_ptr::Ptr{Cint}, jcol_ptr::Ptr{Cint},
    jval_ptr::Ptr{Float64}, jnnz_ptr::Ptr{Cint}, lim::Cint,
    lmem_ptr::Ptr{UInt8}, flag_ptr::Ptr{Cint})

    # Compute gradient of the objective
    x = unsafe_wrap(Array, x_ptr, Int(n))
    f_grad = unsafe_wrap(Array, f_grad_ptr, Int(n))
    grad!(current_algencan_nlpmodel, x, f_grad)

    # I am assuming lmem is always big enough (there is enough memory)
    # Compute Jacobian
    rows, cols, vals= jac_coord(current_algencan_nlpmodel, x)
    nnz = length(rows)
    jcol_ind = unsafe_wrap(Array, jcol_ptr, Int(lim))
    jrow_ind = unsafe_wrap(Array, jrow_ptr, Int(lim))
    J = unsafe_wrap(Array, jval_ptr, Int(lim))
    jrow_ind[1:nnz] .= rows - 1
    jcol_ind[1:nnz] .= cols - 1
    J[1:nnz] .= vals
    unsafe_store!(jnnz_ptr, nnz)
    unsafe_store!(flag_ptr, 0)
    unsafe_store!(lmem_ptr, 0)
    nothing
end

" Compute the Hessian of the Lagrangian as required by Algencan"
function nlpm_hl(n::Cint, x_ptr::Ptr{Float64}, m::Cint,
    mult_ptr::Ptr{Float64}, scale_f::Float64, scale_g_ptr::Ptr{Float64},
    hrow_ptr::Ptr{Cint}, hcol_ptr::Ptr{Cint},
    hval_ptr::Ptr{Float64}, hnnz_ptr::Ptr{Cint}, lim::Cint,
    lmem_ptr::Ptr{UInt8}, flag_ptr::Ptr{Cint})

    σ = scale_f
    μ = unsafe_wrap(Array, mult_ptr, Int(m))
    scale_g = unsafe_wrap(Array, scale_g_ptr, Int(m))
    μ .*= scale_g
    x = unsafe_wrap(Array, x_ptr, Int(n))
    rows, cols, vals = hess_coord(current_algencan_nlpmodel, x; obj_weight=σ, y=μ)
    nnz = length(rows)

    hcol_ind = unsafe_wrap(Array, hcol_ptr, Int(lim))
    hrow_ind = unsafe_wrap(Array, hrow_ptr, Int(lim))
    H = unsafe_wrap(Array, hval_ptr, Int(lim))
    hrow_ind[1:nnz] .= rows - 1
    hcol_ind[1:nnz] .= cols - 1
    H[1:nnz] .= vals
    unsafe_store!(hnnz_ptr, nnz)

    unsafe_store!(flag_ptr, 0)
    unsafe_store!(lmem_ptr, 0)
    nothing
end

function optimize(model::AbstractNLPModel)
    global current_algencan_nlpmodel = model

    ###########################################################################
    # Algencan callback function wrappers
    ###########################################################################

    const c_nlpm_fc = cfunction(nlpm_fc, Void, (Cint, Ptr{Float64},
        Ptr{Float64}, Cint, Ptr{Float64}, Ptr{Cint}))

    const c_nlpm_gjac = cfunction(nlpm_gjac, Void, (Cint, Ptr{Float64},
        Ptr{Float64}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Cint,
        Ptr{UInt8}, Ptr{Cint}))

    const c_nlpm_hl = cfunction(nlpm_hl, Void, (Cint, Ptr{Float64}, Cint,
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
    myevalfc = c_nlpm_fc
    myevalgjac = c_nlpm_gjac
    myevalgjacp = C_NULL
    myevalhl = c_nlpm_hl
    myevalhlp = C_NULL
    jcnnzmax = model.meta.nvar*model.meta.ncon
    hnnzmax = model.meta.nvar*model.meta.nvar
    coded = zeros(UInt8, 11)
    coded[7] = 1
    coded[8] = 1
    coded[10] = 1
    checkder = UInt8(0)

    # Parameters controling precision
    epsfeas = [1.0e-8]
    epsopt = [1.0e-8]
    efstin = [sqrt.(1.0e-8)]
    eostin = [(1.0e-8)^1.5]
    efacc  = [sqrt.(1.0e-8)]
    eoacc  = [sqrt.(1.0e-8)]

    # Extra parameters
    outputfnm = ""
    specfnm   = ""
    vparam = Vector{String}(0)
    nvparam = length(vparam)

    # Return information
    f = [0.0]
    cnorm = [0.0]
    snorm = [0.0]
    nlpsupn = [0.0]
    inform = Vector{Cint}([0])

    mult = zeros(model.meta.ncon)
    equality = zeros(UInt8, model.meta.ncon)
    equality[model.meta.jfix] = 1
    linear  = zeros(UInt8, model.meta.ncon)
    linear[model.meta.lin] = 1

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
        nvparam, vparam, model.meta.nvar, model.meta.x0, model.meta.lvar,
        model.meta.uvar, model.meta.ncon, mult, equality, linear,
        coded, checkder, f, cnorm, snorm, nlpsupn, inform
    )

    obj_val = f[1]
    status = :Optimal #find_status(model, cnorm[1], snorm[1], nlpsupn[1], Int(inform[1]))
    println("Final objetive = $obj_val, final status = $status")
    return Int(inform[1])
end

end
