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
    inner::AlgencanProblem
    options

    function AlgencanMathProgModel(;options...)
        model = new()
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
    m.inner = createProblem(numVar, float(x_l), float(x_u), numConstr,
        float(g_lb), float(g_ub), sense, d)
end

# Simple functions
getsense(m::AlgencanMathProgModel) = m.inner.sense

numvar(m::AlgencanMathProgModel) = m.inner.n

numconstr(m::AlgencanMathProgModel) = m.inner.m

# TODO: Allow this
numlinconstr(m::AlgencanMathProgModel) = 0

numquadconstr(m::AlgencanMathProgModel) = 0

function optimize!(m::AlgencanMathProgModel)
    #copy!(m.inner.x, m.warmstart) # set warmstart
    # TODO: No options for now
    # for (name,value) in m.options
    #     sname = string(name)
    #     if match(r"(^resto_)", sname) != nothing
    #         sname = replace(sname, r"(^resto_)", "resto.")
    #     end
    #     addOption(m.inner, sname, value)
    # end
    solveProblem(m.inner)
end

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

getobjval(m::AlgencanMathProgModel) = m.inner.obj_val * (m.inner.sense == :Max ? -1 : +1)

getsolution(m::AlgencanMathProgModel) = m.inner.x

function getreducedcosts(m::AlgencanMathProgModel)
    # TODO: Verify, I am not thinking about constraints yet.
    # sense = m.inner.sense
    # redcost = m.inner.mult_x_U - m.inner.mult_x_L
    # return sense == :Max ? redcost : -redcost
    return zeros(m.inner.m)
end

function getconstrduals(m::AlgencanMathProgModel)
    # TODO: Verify, I am not thinking about constraints yet.

    v = m.inner.mult_g # return multipliers for all constraints
    return m.inner.sense == :Max ? copy(v) : -v
end

getrawsolver(m::AlgencanMathProgModel) = m.inner

setwarmstart!(m::AlgencanMathProgModel, x) = (m.inner.x = x)
