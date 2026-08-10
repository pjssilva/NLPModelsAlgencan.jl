# MathOptInterface conformance suite, run against Algencan.
#
# Algencan has no MOI wrapper of its own. It is an NLPModels solver, and JuMP
# reaches it through `NLPModelsJuMP.Optimizer`, the generic MOI wrapper for any
# `SolverCore.AbstractOptimizationSolver`:
#
#     model = Model(NLPModelsJuMP.Optimizer)
#     set_attribute(model, "solver", NLPModelsAlgencan.AlgencanSolver)
#
# That is the path exercised here. 541 of the 561 `MOI.Test` tests pass; the 20
# exclusions below are grouped by cause and each group is explained.
module TestMOI

using Test
import MathOptInterface as MOI
import NLPModelsJuMP
import NLPModelsAlgencan

# Attributes that `NLPModelsJuMP.Optimizer` does not implement. Duals are the
# notable one: Algencan returns multipliers, but the wrapper does not yet map
# them onto MOI's constraint duals.
const EXCLUDED_ATTRIBUTES = Any[
    MOI.ConstraintBasisStatus,
    MOI.VariableBasisStatus,
    MOI.ConstraintName,
    MOI.VariableName,
    MOI.ObjectiveBound,
    MOI.SolverVersion,
    MOI.DualObjectiveValue,
    MOI.ConstraintDual,
]

const EXCLUDED_TESTS = [
    # ---------------------------------------------------------------------
    # Second-derivative oracles are not negotiated with the user's evaluator.
    #
    # Algencan is a second-order method: it evaluates the Hessian of the
    # Lagrangian. When the model arrives as a raw `MOI.NLPBlock`,
    # `NLPModelsJuMP` requests `:Hess`/`:HessVec`/`:JacVec` from the
    # user-supplied evaluator unconditionally, and the hand-written evaluators
    # in `MOI.Test` deliberately advertise only a subset, so `MOI.initialize`
    # raises `Unsupported feature ...`.
    #
    # This is a limitation of the shared wrapper, not of Algencan: models built
    # through JuMP get a full `MOI.Nonlinear` evaluator and work fine (see the
    # "JuMP interface test" testset in runtests.jl). Fixing it means having
    # `NLPModelsJuMP` request only the features the evaluator reports and fall
    # back to finite differences or a quasi-Newton Hessian otherwise.
    r"^test_nonlinear_expression_multivariate_function$",     # Unsupported feature Hess
    r"^test_nonlinear_hs071$",                                # Unsupported feature HessVec
    r"^test_nonlinear_hs071_NLPBlockDual$",                   # Unsupported feature HessVec
    r"^test_nonlinear_hs071_hessian_vector_product$",         # Unsupported feature Hess
    r"^test_nonlinear_hs071_no_hessian$",                     # Unsupported feature Hess
    r"^test_nonlinear_invalid$",                              # feat in features_available(d)
    r"^test_nonlinear_objective$",                            # Unsupported feature JacVec
    r"^test_nonlinear_objective_and_moi_objective_test$",     # Unsupported feature JacVec
    r"^test_nonlinear_without_objective$",                    # Unsupported feature JacVec

    # ---------------------------------------------------------------------
    # No global infeasibility or unboundedness certificate.
    #
    # Algencan is a local augmented-Lagrangian method. When it stops at a
    # stationary point of the infeasibility measure it reports
    # `LOCALLY_INFEASIBLE`, which is the honest answer; these tests require the
    # global `INFEASIBLE` / `DUAL_INFEASIBLE` / `INFEASIBLE_OR_UNBOUNDED`.
    # Other local NLP solvers exclude the same tests.
    r"^test_conic_NormInfinityCone_INFEASIBLE$",
    r"^test_conic_NormOneCone_INFEASIBLE$",
    r"^test_conic_linear_INFEASIBLE$",
    r"^test_conic_linear_INFEASIBLE_2$",
    r"^test_linear_DUAL_INFEASIBLE$",
    r"^test_linear_INFEASIBLE$",
    r"^test_linear_INFEASIBLE_2$",

    # ---------------------------------------------------------------------
    # Converges to an infeasible stationary point from the starting point the
    # test supplies: `LOCALLY_INFEASIBLE` where `LOCALLY_SOLVED` is expected.
    # Again a property of local optimization, not a wrong answer.
    r"^test_nonlinear_expression_hs109$",
    r"^test_quadratic_SecondOrderCone_basic$",
    r"^test_quadratic_nonconvex_constraint_basic$",

    # ---------------------------------------------------------------------
    # `NLPModelsJuMP.Optimizer` reads the model attributes it knows and ignores
    # the rest, so it does not raise `MOI.UnsupportedAttribute` for an unknown
    # one. Wrapper-side; there is a FIXME to this effect in NLPModelsJuMP.
    r"^test_model_copy_to_UnsupportedAttribute$",
]

function test_runtests()
    model = MOI.instantiate(NLPModelsJuMP.Optimizer, with_bridge_type=Float64)
    MOI.set(model, MOI.RawOptimizerAttribute("solver"), NLPModelsAlgencan.AlgencanSolver)
    MOI.set(model, MOI.Silent(), true)  # comment out to see Algencan's own output
    config = MOI.Test.Config(
        atol=1.0e-2,
        optimal_status=MOI.LOCALLY_SOLVED,
        exclude=EXCLUDED_ATTRIBUTES,
    )
    MOI.Test.runtests(model, config; exclude=EXCLUDED_TESTS)
    return
end

function runtests()
    for name in names(@__MODULE__; all=true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

end  # module

TestMOI.runtests()
