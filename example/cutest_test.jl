push!(LOAD_PATH, "../src/")
using Ipopt
using MathProgBase
using NLPModels
using CUTEst
using Algencan
using BenchmarkProfiles
using Plots

function cutest_bench(name)
    times, status, values = zeros(1, 2), Array{Symbol}(1, 2), zeros(1, 2)

    # DISCLAIMER: I have chosen the values below to make the stopping criteria
    # closer to each other. However both methods have different stopping
    # criteria as Algencan uses always unscaled versions while Ipopt
    # uses scaled version in the tol criteria. Therefore this benchmarks can
    # only give an idea on how the methods compare and it is not intended to
    # be a definite (note even close) comparison.

    solvers = [
        IpoptSolver(tol=1.0e-7, constr_viol_tol=1.0e-5, compl_inf_tol=1.0e-5,
            print_level=2),
        AlgencanSolver(epsfeas=1.0e-5, epsopt=1.0e-5,
            INNER_ITERATIONS_LIMIT=1000,
            ITERATIONS_OUTPUT_DETAIL=0,
            NUMBER_OF_ARRAYS_COMPONENTS_IN_OUTPUT=0)
    ]

    for i in 1:length(solvers)
        solver = solvers[i]
        nlp = CUTEstModel(name)
        model = NLPtoMPB(nlp, solver)
        bench_data = @timed MathProgBase.optimize!(model)
        times[1, i] = bench_data[2]
        status[1, i] = MathProgBase.status(model)
        values[1, i] = MathProgBase.getobjval(model)
        finalize(nlp)
    end
    return times, status, values
end

function has_lb_const(lb, ub)
    has_lower = (lb .!= -Inf)
    not_equal = .!(lb .== ub)
    sum(has_lower .& not_equal) > 0
end

# First run to compile
cutest_bench("HS6")

# Grab a list of CUTEst tests
test_problems = CUTEst.select(;min_var=10, max_var=1000, min_con=10, max_con=1000)
n_tests = length(test_problems)

# Alocate matrices to store test results
n_solvers = 2
times = Array{Float64}(0, n_solvers)
status = Array{Symbol}(0, n_solvers)
values = Array{Float64}(0, n_solvers)

# Run benchmarks
for i = 1:n_tests
    name = test_problems[i]
    println("\nSolving Problem $name - $i of $n_tests.\n")

    nlp = CUTEstModel(name)
    skip = has_lb_const(nlp.meta.lcon, nlp.meta.ucon)
    finalize(nlp)
    if !skip
        println("\n*************************************************************")
        println("Problem $name ")
        println("it has lower bound constraints check result")
        println("*************************************************************\n")
    end
    t, s, v = cutest_bench(name)
    times = [times; t]
    status = [status; s]
    values = [values; v]
    println("\n*************************************************************")
    println("Times = ", t)
    println("Status = ", s)
    println("Obj values = ", v)
    println("*************************************************************\n")
end
println("Solved ", size(times)[1], " problems.")

# Generate a simple performance profile
T = copy(times)
T[status[:,1] .!= :Optimal, 1] *= -1
T[status[:,2] .!= :Optimal, 2] *= -1
performance_profile(T, ["Ipopt", "Algencan"], title="Simple test")
savefig("profile.png")
