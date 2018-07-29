# Run a series of CUTEst model to test a solver.
#
# Generates the report solvername_cutest.txt that can be used to create a simple
# performance profile comparing solvers. However it should be taken into
# account that the stopping criterion need to be somewhat similar for a
# meaninful comparison.
#
# Obs: It is set to run Algencan now, but you should easily adapt to other
# solvers.

using MathProgBase
using NLPModels
using CUTEst

using Algencan
# Algencan tolerances
const solver = AlgencanSolver(epsfeas=1.0e-5, epsopt=1.0e-5,
    efstain=3.162278e-03, eostain=3.162278e-08, efacc=3.162278e-3,
    eoacc=3.162278e-3,
    ITERATIONS_OUTPUT_DETAIL=0, NUMBER_OF_ARRAYS_COMPONENTS_IN_OUTPUT=0)
    # SKIP_ACCELERATION_PROCESS=1)
const solver_name = "algencan_hsl_accel"

# # Alternative definition to run Ipopt
# using Ipopt
# const solver = IpoptSolver(tol=1.0e-7, constr_viol_tol=1.0e-5,
#     compl_inf_tol=1.0e-5, print_level=5, print_frequency_iter=100,
#     max_iter=10000, max_cpu_time=3600.0)
# const solver_name = "ipopt"

function cutest_bench(name)
    nlp = CUTEstModel(name)
    model = NLPtoMPB(nlp, solver)
    bench_data = @timed MathProgBase.optimize!(model)
    time = bench_data[2]
    status = MathProgBase.status(model)
    objval = MathProgBase.getobjval(model)
    finalize(nlp)
    return time, status, objval
end

function has_lb_const(lb, ub)
    has_lower = (lb .!= -Inf)
    not_equal = .!(lb .== ub)
    sum(has_lower .& not_equal) > 0
end

# First run to compile
cutest_bench("HS6")

# Grab a list of CUTEst tests

# Tests that generete error in Algencan, probably it tries to compute values
# outside the functions domains.
tests_with_error = ["SPINOP", "DITTERT", "LEUVEN4", "KTMODEL"]
test_problems = CUTEst.select(;min_var=200, max_var=2000, min_con=10)
# Exclude tests that are known to take too long
test_problems = filter(name -> name ∉ tests_with_error, test_problems)
n_tests = length(test_problems)

# Create vector to store result
times = Array{Float64}(0)
status = Array{Symbol}(0)
values = Array{Float64}(0)

# Run benchmarks
report = open(string(solver_name, "_cutest.txt"), "w")
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
    push!(times, t)
    push!(status, s)
    push!(values, v)
    println("\n*************************************************************")
    println("Problem name = ", name)
    println("Times = ", t)
    println("Status = ", s)
    println("Obj values = ", v)
    println("*************************************************************\n")
    line = @sprintf("%-14s%-14s%12.4e\t%12.4e\t%12.4e\n", name, s, t,
        (s == :Optimal ? t : -t), v)
    write(report, line)
    flush(report)
end
close(report)

println("Solved ", size(times)[1], " problems.")
