# Run a series of CUTEst models to test a solver.
#
# Generates the report solvername_cutest.txt that can be used to create a simple
# performance profile comparing solvers. However it should be taken into
# account that the stopping criterion need to be somewhat similar for a
# meaninful comparison.
#
# Obs: It is set to run Algencan now, but you should easily adapt to other
# solvers.

using Printf
using MathProgBase
MPB = MathProgBase
using NLPModels
using NLPModelsJuMP
using CUTEst

# # Algencan tolerances
# using Algencan
# const solver = AlgencanSolver(epsfeas=1.0e-5, epsopt=1.0e-5, specfnm="algencan.dat")
# const solver_name = "algencan_hsl_accel"

# Alternative definition to run Ipopt
using Ipopt
const solver = IpoptSolver(tol=1.0e-7, constr_viol_tol=1.0e-5,
    compl_inf_tol=1.0e-5, print_level=5, print_frequency_iter=100,
    max_iter=10000, max_cpu_time=3600.0)
const solver_name = "ipopt"

function cutest_bench(name)
    nlp = CUTEstModel(name)
    model = NLPtoMPB(nlp, solver)
    bench_data = @timed MPB.optimize!(model)
    finalize(nlp)
    etime = bench_data[2]
    flag = MPB.status(model)
    objval = MPB.getobjval(model)
    # n_fc, n_ggrad, n_hl, n_hlp = MPB.getnfevals(model)
    # return flag, etime, n_fc, n_ggrad, n_hl, n_hlp, objval
    return flag, etime, 0, 0, 0, 0, objval
end

function has_lb_const(lb, ub)
    has_lower = (lb .!= -Inf)
    not_equal = .!(lb .== ub)
    sum(has_lower .& not_equal) > 0
end

# First run to compile
cutest_bench("HS6")

# Grab a list of CUTEst tests
test_problems = CUTEst.select(;min_var=200, max_var=2000, min_con=10)

# Avoid tests that generete error in Algencan, probably it tries to compute
# values outside the functions domains and tests that take too long (Algencan
# does not have a timeout option).
avoid = ["SPINOP", "DITTERT", "LEUVEN4", "KTMODEL", "TRO21X5", "NUFFIELD", "SPIN"]
# Exclude problematic tests
test_problems = filter(name -> name ∉ avoid, test_problems)
n_tests = length(test_problems)

# Save the list of test_problems
selection = open("cutest_selection.txt", "w")
for p in test_problems
    write(selection, string(p, "\n"))
end
close(selection)

# Run benchmarks
report = open(string(solver_name, "_cutest.txt"), "w")
for i = 1:n_tests
    name = test_problems[i]
    println("\nSolving Problem $name - $i of $n_tests.\n")

    nlp = CUTEstModel(name)
    skip = has_lb_const(nlp.meta.lcon, nlp.meta.ucon)
    finalize(nlp)
    if skip
        println("\n*************************************************************")
        println("Problem $name ")
        println("it has lower bound constraints check result")
        println("*************************************************************\n")
    end
    s, t, fc, ggrad, hl, hlp, v = cutest_bench(name)
    println("\n*************************************************************")
    println("Problem name = ", name)
    println("Performance = ", t, fc, ggrad, hl, hlp)
    println("Status = ", s)
    println("Obj value = ", v)
    println("*************************************************************\n")
    line = @sprintf("%-14s%-14s%12.4e\t%10.4d\t%10.4d\t%10.4d\t%10.4d\t%12.4e\n",
        name, s, t, fc, ggrad, hl, hlp, v)
    write(report, line)
    flush(report)
end
close(report)

println("Solved ", n_tests, " problems.")
