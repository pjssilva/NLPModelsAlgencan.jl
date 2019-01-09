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

function cutest_bench(name)
    nlp = CUTEstModel(name)
    model = NLPtoMPB(nlp, solver)
    bench_data = @timed MPB.optimize!(model)
    finalize(nlp)
    etime = bench_data[2]
    flag = MPB.status(model)
    objval = MPB.getobjval(model)
    n_fc, n_ggrad, n_hl, n_hlp = MPB.getnfevals(model)
    return flag, etime, n_fc, n_ggrad, n_hl, n_hlp, objval
end

function has_lb_const(lb, ub)
    has_lower = (lb .!= -Inf)
    not_equal = .!(lb .== ub)
    sum(has_lower .& not_equal) > 0
end

# First run to compile
cutest_bench("HS6")

# Grab a list of CUTEst tests
test_problems = readlines(open("cutest_selection.txt"))
n_tests = length(test_problems)

# Run tests
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
