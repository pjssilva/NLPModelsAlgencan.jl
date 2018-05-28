push!(LOAD_PATH, "../src/")
using Ipopt
using MathProgBase
using NLPModels
using CUTEst
using Algencan

function cutest_bench(name)
    times, status, values = zeros(2), Vector{Symbol}(2), zeros(2)

    solvers = [IpoptSolver(), AlgencanSolver(epsfeas=1.0e-5, epsopt=1.0e-5, skip_acceleration_process=1)]
    for i in 1:length(solvers)
        solver = solvers[i]
        nlp = CUTEstModel(name)
        model = NLPtoMPB(nlp, solver)
        bench_data = @timed MathProgBase.optimize!(model)
        times[i] = bench_data[2]
        status[i] = MathProgBase.status(model)
        values[i] = MathProgBase.getobjval(model)
        finalize(nlp)
    end
    return times, status, values
end

# First run to compile
cutest_bench("HS6")

# Grab a list of CUTEst tests
test_problems = CUTEst.select(;min_var=10, max_var=1000, min_con=10, max_con=1000)
n_tests = length(test_problems)

# Alocate matrices to store test results
n_solvers = 2
times = Array{Float64}(n_tests, n_solvers)
status = Array{Symbol}(n_tests, n_solvers)
values = Array{Float64}(n_tests, n_solvers)

# Run benchmarks
for i = 1:n_tests
    println("Solving test ", test_problems[i], "  - $i of $n_tests.")
    times[i,:], status[i, :], values[i, :] = cutest_bench(test_problems[i])
    println("Times = ", times[i,:])
    println("Status = ", status[i,:])
    println("Obj values = ", values[i,:])
end
println("Times =")
print(times)
println("Status =")
println(status)
println("Obj values =")
println(values)
