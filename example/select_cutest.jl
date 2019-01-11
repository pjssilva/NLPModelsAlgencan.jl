"""
Select some CUTEst problem for tests.
"""

using CUTEst

# Grab a list of CUTEst tests
test_problems = CUTEst.select(;min_var=5, max_var=2000, min_con=2)

# Avoid tests that generete error in Algencan, probably it tries to compute
# values outside the functions domains and tests that take too long (Algencan
# does not have a timeout option).
avoid = ["SPINOP", "DITTERT", "LEUVEN4", "KTMODEL", "TRO21X5", "NUFFIELD", "SPIN"]
test_problems = filter(name -> name ∉ avoid, test_problems)

# Save the list of test_problems
selection = open("cutest_selection.txt", "w")
for p in test_problems
    write(selection, string(p, "\n"))
end
close(selection)