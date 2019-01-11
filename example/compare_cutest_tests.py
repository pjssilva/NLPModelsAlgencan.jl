'''
Simple comparasion of two runs of cutest_test.jl.
'''

import sys


def read_results(name):
    res = {}
    for line in open(name):
        cols = line.split()
        res[cols[0]] = cols[1]
    return res


def compare_two(name1, name2):
    res1 = read_results(name1)
    res2 = read_results(name2)
    solved1 = set([i for i in res1 if res1[i] == "Optimal"])
    solved2 = set([i for i in res2 if res2[i] == "Optimal"])
    print("Solved by %s but not by %s =" % (name1, name2), solved1 - solved2)
    print("Solved by %s but not by %s =" % (name2, name1), solved2 - solved1)


if __name__ == "__main__":
    compare_two(sys.argv[1], sys.argv[2])