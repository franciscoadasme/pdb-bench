import sys
from shared import timeit

from chemfiles import Selection, Trajectory


def count(f):
    len(Selection("resname ALA and name CA").evaluate(f))


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

f = Trajectory(pdbfile).read()

print(timeit(count, f, repeats=repeats))
