import sys
from shared import timeit

from chemfiles import Trajectory


def parse(filepath):
    t = Trajectory(filepath)
    for _ in range(t.nsteps):
        t.read()


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

print(timeit(parse, pdbfile, repeats=repeats))
