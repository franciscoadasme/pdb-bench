import sys
from shared import timeit

from chemfiles import Selection, Trajectory


def distance(f):
    ai = Selection("resid 50").evaluate(f)
    bi = Selection("resid 60").evaluate(f)
    min_dist = float("inf")
    for i in ai:
        for j in bi:
            d = f.distance(i, j)
            if d < min_dist:
                min_dist = d


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

f = Trajectory(pdbfile).read()

print(timeit(distance, f, repeats=repeats))
