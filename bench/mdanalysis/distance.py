import sys
from shared import timeit

import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array


def distance(u):
    segA = u.segments[0]
    r50 = segA.atoms.select_atoms("resid 50")
    r60 = segA.atoms.select_atoms("resid 60")
    da = distance_array(r50.positions, r60.positions)
    return da.min()


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

u = mda.Universe(pdbfile)

print(timeit(distance, u, repeats=repeats))
