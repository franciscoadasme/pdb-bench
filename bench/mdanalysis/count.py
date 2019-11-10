import sys
from shared import timeit

import MDAnalysis as mda


def count(u):
    return (u.residues.resnames == "ALA").sum()


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

u = mda.Universe(pdbfile)

print(timeit(count, u, repeats=repeats))
