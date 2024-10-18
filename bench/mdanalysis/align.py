import sys
from shared import timeit

import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj


def align(u):
    AlignTraj(u, u).run()


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

u = mda.Universe(pdbfile)

print(timeit(align, u, repeats=repeats))
