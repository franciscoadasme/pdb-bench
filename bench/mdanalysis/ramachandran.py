import sys
from shared import timeit

import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran


def ramachandran(u):
    Ramachandran(u.select_atoms("protein")).run()


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

u = mda.Universe(pdbfile)

print(timeit(ramachandran, u, repeats=repeats))
