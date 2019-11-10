import sys
from shared import timeit

import MDAnalysis as mda


def parse(filepath):
    mda.Universe(filepath)


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

print(timeit(parse, pdbfile, repeats=repeats))
