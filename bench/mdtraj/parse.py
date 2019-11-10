import sys
from shared import timeit

import mdtraj as md


def parse(filepath):
    md.load(filepath)


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

print(timeit(parse, pdbfile, repeats=repeats))
