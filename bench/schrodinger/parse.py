import sys
from shared import timeit

from schrodinger.structure import StructureReader


def parse(pdbfile):
    list(StructureReader(pdbfile))


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

print(timeit(parse, pdbfile, repeats=repeats))
