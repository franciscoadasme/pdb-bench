import sys
from shared import timeit

from schrodinger.structure import Structure


def count(struc):
    return sum(1 for r in struc.residue if r.pdbres == "ALA")


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

struc = Structure.read(pdbfile)

print(timeit(count, struc, repeats=repeats))
