import sys
from shared import timeit

from schrodinger.structure import Structure
from schrodinger.structutils.rmsd import superimpose


def align(struc, other):
    idxs = list(range(1, len(struc.atom) + 1))
    superimpose(struc, idxs, other, idxs, use_symmetry=False)


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

struc = Structure.read(pdbfile)
other = struc.copy()

print(timeit(align, struc, other, repeats=repeats))
