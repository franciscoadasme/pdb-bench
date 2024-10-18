import sys
from shared import timeit

from Bio.PDB import PDBParser


def distance(struc):
    min_dist = float("inf")
    for atom_a in struc[0]["A"][50]:
        for atom_b in struc[0]["A"][60]:
            if atom_a - atom_b < min_dist:
                min_dist = atom_a - atom_b
    return min_dist


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

parser = PDBParser()
struc = parser.get_structure("", sys.argv[1])

print(timeit(distance, struc, repeats=repeats))
