import sys
from shared import timeit

from Bio.PDB import PDBParser


def distance(st):
    min_dist = float("inf")
    for atom_a in st[0]["A"][50]:
        for atom_b in st[0]["A"][60]:
            if atom_a - atom_b < min_dist:
                min_dist = atom_a - atom_b
    return min_dist


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

parser = PDBParser()
st = parser.get_structure("", sys.argv[1])

print(timeit(distance, st, repeats=repeats))
