import sys
from shared import timeit

from Bio.PDB import PDBParser


def count(struc):
    count = 0
    for res in struc.get_residues():
        if res.get_resname() == "ALA":
            count += 1


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

parser = PDBParser()
struc = parser.get_structure("", sys.argv[1])

print(timeit(count, struc, repeats=repeats))
