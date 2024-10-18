import sys
from shared import timeit

from Bio.PDB import PDBParser


def parse(filepath):
    parser = PDBParser()
    parser.get_structure("", filepath)


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

print(timeit(parse, pdbfile, repeats=repeats))
