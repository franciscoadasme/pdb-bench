import sys
from shared import timeit

from Bio.PDB import PDBParser, qcprot


def align(struc):
    qcp = qcprot.QCPSuperimposer()
    qcp.set_atoms(list(struc.get_atoms()), list(struc.get_atoms()))
    qcp.run()


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

parser = PDBParser()
struc = parser.get_structure("", sys.argv[1])

print(timeit(align, struc, repeats=repeats))
