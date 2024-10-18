import sys
from shared import timeit

from schrodinger.structure import Structure
from schrodinger.structutils.measure import measure_dihedral_angle as dihedral


def ramachandran(struc):
    torsions = []
    for residue in struc.residue:
        try:
            atoms = residue.getDihedralAtoms("Phi")
            phi = dihedral(*atoms)
        except ValueError:
            phi = float("nan")

        try:
            atoms = residue.getDihedralAtoms("Psi")
            psi = dihedral(*atoms)
        except ValueError:
            psi = float("nan")

        torsions.append((phi, psi))


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

struc = Structure.read(pdbfile)

print(timeit(ramachandran, struc, repeats=repeats))
