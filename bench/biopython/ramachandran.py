import sys
from shared import timeit

from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral


def ramachandran(struc):
    phi_angles = []
    psi_angles = []
    residues = list(struc.get_residues())
    for i in range(1, len(residues) - 1):
        res = residues[i]
        res_prev = residues[i - 1]
        res_next = residues[i + 1]
        # Check residues have sequential residue numbers
        if (
            res.get_id()[1] == res_prev.get_id()[1] + 1
            and res_next.get_id()[1] == res.get_id()[1] + 1
        ):
            try:
                phi_angle = calc_dihedral(
                    res_prev["C"].get_vector(),
                    res["N"].get_vector(),
                    res["CA"].get_vector(),
                    res["C"].get_vector(),
                )
                psi_angle = calc_dihedral(
                    res["N"].get_vector(),
                    res["CA"].get_vector(),
                    res["C"].get_vector(),
                    res_next["N"].get_vector(),
                )
                phi_angles.append(phi_angle)
                psi_angles.append(psi_angle)
            except:
                pass
    return phi_angles, psi_angles


pdbfile = sys.argv[1]
repeats = int(sys.argv[2])

parser = PDBParser()
struc = parser.get_structure("", sys.argv[1])

print(timeit(ramachandran, struc, repeats=repeats))
